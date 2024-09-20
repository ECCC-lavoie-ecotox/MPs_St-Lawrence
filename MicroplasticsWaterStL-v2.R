################### -###
### Spatial trend of microplastics in the St. Lawrence
### Script created by Raphael Lavoie
### 05-11-2023
################## -###

# Dependancies installation
libs <- c("readxl","tidyr","plyr","dplyr","ggplot2","forcats","gridExtra","png","grid","ggsignif","ggpubr","RColorBrewer","lme4","lmerTest","emmeans","MuMIn","multcomp","vegan","leaflet")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# data handling
library(readxl)
library(tidyr)
library(plyr)
library(dplyr)
# graphs
library(ggplot2)
library(forcats)
library(gridExtra)
library(png)
library(grid)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)

# models
library(lme4)
library(lmerTest) # forces lmer() to give p-values. The package provides anova function, that gives data frame similar to what gives lme4 package but with p-values calculated from F statistics of types I - III hypotheses.
library(emmeans) # mutliple comparison after lmer() requires latest version
library(MuMIn)
library(multcomp)
library(vegan)

# maps
library(leaflet)


### Data extraction ####
read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames but if you like tidyverse tibbles (the default with read_excel) then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if (!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

mysheets <- read_excel_allsheets("Z:/01-Projets et suivis/Langlois/Microplastics/10. Microplastics/2. Counting/Results2.xlsx")

Micro <- mysheets$TransposeR
colnames(Micro)
remove(mysheets)

## Use transposed DF and pivot it to have 1 column for count
Microp <- Micro %>%
  pivot_longer(
    cols = c(7:23),
    names_to = "Polymer",
    values_to = "Count"
  )


Microp <- data.frame(Microp)
Microp$PPML1 <- (Microp$Count / Microp$DebitL) * 1000000 * 2 # calculate number of particles per million of litres

# write.csv(Microp, "MicroplasticsLong.csv")

plot(Microp$Count ~ Microp$PPML1)
M <- lm(Microp$Count ~ Microp$PPML1)
plot(M)
cooks.distance(M)[which.max(cooks.distance(M))]
plot(M, which = 4)
remove(M)
# an outliers is obvious with value of ~1400 PPML
# Site	NetType	NetTypeRep	Rep	DebitL	PartType	Polymer	Count	PPML1
# TroisPistoles	Polym	Polymer-2	2	38090	Sphere	Polypropylene	26	1365.187713
# this is due to low volume filtered and high count
# Coock's distance is 44, so justifiaqble to remove!

Microp2 <- Microp[-1006, ] # the outlier was removed for the rest of analyses
plot(Microp2$Count ~ Microp2$PPML1)

MP <- Microp2
remove(Microp, Microp2)


## Stacked graph ####

# Mean of replicates (group_by(Site, PartType, Polymer, NetType)), then, sum of all PartType group_by(Site, Polymer, NetType)
AvRepSumPart <- MP %>%
  group_by(Site, PartType, Polymer, NetType) %>%
  summarise(
    Count = mean(Count, na.rm = T),
    PPML1 = mean(PPML1, na.rm = T)
  ) %>%
  group_by(Site, Polymer, NetType) %>%
  summarise(
    Count = sum(Count, na.rm = T),
    PPML1 = sum(PPML1, na.rm = T)
  )
View(AvRepSumPart)


# ggplot(AvRepSumPart, aes(fill=Polymer, y=PPML1, x=Site)) +
#  geom_bar(position="stack", stat="identity")

# New facet label names for net size variable
Netlabs <- c("A. 300 um", "B. 100 um")
names(Netlabs) <- c("Manta", "Polym")

pl1 <- AvRepSumPart %>%
  mutate(Polymer = fct_reorder(Polymer, PPML1)) %>%
  ggplot(aes(fill = Polymer, y = PPML1, x = Site)) +
  theme_classic() +
  geom_bar(position = "fill", stat = "identity") +
  facet_wrap(~NetType) + # panels by net type
  scale_x_discrete(labels = c("Varennes", "Contrecoeur", "Sorel", "Trois-Rivières", "Portneuf", "Québec", "Montmagny", "Baie-Saint-Paul", "La Pocatière", "Malbaie", "Trois-Pistoles")) +
  labs(x = "Site", y = "Plastic particles per million litres (PPML)", fill = "Polymer name") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.99)) +
  # theme(strip.text.x = element_blank())+
  theme(strip.background = element_rect(color = "white", fill = "white", size = 1.5, linetype = "solid")) +
  facet_grid(~NetType, labeller = labeller(NetType = Netlabs))

pl1

ggsave(
  filename = "Figs/Stack-ByNet_bySite_byPoly.png", # Nommez le fichier dans lequel vous voulez enregistrer, ajoutez l'extension du format de fichier que vous voulez utiliser (ex. pdf).
  plot = pl1, # Fournir le nom de l'objet plot dans R
  height = 8.5 * 0.75, # Fournir les dimensions voulues
  width = 11 * 0.75,
  units = "in",
  dpi = 1500
)


# Mean of replicates (group_by(Site, PartType, Polymer, NetType)), then, sum of all Polymer group_by(Site, PartType, NetType)
AvRepSumPoly <- MP %>%
  group_by(Site, PartType, Polymer, NetType) %>%
  summarise(
    Count = mean(Count, na.rm = T),
    PPML1 = mean(PPML1, na.rm = T)
  ) %>%
  group_by(Site, PartType, NetType) %>%
  summarise(
    Count = sum(Count, na.rm = T),
    PPML1 = sum(PPML1, na.rm = T)
  )
# View(AvRepSumPoly)


# New facet label names for net size variable
Netlabs <- c("A. 300 um", "B. 100 um")
names(Netlabs) <- c("Manta", "Polym")

pl2 <- AvRepSumPoly %>%
  mutate(PartType = fct_reorder(PartType, PPML1)) %>%
  ggplot(aes(fill = PartType, y = PPML1, x = Site)) +
  theme_classic() +
  geom_bar(position = "fill", stat = "identity") +
  facet_wrap(~NetType) + # panels by net type
  scale_x_discrete(labels = c("Varennes", "Contrecoeur", "Sorel", "Trois-Rivières", "Portneuf", "Québec", "Montmagny", "Baie-Saint-Paul", "La Pocatière", "Malbaie", "Trois-Pistoles")) +
  labs(x = "Site", y = "Plastic particles per million litres (PPML)", fill = "Shape") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.99)) +
  theme(strip.background = element_rect(color = "white", fill = "white", size = 1.5, linetype = "solid")) +
  facet_grid(~NetType, labeller = labeller(NetType = Netlabs))

pl2

ggsave(
  filename = "Figs/Stack-ByNet_bySite_byFiber.png", # Nommez le fichier dans lequel vous voulez enregistrer, ajoutez l'extension du format de fichier que vous voulez utiliser (ex. pdf).
  plot = pl2, # Fournir le nom de l'objet plot dans R
  height = 8.5 * 0.5, # Fournir les dimensions voulues
  width = 11 * 0.75,
  units = "in",
  dpi = 1500
)
#


#### Models ####

SumPerRep <- MP %>%
  group_by(Site, NetType, NetTypeRep) %>%
  summarise(
    Count = sum(Count, na.rm = T),
    PPML1 = sum(PPML1, na.rm = T)
  )
SumPerRep$SiteRepID <- paste0(SumPerRep$Site, SumPerRep$NetType)

SumPerRep <- SumPerRep %>%
  mutate(UD = case_when(
    Site == "AVarennes" ~ "UP",
    Site == "BContrecœur" ~ "UP",
    Site == "CSorel" ~ "UP",
    Site == "DTroisRivieres" ~ "UP",
    Site == "EPortneuf" ~ "UP",
    Site == "FQuebec" ~ "UP",
    Site == "GMontmagny" ~ "DOWN",
    Site == "HBaieSaintPaul" ~ "DOWN",
    Site == "ILaPocatiere" ~ "DOWN",
    Site == "JMalbaie" ~ "DOWN",
    Site == "KTroisPistoles" ~ "DOWN"
  ))

SumPerRep$Site <- as.factor(SumPerRep$Site)
SumPerRep$NetType <- as.factor(SumPerRep$NetType)
SumPerRep$NetTypeRep <- as.factor(SumPerRep$NetTypeRep)
SumPerRep$SiteRepID <- as.factor(SumPerRep$SiteRepID)
SumPerRep$UD <- as.factor(SumPerRep$UD)

SumPerRep$LPPML1 <- log10(SumPerRep$PPML1)

hist(sqrt(SumPerRep$PPML1), xlab = "PPML", main = "")

#

#### //##// Sum polymer - replicates without lmer() ####

# justification: replicates can be seen as true replicates (i.e., several birds sampled at a site). For mixed effect models, see lmer() in other script

dat <- SumPerRep

rm(M0, M1, M2, M3, M4, M5, M6, M7)

M0 <- lm(LPPML1 ~ 1, data = dat)
M1 <- lm(LPPML1 ~ Site, data = dat)
M2 <- lm(LPPML1 ~ NetType, data = dat)
M3 <- lm(LPPML1 ~ Site + NetType, data = dat)
M4 <- lm(LPPML1 ~ UD, data = dat)
M5 <- lm(LPPML1 ~ UD + Site, data = dat)
M6 <- lm(LPPML1 ~ UD + Site + NetType, data = dat)

AIC.table <- MuMIn::model.sel(M0, M1, M2, M3, M4, M5, M6)
# (AIC.table <- AIC.table[ , c("df", "logLik", "AICc", "delta")])
AICdf <- data.frame(AIC.table)
AICdf <- AICdf[-c(1, 6)]
write.csv(AICdf, "Figs/AIC_lm_bothNet.csv")

# M4, 0, 2,1
anova(M4) # 1
summary(M4)
t.test(LPPML1 ~ UD, data = dat)
plot(dat$LPPML1 ~ dat$UD)
anova(M2) # 3 - net
summary(M2)

anova(M1) # 3 - site
summary(M1)

# site_comparisons <- glht(FM, mcp(Site = "Tukey"))
# summary(site_comparisons)

mod <- M1

tuk <- glht(mod, linfct = mcp(Site = "Tukey"))
summary(tuk)
tuk.cld <- cld(tuk, decreasing = F) # decreasing=T to get higher values as "a"
tuk.cld
plot(tuk.cld)
TUK <- data.frame(tuk.cld[["mcletters"]][["monospacedLetters"]]) # to include letters in plot
TUK$tuk <- TUK$tuk.cld...mcletters......monospacedLetters...
TUK <- TUK[-1]


shapiro.test(resid(mod))
hist(resid(mod), breaks = 10)
plot(abs(resid(mod)) ~ fitted(mod)) # What about independence between residual and fitted values?
cor.test(abs(resid(mod)), fitted(mod))
summary(lm(abs(resid(mod)) ~ dat$Site)) # Levene test
boxplot(abs(resid(mod)) ~ Site, data = dat)

# Site not in best model, but significant when lm(~site). Respects condition of parametric test when log10

pb1 <- ggplot(dat, aes(x = Site, y = PPML1)) +
  theme_classic() +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.99)) +
  scale_x_discrete(labels = c("Varennes", "Contrecoeur", "Sorel", "Trois-Rivières", "Portneuf", "Québec", "Montmagny", "Baie-Saint-Paul", "La Pocatière", "Malbaie", "Trois-Pistoles")) +
  labs(x = "Site", y = "Plastic particles per million litres (PPML)", fill = "Polymer name") +
  annotate("text", y = c(1500:1500), x = c(1:11), label = TUK$tuk) +
  geom_vline(xintercept = 6.5, linetype = 2)

pb1

ggsave(
  filename = "Figs/BoxTukey-bySite.png", # Nommez le fichier dans lequel vous voulez enregistrer, ajoutez l'extension du format de fichier que vous voulez utiliser (ex. pdf).
  plot = pb1, # Fournir le nom de l'objet plot dans R
  height = 8.5 * 0.75, # Fournir les dimensions voulues
  width = 11 * 0.75,
  units = "in",
  dpi = 1500
)

####  //##// Sum Fiber ####

# keep replicates to include as a random variable (values comparable to manual technique on Excel)
SumPerRepPart <- MP %>%
  group_by(Site, NetType, NetTypeRep, PartType) %>%
  summarise(
    Count = sum(Count, na.rm = T),
    PPML1 = sum(PPML1, na.rm = T)
  )
SumPerRepPart$SiteRepID <- paste0(SumPerRepPart$Site, SumPerRepPart$NetType)

SumPerRepPart <- SumPerRepPart %>%
  mutate(UD = case_when(
    Site == "AVarennes" ~ "UP",
    Site == "BContrecœur" ~ "UP",
    Site == "CSorel" ~ "UP",
    Site == "DTroisRivieres" ~ "UP",
    Site == "EPortneuf" ~ "UP",
    Site == "FQuebec" ~ "UP",
    Site == "GMontmagny" ~ "DOWN",
    Site == "HBaieSaintPaul" ~ "DOWN",
    Site == "ILaPocatiere" ~ "DOWN",
    Site == "JMalbaie" ~ "DOWN",
    Site == "KTroisPistoles" ~ "DOWN"
  ))

SumPerRepPart$Site <- as.factor(SumPerRepPart$Site)
SumPerRepPart$NetType <- as.factor(SumPerRepPart$NetType)
SumPerRepPart$NetTypeRep <- as.factor(SumPerRepPart$NetTypeRep)
SumPerRepPart$SiteRepID <- as.factor(SumPerRepPart$SiteRepID)
SumPerRepPart$LPPML1 <- log10(SumPerRepPart$PPML1)
SumPerRepPart$UD <- as.factor(SumPerRepPart$UD)

dat <- filter(SumPerRepPart, PartType == "Fiber")
dat <- filter(dat, PPML1 > 0)

rm(M0, M1, M2, M3, M4, M5, M6)

M0 <- lm(LPPML1 ~ 1, data = dat)
M1 <- lm(LPPML1 ~ Site, data = dat)
M2 <- lm(LPPML1 ~ NetType, data = dat)
M3 <- lm(LPPML1 ~ Site + NetType, data = dat)
M4 <- lm(LPPML1 ~ UD, data = dat)
M5 <- lm(LPPML1 ~ UD + Site, data = dat)
M6 <- lm(LPPML1 ~ UD + Site + NetType, data = dat)

AIC.table <- MuMIn::model.sel(M0, M1, M2, M3, M4, M5, M6)
# (AIC.table <- AIC.table[ , c("df", "logLik", "AICc", "delta")])
AICdf <- data.frame(AIC.table)
AICdf <- AICdf[-c(1, 6)]
write.csv(AICdf, "Figs/AIC_lm_Fiber.csv")

anova(M0) # 1
summary(M0)
anova(M2) # 2
summary(M2)
anova(M1) # 3 - site
summary(M1)

# nothing significant

tuk <- glht(M2, linfct = mcp(NetType = "Tukey"))
summary(tuk)
tuk.cld <- cld(tuk, decreasing = T) # decreasing=T to get higher values as "a"
tuk.cld
plot(tuk.cld)

shapiro.test(resid(M2))
hist(resid(M2), breaks = 10)
plot(abs(resid(M2)) ~ fitted(M2)) # What about independence between residual and fitted values?
cor.test(abs(resid(M2)), fitted(M2))
X <- terms(M2) # this gets the variables in the parametric model
summary(lm(abs(resid(M2)) ~ dat$Site)) # Levene test
boxplot(abs(resid(M2)) ~ Site, data = dat)


pb3 <- ggplot(dat, aes(x = Site, y = PPML1)) +
  theme_classic() +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.99)) +
  scale_x_discrete(labels = c("Varennes", "Contrecoeur", "Sorel", "Trois-Rivières", "Portneuf", "Québec", "Montmagny", "Baie-Saint-Paul", "La Pocatière", "Malbaie", "Trois-Pistoles")) +
  labs(x = "Site", y = "Plastic particles per million litres (PPML)", fill = "Polymer name") +
  geom_vline(xintercept = 6.5, linetype = 2)

pb3

ggsave(
  filename = "Figs/BoxTukey-bySite-Fibre.png", # Nommez le fichier dans lequel vous voulez enregistrer, ajoutez l'extension du format de fichier que vous voulez utiliser (ex. pdf).
  plot = pb3, # Fournir le nom de l'objet plot dans R
  height = 8.5 * 0.75, # Fournir les dimensions voulues
  width = 11 * 0.75,
  units = "in",
  dpi = 1500
)


####  //##// Sum Fragments ####

# see previous dataset SumPerRepPart to subset for Fragment
dat <- filter(SumPerRepPart, PartType == "Fragment")
dat <- filter(dat, PPML1 > 0)

rm(M0, M1, M2, M3, M4, M5, M6)

M0 <- lm(LPPML1 ~ 1, data = dat)
M1 <- lm(LPPML1 ~ Site, data = dat)
M2 <- lm(LPPML1 ~ NetType, data = dat)
M3 <- lm(LPPML1 ~ Site + NetType, data = dat)
M4 <- lm(LPPML1 ~ UD, data = dat)
M5 <- lm(LPPML1 ~ UD + Site, data = dat)
M6 <- lm(LPPML1 ~ UD + Site + NetType, data = dat)

AIC.table <- MuMIn::model.sel(M0, M1, M2, M3, M4, M5, M6)
# (AIC.table <- AIC.table[ , c("df", "logLik", "AICc", "delta")])
AICdf <- data.frame(AIC.table)
AICdf <- AICdf[-c(1, 6)]
write.csv(AICdf, "Figs/AIC_lm_Fragment.csv")

anova(M4) # 1
summary(M4)
Tt <- t.test(PPML1 ~ UD, data = dat)
plot(dat$LPPML1 ~ dat$UD)

anova(M5) # 2
summary(M5)
anova(M1) # 3
summary(M1)



# site_comparisons <- glht(FM, mcp(Site = "Tukey"))
# summary(site_comparisons)

mod <- M1

tuk <- glht(mod, linfct = mcp(Site = "Tukey"))
summary(tuk)
tuk.cld <- cld(tuk, level = 0.05, decreasing = T) # decreasing=T to get higher values as "a"
tuk.cld
plot(tuk.cld)
TUK <- data.frame(tuk.cld[["mcletters"]][["monospacedLetters"]]) # to include letters in plot
TUK$tuk <- TUK$tuk.cld...mcletters......monospacedLetters...
TUK <- TUK[-1]

shapiro.test(resid(mod))
hist(resid(mod), breaks = 10)
plot(abs(resid(mod)) ~ fitted(mod)) # What about independence between residual and fitted values?
cor.test(abs(resid(mod)), fitted(mod))
X <- terms(mod) # this gets the variables in the parametric model
summary(lm(abs(resid(mod)) ~ dat$Site)) # Levene test
boxplot(abs(resid(mod)) ~ Site, data = dat)

# Site IS best model (and dAIC is >2) and IS significant when lm(~site). Respects condition of parametric test when log10 (but not if not log)


pb2 <- ggplot(dat, aes(x = Site, y = PPML1)) +
  theme_classic() +
  geom_boxplot() +
  scale_x_discrete(labels = c("Varennes", "Contrecoeur", "Sorel", "Trois-Rivières", "Portneuf", "Québec", "Montmagny", "Baie-Saint-Paul", "La Pocatière", "Malbaie", "Trois-Pistoles")) +
  labs(x = "Site", y = "Plastic particles per million litres (PPML)", fill = "Polymer name") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.99)) +
  annotate("text", y = c(800:800), x = c(1:11), label = TUK$tuk) +
  annotate("text", y = c(850:850), x = c(3.5, 9), label = c("A", "B")) +
  geom_vline(xintercept = 6.5, linetype = 2)

pb2

ggsave(
  filename = "Figs/BoxTukey-bySite-Fragment.png", # Nommez le fichier dans lequel vous voulez enregistrer, ajoutez l'extension du format de fichier que vous voulez utiliser (ex. pdf).
  plot = pb2, # Fournir le nom de l'objet plot dans R
  height = 8.5 * 0.75, # Fournir les dimensions voulues
  width = 11 * 0.75,
  units = "in",
  dpi = 1500
)

# Combine 3 figures together

F1 <- grid.arrange(arrangeGrob(pb1, left = textGrob("A.", x = unit(1, "npc"), y = unit(.95, "npc"))),
  arrangeGrob(pb3, left = textGrob("B.", x = unit(1, "npc"), y = unit(.95, "npc"))),
  arrangeGrob(pb2, left = textGrob("C.", x = unit(1, "npc"), y = unit(.95, "npc"))),
  layout_matrix = rbind(
    c(1, 1, 2),
    c(1, 1, 3)
  )
)
F1
ggsave(
  filename = "Figs/BoxTukey-3plots.png", # Nommez le fichier dans lequel vous voulez enregistrer, ajoutez l'extension du format de fichier que vous voulez utiliser (ex. pdf).
  plot = F1, # Fournir le nom de l'objet plot dans R
  height = 8.5, # Fournir les dimensions voulues
  width = 11,
  units = "in",
  dpi = 1500
)


## Non-metric-multidimensional-scaling (nMDR) ####

# subset for each net type and restructure for NMDS


NMDSdf <- MP %>%
  group_by(Site, NetType, NetTypeRep, Polymer) %>%
  summarise(
    Count = sum(Count, na.rm = T),
    PPML1 = sum(PPML1, na.rm = T)
  )

NMDSdf$SiteRepID <- paste0(NMDSdf$Site, NMDSdf$NetTypeRep)
NMDSdf <- NMDSdf %>%
  mutate(UD = case_when(
    Site == "AVarennes" ~ "UP",
    Site == "BContrecœur" ~ "UP",
    Site == "CSorel" ~ "UP",
    Site == "DTroisRivieres" ~ "UP",
    Site == "EPortneuf" ~ "UP",
    Site == "FQuebec" ~ "UP",
    Site == "GMontmagny" ~ "DOWN",
    Site == "HBaieSaintPaul" ~ "DOWN",
    Site == "ILaPocatiere" ~ "DOWN",
    Site == "JMalbaie" ~ "DOWN",
    Site == "KTroisPistoles" ~ "DOWN"
  ))

allNMDS <- NMDSdf
ManNMDS <- filter(NMDSdf, NetType == "Manta")
PolNMDS <- filter(NMDSdf, NetType == "Polym")

## //##// All (Manta + Poly-mer) ####

# Pool Manta and Poly-mer because there is no difference with univariate stats

# use "allNMDS" dataset

# add matching abbreviations of polymers
PolNam <- read.csv("Z:/01-Projets et suivis/Langlois/Microplastics/10. Microplastics/Stats/PolymerAbbreviationNames.csv", sep = ";")
allNMDS$Polymer <- as.factor(allNMDS$Polymer)
levels(allNMDS$Polymer)

allNMDS <- merge(allNMDS, PolNam)
allNMDS <- arrange(allNMDS, Site, NetTypeRep, Polymer)

allNMDS$PPML1 <- as.integer(allNMDS$PPML1)
allNMDS <- subset(allNMDS, select = -c(5)) # remove count
allNMDS <- subset(allNMDS, select = -c(1)) # remove Polymer bcs it interacts with Abbrev

allNMDSw <- allNMDS %>%
  pivot_wider(
    names_from = "Abbrev",
    values_from = "PPML1",
    values_fill = 0
  )
allNMDSwFull <- allNMDSw
allNMDSw <- subset(allNMDSw, select = -c(1:5)) # remove all categorical variables


# Diversity indices

simpson <- diversity(allNMDSw, "simpson") # calculate Simpson's 1-D Index of Diversity for each site. # closer to 1 = greater diversity
simpson
shannon <- diversity(allNMDSw) # note that Shannon's is default
shannon # Typically ranges from 1.5 - 3.4, higher = more diverse

DATAA <- data.frame(Site = factor())
DATAA <- unique(allNMDS$SiteRepID)
DATAA <- data.frame(DATAA)
DATAA$Site <- substr(DATAA$DATAA, 1, 7)
DATAA <- DATAA %>%
  mutate(UD = case_when(
    Site == "AVarenn" ~ "UP",
    Site == "BContre" ~ "UP",
    Site == "CSorelM" ~ "UP",
    Site == "CSorelP" ~ "UP",
    Site == "DTroisR" ~ "UP",
    Site == "EPortne" ~ "UP",
    Site == "FQuebec" ~ "UP",
    Site == "GMontma" ~ "DOWN",
    Site == "HBaieSa" ~ "DOWN",
    Site == "ILaPoca" ~ "DOWN",
    Site == "JMalbai" ~ "DOWN",
    Site == "KTroisP" ~ "DOWN"
  ))

DATAA$UD <- as.factor(DATAA$UD)
DATAA$Site <- as.factor(DATAA$Site)

resultSh <- adonis(shannon ~ Site, data = DATAA)
resultShUD <- adonis(shannon ~ UD, data = DATAA)
resultSi <- adonis(simpson ~ Site, data = DATAA)
resultSiUD <- adonis(simpson ~ UD, data = DATAA)


# Run the NMDS
allNMDS <- metaMDS(allNMDSw, distance = "bray", k = 2, trymax = 50)

# Extract the results
allNMDS

# Assess the goodness of fit and draw a Shepard plot
allNMDS$stress

# A good rule of thumb: stress < 0.05 provides an excellent representation in reduced dimensions, < 0.1 is good, < 0.2 is fair, and stress > 0.2 provides a poor representation.** To reiterate: high stress is bad, low stress is good!

stressplot(allNMDS, main = "Shepard plot")

DFFF <- data.frame(allNMDS[["points"]])
NMDSall <- cbind(DATAA, DFFF)
NMDSall$Nb <- c(1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11)

# Construct the biplot

# with all replicates

png(
  file = "Figs/nMDS_Both-rep.png",
  height = 5,
  width = 5, units = "in", res = 1500
)
treat <- c(rep("Flu", 32), rep("Est", 29)) # for each replicates
ordiplot(allNMDS, type = "n", main = paste("NMDS/Bray - Stress=", round(allNMDS$stress, 3)))
orditorp(allNMDS, display = "species", col = "darkorange", air = 0.01, cex = 1)
# ordiellipse(allNMDS,groups=treat,draw="polygon",col=c("steelblue","darkgreen"),label=F)
ordihull(allNMDS, groups = treat, draw = "polygon", col = c("steelblue", "darkgreen"), label = F, alpha = 0.25)
text(y = c(NMDSall$MDS2[which(NMDSall$UD == "UP")]), x = c(NMDSall$MDS1[which(NMDSall$UD == "UP")]), labels = NMDSall$Nb[which(NMDSall$UD == "UP")], col = "darkgreen")
text(y = c(NMDSall$MDS2[which(NMDSall$UD == "DOWN")]), x = c(NMDSall$MDS1[which(NMDSall$UD == "DOWN")]), labels = NMDSall$Nb[which(NMDSall$UD == "DOWN")], col = "darkblue")
dev.off()

## PERMANOVA ###

# with all replicates
MPPerm <- adonis2(allNMDSw ~ Site, data = NMDSall, permutation = 999, method = "bray")
MPPerm
# groups are significantly different

MPdist <- vegdist(allNMDSw, method = "bray")
dispersion <- betadisper(MPdist, group = NMDSall$Site)
permutest(dispersion)
plot(dispersion, hull = FALSE, ellipse = TRUE) ## sd ellipse

MPPerm <- adonis2(allNMDSw ~ UD, data = NMDSall, permutation = 999, method = "bray")
MPPerm

MPdist <- vegdist(allNMDSw, method = "bray")
dispersion <- betadisper(MPdist, group = NMDSall$UD)
permutest(dispersion)
plot(dispersion, hull = FALSE, ellipse = TRUE) ## sd ellipse


## //##// Manta ####

# use "ManNMDS" dataset

# add matching abbreviations of polymers
PolNam <- read.csv("Z:/01-Projets et suivis/Langlois/Microplastics/10. Microplastics/Stats/PolymerAbbreviationNames.csv", sep = ";")
ManNMDS$Polymer <- as.factor(ManNMDS$Polymer)
levels(ManNMDS$Polymer)

ManNMDS <- merge(ManNMDS, PolNam)
ManNMDS <- arrange(ManNMDS, Site, NetTypeRep, Polymer)

ManNMDS$PPML1 <- as.integer(ManNMDS$PPML1)
ManNMDS <- subset(ManNMDS, select = -c(5)) # remove count
ManNMDS <- subset(ManNMDS, select = -c(1)) # remove Polymer bcs it interacts with Abbrev

ManNMDSw <- ManNMDS %>%
  pivot_wider(
    names_from = "Abbrev",
    values_from = "PPML1",
    values_fill = 0
  )

ManNMDSwFull <- ManNMDSw
ManNMDSw <- subset(ManNMDSw, select = -c(1:5)) # remove all categorical variables

# Diversity indices

simpson <- diversity(ManNMDSw, "simpson") # calculate Simpson's 1-D Index of Diversity for each site. # closer to 1 = greater diversity
simpson
shannon <- diversity(ManNMDSw) # note that Shannon's is default
shannon # Typically ranges from 1.5 - 3.4, higher = more diverse


DATAA <- data.frame(Site = factor())
DATAA <- unique(ManNMDS$SiteRepID)
DATAA <- data.frame(DATAA)
DATAA$Site <- substr(DATAA$DATAA, 1, 7)
DATAA <- DATAA %>%
  mutate(UD = case_when(
    Site == "AVarenn" ~ "UP",
    Site == "BContre" ~ "UP",
    Site == "CSorelM" ~ "UP",
    Site == "DTroisR" ~ "UP",
    Site == "EPortne" ~ "UP",
    Site == "FQuebec" ~ "UP",
    Site == "GMontma" ~ "DOWN",
    Site == "HBaieSa" ~ "DOWN",
    Site == "ILaPoca" ~ "DOWN",
    Site == "JMalbai" ~ "DOWN",
    Site == "KTroisP" ~ "DOWN"
  ))

DATAA$UD <- as.factor(DATAA$UD)
DATAA$Site <- as.factor(DATAA$Site)

resultSh <- adonis(shannon ~ Site, data = DATAA)
resultShUD <- adonis(shannon ~ UD, data = DATAA)
resultSi <- adonis(simpson ~ Site, data = DATAA)
resultSiUD <- adonis(simpson ~ UD, data = DATAA)


# Run the NMDS model
MANnmds <- metaMDS(ManNMDSw, distance = "bray", k = 2, trymax = 50)

# Extract the results
MANnmds

# Assess the goodness of fit and draw a Shepard plot
MANnmds$stress

# A good rule of thumb: stress < 0.05 provides an excellent representation in reduced dimensions, < 0.1 is good, < 0.2 is fair, and stress > 0.2 provides a poor representation.** To reiterate: high stress is bad, low stress is good!

stressplot(MANnmds, main = "Shepard plot")

DFFF <- data.frame(MANnmds[["points"]])
NMDSmanta <- cbind(DATAA, DFFF)
NMDSmanta$Nb <- c(1, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10, 11, 11, 11)

# Construct the biplot

# with all replicates (for average per site, see other script)

png(
  file = "Z:/01-Projets et suivis/Langlois/Microplastics/10. Microplastics/Stats/Figs/nMDS_Manta-rep.png",
  height = 5,
  width = 5.5, units = "in", res = 1500
)
treat <- c(rep("Flu", 18), rep("Est", 14)) # for each replicates
ordiplot(MANnmds, type = "n", main = paste("NMDS/Bray - Stress=", round(MANnmds$stress, 3)))
orditorp(MANnmds, display = "species", col = "darkorange", air = 0.01, cex = 1)
# ordiellipse(MANnmds,groups=treat,draw="polygon",col=c("steelblue","darkgreen"),label=F)
ordihull(MANnmds, groups = treat, draw = "polygon", col = c("steelblue", "darkgreen"), alpha = 0.25, label = F)
text(y = c(NMDSmanta$MDS2[which(NMDSmanta$UD == "UP")]), x = c(NMDSmanta$MDS1[which(NMDSmanta$UD == "UP")]), labels = NMDSmanta$Nb[which(NMDSmanta$UD == "UP")], col = "darkgreen")
text(y = c(NMDSmanta$MDS2[which(NMDSmanta$UD == "DOWN")]), x = c(NMDSmanta$MDS1[which(NMDSmanta$UD == "DOWN")]), labels = NMDSmanta$Nb[which(NMDSmanta$UD == "DOWN")], col = "darkblue")

dev.off()

#


## PERMANOVA ###

# with all replicates
MPPerm <- adonis2(ManNMDSw ~ Site, data = NMDSmanta, permutation = 999, method = "bray")
MPPerm

MPdist <- vegdist(ManNMDSw, method = "bray")
dispersion <- betadisper(MPdist, group = NMDSmanta$Site)
permutest(dispersion)
plot(dispersion, hull = FALSE, ellipse = TRUE) ## sd ellipse

MPPerm <- adonis2(ManNMDSw ~ UD, data = NMDSmanta, permutation = 999, method = "bray")
MPPerm

MPdist <- vegdist(ManNMDSw, method = "bray")
dispersion <- betadisper(MPdist, group = NMDSmanta$UD)
permutest(dispersion)
plot(dispersion, hull = FALSE, ellipse = TRUE) ## sd ellipse


## //##// Poly-mer ####

# use "PolNMDS" dataset

# add matching abbreviations of polymers
PolNam <- read.csv("Z:/01-Projets et suivis/Langlois/Microplastics/10. Microplastics/Stats/PolymerAbbreviationNames.csv", sep = ";")
PolNMDS$Polymer <- as.factor(PolNMDS$Polymer)
levels(PolNMDS$Polymer)

PolNMDS <- merge(PolNMDS, PolNam)
PolNMDS <- arrange(PolNMDS, Site, NetTypeRep, Polymer)

PolNMDS$PPML1 <- as.integer(PolNMDS$PPML1)
PolNMDS <- subset(PolNMDS, select = -c(5)) # remove count
PolNMDS <- subset(PolNMDS, select = -c(1)) # remove Polymer bcs it interacts with Abbrev

PolNMDSw <- PolNMDS %>%
  pivot_wider(
    names_from = "Abbrev",
    values_from = "PPML1",
    values_fill = 0
  )

PolNMDSwFull <- PolNMDSw
PolNMDSw <- subset(PolNMDSw, select = -c(1:5)) # remove all categorical variables

# Diversity indices

simpson <- diversity(PolNMDSw, "simpson") # calculate Simpson's 1-D Index of Diversity for each site. # closer to 1 = greater diversity
simpson
shannon <- diversity(PolNMDSw) # note that Shannon's is default
shannon # Typically ranges from 1.5 - 3.4, higher = more diverse


DATAA <- data.frame(Site = factor())
DATAA <- unique(PolNMDS$SiteRepID)
DATAA <- data.frame(DATAA)
DATAA$Site <- substr(DATAA$DATAA, 1, 7)
DATAA <- DATAA %>%
  mutate(UD = case_when(
    Site == "AVarenn" ~ "UP",
    Site == "BContre" ~ "UP",
    Site == "CSorelM" ~ "UP",
    Site == "CSorelP" ~ "UP",
    Site == "DTroisR" ~ "UP",
    Site == "EPortne" ~ "UP",
    Site == "FQuebec" ~ "UP",
    Site == "GMontma" ~ "DOWN",
    Site == "HBaieSa" ~ "DOWN",
    Site == "ILaPoca" ~ "DOWN",
    Site == "JMalbai" ~ "DOWN",
    Site == "KTroisP" ~ "DOWN"
  ))

DATAA$UD <- as.factor(DATAA$UD)
DATAA$Site <- as.factor(DATAA$Site)

resultSh <- adonis(shannon ~ Site, data = DATAA)
resultShUD <- adonis(shannon ~ UD, data = DATAA)
resultSi <- adonis(simpson ~ Site, data = DATAA)
resultSiUD <- adonis(simpson ~ UD, data = DATAA)


# Run the NMDS model
PolNMDS <- metaMDS(PolNMDSw, distance = "bray", k = 2, trymax = 50)

# Extract the results
PolNMDS

# Assess the goodness of fit and draw a Shepard plot
PolNMDS$stress

# A good rule of thumb: stress < 0.05 provides an excellent representation in reduced dimensions, < 0.1 is good, < 0.2 is fair, and stress > 0.2 provides a poor representation.** To reiterate: high stress is bad, low stress is good!

stressplot(PolNMDS, main = "Shepard plot")

DFFF <- data.frame(PolNMDS[["points"]])
NMDSpoly <- cbind(DATAA, DFFF)
NMDSpoly$Nb <- c(2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 7, 7, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 11, 11, 11)

# Construct the biplot

# with all replicates (for average per site, see other script)

png(
  file = "Z:/01-Projets et suivis/Langlois/Microplastics/10. Microplastics/Stats/Figs/nMDS_Poly-rep.png",
  height = 5,
  width = 5, units = "in", res = 1500
)
treat <- c(rep("Flu", 14), rep("Est", 15)) # for each replicates
ordiplot(PolNMDS, type = "n", main = paste("NMDS/Bray - Stress=", round(PolNMDS$stress, 3)))
orditorp(PolNMDS, display = "species", col = "darkorange", air = 0.01, cex = 1)
# ordiellipse(PolNMDS,groups=treat,draw="polygon",col=c("steelblue","darkgreen"),label=F)
ordihull(PolNMDS, groups = treat, draw = "polygon", col = c("steelblue", "darkgreen"), alpha = 0.25, label = F)
text(y = c(NMDSpoly$MDS2[which(NMDSpoly$UD == "UP")]), x = c(NMDSpoly$MDS1[which(NMDSpoly$UD == "UP")]), labels = NMDSpoly$Nb[which(NMDSpoly$UD == "UP")], col = "darkgreen")
text(y = c(NMDSpoly$MDS2[which(NMDSpoly$UD == "DOWN")]), x = c(NMDSpoly$MDS1[which(NMDSpoly$UD == "DOWN")]), labels = NMDSpoly$Nb[which(NMDSpoly$UD == "DOWN")], col = "darkblue")

dev.off()
#

## PERMANOVA ###

# with all replicates
MPPerm <- adonis2(PolNMDSw ~ Site, data = NMDSpoly, permutation = 999, method = "bray")
MPPerm

MPdist <- vegdist(PolNMDSw, method = "bray")
dispersion <- betadisper(MPdist, group = NMDSpoly$Site)
permutest(dispersion)
plot(dispersion, hull = FALSE, ellipse = TRUE) ## sd ellipse

MPPerm <- adonis2(PolNMDSw ~ UD, data = NMDSpoly, permutation = 999, method = "bray")
MPPerm

MPdist <- vegdist(PolNMDSw, method = "bray")
dispersion <- betadisper(MPdist, group = NMDSpoly$UD)
permutest(dispersion)
plot(dispersion, hull = FALSE, ellipse = TRUE) ## sd ellipse
#

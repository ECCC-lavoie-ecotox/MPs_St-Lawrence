###########################################
### Spatial trend of microplastics in the St. Lawrence
### Script created by Raphael Lavoie
### 05-11-2023
##########################################

#####################
# Libraries management
#####################

# Install packages if missing
libs <- c("readxl","tidyr","dplyr","ggplot2","forcats","gridExtra","RColorBrewer","dichromat","lme4","multcomp","vegan")
new_libs <- libs[!(libs %in% installed.packages()[,"Package"])]
if(length(new_libs)) install.packages(new_libs)

# data handling
library(readxl)
library(tidyr)
library(dplyr)
# graphs
library(ggplot2)
library(forcats)
library(gridExtra)
library(RColorBrewer)
library(dichromat)
# models
library(lme4)
library(multcomp)
library(vegan)

#####################
### Data ####
#####################

# Load data with microplastic count
Micro <- readxl::read_xlsx("data/Results_2.xlsx", sheet = "TransposeR")

## Use transposed DF and pivot it to have 1 column for count
Microp <- Micro |>
  tidyr::pivot_longer(
    cols = c(Polyethylene:Silicone),
    names_to = "Polymer",
    values_to = "Count"
  ) |> 
  # calculate number of particles per million of litres  
  dplyr::mutate(PPML1 = (Count / DebitL) * 1000000 * 2) |> 
  as.data.frame()


# Outliers detection
#####################

plot(Microp$Count ~ Microp$PPML1)
M <- lm(Microp$Count ~ Microp$PPML1)
plot(M)
cooks.distance(M)[which.max(cooks.distance(M))]
plot(M, which = 4)
rm(M)


# an outliers is obvious with value of ~1400 PPML
# Site	NetType	NetTypeRep	Rep	DebitL	PartType	Polymer	Count	PPML1
# TroisPistoles	Polym	Polymer-2	2	38090	Sphere	Polypropylene	26	1365.187713
# this is due to low volume filtered and high count
# Coock's distance is 44, so justifiaqble to remove!
Microp2 <- Microp[-1006, ] # the outlier was removed for the rest of analyses
plot(Microp2$Count ~ Microp2$PPML1)

MP <- Microp2
rm(Microp, Microp2)

######################
### Stacked graph ####
######################

# Mean of replicates (group_by(Site, PartType, Polymer, NetType)), then, sum of all PartType group_by(Site, Polymer, NetType)
AvRepSumPart <- MP |>
  dplyr::group_by(Site, PartType, Polymer, NetType) |>
  dplyr::summarise(
    Count = mean(Count, na.rm = T),
    PPML1 = mean(PPML1, na.rm = T)
  ) |>
  dplyr::group_by(Site, Polymer, NetType) |>
  dplyr::summarise(
    Count = sum(Count, na.rm = T),
    PPML1 = sum(PPML1, na.rm = T)
  )

View(AvRepSumPart)

# Descriptives stats by sites
descSites <- MP |> dplyr::group_by(Site) |>
  dplyr::summarise(
    mean = mean(PPML1, na.rm = TRUE),
    sd = sd(PPML1, na.rm = TRUE),
    n = n()
  ) |>
  dplyr::arrange(desc(mean))

# New facet label names for net size variable
Netlabs <- c("A. 300 um", "B. 100 um")
names(Netlabs) <- c("Manta", "Polym")

pl1 <- AvRepSumPart |>
  dplyr::mutate(Polymer = forcats::fct_reorder(Polymer, PPML1)) |>
  ggplot(aes(fill = Polymer, y = PPML1, x = Site)) +
  scale_fill_manual(values=dichromat::colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(17)) +
  theme_classic() +
  geom_bar(position = "fill", stat = "identity", colour="white") +
  facet_wrap(~NetType) + # panels by net type
  scale_x_discrete(labels = c("Varennes", "Contrecoeur", "Sorel", "Trois-Rivières", "Portneuf", "Québec", "Montmagny", "Baie-Saint-Paul", "La Pocatière", "Malbaie", "Trois-Pistoles")) +
  labs(x = "Site", y = "Plastic particles per million litres (PPML)", fill = "Polymer name") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.99)) +
  # theme(strip.text.x = element_blank())+
  theme(strip.background = element_blank()) +
  facet_grid(~NetType, labeller = labeller(NetType = Netlabs))

pl1



ggsave(
  filename = "figs/Stack-ByNet_bySite_byPoly.png", # Nommez le fichier dans lequel vous voulez enregistrer, ajoutez l'extension du format de fichier que vous voulez utiliser (ex. pdf).
  plot = pl1, # Fournir le nom de l'objet plot dans R
  height = 8.5 * 0.75, # Fournir les dimensions voulues
  width = 11 * 0.75,
  units = "in",
  dpi = 300
)

# Mean of replicates (group_by(Site, PartType, Polymer, NetType)), then, sum of all Polymer group_by(Site, PartType, NetType)
AvRepSumPoly <- MP |>
  dplyr::group_by(Site, PartType, Polymer, NetType) |>
  dplyr::summarise(
    Count = mean(Count, na.rm = T),
    PPML1 = mean(PPML1, na.rm = T)
  ) |>
  dplyr::group_by(Site, PartType, NetType) |>
  dplyr::summarise(
    Count = sum(Count, na.rm = T),
    PPML1 = sum(PPML1, na.rm = T)
  )

View(AvRepSumPoly)

# New facet label names for net size variable
Netlabs <- c("A. 300 um", "B. 100 um")
names(Netlabs) <- c("Manta", "Polym")

pl2 <- AvRepSumPoly |>
  dplyr::mutate(PartType = forcats::fct_reorder(PartType, PPML1)) |>
  ggplot(aes(fill = PartType, y = PPML1, x = Site)) +
  theme_classic() +
  scale_fill_manual(values=RColorBrewer::brewer.pal(3,"Paired")) +
  geom_bar(position = "fill", stat = "identity", colour="white") +
  facet_wrap(~NetType) + # panels by net type
  scale_x_discrete(labels = c("Varennes", "Contrecoeur", "Sorel", "Trois-Rivières", "Portneuf", "Québec", "Montmagny", "Baie-Saint-Paul", "La Pocatière", "Malbaie", "Trois-Pistoles")) +
  labs(x = "Site", y = "Plastic particles per million litres (PPML)", fill = "Shape") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.99)) +
  theme(strip.background = element_blank()) +
  facet_grid(~NetType, labeller = labeller(NetType = Netlabs))

pl2

ggsave(
  filename = "figs/Stack-ByNet_bySite_byFiber.png", # Nommez le fichier dans lequel vous voulez enregistrer, ajoutez l'extension du format de fichier que vous voulez utiliser (ex. pdf).
  plot = pl2, # Fournir le nom de l'objet plot dans R
  height = 8.5 * 0.5, # Fournir les dimensions voulues
  width = 11 * 0.75,
  units = "in",
  dpi = 300
)

######################
#### Models ####
######################

SumPerRep <- MP |>
  dplyr::group_by(Site, NetType, NetTypeRep) |>
  dplyr::summarise(
    Count = sum(Count, na.rm = T),
    PPML1 = sum(PPML1, na.rm = T)
  ) |> 
  dplyr::mutate(SiteRepID = paste0(Site, NetType))

SumPerRep <- SumPerRep |>
  dplyr::mutate(UD = dplyr::case_when(
    Site %in% c("AVarennes","BContrecœur","CSorel","DTroisRivieres","EPortneuf","FQuebec") ~ "UP",
    Site %in% c("GMontmagny", "HBaieSaintPaul", "ILaPocatiere", "JMalbaie", "KTroisPistoles") ~ "DOWN"
  )) |> dplyr::mutate_at(vars("Site", "NetType", "NetTypeRep", "SiteRepID", "UD"), as.factor) |>
  dplyr::mutate(LPPML1 = log10(PPML1))

hist(sqrt(SumPerRep$PPML1), xlab = "PPML", main = "")

##############################################################
#### //##// Sum polymer - replicates without lmer() ####
##############################################################

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
AICdf <- data.frame(AIC.table) |> dplyr::select(-X.Intercept., -weight)
write.csv(AICdf, "figs/AIC_lm_bothNet.csv")

# M4, 0, 2,1
anova(M4) # 1
summary(M4)
t.test(LPPML1 ~ UD, data = dat)
plot(dat$LPPML1 ~ dat$UD)
anova(M2) # 3 - net
summary(M2)
anova(M1) # 3 - site
summary(M1)

# Selected model
###################
mod <- M1
tuk <- multcomp::glht(mod, linfct = multcomp::mcp(Site = "Tukey"))
summary(tuk)
tuk.cld <- multcomp::cld(tuk, decreasing = F) # decreasing=T to get higher values as "a"
tuk.cld
plot(tuk.cld)
TUK <- data.frame(tuk.cld[["mcletters"]][["monospacedLetters"]]) # to include letters in plot
colnames(TUK) <- c("tuk")

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
  filename = "figs/BoxTukey-bySite.png", # Nommez le fichier dans lequel vous voulez enregistrer, ajoutez l'extension du format de fichier que vous voulez utiliser (ex. pdf).
  plot = pb1, # Fournir le nom de l'objet plot dans R
  height = 8.5 * 0.75, # Fournir les dimensions voulues
  width = 11 * 0.75,
  units = "in",
  dpi = 300
)

###############################
####  //##// Sum Fiber ####
###############################

# keep replicates to include as a random variable (values comparable to manual technique on Excel)
SumPerRepPart <- MP |>
  dplyr::group_by(Site, NetType, NetTypeRep, PartType) |>
  dplyr::summarise(
    Count = sum(Count, na.rm = T),
    PPML1 = sum(PPML1, na.rm = T)
  ) |> dplyr::mutate(SiteRepID = paste0(Site, NetType))

SumPerRepPart <- SumPerRepPart |>
  dplyr::mutate(UD = dplyr::case_when(
    Site %in% c("AVarennes","BContrecœur","CSorel","DTroisRivieres","EPortneuf","FQuebec") ~ "UP",
    Site %in% c("GMontmagny", "HBaieSaintPaul", "ILaPocatiere", "JMalbaie", "KTroisPistoles") ~ "DOWN"
  ))|> dplyr::mutate_at(vars("Site", "NetType", "NetTypeRep", "SiteRepID", "UD"), as.factor) |>
  dplyr::mutate(LPPML1 = log10(PPML1))

dat <- dplyr::filter(SumPerRepPart, PartType == "Fiber")
dat <- dplyr::filter(dat, PPML1 > 0)

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
AICdf <- data.frame(AIC.table) |> dplyr::select(-X.Intercept., -weight)
write.csv(AICdf, "figs/AIC_lm_Fiber.csv")

anova(M0) # 1
summary(M0)
anova(M2) # 2
summary(M2)
anova(M1) # 3 - site
summary(M1)

# nothing significant

tuk <- multcomp::glht(M2, linfct = multcomp::mcp(NetType = "Tukey"))
summary(tuk)
tuk.cld <- multcomp::cld(tuk, decreasing = T) # decreasing=T to get higher values as "a"
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
  filename = "figs/BoxTukey-bySite-Fibre.png", # Nommez le fichier dans lequel vous voulez enregistrer, ajoutez l'extension du format de fichier que vous voulez utiliser (ex. pdf).
  plot = pb3, # Fournir le nom de l'objet plot dans R
  height = 8.5 * 0.75, # Fournir les dimensions voulues
  width = 11 * 0.75,
  units = "in",
  dpi = 300
)

###############################
####  //##// Sum Fragments ####
###############################

# see previous dataset SumPerRepPart to subset for Fragment
dat <- dplyr::filter(SumPerRepPart, PartType == "Fragment")
dat <- dplyr::filter(dat, PPML1 > 0)

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
AICdf <- data.frame(AIC.table) |> dplyr::select(-X.Intercept., -weight)
write.csv(AICdf, "figs/AIC_lm_Fragment.csv")

anova(M4) # 1
summary(M4)
Tt <- t.test(PPML1 ~ UD, data = dat)
plot(dat$LPPML1 ~ dat$UD)

anova(M5) # 2
summary(M5)
anova(M1) # 3
summary(M1)

# Selected model
###################
mod <- M1

tuk <- multcomp::glht(mod, linfct = multcomp::mcp(Site = "Tukey"))
summary(tuk)
tuk.cld <- multcomp::cld(tuk, level = 0.05, decreasing = T) # decreasing=T to get higher values as "a"
tuk.cld
plot(tuk.cld)
TUK <- data.frame(tuk.cld[["mcletters"]][["monospacedLetters"]]) # to include letters in plot
colnames(TUK) <- c("tuk")

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
  filename = "figs/BoxTukey-bySite-Fragment.png", # Nommez le fichier dans lequel vous voulez enregistrer, ajoutez l'extension du format de fichier que vous voulez utiliser (ex. pdf).
  plot = pb2, # Fournir le nom de l'objet plot dans R
  height = 8.5 * 0.75, # Fournir les dimensions voulues
  width = 11 * 0.75,
  units = "in",
  dpi = 300
)

# Combine 3 figures together
F1 <- gridExtra::grid.arrange(gridExtra::arrangeGrob(pb1, left = grid::textGrob("A.", x = unit(1, "npc"), y = unit(.95, "npc"))),
  gridExtra::arrangeGrob(pb3, left = grid::textGrob("B.", x = unit(1, "npc"), y = unit(.95, "npc"))),
  gridExtra::arrangeGrob(pb2, left = grid::textGrob("C.", x = unit(1, "npc"), y = unit(.95, "npc"))),
  layout_matrix = rbind(
    c(1, 1, 2),
    c(1, 1, 3)
  )
)

F1
ggsave(
  filename = "figs/BoxTukey-3plots.png", # Nommez le fichier dans lequel vous voulez enregistrer, ajoutez l'extension du format de fichier que vous voulez utiliser (ex. pdf).
  plot = F1, # Fournir le nom de l'objet plot dans R
  height = 8.5, # Fournir les dimensions voulues
  width = 11,
  units = "in",
  dpi = 300
)

## Non-metric-multidimensional-scaling (nMDR) ####
# subset for each net type and restructure for NMDS
NMDSdf <- MP |>
  dplyr::group_by(Site, NetType, NetTypeRep, Polymer) |>
  dplyr::summarise(
    Count = sum(Count, na.rm = T),
    PPML1 = sum(PPML1, na.rm = T)
  ) |> dplyr::mutate(SiteRepID = paste0(Site, NetTypeRep))

NMDSdf <- NMDSdf |>
  dplyr::mutate(UD = dplyr::case_when(
    Site %in% c("AVarennes","BContrecœur","CSorel","DTroisRivieres","EPortneuf","FQuebec") ~ "UP",
    Site %in% c("GMontmagny", "HBaieSaintPaul", "ILaPocatiere", "JMalbaie", "KTroisPistoles") ~ "DOWN"
  ))

allNMDS <- NMDSdf
ManNMDS <- dplyr::filter(NMDSdf, NetType == "Manta")
PolNMDS <- dplyr::filter(NMDSdf, NetType == "Polym")

#####################################
## //##// All (Manta + Poly-mer) ####
#####################################

# Pool Manta and Poly-mer because there is no difference with univariate stats
# use "allNMDS" dataset

# add matching abbreviations of polymers
PolNam <- read.csv("data/PolymerAbbreviationNames.csv", sep = ";")
allNMDS$Polymer <- as.factor(allNMDS$Polymer)
levels(allNMDS$Polymer)

allNMDS <- merge(allNMDS, PolNam) |>
  dplyr::arrange(Site, NetTypeRep, Polymer) |>
  dplyr::mutate(PPML1 = as.integer(PPML1)) |>
  dplyr::select(-Count, -Polymer)

allNMDSw <- allNMDS |>
  tidyr::pivot_wider(
    names_from = "Abbrev",
    values_from = "PPML1",
    values_fill = 0
  )

allNMDSwFull <- allNMDSw
allNMDSw <- dplyr::select(allNMDSw, -c(Site:UD)) # remove all categorical variables

# Diversity indices
simpson <- vegan::diversity(allNMDSw, "simpson") # calculate Simpson's 1-D Index of Diversity for each site. # closer to 1 = greater diversity
simpson
shannon <- vegan::diversity(allNMDSw) # note that Shannon's is default
shannon # Typically ranges from 1.5 - 3.4, higher = more diverse

DATAA <- data.frame(Site = unique(allNMDS$SiteRepID)) |>
  dplyr::mutate(UD = dplyr::case_when(
    stringr::str_detect(Site, "AVarennes|BContrecœur|CSorel|DTroisRivieres|EPortneuf|FQuebec") ~ "UP",
    stringr::str_detect(Site, "GMontmagny|HBaieSaintPaul|ILaPocatiere|JMalbaie|KTroisPistoles") ~ "DOWN"
  )) |> dplyr::mutate_at(vars("UD", "Site"), as.factor)

resultSh <- adonis2(shannon ~ Site, data = DATAA)
resultShUD <- adonis2(shannon ~ UD, data = DATAA)
resultSi <- adonis2(simpson ~ Site, data = DATAA)
resultSiUD <- adonis2(simpson ~ UD, data = DATAA)


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
  file = "figs/nMDS_Both-rep.png",
  height = 5,
  width = 5, units = "in", res = 300
)

treat <- c(rep("Flu", 32), rep("Est", 29)) # for each replicates
ordiplot(allNMDS, type = "n", main = paste("NMDS/Bray - Stress=", round(allNMDS$stress, 3)))
orditorp(allNMDS, display = "species", col = "darkorange", air = 0.01, cex = 1)
# ordiellipse(allNMDS,groups=treat,draw="polygon",col=c("steelblue","darkgreen"),label=F)
ordihull(allNMDS, groups = treat, draw = "polygon", col = c("steelblue", "darkgreen"), label = F, alpha = 0.25)
text(y = NMDSall$MDS2[which(NMDSall$UD == "UP")], x = NMDSall$MDS1[which(NMDSall$UD == "UP")], labels = NMDSall$Nb[which(NMDSall$UD == "UP")], col = "darkgreen")
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
PolNam <- read.csv("data/PolymerAbbreviationNames.csv", sep = ";")
head(ManNMDS)
ManNMDS$Polymer <- as.factor(ManNMDS$Polymer)
levels(ManNMDS$Polymer)

ManNMDS <- merge(ManNMDS, PolNam) |>
  dplyr::arrange(Site, NetTypeRep, Polymer) |>
  dplyr::mutate(PPML1 = as.integer(PPML1)) |>
  dplyr::select(-Count, -Polymer)

ManNMDSw <- ManNMDS |>
  tidyr::pivot_wider(
    names_from = "Abbrev",
    values_from = "PPML1",
    values_fill = 0
  )

ManNMDSwFull <- ManNMDSw
ManNMDSw <- dplyr::select(ManNMDSw, -c(Site:UD)) # remove all categorical variables

# Diversity indices
simpson <- diversity(ManNMDSw, "simpson") # calculate Simpson's 1-D Index of Diversity for each site. # closer to 1 = greater diversity
simpson
shannon <- diversity(ManNMDSw) # note that Shannon's is default
shannon # Typically ranges from 1.5 - 3.4, higher = more diverse

DATAA <- data.frame(Site = unique(ManNMDS$SiteRepID)) |>
  dplyr::mutate(UD = dplyr::case_when(
    stringr::str_detect(Site, "AVarennes|BContrecœur|CSorel|DTroisRivieres|EPortneuf|FQuebec") ~ "UP",
    stringr::str_detect(Site, "GMontmagny|HBaieSaintPaul|ILaPocatiere|JMalbaie|KTroisPistoles") ~ "DOWN"
  )) |> dplyr::mutate_at(vars("UD", "Site"), as.factor)

resultSh <- adonis2(shannon ~ Site, data = DATAA)
resultShUD <- adonis2(shannon ~ UD, data = DATAA)
resultSi <- adonis2(simpson ~ Site, data = DATAA)
resultSiUD <- adonis2(simpson ~ UD, data = DATAA)


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
  file = "figs/nMDS_Manta-rep.png",
  height = 5,
  width = 5.5, units = "in", res = 300
)
treat <- c(rep("Flu", 18), rep("Est", 14)) # for each replicates
ordiplot(MANnmds, type = "n", main = paste("NMDS/Bray - Stress=", round(MANnmds$stress, 3)))
orditorp(MANnmds, display = "species", col = "darkorange", air = 0.01, cex = 1)
# ordiellipse(MANnmds,groups=treat,draw="polygon",col=c("steelblue","darkgreen"),label=F)
ordihull(MANnmds, groups = treat, draw = "polygon", col = c("steelblue", "darkgreen"), alpha = 0.25, label = F)
text(y = c(NMDSmanta$MDS2[which(NMDSmanta$UD == "UP")]), x = c(NMDSmanta$MDS1[which(NMDSmanta$UD == "UP")]), labels = NMDSmanta$Nb[which(NMDSmanta$UD == "UP")], col = "darkgreen")
text(y = c(NMDSmanta$MDS2[which(NMDSmanta$UD == "DOWN")]), x = c(NMDSmanta$MDS1[which(NMDSmanta$UD == "DOWN")]), labels = NMDSmanta$Nb[which(NMDSmanta$UD == "DOWN")], col = "darkblue")

dev.off()

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
PolNam <- read.csv("data/PolymerAbbreviationNames.csv", sep = ";")
PolNMDS$Polymer <- as.factor(PolNMDS$Polymer)
levels(PolNMDS$Polymer)

PolNMDS <- merge(PolNMDS, PolNam) |>
  dplyr::arrange(Site, NetTypeRep, Polymer) |>
  dplyr::mutate(PPML1 = as.integer(PPML1)) |>
  dplyr::select(-Count, -Polymer)

PolNMDSw <- PolNMDS |>
  tidyr::pivot_wider(
    names_from = "Abbrev",
    values_from = "PPML1",
    values_fill = 0
  )

PolNMDSwFull <- PolNMDSw
PolNMDSw <- dplyr::select(PolNMDSw, -c(Site:UD)) # remove all categorical variables

# Diversity indices
simpson <- diversity(PolNMDSw, "simpson") # calculate Simpson's 1-D Index of Diversity for each site. # closer to 1 = greater diversity
simpson
shannon <- diversity(PolNMDSw) # note that Shannon's is default
shannon # Typically ranges from 1.5 - 3.4, higher = more diverse

DATAA <- data.frame(Site = unique(PolNMDS$SiteRepID)) |>
  dplyr::mutate(UD = dplyr::case_when(
    stringr::str_detect(Site, "AVarennes|BContrecœur|CSorel|DTroisRivieres|EPortneuf|FQuebec") ~ "UP",
    stringr::str_detect(Site, "GMontmagny|HBaieSaintPaul|ILaPocatiere|JMalbaie|KTroisPistoles") ~ "DOWN"
  )) |> dplyr::mutate_at(vars("UD", "Site"), as.factor)


resultSh <- adonis2(shannon ~ Site, data = DATAA)
resultShUD <- adonis2(shannon ~ UD, data = DATAA)
resultSi <- adonis2(simpson ~ Site, data = DATAA)
resultSiUD <- adonis2(simpson ~ UD, data = DATAA)

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
  file = "figs/nMDS_Poly-rep.png",
  height = 5,
  width = 5, units = "in", res = 300
)
treat <- c(rep("Flu", 14), rep("Est", 15)) # for each replicates
ordiplot(PolNMDS, type = "n", main = paste("NMDS/Bray - Stress=", round(PolNMDS$stress, 3)))
orditorp(PolNMDS, display = "species", col = "darkorange", air = 0.01, cex = 1)
# ordiellipse(PolNMDS,groups=treat,draw="polygon",col=c("steelblue","darkgreen"),label=F)
ordihull(PolNMDS, groups = treat, draw = "polygon", col = c("steelblue", "darkgreen"), alpha = 0.25, label = F)
text(y = c(NMDSpoly$MDS2[which(NMDSpoly$UD == "UP")]), x = c(NMDSpoly$MDS1[which(NMDSpoly$UD == "UP")]), labels = NMDSpoly$Nb[which(NMDSpoly$UD == "UP")], col = "darkgreen")
text(y = c(NMDSpoly$MDS2[which(NMDSpoly$UD == "DOWN")]), x = c(NMDSpoly$MDS1[which(NMDSpoly$UD == "DOWN")]), labels = NMDSpoly$Nb[which(NMDSpoly$UD == "DOWN")], col = "darkblue")

dev.off()

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


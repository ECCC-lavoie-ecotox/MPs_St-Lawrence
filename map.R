# Map of the St. Lawrence River
# Written by Raphael Lavoie, Sep 2024
library(prettymapr)
library(sf) 

### Load sites location ####
sites <- readxl::read_xlsx("data/sites_locations.xlsx") |>
    sf::st_as_sf(coords = c("long", "lat"), crs = 4326)

# waters - CANVEC
waters <- sf::st_read("data/waterbody_2.shp") |>
  dplyr::mutate(area = sf::st_area(geometry)) |>
  sf::st_transform(sf::st_crs(4326)) |>
  dplyr::filter(area > units::set_units(10000000, m^2))

# Load USA - naturalearth (CC-0 Licence)
if(!require("rnaturalearth")) install.packages("rnaturalearth")
usa <- rnaturalearth::ne_countries(country = c("united states of america"), scale = 10)

# Load world shapfiles
sf_use_s2(FALSE)
world <- rnaturalearth::ne_countries(type = "countries", scale = "small")
north_america <- sf::st_crop(sf::st_make_valid(world), xmin = -84, ymin = 24, xmax = -50, ymax = 64)

# ### Set up the plot parameters ####
png(file = "StLawrenceRiverMap.png", height = 5, width = 5.3, units="in", res=300)

plot(1, type="n", xlab="Longitude", ylab="Latitude", 
     xlim=c(sf::st_bbox(sites)[1], sf::st_bbox(sites)[3]), ylim=c(sf::st_bbox(sites)[2], sf::st_bbox(sites)[4] + 0.05))

# Add layers to the map
plot(sf::st_geometry(usa), add=TRUE)
plot(sf::st_geometry(waters), col="lightblue", lwd=0.2, add=TRUE)

points(sf::st_coordinates(sites), col = "black", bg="blue", cex = 0.8, pch=21)
points(-73.56, 45.508, col = "black", bg="black", cex = 1, pch=21)
points(-71.254028, 46.829853, col = "black", bg="black", cex = 1, pch=21)
#rect(min(site$long)-0.45,min(site$lat)-0.45,max(site$long)+0.45,max(site$lat)+0.45)# xleft, ybottom, xright, ytop
text(sf::st_coordinates(sites)[,"X"]-0.06, sf::st_coordinates(sites)[,"Y"]+0.08, sites$Number, cex = 0.8)

# Add annotations
text(-71.55, 46.87, "Québec", cex=0.8)
text(-70.9, 46.1, "CANADA", cex=0.8)
text(-69.8, 46.1, "USA", cex=0.8)

addnortharrow(pos="bottomright", scale=0.6, padin=c(0.05, 0.25))
addscalebar(style="bar", label.cex=0.8, pos="bottomright", padin=c(0.05, 0.05))

par(fig = c(0.18, 0.4, 0.60, 0.90), mar = c(0, 0, 2, 2), cex.axis = 0.75, new = TRUE, asp = 1)
par(xaxs = "i", yaxs = "i", mgp = c(0, 0.25, 0), tck = -0.02)

plot(1, type="n", xlab="", ylab="", xlim=c(sf::st_bbox(north_america)[1], sf::st_bbox(north_america)[3]), ylim=c(sf::st_bbox(north_america)[2], sf::st_bbox(north_america)[4]), xaxt="n", yaxt="n") 

rect(sf::st_bbox(north_america)[1], sf::st_bbox(north_america)[3], sf::st_bbox(north_america)[2], sf::st_bbox(north_america)[4], col = "lightblue")
plot(north_america, col="white", add=T, lwd = 0.4)
plot(sf::st_geometry(st_as_sfc(sf::st_bbox(sites), crs = 4326)), add = TRUE, border = "black")
dev.off()


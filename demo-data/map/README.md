

From https://github.com/mtennekes/tmap/blob/master/demo/LondonCrimes/crimes_in_Greater_London.R

## London shape

`library(tmap)
library(tmaptools)
library(rnaturalearth)
library(sp)
library(sf)
regions <- ne_load(scale = "large", type = "states", category = "cultural", destdir = "demo-data/ne")
london <- regions[which(regions$region == "Greater London"),]
london <- set_projection(london, projection = 27700)
writeOGR(london, "demo-data/map/", "London", driver="ESRI Shapefile")
`

## River Thames

`rivers <- ne_load(scale="large", type="rivers_lake_centerlines", category="physical", destdir = "demo-data/ne")
thames <- crop_shape(rivers, london)
writeOGR(thames, "demo-data/map/", "Thames", driver="ESRI Shapefile")
`

## Crime densities

`# download http://www.von-tijn.nl/tijn/research/tmap/crimes_in_Greater_London_2015-10.zip
crimes <- rbind(read.csv("2015-10-city-of-london-street.csv"),
          read.csv("2015-10-metropolitan-street.csv"))

crimes <- crimes[!is.na(crimes$Longitude) & !is.na(crimes$Latitude), ]
coordinates(crimes) <- ~ Longitude + Latitude
crimes <- set_projection(crimes, current.projection = "longlat", projection = 27700)

crimes_london <- crop_shape(crimes, london, polygon =  TRUE)

crime_densities <- my_smooth_map_hack(crimes_london, bandwidth = 0.5, breaks = c(0, 50, 100, 250, 500, 1000), cover = london)

# write out
crime_densities$polygons@data$level <- as.character(crime_densities$polygons@data$level) # cannot save as factors
writeOGR(crime_densities$polygons, "demo-data/map/", "Crimes", driver = "ESRI Shapefile")

`


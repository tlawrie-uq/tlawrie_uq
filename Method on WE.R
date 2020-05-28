
############################################# Method on WE

rm(list = ls(all.names = TRUE)) 

require(sp)
require(raster)
require(rgeos)
require(simba)
require(rgdal)
require(CoordinateCleaner)
require(alphahull)
require(ggplot2)
require(scales)
require(stringr)
require(rasterVis)
require(classInt)
require(dismo)
require(dplyr)
require(methods)

############################################# Retrieving and creating all neccesary arguments for WE calculation

## Retreiving Eucalyptus data (cleaned using the 95% centroid method)
Eucalyptus_cleaned<-read.csv("/Users/Tayla Lawrie/Documents/project/Data/Eucalyptus_cleaned.csv", row.names=NULL)

## Creating a raster (grid map) for WE calaculation
myraster<-raster()
my_crs<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
projection(myraster) <- CRS(my_crs)
extent(myraster)<-c(113.09, 153.38, -43.39, -9.81)
res(myraster)<-1
myraster
mean(values(area(myraster)))

## Creating the outline of Australia for mapping
shapeImport <- readOGR(dsn = "C:/Users/Tayla Lawrie/Documents/project2/Shapefiles", 
                       layer = "GSHHS_h_L1" )
testBB <- bbox(shapeImport)
testBB[1 , ] <- c(113.09, 153.38)
testBB[2 , ] <- c(-43.39, -9.81)
gClip <- function(shp, bb){
  if(class(bb) == "matrix") b_poly <- as(extent(as.vector(t(bb))), "SpatialPolygons")
  else b_poly <- as(extent(bb), "SpatialPolygons")
  gIntersection(shp, b_poly, byid = T)
}
aus_coastline <- gClip(shapeImport, testBB)
outline<-as(aus_coastline, 'SpatialLines')

############################################# Calculating and mapping AOO WE (on both normalised and logarithmic scales)

## Calculate AOO weighted endemism 
aoo<-weighted_endemism(species_records=Eucalyptus_cleaned, species="species", longitude="longitude", 
                       latitude="latitude", frame.raster=myraster, weight.type="cell", 
                       two_method="span", three_plus_method="span", assigned_area=1)
aoo_raster<-aoo$WE_raster
aoo_ranges<-aoo$weights

## Scaling the WE values (normalised and logarithmic)
normal_aoo_raster<-aoo_raster
log_aoo_raster<-aoo_raster
values(normal_aoo_raster)<-rescale(values(normal_aoo_raster), to = c(0, 1))
values(log_aoo_raster)<-log(values(log_aoo_raster))
normal_aoo_values<-values(normal_aoo_raster)
log_aoo_values<-values(log_aoo_raster)

## Plotting the WE maps 
plot(aoo_raster, colNA="NA")
plot(outline, add=TRUE)

plot(normal_aoo_raster, colNA="NA", xlab="Longitude (degrees)", ylab="Latitude (degrees)",
     legend.args=list(text='WE', side=3, font=2, line=0.6, cex=0.8))
plot(outline, add=TRUE)

plot(log_aoo_raster, colNA="NA",  xlab="Longitude (degrees)", ylab="Latitude (degrees)",
     legend.args=list(text=' WE (log)', side=3, font=2, line=0.6, cex=0.8))
plot(outline, add=TRUE)

############################################# Calculating and mapping EOO WE (on both normalised and logarithmic scales)

## Calculate EOO weighted endemism 
eoo<-weighted_endemism(species_records=Eucalyptus_cleaned, species="species", longitude="longitude", 
                       latitude="latitude", frame.raster=myraster, weight.type="geo", 
                       two_method="span", three_plus_method="span", assigned_area=1)

eoo_raster<-eoo$WE_raster
eoo_ranges<-eoo$weights

## Scaling the WE values (normalised and logarithmic)
normal_eoo_raster<-eoo_raster
log_eoo_raster<-eoo_raster
values(normal_eoo_raster)<-rescale(values(normal_eoo_raster), to = c(0, 1))
values(log_eoo_raster)<-log(values(log_eoo_raster))
normal_eoo_values<-values(normal_eoo_raster)
log_eoo_values<-values(log_eoo_raster)

## Plotting the WE maps
plot(eoo_raster, colNA="NA")
plot(outline, add=TRUE)

plot(normal_eoo_raster, colNA="NA", xlab="Longitude (degrees)", ylab="Latitude (degrees)",
     legend.args=list(text='WE', side=3, font=2, line=0.6, cex=0.8))
plot(outline, add=TRUE)

plot(log_eoo_raster, colNA="NA", xlab="Longitude (degrees)", ylab="Latitude (degrees)",
     legend.args=list(text=' WE (log)', side=3, font=2, line=0.6, cex=0.8))
plot(outline, add=TRUE)

############################################# AOO vs. EOO WE cell value scatterplot

## plot AOO versus EOO WE cell values
plot(log_eoo_values, log_aoo_values, xlab="EOO weighted endemism (log)", ylab="AOO weighted endemism (log)")
plot(normal_eoo_values, normal_aoo_values,xlab="EOO weighted endemism", ylab="AOO weighted endemism")

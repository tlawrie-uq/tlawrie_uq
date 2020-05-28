
############################################# Scale on method

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

## Retreiving all taxa data (cleaned using the 95% centroid method)
Eucalyptus_cleaned<-read.csv("/Users/Tayla Lawrie/Documents/project/Data/Eucalyptus_cleaned.csv", row.names=NULL)
Melaleuca_cleaned<-read.csv("/Users/Tayla Lawrie/Documents/project/Data/Melaleuca_cleaned.csv", row.names=NULL)
Corymbia_cleaned<-read.csv("/Users/Tayla Lawrie/Documents/project/Data/Corymbia_cleaned.csv", row.names=NULL)
Leptospermum_cleaned<-read.csv("/Users/Tayla Lawrie/Documents/project/Data/Leptospermum_cleaned.csv", row.names=NULL)
Callistemon_cleaned<-read.csv("/Users/Tayla Lawrie/Documents/project/Data/Callistemon_cleaned.csv", row.names=NULL)

## Creating a raster (grid map) for WE calaculation
myraster<-raster()
my_crs<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
projection(myraster) <- CRS(my_crs)
extent(myraster)<-c(113.09, 153.38, -43.39, -9.81)
res(myraster)<-1
myraster
mean(values(area(myraster)))

## Creating the outline of Australia for mapping
shapeImport <- readOGR(dsn = "C:/Users/Tayla Lawrie/Documents/project/Shapefiles", 
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

############################################# Choosing the species 

## Calculate AOO weighted endemism 
aoo<-weighted_endemism(species_records=Corymbia_cleaned, species="species", longitude="longitude", 
                       latitude="latitude", frame.raster=myraster, weight.type="cell", 
                       two_method="span", three_plus_method="span", assigned_area=1)
aoo_ranges<-aoo$weights

## Calculate EOO weighted endemism 
eoo<-weighted_endemism(species_records=Corymbia_cleaned, species="species", longitude="longitude", 
                       latitude="latitude", frame.raster=myraster, weight.type="geo", 
                       two_method="span", three_plus_method="span", assigned_area=1)
eoo_ranges<-eoo$weights

## Compare EOO and AOO range size values
comparison<-vector("numeric", length(aoo_ranges))
names<-names(aoo_ranges)
aoo_r<-as.numeric(unname(aoo_ranges))
eoo_r<-as.numeric(unname(eoo_ranges))
AOOvEOO<-cbind.data.frame(names, eoo_r, aoo_r, comparison)
AOOvEOO
for (i in 1:nrow(AOOvEOO)){
  AOOvEOO$comparison[i]<-AOOvEOO$eoo_r[i]-AOOvEOO$aoo_r[i]
}
AOOvEOO[,-1]<-round(AOOvEOO[,-1],2) 

## Retrieving species specific occurence records (termed equal and large based on their range size patterns)
equal<-Corymbia_cleaned[Corymbia_cleaned$species == "Corymbia blakei",]
large<-Corymbia_cleaned[Corymbia_cleaned$species == "Corymbia papuana",]
rownames(equal)<-NULL
rownames(large)<-NULL

## Plot species occurence records
plot(equal$longitude, equal$latitude)
plot(large$longitude, large$latitude)

## Create Spatial Points objects for each species 
equal_new<-data.frame(longitude=equal$longitude, latitude=equal$latitude)
equal_mat<-as.matrix(equal_new)
equal_sp<-SpatialPoints(equal_mat, proj4string=crs(aoo_raster))
large_new<-data.frame(longitude=large$longitude, latitude=large$latitude)
large_mat<-as.matrix(large_new)
large_sp<-SpatialPoints(large_mat, proj4string=crs(aoo_raster))

## Calculate the convex hull polygons for each species (independent of scale)
e_poly<-gConvexHull(equal_sp)
l_poly<-gConvexHull(large_sp)

## Create an empty raster for plotting
empty_ras<-myraster
values(empty_ras)<-NA

## Create an extent vector for plotting 
ext<-c(123, 150, -31.0, -10.2)

############################################# AOO vs. EOO at SMALL RESOLUTION

## Create empty raster for both species at a small resolution 
e_raster_small<-myraster
l_raster_small<-myraster
res(e_raster_small)<-0.7
res(l_raster_small)<-0.7
values(e_raster_small)<-NA
values(l_raster_small)<-NA

## Measure AOO at this small scale for each species 
e_cells_small<-cellFromXY(e_raster_small, equal_sp)
e_cells_small<-unique(e_cells_small)
e_raster_small[e_cells_small] <- 1
l_cells_small<-cellFromXY(l_raster_small, large_sp)
l_cells_small<-unique(l_cells_small)
l_raster_small[l_cells_small] <- 1

## Create a small resolution grid for mapping 
grid_small <- raster(extent(l_raster_small), res=c(0.7,0.7))
values(grid_small)<-0
proj4string(grid_small)<-proj4string(l_raster_small)
gridpolygon_small <- rasterToPolygons(grid_small)
grid_small<- spTransform(gridpolygon_small, "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ")

## Plotting the EOO and AOO of each species at this small resolution 
plot(empty_ras, col="white", xlab="Longitude (degrees)", ylab="Latitude (degrees)", legend=FALSE, ext=ext)
plot(geometry(grid_small), add=TRUE)
plot(outline, add=TRUE)
plot(e_poly, border=NA, col="royalblue", lwd=2, add=TRUE)
plot(e_raster_small, add=TRUE, col="tomato", legend=FALSE)
plot(e_poly, border="royalblue", col=NA, lwd=3, add=TRUE)
points(equal_sp, pch=20, col="green3", cex=1)

plot(empty_ras, col="coral", xlab="Longitude (degrees)", ylab="Latitude (degrees)", legend=FALSE, ext=ext)
plot(geometry(grid_small), add=TRUE)
plot(outline, add=TRUE)
plot(l_poly, border=NA, col="royalblue", lwd=2, add=TRUE)
plot(l_raster, col="tomato", add=TRUE, legend=FALSE)
plot(l_poly, border="royalblue", col=NA, lwd=3, add=TRUE)
points(large_sp, pch=20, col="green3", cex=1)

############################################# AOO vs. EOO at LARGE RESOLUTION

## Create empty raster for both species at a large resolution 
e_raster_large<-myraster
l_raster_large<-myraster
res(e_raster_large)<-1.3
res(l_raster_large)<-1.3
values(e_raster_large)<-NA
values(l_raster_large)<-NA

## Measure AOO at this large scale for each species 
e_cells_large<-cellFromXY(e_raster_large, equal_sp)
e_cells_large<-unique(e_cells_large)
e_raster_large[e_cells_large] <- 1
l_cells_large<-cellFromXY(l_raster_large, large_sp)
l_cells_large<-unique(l_cells_large)
l_raster_large[l_cells_large] <- 1

## Create a large resolution grid for mapping 
grid_large <- raster(extent(l_raster_large), res=c(1.3,1.3))
values(grid_large)<-0
proj4string(grid_large)<-proj4string(l_raster_large)
gridpolygon_large <- rasterToPolygons(grid_large)
grid_large <- spTransform(gridpolygon_large, "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ")

## Plotting the EOO and AOO of each species at this small resolution
plot(empty_ras, col="white", xlab="Longitude (degrees)", ylab="Latitude (degrees)", legend=FALSE, ext=ext)
plot(geometry(grid_large), add=TRUE)
plot(outline, add=TRUE)
plot(e_poly, border=NA, col="royalblue", lwd=2, add=TRUE)
plot(e_raster_large, add=TRUE, col="tomato", legend=FALSE)
plot(e_poly, border="royalblue", col=NA, lwd=3, add=TRUE)
points(equal_sp, pch=20, col="green3", cex=1)

plot(empty_ras, col="coral", xlab="Longitude (degrees)", ylab="Latitude (degrees)", legend=FALSE, ext=ext)
plot(geometry(grid_large), add=TRUE)
plot(outline, add=TRUE)
plot(l_poly, border=NA, col="royalblue", lwd=2, add=TRUE)
plot(l_raster_large, col="tomato", add=TRUE, legend=FALSE)
plot(l_poly, border="royalblue", col=NA, lwd=3, add=TRUE)
points(large_sp, pch=20, col="green3", cex=1)

############################################# Measuring the effect of scale on EOO vs AOO for all taxa

## Create and set raster of four different grid cell sizes 
small_ras<-myraster
medium_ras<-myraster
large_ras<-myraster
elarge_ras<-myraster
res(small_ras)<-0.01
res(medium_ras)<-0.1
res(large_ras)<-1
res(elarge_ras)<-1.35

## Retrive taxa occurence data
species_records<-list(Eucalyptus_cleaned, Callistemon_cleaned, Corymbia_cleaned, Melaleuca_cleaned, Leptospermum_cleaned)
names(species_records)<-c("Eucalyptus", "Callistemon", "Corymbia", "Melaleuca", "Leptospermum")

## Create empty data frame to be used later 
summary_table<-data.frame(Taxa=numeric(), small_AE1=numeric(), small_AE2=numeric(), small_AE3=numeric(), 
                          medium_AE1=numeric(), medium_AE2=numeric(), medium_AE3=numeric(),
                          large_AE1=numeric(), large_AE2=numeric(), large_AE3=numeric(),
                          elarge_AE1=numeric(), elarge_AE2=numeric(), elarge_AE3=numeric())

## Calaculating the percent of species whose EOO is larger than AOO for each taxa across the 4 grid cell sizes
for (s in 1:length(species_records)){
  records<-species_records[[s]]
  taxa_name<-names(species_records)[s]
  aoo_small<-weighted_endemism(species_records=records, species="species", longitude="longitude", 
                               latitude="latitude", frame.raster=small_ras, weight.type="cell", 
                               two_method="span", three_plus_method="span", assigned_area=1)
  aoo_small_ranges<-aoo_small$weights
  aoo_medium<-weighted_endemism(species_records=records, species="species", longitude="longitude", 
                                latitude="latitude", frame.raster=medium_ras, weight.type="cell", 
                                two_method="span", three_plus_method="span", assigned_area=1)
  aoo_medium_ranges<-aoo_medium$weights
  aoo_large<-weighted_endemism(species_records=records, species="species", longitude="longitude", 
                               latitude="latitude", frame.raster=large_ras, weight.type="cell", 
                               two_method="span", three_plus_method="span", assigned_area=1)
  aoo_large_ranges<-aoo_large$weights
  aoo_elarge<-weighted_endemism(species_records=records, species="species", longitude="longitude", 
                                latitude="latitude", frame.raster=elarge_ras, weight.type="cell", 
                                two_method="span", three_plus_method="span", assigned_area=1)
  aoo_elarge_ranges<-aoo_elarge$weights
  eoo<-weighted_endemism(species_records=records, species="species", longitude="longitude", 
                         latitude="latitude", frame.raster=small_ras, weight.type="geo", 
                         two_method="span", three_plus_method="span", assigned_area=1)
  eoo_ranges<-eoo$weights
  
  taxa_diff<-as.data.frame(cbind(aoo_small_ranges, aoo_medium_ranges, aoo_large_ranges, 
                                 aoo_elarge_ranges, eoo_ranges))
  
  small_AE<-vector("numeric", length=nrow(taxa_diff))
  medium_AE<-vector("numeric", length=nrow(taxa_diff))
  large_AE<-vector("numeric", length=nrow(taxa_diff))
  elarge_AE<-vector("numeric", length=nrow(taxa_diff))
  
  data_comparison<-cbind.data.frame(taxa_diff, small_AE, medium_AE, large_AE, elarge_AE)
  data_comparison
  
  data_comparison<-round(data_comparison, digits = 2)
  
  for (i in 1:nrow(data_comparison)){
    data_comparison$small_AE[i]<-data_comparison$eoo_ranges[i]-data_comparison$aoo_small_ranges[i]
    data_comparison$medium_AE[i]<-data_comparison$eoo_ranges[i]-data_comparison$aoo_medium_ranges[i]
    data_comparison$large_AE[i]<-data_comparison$eoo_ranges[i]-data_comparison$aoo_large_ranges[i]
    data_comparison$elarge_AE[i]<-data_comparison$eoo_ranges[i]-data_comparison$aoo_elarge_ranges[i]
    if ((data_comparison$small_AE[i])==0){
      data_comparison$small_AE[i]<-"equal" # if EOO = AOO
    }
    if ((data_comparison$small_AE[i])>0){
      data_comparison$small_AE[i]<-1 # if EOO > AOO
    }
    if ((data_comparison$small_AE[i])<0){
      data_comparison$small_AE[i]<-0 # if EOO < AOO
    }
    if ((data_comparison$medium_AE[i])==0){
      data_comparison$medium_AE[i]<-"equal" # if EOO = AOO
    }
    if ((data_comparison$medium_AE[i])>0){
      data_comparison$medium_AE[i]<-1 # if EOO > AOO
    }
    if ((data_comparison$medium_AE[i])<0){
      data_comparison$medium_AE[i]<-0 # if EOO < AOO
    }
    if ((data_comparison$large_AE[i])==0){
      data_comparison$large_AE[i]<-"equal" # if EOO = AOO
    }
    if ((data_comparison$large_AE[i])>0){
      data_comparison$large_AE[i]<-1 # if EOO > AOO
    }
    if ((data_comparison$large_AE[i])<0){
      data_comparison$large_AE[i]<-0 # if EOO < AOO
    }
    if ((data_comparison$elarge_AE[i])==0){
      data_comparison$elarge_AE[i]<-"equal" # if EOO = AOO
    }
    if ((data_comparison$elarge_AE[i])>0){
      data_comparison$elarge_AE[i]<-1 # if EOO > AOO
    }
    if ((data_comparison$elarge_AE[i])<0){
      data_comparison$elarge_AE[i]<-0 # if EOO < AOO
    }
  }
  
  Taxa<-taxa_name
  small_AE1<-length(which(data_comparison$small_AE == "equal"))
  small_AE2<-length(which(data_comparison$small_AE == 1))
  small_AE3<-length(which(data_comparison$small_AE == 0))
  medium_AE1<-length(which(data_comparison$medium_AE == "equal"))
  medium_AE2<-length(which(data_comparison$medium_AE == 1))
  medium_AE3<-length(which(data_comparison$medium_AE == 0))
  large_AE1<-length(which(data_comparison$large_AE == "equal"))
  large_AE2<-length(which(data_comparison$large_AE == 1))
  large_AE3<-length(which(data_comparison$large_AE == 0))
  elarge_AE1<-length(which(data_comparison$elarge_AE == "equal"))
  elarge_AE2<-length(which(data_comparison$elarge_AE == 1))
  elarge_AE3<-length(which(data_comparison$elarge_AE == 0))
  
  taxa_summary_data<-as.data.frame(cbind(Taxa, small_AE1, small_AE2, small_AE3, 
                                         medium_AE1, medium_AE2, medium_AE3, 
                                         large_AE1, large_AE2, large_AE3,
                                         elarge_AE1, elarge_AE2, elarge_AE3))
  taxa_summary<-rbind(summary_table, taxa_summary_data)
  summary_table<-taxa_summary
}

## Saving the data
summary_table
write.csv(summary_table, "scale_on_method.csv")

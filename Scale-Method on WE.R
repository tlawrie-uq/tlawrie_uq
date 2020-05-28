
############################################# Scale + Method on WE

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

## Retreiving all taxa data
Eucalyptus<-read.csv("/Users/Tayla Lawrie/Documents/project/Data/Eucalyptus.csv", row.names=NULL)
Melaleuca<-read.csv("/Users/Tayla Lawrie/Documents/project/Data//Melaleuca.csv", row.names=NULL)
Corymbia<-read.csv("/Users/Tayla Lawrie/Documents/project/Data/Corymbia.csv", row.names=NULL)
Leptospermum<-read.csv("/Users/Tayla Lawrie/Documents/project/Data/Leptospermum.csv", row.names=NULL)
Callistemon<-read.csv("/Users/Tayla Lawrie/Documents/project/Data/Callistemon.csv", row.names=NULL)

## Create a list of species occurences and a vector of their names 
records<-list(Eucalyptus, Corymbia, Callistemon, Melaleuca, Leptospermum) 
name_vec<-c("Eucalyptus", "Corymbia", "Callistemon", "Melaleuca", "Leptospermum")

## Create a vector of 31 different resolutions to be used later 
res_vec<-c(1:31)
res_vec[1]<-0.01
for (i in 2:length(res_vec)){
  if ( i == 1){
    res_vec[1]<-0.01
  }
  if (i == 2){
    res_vec[2]<-0.05
  }
  else {res_vec[i]<-res_vec[i-1]+0.05}
}

## Creating a raster (grid map) for WE calaculation
myraster<-raster()
my_crs<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
projection(myraster) <- CRS(my_crs)
extent(myraster)<-c(113.09, 153.38, -43.39, -9.81)
res(myraster)<-0.01
myraster

## Calculating Endemism Difference (ED) and Hotspot Difference (HD) for each taxa using the defined resolutions
SM_WE<-automated_endemism(resolutions=res_vec, species_records=records, 
                               species="species",longitude="longitude", latitude="latitude", 
                               frame.raster=myraster, names=name_vec, minimal_method="keep",
                               two_method="span", three_plus_method="span",assigned_area=1,
                               difference_method="absolute", outlier_method="centroid", 
                               outlier_crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", 
                               outlier_percent=5, min_occs=6,scaling="yes", 
                               hotspot_percent=5, location_area=7688287) 

## Retreving relevent function outputs 
ED<-SM_WE$ED_data
HD<-SM_WE$HD_data
ED.plot<-SM_WE$ED_plot
HD.plot<-SM_WE$HD_plot

## Saving the data
write.csv(ED, 'ED.csv')
write.csv(HD,'HD.csv')

############################################# Plotting HD for Eucalyptus

## Retreiving Eucalyptus data (cleaned using the 95% centroid method)
Eucalyptus_cleaned<-read.csv("/Users/Tayla Lawrie/Documents/project/Data/Eucalyptus_cleaned.csv", row.names=NULL)

## Cretae raster with small and large cell sizes 
small_ras<-myraster
res(small_ras)<-0.7
mean(values(area(small_ras)))
large_ras<-myraster
res(large_ras)<-1.35
mean(values(area(large_ras)))

## Calculate AOO and EOO weighted endemism at both small and large scales 
aoo_small<-weighted_endemism(species_records=Eucalyptus_cleaned, species="species", longitude="longitude", 
                             latitude="latitude", frame.raster=small_ras, weight.type="cell", 
                             two_method="span", three_plus_method="span", assigned_area=1)
aoo_small_raster<-aoo_small$WE_raster
aoo_large<-weighted_endemism(species_records=Eucalyptus_cleaned, species="species", longitude="longitude", 
                             latitude="latitude", frame.raster=large_ras, weight.type="cell", 
                             two_method="span", three_plus_method="span", assigned_area=1)
aoo_large_raster<-aoo_large$WE_raster
eoo_small<-weighted_endemism(species_records=Eucalyptus_cleaned, species="species", longitude="longitude", 
                             latitude="latitude", frame.raster=small_ras, weight.type="geo", 
                             two_method="span", three_plus_method="span", assigned_area=1)
eoo_small_raster<-eoo_small$WE_raster
eoo_large<-weighted_endemism(species_records=Eucalyptus_cleaned, species="species", longitude="longitude", 
                             latitude="latitude", frame.raster=large_ras, weight.type="geo", 
                             two_method="span", three_plus_method="span", assigned_area=1)
eoo_large_raster<-eoo_large$WE_raster

## Define the total land area and percent of that land area to be used in identifying hotspots 
loc_percent<-5
loc_area<-7688287
p<-loc_percent/100
a<-loc_area
hot_area<-p*a

#################### Small grid cell size

## Determine which WE cell numbers of chosen by both AOO (v1) and EOO (v2) at this small scale 
v1_small<-values(aoo_small_raster)
names(v1_small)<-1:length(v1_small)
v2_small<-values(eoo_small_raster)
names(v2_small)<-1:length(v2_small)
v1_small<-v1_small[!is.na(v1_small)]
v2_small<-v2_small[!is.na(v2_small)]
c_size_small<-mean(values(area(aoo_small_raster)))
num_cells_small<-hot_area/c_size_small
num_cells_small<-round(num_cells_small)
v1_small<-sort(v1_small, decreasing=TRUE)
v1_small<-v1_small[1:num_cells_small]
v1_small<-as.numeric(names(v1_small))
v2_small<-sort(v2_small, decreasing=TRUE)
v2_small<-v2_small[1:num_cells_small]
v2_small<-as.numeric(names(v2_small))
sam_cells_small<-match(v1_small, v2_small)
sam_cells_small<-sam_cells_small[which(!is.na(sam_cells_small))]
sam_cells_small<-v2[sam_cells_small]

## Calculate the percent difference in chosen cells for AOO and EOO (as termed Hotspot Difference (HD)) at this small scale
select_small<-rep(0, length(v2_small))
for (i in 1:length(v2_small)){
  check_small<-match(v1_small[i],v2_small)
  if (is.na(check_small)==TRUE) {
    select_small[i]<-0
  }
  if (is.na(check_small)==FALSE) {
    select_small[i]<-1
  }
}
num_sam_cells_small<-length(which(select_small==1))
num_cells_small<-length(select_small)
small_percent_sam<-100-((num_sam_cells_small/num_cells_small)*100)

## Mapping hotspots at this small scale (both EOO and AOO cells, EOO only cells, and AOO only cells)
small<-aoo_small_raster
values(small)<-NA
small[v1_small]<--0.5
small[v2_small]<-1.5
small[sam_cells_small]<-0.5
plot(small, breaks = c(-1, 0, 1, 2), col = rainbow(3), xlab="Longitude (degrees)", ylab="Latitude (degrees)")
plot(outline, add=TRUE)

#################### Large grid cell size

## Determine which WE cell numbers of chosen by both AOO (v1) and EOO (v2) at this large scale
v1_large<-values(aoo_large_raster)
names(v1_large)<-1:length(v1_large)
v2_large<-values(eoo_large_raster)
names(v2_large)<-1:length(v2_large)
v1_large<-v1_large[!is.na(v1_large)]
v2_large<-v2_large[!is.na(v2_large)]
c_size_large<-mean(values(area(aoo_large_raster)))
num_cells_large<-hot_area/c_size_large
num_cells_large<-round(num_cells_large)
v1_large<-sort(v1_large, decreasing=TRUE)
v1_large<-v1_large[1:num_cells_large]
v1_large<-as.numeric(names(v1_large))
v2_large<-sort(v2_large, decreasing=TRUE)
v2_large<-v2_large[1:num_cells_large]
v2_large<-as.numeric(names(v2_large))
sam_cells_large<-match(v1_large, v2_large)
sam_cells_large<-sam_cells_large[which(!is.na(sam_cells_large))]
sam_cells_large<-v2[sam_cells_large]

## Calculate the percent difference in chosen cells for AOO and EOO (as termed Hotspot Difference (HD)) at this large scale
select_large<-rep(0, length(v2_large))
for (i in 1:length(v2_large)){
  check_large<-match(v1_large[i],v2_large)
  if (is.na(check_large)==TRUE) {
    select_large[i]<-0
  }
  if (is.na(check_large)==FALSE) {
    select_large[i]<-1
  }
}
num_sam_cells_large<-length(which(select_large==1))
num_cells_large<-length(select_large)
large_percent_sam<-100-((num_sam_cells_large/num_cells_large)*100)

## Mapping hotspots at this large scale (both EOO and AOO cells, EOO only cells, and AOO only cells)
large<-aoo_large_raster
values(large)<-NA
large[v1_large]<--0.5
large[v2_large]<-1.5
large[sam_cells_large]<-0.5
plot(large, breaks = c(-1, 0, 1, 2), col = rainbow(3), xlab="Longitude (degrees)", ylab="Latitude (degrees)")
plot(outline, add=TRUE)


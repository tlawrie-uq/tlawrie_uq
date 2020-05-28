
############################################# Minimal method on WE

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

## Retreiving all relevant taxa data
Eucalyptus<-read.csv("/Users/Tayla Lawrie/Documents/project/Data/Eucalyptus.csv", row.names=NULL)
Melaleuca<-read.csv("/Users/Tayla Lawrie/Documents/project/Data//Melaleuca.csv", row.names=NULL)
Leptospermum<-read.csv("/Users/Tayla Lawrie/Documents/project/Data/Leptospermum.csv", row.names=NULL)

## Create a list of species occurences, a vector of their names and a vector of the minimal methods and areas to be anlaysed
records<-list(Eucalyptus, Melaleuca, Leptospermum) 
name_vec<-c("Eucalyptus", "Melaleuca", "Leptospermum")
min.methods<-c("remove", "assign", "span", "jitter")
assign.areas<-c(0.1, 1, 10)

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

## Calculating Endemism Difference (ED) and Hotspot Difference (HD) for each defined taxa, resolution and minimal method
minimal_WE<-outlier_method(resolutions=res_vec, species_records=records, 
                           species="species", longitude="longitude", latitude="latitude", 
                           frame.raster=myraster,names=name_vec, minimal_method="keep",
                           difference_method="absolute", minimal_methods=min.methods, assign_values=assign.areas, 
                           outlier_crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", 
                           outlier_percent=5, min_occs=6, outlier_threshold=5, outlier_raster.res=0.1,
                           package_mltpl=5, package_tdi=1000, hotspot_percent=5, location_area=7688287)

## Retreving relevent function outputs 
min.ED<-minimal_WE$ED_data
min.HD<-minimal_WE$HD_data 

## Saving the data
write.csv(min.ED, 'min.ED.csv')
write.csv(min.HD,'min.HD.csv')

############################################# Plotting ED

min.ED.euc<-min.ED_final[which(min.ED$Taxa=="Eucalyptus"),]
min.ED.lep<-min.ED_final[which(min.ED$Taxa=="Leptospermum"),]
min.ED.mel<-min.ED_final[which(min.ED$Taxa=="Melaleuca"),]

min.ED.euc_plot<- ggplot(min.ED.euc, aes(Cell_size, ED)) + geom_line(aes(group = Outlier_method), colour = "grey50") + 
  geom_point(aes(colour = Minimal_method)) + ylab("Endemism Difference") + xlab("Cell size (km2)") + theme_bw() +
  scale_x_continuous(breaks=seq(0,20000,2500)) +  scale_y_continuous(breaks=seq(0,80,10))
min.ED.euc_plot

min.ED.mel_plot<- ggplot(min.ED.mel, aes(Cell_size, ED)) + geom_line(aes(group = Outlier_method), colour = "grey50") + 
  geom_point(aes(colour = Minimal_method)) + ylab("Endemism Difference") + xlab("Cell size (km2)") + theme_bw()+
  scale_x_continuous(breaks=seq(0,20000,2500)) +  scale_y_continuous(breaks=seq(0,50,5))
min.ED.mel_plot

min.ED.lep_plot<- ggplot(min.ED.lep, aes(Cell_size, ED)) + geom_line(aes(group = Outlier_method), colour = "grey50") + 
  geom_point(aes(colour = Minimal_method)) + ylab("Endemism Difference") + xlab("Cell size (km2)") + theme_bw()+
  scale_x_continuous(breaks=seq(0,20000,2500)) +  scale_y_continuous(breaks=seq(0,40,5))
min.ED.lep_plot

############################################# Plotting HD

min.HD.euc<-min.HD_final[which(min.HD$Taxa=="Eucalyptus"),]
min.HD.lep<-min.HD_final[which(min.HD$Taxa=="Leptospermum"),]
min.HD.mel<-min.HD_final[which(min.HD$Taxa=="Melaleuca"),]

min.HD.euc_plot<- ggplot(min.HD.euc, aes(Cell_size, HD)) + geom_line(aes(group = Outlier_method), colour = "grey50") + 
  geom_point(aes(colour = Minimal_method)) + ylab("Hotspot Difference (%)") + xlab("Cell size (km2)") + theme_bw() +
  scale_x_continuous(breaks=seq(0,20000,2500)) +  scale_y_continuous(breaks=seq(0,80,10))
min.HD.euc_plot

min.HD.mel_plot<- ggplot(min.HD.mel, aes(Cell_size, HD)) + geom_line(aes(group = Outlier_method), colour = "grey50") + 
  geom_point(aes(colour = Minimal_method)) + ylab("Hotspot Difference (%)") + xlab("Cell size (km2)") + theme_bw()+
  scale_x_continuous(breaks=seq(0,20000,2500)) +  scale_y_continuous(breaks=seq(0,50,5))
min.HD.mel_plot

min.HD.lep_plot<- ggplot(min.HD.lep, aes(Cell_size, HD)) + geom_line(aes(group = Outlier_method), colour = "grey50") + 
  geom_point(aes(colour = Minimal_method)) + ylab("Hotspot Difference (%)") + xlab("Cell size (km2)") + theme_bw()+
  scale_x_continuous(breaks=seq(0,20000,2500)) +  scale_y_continuous(breaks=seq(0,40,5))
min.HD.lep_plot

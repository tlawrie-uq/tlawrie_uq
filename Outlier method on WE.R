
############################################# Outlier method on WE

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

## Create a list of species occurences, a vector of their names and a vector of the outlier methods to be anlaysed
records<-list(Eucalyptus, Corymbia, Callistemon, Melaleuca, Leptospermum) 
name_vec<-c("Eucalyptus", "Corymbia", "Callistemon", "Melaleuca", "Leptospermum")
out.methods<-c("centroid", "relative.impact", "distance", "quantile", "mad")

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

## Calculating Endemism Difference (ED) and Hotspot Difference (HD) for each defined taxa, resolution and outlier method
outlier_WE<-outlier_method(resolutions=res_vec, species_records=records, 
                              species="species", longitude="longitude", latitude="latitude", 
                              frame.raster=myraster,names=name_vec, minimal_method="keep",
                              difference_method="absolute",two_method="span", three_plus_method="span", 
                              outlier_methods=out.methods, assigned_area=1, 
                              outlier_crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", 
                              outlier_percent=5, min_occs=6, outlier_threshold=5, outlier_raster.res=0.1,
                              package_mltpl=5, package_tdi=1000, hotspot_percent=5, location_area=7688287)

## Retreving relevent function outputs 
out.ED<-outlier_WE$ED_data
out.HD<-outlier_WE$HD_data 

## Saving the data
write.csv(out.ED, 'out.ED.csv')
write.csv(out.HD,'out.HD.csv')

############################################# Plotting ED

out.ED.euc<-out.ED_final[which(out.ED$Taxa=="Eucalyptus"),]
out.ED.cor<-out.ED_final[which(out.ED$Taxa=="Corymbia"),]
out.ED.call<-out.ED_final[which(out.ED$Taxa=="Callistemon"),]
out.ED.lep<-out.ED_final[which(out.ED$Taxa=="Leptospermum"),]
out.ED.mel<-out.ED_final[which(out.ED$Taxa=="Melaleuca"),]

out.ED.euc_plot<- ggplot(out.ED.euc, aes(Cell_size, ED)) + geom_line(aes(group = Outlier_method), colour = "grey50") + 
  geom_point(aes(colour = Outlier_method)) + ylab("Endemism Difference") + xlab("Cell size (km2)") + theme_bw() +
  scale_x_continuous(breaks=seq(0,20000,2500)) +  scale_y_continuous(breaks=seq(0,80,10))
out.ED.euc_plot

out.ED.cor_plot<- ggplot(out.ED.cor, aes(Cell_size, ED)) + geom_line(aes(group = Outlier_method), colour = "grey50") + 
  geom_point(aes(colour = Outlier_method)) + ylab("Endemism Difference") + xlab("Cell size (km2)") + theme_bw()+
  scale_x_continuous(breaks=seq(0,20000,2500)) +  scale_y_continuous(breaks=seq(0,45,5))
out.ED.cor_plot

out.ED.call_plot<- ggplot(out.ED.call, aes(Cell_size, ED)) + geom_line(aes(group = Outlier_method), colour = "grey50") + 
  geom_point(aes(colour = Outlier_method)) + ylab("Endemism Difference") + xlab("Cell size (km2)") + theme_bw()+
  scale_x_continuous(breaks=seq(0,20000,2500)) +  scale_y_continuous(breaks=seq(0,40,5))
out.ED.call_plot

out.ED.mel_plot<- ggplot(out.ED.mel, aes(Cell_size, ED)) + geom_line(aes(group = Outlier_method), colour = "grey50") + 
  geom_point(aes(colour = Outlier_method)) + ylab("Endemism Difference") + xlab("Cell size (km2)") + theme_bw()+
  scale_x_continuous(breaks=seq(0,20000,2500)) +  scale_y_continuous(breaks=seq(0,50,5))
out.ED.mel_plot

out.ED.lep_plot<- ggplot(out.ED.lep, aes(Cell_size, ED)) + geom_line(aes(group = Outlier_method), colour = "grey50") + 
  geom_point(aes(colour = Outlier_method)) + ylab("Endemism Difference") + xlab("Cell size (km2)") + theme_bw()+
  scale_x_continuous(breaks=seq(0,20000,2500)) +  scale_y_continuous(breaks=seq(0,40,5))
out.ED.lep_plot


############################################# Plotting HD

out.HD.euc<-out.HD_final[which(out.HD$Taxa=="Eucalyptus"),]
out.HD.cor<-out.HD_final[which(out.HD$Taxa=="Corymbia"),]
out.HD.call<-out.HD_final[which(out.HD$Taxa=="Callistemon"),]
out.HD.lep<-out.HD_final[which(out.HD$Taxa=="Leptospermum"),]
out.HD.mel<-out.HD_final[which(out.HD$Taxa=="Melaleuca"),]

out.HD.euc_plot<- ggplot(out.HD.euc, aes(Cell_size, HD)) + geom_line(aes(group = Outlier_method), colour = "grey50") + 
  geom_point(aes(colour = Outlier_method)) + ylab("Hotspot Difference (%)") + xlab("Cell size (km2)") + theme_bw() +
  scale_x_continuous(breaks=seq(0,20000,2500)) +  scale_y_continuous(breaks=seq(0,80,10))
out.HD.euc_plot

out.HD.cor_plot<- ggplot(out.HD.cor, aes(Cell_size, HD)) + geom_line(aes(group = Outlier_method), colour = "grey50") + 
  geom_point(aes(colour = Outlier_method)) + ylab("Hotspot Difference (%)") + xlab("Cell size (km2)") + theme_bw()+
  scale_x_continuous(breaks=seq(0,20000,2500)) +  scale_y_continuous(breaks=seq(0,45,5))
out.HD.cor_plot

out.HD.call_plot<- ggplot(out.HD.call, aes(Cell_size, HD)) + geom_line(aes(group = Outlier_method), colour = "grey50") + 
  geom_point(aes(colour = Outlier_method)) + ylab("Hotspot Difference (%)") + xlab("Cell size (km2)") + theme_bw()+
  scale_x_continuous(breaks=seq(0,20000,2500)) +  scale_y_continuous(breaks=seq(0,40,5))
out.HD.call_plot

out.HD.mel_plot<- ggplot(out.HD.mel, aes(Cell_size, HD)) + geom_line(aes(group = Outlier_method), colour = "grey50") + 
  geom_point(aes(colour = Outlier_method)) + ylab("Hotspot Difference (%)") + xlab("Cell size (km2)") + theme_bw()+
  scale_x_continuous(breaks=seq(0,20000,2500)) +  scale_y_continuous(breaks=seq(0,50,5))
out.HD.mel_plot

out.HD.lep_plot<- ggplot(out.HD.lep, aes(Cell_size, HD)) + geom_line(aes(group = Outlier_method), colour = "grey50") + 
  geom_point(aes(colour = Outlier_method)) + ylab("Hotspot Difference (%)") + xlab("Cell size (km2)") + theme_bw()+
  scale_x_continuous(breaks=seq(0,20000,2500)) +  scale_y_continuous(breaks=seq(0,40,5))
out.HD.lep_plot

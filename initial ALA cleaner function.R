initial_cleaner<-function(data=occ_dat, accepted=names, ext_raster=aus)
#Description --
  #This function removes missing data, duplicates, ocean occuring records and records outside mainland Australia
#Usage -- 
  #For example: initial_cleaner(data=occ_dat, accepted=names, ext_raster=aus)
#Arguments --
  #data
    #The occurence data set that you want cleaned (with "species", "latitude" and "longitude" column names - in that order)
  #accepted
    #A data frame downloaded from ALA containing infromation regarding the accpted status of each species names (column names should be "species" and "status")
  #ext_raster:
    #A raster with the extent of the region of mainland Australia of interest 
#Value --
  #Returns a data frame with 2 columns (e.g. Resolution and ED) with 1 row if the input data frame is empty and n+1 rows if the input data frame contains n rows) 
#Required packages --
  #raster, CoordinateCleaner
{
  require(raster)
  require(CoordinateCleaner)
  
  data<-data[complete.cases(data), ]
  two<-vector("numeric", length=nrow(data))
  check<-vector("numeric", length=nrow(data))
  data<-cbind(data, two, check)
  data$two<-word(data$species, 1,2)
  rownames(data)<-NULL
  for (i in 1:nrow(data)){
    if ((data$two[i]==data$species[i])==FALSE){
      data$check[i]<-1
    }
  }
  remove<-which(data$check==1)
  if (length(remove)==0){
    data<-data[,1:3]
  }
  if (length(remove)>0){
    data<-data[-which(data$check==1),1:3]
  }
  accepted<-as.vector(accepted[which(accepted$status=="accepted"),1])
  check<-vector("numeric", length=nrow(data))
  data<-cbind(data, check)
  rownames(data)<-NULL
  for (i in 1:nrow(data)) {
    if ((is.na(match(data$species[i], accepted)))==TRUE){
      data$check[i]<-1
    }
  }
  remove<-which(data$check==1)
  if (length(remove)==0){
    data<-data[,1:3]
  }
  if (length(remove)>0){
    data<-data[-which(data$check==1),1:3]
  }
  rownames(data)<-NULL
  data<-as.data.frame(cc_dupl(data, value = "clean", lon=longitude, lat=latitude, species=species))
  data<-as.data.frame(cc_sea(data, lon="longitude", lat="latitude"))
  species_records<-data
  frame.raster<-ext_raster
  colnames(species_records)<-c("SPECIES", "LATITUDE", "LONGITUDE")
  coordinates(species_records) <- c("LONGITUDE", "LATITUDE")
  if(!(extent(species_records)@xmin >= extent(frame.raster)@xmin & extent(species_records)@xmax <= extent(frame.raster)@xmax & extent(species_records)@ymin >= extent(frame.raster)@ymin & extent(species_records)@ymax <= extent(frame.raster)@ymax)) {
    species_record_COORDS <- as.data.frame(coordinates(species_records))
    species_records <- species_records[-which(species_record_COORDS$LONGITUDE < extent(frame.raster)@xmin | species_record_COORDS$LONGITUDE > extent(frame.raster)@xmax | species_record_COORDS$LATITUDE < extent(frame.raster)@ymin | species_record_COORDS$LATITUDE > extent(frame.raster)@ymax),]
  } 
  species_records<-as.data.frame(species_records)
  colnames(species_records)<-c("species", "latitude", "longitude")
  if((length(which(species_records$longitude>145 & species_records$latitude>-12)))>=1){
    species_records<-species_records[-which(species_records$longitude>145 & species_records$latitude>-12),]
  }
  if ((length(which(species_records$longitude<127.5 & species_records$latitude>-12)))>=1){
    species_records<-species_records[-which(species_records$longitude<127.5 & species_records$latitude>-12),]
  }
  row.names(species_records)<-NULL
  return(species_records)
}

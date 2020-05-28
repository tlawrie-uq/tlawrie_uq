hotspot_difference<-function(raster1, raster2,  hot_data, loc_percent, loc_area) 
#Description --
  #This function calculates hotspot difference (HD) based on two WE raster (e.g. ones calculated using area of occupancy (AOO) and extent of occurrence (EOO) range weight methods) 
#Usage -- 
  #For example: hotspot_difference(raster1=AOO_ras, raster2=EOO, hot_data=)
#Arguments --
  #raster1 and raster2
    #The 2 endemism rasters you wish to compare 
  #hot_data
    #A data frame with 2 numeric columns (e.g. data.frame(Cell_size=numeric(0), HD=numeric(0)) - it can be empty or already contain values)
  #loc_percent
    #The percent of total area to be identified as endemic hotspots 
  #loc_area
    #The total area of the region of study
#Value --
  #Returns a data frame with 2 columns (e.g. Resolution and HD) with 1 row if the input data frame is empty and n+1 rows if the input data frame contains n rows) 
#Required packages --
  #raster
{
  require(raster)
  
  v1<-values(raster1)
  names(v1)<-1:length(v1)
  v2<-values(raster2)
  names(v2)<-1:length(v2)
  v1<-v1[!is.na(v1)]
  v2<-v2[!is.na(v2)]
  p<-loc_percent/100
  a<-loc_area
  hot_area<-p*a
  c_size<-mean(values(area(raster1)))
  num_cells<-hot_area/c_size
  num_cells<-round(num_cells)
  v1<-sort(v1, decreasing=TRUE)
  v1<-v1[1:num_cells]
  v1<-as.numeric(names(v1))
  v2<-sort(v2, decreasing=TRUE)
  v2<-v2[1:num_cells]
  v2<-as.numeric(names(v2))
  select<-rep(0, length(v2))
  for (i in 1:length(v2)){
    check<-match(v1[i],v2)
    if (is.na(check)==TRUE) {
      select[i]<-0
    }
    if (is.na(check)==FALSE) {
      select[i]<-1
    }
  }
  num_sam_cells<-length(which(select==1))
  num_cells<-length(select)
  percent_sam<-100-((num_sam_cells/num_cells)*100)
  hot_data<-rbind(hot_data,setNames(as.list(c(c_size,percent_sam)), names(hot_data)))
  return(hot_data)
}
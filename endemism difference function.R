endemism_difference<-function(raster1,raster2,diff_data,method)
#Description --
  #This function calculates the overall difference in endemism scores of cells across two rasters (e.g. ones calculated using area of occupancy (AOO) and extent of occurrence (EOO) range weight methods) 
#Usage -- 
  #For example: endemism_difference(raster1=AOO_ras, raster2=EOO, diff_data=dat, method="absolute")
#Arguments --
  #raster1 and raster2
    #The 2 endemism rasters you wish to compare 
  #diff_data
    #A data frame with 2 numeric columns (e.g. data.frame(Cell_size=numeric(0), ED=numeric(0)) - it can be empty or already contain values)
  #method:
    # = "absolute" will take the absolute difference between cell values 
    # = "direction" will take the difference between cell values (could be positive or negative > is always raster1-raster2)
#Value --
  #Returns a data frame with 2 columns (e.g. Resolution and ED) with 1 row if the input data frame is empty and n+1 rows if the input data frame contains n rows) 
#Required packages --
    #raster
{
  require(raster)
  
  v1<-values(raster1)
  v2<-values(raster2)
  v1<-v1[!is.na(v1)]
  v2<-v2[!is.na(v2)]
  resv<-vector(length=length(v1))
  if (method == "absolute") {
    for (x in 1:length(resv)){
      resv[x]=abs(v1[x]-v2[x])
    }
  }
  if (method == "direction") {
    for (x in 1:length(resv)){
      resv[x]=v1[x]-v2[x]
    }
  }
  finalresult<-(sum(resv))/length(resv)
  size<-mean(values(area(raster1)))
  diff_data<-rbind(diff_data,setNames(as.list(c(size,finalresult)), names(diff_data)))
  return(diff_data)
}

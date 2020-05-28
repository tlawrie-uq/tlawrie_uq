coordinate_cleaner<- function(occurence_data, method, crs, percent, threshold, raster.res,
                              package_method, package_mltpl, package_tdi, min_occs)
#Description --
  #This function removes duplicate and outlier coordinates from a data set based on a defined method 
# Usage --
  #For example: coordinate_cleaner(occurence_data=mydata, method="centroid", crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", percent=5, min_occs=6)
#Arguments --
  #occurence_data
    #The occurence data set that you want cleaned (with "species", "latitude" and "longitude" column names - in that order)
  #method:
    # = "centroid" will use the centroid method to remove outliers (remove "percent" number coordinates furthest from polygon centre)
    # = "relative.impact" will remove all occurence points that add amounts equal to/larger than the percentage "threshold" to the polygon area relative to all other points 
    # = "package.cleaner" will use the CoordinateCleaner package to clean the data set (see ?cc_outl)
  #crs
    #The cordinate reference system to be used 
  #percent
    #The percent of all data points that should be removed when using the centroid method 
  #threshold
    #The threshold area to be used for the relative.impact method 
  #raster.res
    #The resolution of the raster to be used with the "distance" method 
  #package_method:
    # = "distance" > records are flagged as outliers, if the minimum distance to the next record of the species is > tdi. For species with records from > 10000 unique locations a random sample of 1000 records is used for the distance matrix calculation. The test skipps species with less than min_occs, geographically unique records     
    # = "quantile" > a boxplot method is used and records are flagged as outliers if their mean distance to all other records of the same species is larger than mltpl * the interquartile range of the mean distance of all records of this species
    # = "mad" > the median absolute deviation is used. In this case a record is flagged as outlier, if the mean distance to all other records of the same species is larger than the median of the mean distance of all points plus/minus the mad of the mean distances of all records of the species * mltpl
  #package_mltpl
    #The multiplier of the interquartile range (method == 'quantile') or median absolute deviation (method == 'mad')to identify outliers
  #package_tdi
    #The minimum absolute distance (method == 'distance') of a record to all other records of a species to be identified as outlier, in km
  #min_occs
      #Minimum number of geographically unique datapoints needed for a species to be tested. This is necessary for reliable outlier estimation. Species wit less than min_occs records will not be tested and the output value will be 'TRUE'. If method == 'distance', consider a lower threshold.
#Value --
  #Returns the cleaned data set   
#Required packages --
  #raster, sp, rgeos, CoordinateCleaner
{
  require(sp)
  require(raster)
  require(rgeos)
  require(CoordinateCleaner)
  
  if(class(occurence_data) != "data.frame") {
    stop("Species data must be in a data.frame")
  } 
  spec<-unique(occurence_data[,1])
  
  if(method == "centroid") {
    crdref<-CRS(crs)
    occur_data<-data.frame(species=character(), latitude=numeric(), longitude=numeric())
    for (i in 1:length(spec)) {
      spec_data<-as.data.frame(occurence_data[which(occurence_data$species == spec[i]),])
      row.names(spec_data)<-NULL
      if (nrow(spec_data) > min_occs) {
        lon.lat<-cbind(spec_data$longitude, spec_data$latitude)
        data.pts<-SpatialPoints(lon.lat, proj4string =crdref)
        hull<-gConvexHull(data.pts)
        centre<-gCentroid(hull)
        centre_pts<-as.vector(extent(centre))[c(2,4)]
        x_c = centre_pts[1]
        y_c = centre_pts[2]
        lon.lat.2<-as.data.frame(lon.lat)
        x<-as.vector(lon.lat.2[,1])
        y<-as.vector(lon.lat.2[,2])
        d = sapply(seq_along(x), function(ind){
          dist(rbind(c(x_c, y_c), c(x[ind], y[ind])))
        }) #Calculate distance of each point from the center
        k = (percent/100)*(length(x))
        x2 = head(x[order(d)], -k)
        y2 = head(y[order(d)], -k)
        new_lon_lat<-as.data.frame(cbind(y2, x2))
        name_vec<-rep(spec[i], each=nrow(new_lon_lat))
        new_spec_data<-as.data.frame(cbind(name_vec,new_lon_lat))
        colnames(new_spec_data)<-c("species", "latitude", "longitude")
        cumul_spec_data<-rbind(occur_data, new_spec_data)
        occur_data<-cumul_spec_data
      }
      else {
        new_spec_data<-spec_data
        cumul_spec_data<-rbind(occur_data, new_spec_data)
        occur_data<-cumul_spec_data
      }
    }
    cleaned_occurence_data<-occur_data
  } #cls if(method == "centroid")

  if(method == "relative.impact") {
    genus_occ_data<-data.frame(species=character(), latitude=numeric(), longitude=numeric(), relative_impact=numeric())
    for (i in 1:length(spec)) {
      data1<-occurence_data[which(occurence_data$species == spec[i]),]
      row.names(data1)<-NULL
      if (nrow(data1) > min_occs) {
        lon.lat<-cbind(data1$longitude, data1$latitude)
        crdref<-CRS(crs)
        data.pts<-SpatialPoints(lon.lat, proj4string =crdref)
        percent2<- vector(length=nrow(data1))
        data2 <- cbind(data1, percent2)
        poly<-gConvexHull(data.pts)
        original<-(area(poly))/10^6
        for (b in 1:nrow(data2)) {
          poly<-gConvexHull(data.pts[-b,])
          new_poly<-(area(poly))/10^6
          if (((new_poly/original)*100)>100) {
            data2[b,4]<-0
          }
          else {
            data2[b,4]<-100-((new_poly/original)*100)
          }
        }
        specdata<-data2 #species data with each points impact on the range size 
        sorted_specdata <- specdata[with(specdata, order(-percent2)), ] #species data with each points impact on the range size (ordered with largest impact first)
        sorted_specdata<-cbind(sorted_specdata, relative_impact=vector(length=nrow(sorted_specdata)))
        for (g in 1:nrow(sorted_specdata)) {
          relative_vec<-as.numeric(sorted_specdata[-g,4])
          sorted_specdata[g,5]<-sorted_specdata[g,4]-mean(relative_vec)
          if(sorted_specdata[g,5]<0){
            sorted_specdata[g,5]<-0
          }
        }
        species_occ_data<-as.data.frame(sorted_specdata[,-4])
        occ_data<-as.data.frame(rbind(genus_occ_data, species_occ_data))
        genus_occ_data<-occ_data
      }
      else {
        minim_dat<-data1
        minim_dat$relative_impact<-rep(0, length=nrow(minim_dat))
        occ_data<-as.data.frame(rbind(genus_occ_data, minim_dat))
        genus_occ_data<-occ_data
      }
    }
    genus_occ_data_new<-genus_occ_data[which(genus_occ_data$relative_impact<=threshold),1:3]
    cleaned_occurence_data<-as.data.frame(genus_occ_data_new)
  } #cls if(method == "relative.impact")
  if(method == "package.cleaner") {
    outl <- cc_outl(occurence_data, method=package_method, mltpl=package_mltpl, tdi=package_tdi, min_occs=min_occs, thinning_res=raster.res, value = "clean", lon="longitude", lat="latitude", species="species")
    cleaned_occurence_data<-outl
  } #cls if(method == "package.cleaner")
  row.names(cleaned_occurence_data)<-NULL
  return(cleaned_occurence_data)
} # cls function 

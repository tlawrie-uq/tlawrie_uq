minimal_method<-function(resolutions=resolution_vec, species_records=cal, species="SPECIES",
                         longitude="LONGITUDE", latitude="LATITUDE", frame.raster=myraster,
                         names=names_vec, minimal_methods=keep_mth,
                         jitter_amount=0.1, difference_method="absolute", 
                         outlier_method="centroid", outlier_crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", 
                         outlier_percent=0.05, min_occs=5, assign_values=assign_vec, scaling="yes",
                         hotspot_percent=5, location_area=45599)
#Description -- 
  #Calculates Endemism Difference (ED) and Hotspot Difference (HD) for each inputted taxa, resolution and minimal method (based on the 'weighted_endemism' function)
#Usage --
  #For example: outlier_method((resolutions=resolution_vec, species_records=mydata, minimal_methods=min_methds, species="SPECIES",longitude="LONGITUDE", latitude="LATITUDE", frame.raster=myraster, names=names_vec, minimal_method="keep", two_method="assign", three_plus_method="assign", assigned_area=0.0001, difference_method="absolute", outlier_method="centroid", outlier_crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", outlier_percent=0.05, min_occs=6,  scaling="yes", hotspot_percent=5, location_area=45599)  
#Arguments --
  #resolutions
    #A vector containing the resolutions to be analysed
    #Must contain resolutions that are whole interger multiples of the resolution of the frame.raster (e.g. if res(frame.raster)=0.1, resolutions could be 0.2 or 1 but not 0.45)
  #species_records
    #A data.frame with rows as individual species records, and columns that include fields for species name, longitude and latitude (see species, longitude, latitude below)
    #Or a list of data.frames
  #minimal_methods
    #A character vector containing the methods wished to be analysed (based on cooridnate_cleaner inputs)
    #For example: methds<-c("centroid", "mad)
  #names
    #A vector containing the taxa names of the species_records (in the correct order)
  #species
    #What colname in the supplied species_records contains species names
  #latitude
    #What colname in the supplied species_records contains latitude values
  #longitude
    #What colname in the supplied species_records contains longitude values
  #frame.raster 
    #An existing raster object used as the basis for testing each new resolution (with a specified res, extent and crs)
  #assigned_area
    #The area (in square kilometers) that will be assigned to the range values of species with minimal occurences, if the assign or assign_span methods are chosen (see below)
  #jitter_amount
    #The amount that the data will be jitterred, if the jitter method is chosen for three occurences (see below)
  #difference_method:
    # ="absolute" will calculate absolute ED
    # ="direction" will calculate directional ED
  #outlier_method:
    # = "centroid" will use the centroid method to remove outliers (remove "percent" number coordinates furthest from polygon centre)
    # = "relative.impact" will remove all occurence points that add the largest amount to the polygon area relative to all other points (based on the area "threshold")  
    # = "package.cleaner" will use the CoordinateCleaner package to clean the data set 
  #outlier_crs
    #The cordinate reference system to be used in the coordinate_cleaner function
    #outlier_percent
  #The percent of all data points that should be removed when using the centroid method
    #outlier_threshold
    #The threshold area to be used for the relative.impact method
  #outlier_raster.res
    #The resolution of the raster to be used with the "distance" method
  #package_method
    # = "distance" > records are flagged as outliers, if the minimum distance to the next record of the species is > tdi. For species with records from > 10000 unique locations a random sample of 1000 records is used for the distance matrix calculation. The test skipps species with less than min_occs, geographically unique records     
    # = "quantile" > a boxplot method is used and records are flagged as outliers if their mean distance to all other records of the same species is larger than mltpl * the interquartile range of the mean distance of all records of this species
    # = "mad" > the median absolute deviation is used. In this case a record is flagged as outlier, if the mean distance to all other records of the same species is larger than the median of the mean distance of all points plus/minus the mad of the mean distances of all records of the species * mltpl
  #package_mltpl
    #The multiplier of the interquartile range (method == 'quantile') or median absolute deviation (method == 'mad')to identify outliers
  #package_tdi
    #The minimum absolute distance (method == 'distance') of a record to all other records of a species to be identified as outlier, in km
  #min_occs
    #Minimum number of geographically unique datapoints needed for a species to be tested. This is necessary for reliable outlier estimation. Species wit less than min_occs records will not be tested and the output value will be 'TRUE'. If method == 'distance', consider a lower threshold.
  #scaling
    #="yes" will normally scale the WE values
    #="no" will leave WE values as those originally produced by the weighted_endemism function
  #hotspot_percent
    #The percent of total area to be identified as endemic hotspots for HD
  #location_area
    #The total area of the region of study (to be using in measuring HD)
#Value --
  #Returns a list of length 2:
    #ED_data: Data frame of ED measures for each method supplied 
    #HD_data: Data frame of HD measures for each method supplied 
#Required packages --
  #sp, raster, regeos, simba, CoordinateCleaner
#Required functions --
  #automated_endemism, coordinate_cleaner, weighted_endemism, hotspot_difference, endemism_difference
{
  require(sp)
  require(raster)
  require(rgeos)
  require(simba)
  require(CoordinateCleaner)
  
  if (class(species_records)=="data.frame") {
    endemdiff<-data.frame(Taxa=character(), Cell_size=numeric(0), ED=numeric(0), Minimal_method=character())
    hotspotdiff<-data.frame(Taxa=character(), Cell_size=numeric(0), HD=numeric(0), Minimal_method=character())
    names<-names
    for (i in 1:length(minimal_methods)){
      method<-minimal_methods[i]
      if ((method == "remove")==TRUE){
        cat("Method: remove ", "\n")
        method_dat<-automated_endemism(resolutions=resolutions, species_records=species_records, species=species,
                                       longitude=longitude, latitude=latitude, frame.raster=frame.raster,
                                       names=names, minimal_method="remove",
                                       difference_method=difference_method, 
                                       outlier_method=outlier_method, outlier_crs=outlier_crs,
                                       outlier_percent=outlier_percent, min_occs=min_occs, scaling=scaling, 
                                       hotspot_percent=hotspot_percent, location_area=location_area)
        method_data1<-method_dat$ED_data
        method_data2<-method_dat$HD_data
        method_vec<-rep("remove", nrow(method_data1))
        data1<-cbind(method_data1, method_vec)
        data2<-cbind(method_data2, method_vec)
        colnames(data1)[4]<-"Minimal_method"
        colnames(data2)[4]<-"Minimal_method"
        data_sum1<-rbind(endemdiff, data1)
        endemdiff<-data_sum1
        data_sum2<-rbind(hotspotdiff, data2)
        hotspotdiff<-data_sum2
      }
      if ((method != "remove")==TRUE){
        cat("Method: keep ", "\n")
        if((method=="assign")==TRUE) {
          cat("Keep method: assign ", "\n")
          if ((length(assign_values)==1)){
            assigned_area<-assign_values
            method_dat<-automated_endemism(resolutions=resolutions, species_records=species_records, species=species,
                                           longitude=longitude, latitude=latitude, frame.raster=frame.raster,
                                           names=names, minimal_method="keep",
                                           two_method=method, three_plus_method=method,
                                           assigned_area=assigned_area, jitter_amount=jitter_amount, difference_method=difference_method, 
                                           outlier_method=outlier_method, outlier_crs=outlier_crs,
                                           outlier_percent=outlier_percent, min_occs=min_occs, scaling=scaling,
                                           hotspot_percent=hotspot_percent, location_area=location_area)
            method_data1<-method_dat$ED_data
            method_data2<-method_dat$HD_data
            method_vec<-rep("assign", nrow(method_data1))
            data1<-cbind(method_data1, method_vec)
            data2<-cbind(method_data2, method_vec)
            colnames(data1)[4]<-"Minimal_method"
            colnames(data2)[4]<-"Minimal_method"
            data_sum1<-rbind(endemdiff, data1)
            endemdiff<-data_sum1
            data_sum2<-rbind(hotspotdiff, data2)
            hotspotdiff<-data_sum2
          }
          if ((length(assign_values)>1)){
            for (a in 1:length(assign_values)) {
              assigned_area<-assign_values[a]
              method_dat<-automated_endemism(resolutions=resolutions, species_records=species_records, species=species,
                                             longitude=longitude, latitude=latitude, frame.raster=frame.raster,
                                             names=names, minimal_method="keep",
                                             two_method=method, three_plus_method=method,
                                             assigned_area=assigned_area, jitter_amount=jitter_amount, difference_method=difference_method, 
                                             outlier_method=outlier_method, outlier_crs=outlier_crs,
                                             outlier_percent=outlier_percent, min_occs=min_occs, scaling=scaling,
                                             hotspot_percent=hotspot_percent, location_area=location_area)
              vec<-paste(method,as.character(assigned_area))
              vec <- sub(" ", "_", vec)
              method_vec<-rep(vec, nrow(method_data1))
              method_data1<-method_dat$ED_data
              method_data2<-method_dat$HD_data
              data1<-cbind(method_data1, method_vec)
              data2<-cbind(method_data2, method_vec)
              colnames(data1)[4]<-"Minimal_method"
              colnames(data2)[4]<-"Minimal_method"
              data_sum1<-rbind(endemdiff, data1)
              endemdiff<-data_sum1
              data_sum2<-rbind(hotspotdiff, data2)
              hotspotdiff<-data_sum2
            }
          }
        }
        if ((method=="jitter")==TRUE){
          cat("Keep method: jitter ", "\n")
          method_dat<-automated_endemism(resolutions=resolutions, species_records=species_records, species=species,
                                         longitude=longitude, latitude=latitude, frame.raster=frame.raster,
                                         names=names, minimal_method="keep",
                                         two_method="span", three_plus_method=method,
                                         assigned_area=1, jitter_amount=jitter_amount, difference_method=difference_method, 
                                         outlier_method=outlier_method, outlier_crs=outlier_crs,
                                         outlier_percent=outlier_percent, min_occs=min_occs, scaling=scaling,
                                         hotspot_percent=hotspot_percent, location_area=location_area)
          method_data1<-method_dat$ED_data
          method_data2<-method_dat$HD_data
          method_vec<-rep("jitter", nrow(method_data1))
          data1<-cbind(method_data1, method_vec)
          data2<-cbind(method_data2, method_vec)
          colnames(data1)[4]<-"Minimal_method"
          colnames(data2)[4]<-"Minimal_method"
          data_sum1<-rbind(endemdiff, data1)
          endemdiff<-data_sum1
          data_sum2<-rbind(hotspotdiff, data2)
          hotspotdiff<-data_sum2
        }
        if ((method=="span")==TRUE){
          cat("Keep method: span ", "\n")
          method_dat<-automated_endemism(resolutions=resolutions, species_records=species_records, species=species,
                                         longitude=longitude, latitude=latitude, frame.raster=frame.raster,
                                         names=names, minimal_method="keep",
                                         two_method=method, three_plus_method=method,
                                         assigned_area=1, jitter_amount=jitter_amount, difference_method=difference_method, 
                                         outlier_method=outlier_method, outlier_crs=outlier_crs,
                                         outlier_percent=outlier_percent, min_occs=min_occs, scaling=scaling,
                                         hotspot_percent=hotspot_percent, location_area=location_area)
          method_data1<-method_dat$ED_data
          method_data2<-method_dat$HD_data
          method_vec<-rep("span", nrow(method_data1))
          data1<-cbind(method_data1, method_vec)
          data2<-cbind(method_data2, method_vec)
          colnames(data1)[4]<-"Minimal_method"
          colnames(data2)[4]<-"Minimal_method"
          data_sum1<-rbind(endemdiff, data1)
          endemdiff<-data_sum1
          data_sum2<-rbind(hotspotdiff, data2)
          hotspotdiff<-data_sum2
        }
      }
    }
    taxa_summary_data1<-endemdiff
    taxa_summary_data2<-hotspotdiff
  }
  if(class(species_records)=="list") {
    total_endemdiff<-data.frame(Taxa=character(), Cell_size=numeric(0), ED=numeric(0), Minimal_method=character()) #create an empty dataframe (used later for the ED function)
    total_hotspotdiff<-data.frame(Taxa=character(), Cell_size=numeric(0), HD=numeric(0), Minimal_method=character())
    for (x in 1:length(species_records)){
      endemdiff<-data.frame(Taxa=character(), Cell_size=numeric(0), ED=numeric(0), Minimal_method=character())
      hotspotdiff<-data.frame(Taxa=character(), Cell_size=numeric(0), HD=numeric(0), Minimal_method=character())
      occurrence_records<-species_records[[x]]
      name<-names[x]
      for (i in 1:length(minimal_methods)){
        method<-minimal_methods[i]
        if ((method == "remove")==TRUE){
          cat("Method: remove ", "\n")
          method_dat<-automated_endemism(resolutions=resolutions, species_records=occurrence_records, species=species,
                                         longitude=longitude, latitude=latitude, frame.raster=frame.raster,
                                         names=name, minimal_method="remove",
                                         difference_method=difference_method, 
                                         outlier_method=outlier_method, outlier_crs=outlier_crs,
                                         outlier_percent=outlier_percent, min_occs=min_occs, scaling=scaling,
                                         hotspot_percent=hotspot_percent, location_area=location_area)
          method_data1<-method_dat$ED_data
          method_data2<-method_dat$HD_data
          method_vec<-rep("remove", nrow(method_data1))
          data1<-cbind(method_data1, method_vec)
          data2<-cbind(method_data2, method_vec)
          colnames(data1)[4]<-"Minimal_method"
          colnames(data2)[4]<-"Minimal_method"
          data_sum1<-rbind(endemdiff, data1)
          endemdiff<-data_sum1
          data_sum2<-rbind(hotspotdiff, data2)
          hotspotdiff<-data_sum2
        }
        if ((method != "remove")==TRUE){
          cat("Method: keep ", "\n")
          if((method=="assign")==TRUE) {
            cat("Keep method: assign ", "\n")
            if ((length(assign_values)==1)){
              assigned_area<-assign_values
              method_dat<-automated_endemism(resolutions=resolutions, species_records=occurrence_records, species=species,
                                             longitude=longitude, latitude=latitude, frame.raster=frame.raster,
                                             names=name, minimal_method="keep",
                                             two_method=method, three_plus_method=method,
                                             assigned_area=assigned_area, jitter_amount=jitter_amount, difference_method=difference_method, 
                                             outlier_method=outlier_method, outlier_crs=outlier_crs,
                                             outlier_percent=outlier_percent, min_occs=min_occs, scaling=scaling,
                                             hotspot_percent=hotspot_percent, location_area=location_area)
              method_data1<-method_dat$ED_data
              method_data2<-method_dat$HD_data
              method_vec<-rep("assign", nrow(method_data1))
              data1<-cbind(method_data1, method_vec)
              data2<-cbind(method_data2, method_vec)
              colnames(data1)[4]<-"Minimal_method"
              colnames(data2)[4]<-"Minimal_method"
              data_sum1<-rbind(endemdiff, data1)
              endemdiff<-data_sum1
              data_sum2<-rbind(hotspotdiff, data2)
              hotspotdiff<-data_sum2
            }
            if ((length(assign_values)>1)){
              for (a in 1:length(assign_values)) {
                assigned_area<-assign_values[a]
                method_dat<-automated_endemism(resolutions=resolutions, species_records=occurrence_records, species=species,
                                               longitude=longitude, latitude=latitude, frame.raster=frame.raster,
                                               names=name, minimal_method="keep",
                                               two_method=method, three_plus_method=method,
                                               assigned_area=assigned_area, jitter_amount=jitter_amount, difference_method=difference_method, 
                                               outlier_method=outlier_method, outlier_crs=outlier_crs,
                                               outlier_percent=outlier_percent, min_occs=min_occs, scaling=scaling,
                                               hotspot_percent=hotspot_percent, location_area=location_area)
                vec<-paste(method,as.character(assigned_area))
                vec <- sub(" ", "_", vec)
                method_vec<-rep(vec, nrow(method_data1))
                method_data1<-method_dat$ED_data
                method_data2<-method_dat$HD_data
                data1<-cbind(method_data1, method_vec)
                data2<-cbind(method_data2, method_vec)
                colnames(data1)[4]<-"Minimal_method"
                colnames(data2)[4]<-"Minimal_method"
                data_sum1<-rbind(endemdiff, data1)
                endemdiff<-data_sum1
                data_sum2<-rbind(hotspotdiff, data2)
                hotspotdiff<-data_sum2
              }
            }
          }
          if ((method=="jitter")==TRUE){
            cat("Keep method: jitter ", "\n")
            method_dat<-automated_endemism(resolutions=resolutions, species_records=occurrence_records, species=species,
                                           longitude=longitude, latitude=latitude, frame.raster=frame.raster,
                                           names=name, minimal_method="keep",
                                           two_method="span", three_plus_method=method,
                                           assigned_area=1, jitter_amount=jitter_amount, difference_method=difference_method, 
                                           outlier_method=outlier_method, outlier_crs=outlier_crs,
                                           outlier_percent=outlier_percent, min_occs=min_occs, scaling=scaling,
                                           hotspot_percent=hotspot_percent, location_area=location_area)
            method_data1<-method_dat$ED_data
            method_data2<-method_dat$HD_data
            method_vec<-rep("jitter", nrow(method_data1))
            data1<-cbind(method_data1, method_vec)
            data2<-cbind(method_data2, method_vec)
            colnames(data1)[4]<-"Minimal_method"
            colnames(data2)[4]<-"Minimal_method"
            data_sum1<-rbind(endemdiff, data1)
            endemdiff<-data_sum1
            data_sum2<-rbind(hotspotdiff, data2)
            hotspotdiff<-data_sum2
          }
          if ((method=="span")==TRUE){
            cat("Keep method: span ", "\n")
            method_dat<-automated_endemism(resolutions=resolutions, species_records=occurrence_records, species=species,
                                           longitude=longitude, latitude=latitude, frame.raster=frame.raster,
                                           names=name, minimal_method="keep",
                                           two_method=method, three_plus_method=method,
                                           assigned_area=1, jitter_amount=jitter_amount, difference_method=difference_method, 
                                           outlier_method=outlier_method, outlier_crs=outlier_crs,
                                           outlier_percent=outlier_percent, min_occs=min_occs, scaling=scaling,
                                           hotspot_percent=hotspot_percent, location_area=location_area)
            method_data1<-method_dat$ED_data
            method_data2<-method_dat$HD_data
            method_vec<-rep("span", nrow(method_data1))
            data1<-cbind(method_data1, method_vec)
            data2<-cbind(method_data2, method_vec)
            colnames(data1)[4]<-"Minimal_method"
            colnames(data2)[4]<-"Minimal_method"
            data_sum1<-rbind(endemdiff, data1)
            endemdiff<-data_sum1
            data_sum2<-rbind(hotspotdiff, data2)
            hotspotdiff<-data_sum2
          }
        }
      }
      taxa_summary_data1<-endemdiff
      taxa_summary_data2<-hotspotdiff
      spec_dat_meth1<-rbind(total_endemdiff, taxa_summary_data1)
      spec_dat_meth2<-rbind(total_hotspotdiff, taxa_summary_data2)
      total_endemdiff<-spec_dat_meth1
      total_hotspotdiff<-spec_dat_meth2
    }
    taxa_summary_data1<-total_endemdiff
    taxa_summary_data2<-total_hotspotdiff
  }
  outputs<- list(ED_data=taxa_summary_data1, HD_data=taxa_summary_data2)
  return(outputs)
}

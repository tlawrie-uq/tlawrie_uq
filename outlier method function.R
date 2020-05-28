outlier_method<-function(resolutions=resolution_vec, species_records=cal, species="SPECIES",
                         longitude="LONGITUDE", latitude="LATITUDE", frame.raster=myraster,
                         names=names_vec, minimal_method="keep",
                         assigned_area=0.1, jitter_amount=0.1, difference_method="absolute", 
                         two_method="span", three_plus_method="span", outlier_methods=methds, 
                         outlier_crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", 
                         outlier_percent=0.05, min_occs=5, outlier_threshold=5, outlier_raster.res=0.1,
                         package_mltpl=5, package_tdi=5, hotspot_percent=5, location_area=45599)
#Description -- 
  #Calculates Endemism Difference (ED) and Hotspot Difference (HD) for each inputted taxa, resolution and outlier method (based on the 'coordinate_cleaner' function)
#Usage --
  #For example: outlier_method(resolutions=resolution_vec, species_records=mydata, outlier_methods=out_methds, species="SPECIES",longitude="LONGITUDE", latitude="LATITUDE", frame.raster=myraster, names=names_vec, minimal_method="keep", two_method="assign", three_plus_method="assign", assigned_area=0.0001, difference_method="absolute", outlier_crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", outlier_percent=0.05, min_occs=6,  scaling="yes", hotspot_percent=5, location_area=45599)  
#Arguments --
  #resolutions
    #A vector containing the resolutions to be analysed
    #Must contain resolutions that are whole interger multiples of the resolution of the frame.raster (e.g. if res(frame.raster)=0.1, resolutions could be 0.2 or 1 but not 0.45)
  #species_records
    #A data.frame with rows as individual species records, and columns that include fields for species name, longitude and latitude (see species, longitude, latitude below)
    #Or a list of data.frames
  #outlier_methods
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
  #minimal_method:
    # ="keep" will keep non-polygon occuring species in the data set
    # ="remove" will remove non-polygon occuring species in the data set 
  #assigned_area
    #The area (in square kilometers) that will be assigned to the range values of species with minimal occurences, if the assign or assign_span methods are chosen (see below)
  #jitter_amount
    #The amount that the data will be jitterred, if the jitter method is chosen for three occurences (see below)
  #two_method:
    # = "assign" will assign a specifiied area per record as the range of the species based on "assigned_area"
    # = "span" will set the range value of the species to the distance between the two occurence points
    # = "assign_span" will set the range of the species as a combination of the assigned areas and distance between the two points 
  #three_plus_method:
    # = "assign" will assign a specifiied area as the range of the species based on "assigned_area"
    # = "span" will set the range value of the species to the maximum distance between the three occurence points
    # = "assign_span" will set the range of the species as a combination of the assigned areas and distance between the three points 
    # = "jitter" will randomly jitter a single point based on the "jitter_emount" specified 
  #difference_method:
    # ="absolute" will calculate absolute ED
    # ="direction" will calculate directional ED
  #outlier_crs
    #The cordinate reference system to be used in the coordinate_cleaner function
  #outlier_percent
    #The percent of all data points that should be removed when using the centroid method
  #outlier_threshold
    #The threshold area to be used for the relative.impact method
  #outlier_raster.res
    #The resolution of the raster to be used with the "distance" method
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
    endemdiff<-data.frame(Taxa=character(), Cell_size=numeric(0), ED=numeric(0), Outlier_method=character())
    hotspotdiff<-data.frame(Taxa=character(), Cell_size=numeric(0), HD=numeric(0), Minimal_method=character())
    names<-names
    for (i in 1:length(outlier_methods)){
      method<-outlier_methods[i]
      if ((method == "quantile")== TRUE) {
        cat ("Quantile method ", "\n")
        method_dat<-automated_endemism(resolutions=resolutions, species_records=species_records, species=species,
                                       longitude=longitude, latitude=latitude, frame.raster=frame.raster,
                                       names=names, minimal_method=minimal_method,
                                       two_method=two_method, three_plus_method=three_plus_method, 
                                       assigned_area=assigned_area, jitter_amount=jitter_amount, 
                                       difference_method=difference_method, 
                                       outlier_method="package.cleaner", outlier_crs=outlier_crs,
                                       package_method=method, min_occs=min_occs, package_mltpl=package_mltpl, 
                                       outlier_raster.res=outlier_raster.res, hotspot_percent=hotspot_percent, location_area=location_area)
        method_data1<-method_dat$ED_data
        method_data2<-method_dat$HD_data
        method_vec<-rep("quantile", nrow(method_data1))
        data1<-cbind(method_data1, method_vec)
        data2<-cbind(method_data2, method_vec)
        colnames(data1)[4]<-"Outlier_method"
        colnames(data2)[4]<-"Outlier_method"
        data_sum1<-rbind(endemdiff, data1)
        endemdiff<-data_sum1
        data_sum2<-rbind(hotspotdiff, data2)
        hotspotdiff<-data_sum2
      }
      if ((method=="mad")==TRUE) {
        cat ("Mad method ", "\n")
        method_dat<-automated_endemism(resolutions=resolutions, species_records=species_records, species=species,
                                       longitude=longitude, latitude=latitude, frame.raster=frame.raster,
                                       names=names, minimal_method=minimal_method,
                                       two_method=two_method, three_plus_method=three_plus_method, 
                                       assigned_area=assigned_area, jitter_amount=jitter_amount, 
                                       difference_method=difference_method, 
                                       outlier_method="package.cleaner", outlier_crs=outlier_crs,
                                       package_method=method, min_occs=min_occs,package_mltpl=package_mltpl, 
                                       outlier_raster.res=outlier_raster.res, hotspot_percent=hotspot_percent, location_area=location_area)
        method_data1<-method_dat$ED_data
        method_data2<-method_dat$HD_data
        method_vec<-rep("mad", nrow(method_data1))
        data1<-cbind(method_data1, method_vec)
        data2<-cbind(method_data2, method_vec)
        colnames(data1)[4]<-"Outlier_method"
        colnames(data2)[4]<-"Outlier_method"
        data_sum1<-rbind(endemdiff, data1)
        endemdiff<-data_sum1
        data_sum2<-rbind(hotspotdiff, data2)
        hotspotdiff<-data_sum2
      }
      if ((method == "distance")==TRUE) {
        cat ("Distance method ", "\n")
        method_dat<-automated_endemism(resolutions=resolutions, species_records=species_records, species=species,
                                       longitude=longitude, latitude=latitude, frame.raster=frame.raster,
                                       names=names, minimal_method=minimal_method,
                                       two_method=two_method, three_plus_method=three_plus_method, 
                                       assigned_area=assigned_area, jitter_amount=jitter_amount, 
                                       difference_method=difference_method,
                                       outlier_method="package.cleaner", outlier_crs=outlier_crs,
                                       package_method=method, min_occs=min_occs,
                                       package_tdi=package_tdi, outlier_raster.res=outlier_raster.res, hotspot_percent=hotspot_percent, location_area=location_area)
        method_data1<-method_dat$ED_data
        method_data2<-method_dat$HD_data
        method_vec<-rep("distance", nrow(method_data1))
        data1<-cbind(method_data1, method_vec)
        data2<-cbind(method_data2, method_vec)
        colnames(data1)[4]<-"Outlier_method"
        colnames(data2)[4]<-"Outlier_method"
        data_sum1<-rbind(endemdiff, data1)
        endemdiff<-data_sum1
        data_sum2<-rbind(hotspotdiff, data2)
        hotspotdiff<-data_sum2
      }
      if ((method == "centroid")==TRUE) {
        cat ("Centroid method ", "\n")
        method_dat<-automated_endemism(resolutions=resolutions, species_records=species_records, species=species,
                                       longitude=longitude, latitude=latitude, frame.raster=frame.raster,
                                       names=names, minimal_method=minimal_method,
                                       two_method=two_method, three_plus_method=three_plus_method, 
                                       assigned_area=assigned_area, jitter_amount=jitter_amount, 
                                       difference_method=difference_method,
                                       outlier_method=method, outlier_crs=outlier_crs,
                                       outlier_percent=outlier_percent, min_occs=min_occs, hotspot_percent=hotspot_percent, location_area=location_area)
        method_data1<-method_dat$ED_data
        method_data2<-method_dat$HD_data
        method_vec<-rep("centroid", nrow(method_data1))
        data1<-cbind(method_data1, method_vec)
        data2<-cbind(method_data2, method_vec)
        colnames(data1)[4]<-"Outlier_method"
        colnames(data2)[4]<-"Outlier_method"
        data_sum1<-rbind(endemdiff, data1)
        endemdiff<-data_sum1
        data_sum2<-rbind(hotspotdiff, data2)
        hotspotdiff<-data_sum2
      }
      if ((method == "relative.impact")==TRUE) {
        cat ("Relative impact method ", "\n")
        method_dat<-automated_endemism(resolutions=resolutions, species_records=species_records, species=species,
                                       longitude=longitude, latitude=latitude, frame.raster=frame.raster,
                                       names=names, minimal_method=minimal_method,
                                       two_method=two_method, three_plus_method=three_plus_method, 
                                       assigned_area=assigned_area, jitter_amount=jitter_amount, 
                                       difference_method=difference_method,
                                       outlier_method=method, outlier_crs=outlier_crs,
                                       outlier_threshold=outlier_threshold, min_occs=min_occs, hotspot_percent=hotspot_percent, location_area=location_area)
        method_data1<-method_dat$ED_data
        method_data2<-method_dat$HD_data
        method_vec<-rep("relative_impact", nrow(method_data1))
        data1<-cbind(method_data1, method_vec)
        data2<-cbind(method_data2, method_vec)
        colnames(data1)[4]<-"Outlier_method"
        colnames(data2)[4]<-"Outlier_method"
        data_sum1<-rbind(endemdiff, data1)
        endemdiff<-data_sum1
        data_sum2<-rbind(hotspotdiff, data2)
        hotspotdiff<-data_sum2
      } 
    }
    taxa_summary_data1<-endemdiff
    taxa_summary_data2<-hotspotdiff
  }
  if(class(species_records)=="list") {
    total_endemdiff<-data.frame(Taxa=character(), Cell_size=numeric(0), ED=numeric(0), Outlier_method=character()) #create an empty dataframe (used later for the ED function)
    total_hotspotdiff<-data.frame(Taxa=character(), Cell_size=numeric(0), HD=numeric(0), Minimal_method=character())
    for (x in 1:length(species_records)){
      endemdiff<-data.frame(Taxa=character(), Cell_size=numeric(0), ED=numeric(0), Outlier_method=character())
      hotspotdiff<-data.frame(Taxa=character(), Cell_size=numeric(0), HD=numeric(0), Minimal_method=character())
      occurrence_records<-species_records[[x]]
      name<-names[x]
      for (i in 1:length(outlier_methods)){
        method<-outlier_methods[i]
        if ((method == "quantile")== TRUE) {
          cat ("Quantile method ", "\n")
          method_dat<-automated_endemism(resolutions=resolutions, species_records=occurrence_records, species=species,
                                         longitude=longitude, latitude=latitude, frame.raster=frame.raster,
                                         names=name, minimal_method=minimal_method,
                                         two_method=two_method, three_plus_method=three_plus_method, 
                                         assigned_area=assigned_area, jitter_amount=jitter_amount, 
                                         difference_method=difference_method,
                                         outlier_method="package.cleaner", outlier_crs=outlier_crs,
                                         package_method=method, min_occs=min_occs, package_mltpl=package_mltpl, 
                                         outlier_raster.res=outlier_raster.res, hotspot_percent=hotspot_percent, location_area=location_area)
          method_data1<-method_dat$ED_data
          method_data2<-method_dat$HD_data
          method_vec<-rep("quantile", nrow(method_data1))
          data1<-cbind(method_data1, method_vec)
          data2<-cbind(method_data2, method_vec)
          colnames(data1)[4]<-"Outlier_method"
          colnames(data2)[4]<-"Outlier_method"
          data_sum1<-rbind(endemdiff, data1)
          endemdiff<-data_sum1
          data_sum2<-rbind(hotspotdiff, data2)
          hotspotdiff<-data_sum2
        }
        if ((method=="mad")==TRUE) {
          cat ("Mad method ", "\n")
          method_dat<-automated_endemism(resolutions=resolutions, species_records=occurrence_records, species=species,
                                         longitude=longitude, latitude=latitude, frame.raster=frame.raster,
                                         names=name, minimal_method=minimal_method,
                                         two_method=two_method, three_plus_method=three_plus_method, 
                                         assigned_area=assigned_area, jitter_amount=jitter_amount, 
                                         difference_method=difference_method,
                                         outlier_method="package.cleaner", outlier_crs=outlier_crs,
                                         package_method=method, min_occs=min_occs,package_mltpl=package_mltpl, 
                                         outlier_raster.res=outlier_raster.res, hotspot_percent=hotspot_percent, location_area=location_area)
          method_data1<-method_dat$ED_data
          method_data2<-method_dat$HD_data
          method_vec<-rep("mad", nrow(method_data1))
          data1<-cbind(method_data1, method_vec)
          data2<-cbind(method_data2, method_vec)
          colnames(data1)[4]<-"Outlier_method"
          colnames(data2)[4]<-"Outlier_method"
          data_sum1<-rbind(endemdiff, data1)
          endemdiff<-data_sum1
          data_sum2<-rbind(hotspotdiff, data2)
          hotspotdiff<-data_sum2
        }
        if ((method == "distance")==TRUE) {
          cat ("Distance method ", "\n")
          method_dat<-automated_endemism(resolutions=resolutions, species_records=occurrence_records, species=species,
                                         longitude=longitude, latitude=latitude, frame.raster=frame.raster,
                                         names=name, minimal_method=minimal_method,
                                         two_method=two_method, three_plus_method=three_plus_method, 
                                         assigned_area=assigned_area, jitter_amount=jitter_amount, 
                                         difference_method=difference_method, 
                                         outlier_method="package.cleaner", outlier_crs=outlier_crs,
                                         package_method=method, min_occs=min_occs,
                                         package_tdi=package_tdi, outlier_raster.res=outlier_raster.res, hotspot_percent=hotspot_percent, location_area=location_area)
          method_data1<-method_dat$ED_data
          method_data2<-method_dat$HD_data
          method_vec<-rep("distance", nrow(method_data1))
          data1<-cbind(method_data1, method_vec)
          data2<-cbind(method_data2, method_vec)
          colnames(data1)[4]<-"Outlier_method"
          colnames(data2)[4]<-"Outlier_method"
          data_sum1<-rbind(endemdiff, data1)
          endemdiff<-data_sum1
          data_sum2<-rbind(hotspotdiff, data2)
          hotspotdiff<-data_sum2
        }
        if ((method == "centroid")==TRUE) {
          cat ("Centroid method ", "\n")
          method_dat<-automated_endemism(resolutions=resolutions, species_records=occurrence_records, species=species,
                                         longitude=longitude, latitude=latitude, frame.raster=frame.raster,
                                         names=name, minimal_method=minimal_method,
                                         two_method=two_method, three_plus_method=three_plus_method, 
                                         assigned_area=assigned_area, jitter_amount=jitter_amount, 
                                         difference_method=difference_method, 
                                         outlier_method=method, outlier_crs=outlier_crs,
                                         outlier_percent=outlier_percent, min_occs=min_occs, hotspot_percent=hotspot_percent, location_area=location_area)
          method_data1<-method_dat$ED_data
          method_data2<-method_dat$HD_data
          method_vec<-rep("centroid", nrow(method_data1))
          data1<-cbind(method_data1, method_vec)
          data2<-cbind(method_data2, method_vec)
          colnames(data1)[4]<-"Outlier_method"
          colnames(data2)[4]<-"Outlier_method"
          data_sum1<-rbind(endemdiff, data1)
          endemdiff<-data_sum1
          data_sum2<-rbind(hotspotdiff, data2)
          hotspotdiff<-data_sum2
        }
        if ((method == "relative.impact")==TRUE) {
          cat ("Relative impact method ", "\n")
          method_dat<-automated_endemism(resolutions=resolutions, species_records=occurrence_records, species=species,
                                         longitude=longitude, latitude=latitude, frame.raster=frame.raster,
                                         names=name, minimal_method=minimal_method,
                                         two_method=two_method, three_plus_method=three_plus_method, 
                                         assigned_area=assigned_area, jitter_amount=jitter_amount, 
                                         difference_method=difference_method, 
                                         outlier_method=method, outlier_crs=outlier_crs,
                                         outlier_threshold=outlier_threshold, min_occs=min_occs, hotspot_percent=hotspot_percent, location_area=location_area)
          method_data1<-method_dat$ED_data
          method_data2<-method_dat$HD_data
          method_vec<-rep("relative_impact", nrow(method_data1))
          data1<-cbind(method_data1, method_vec)
          data2<-cbind(method_data2, method_vec)
          colnames(data1)[4]<-"Outlier_method"
          colnames(data2)[4]<-"Outlier_method"
          data_sum1<-rbind(endemdiff, data1)
          endemdiff<-data_sum1
          data_sum2<-rbind(hotspotdiff, data2)
          hotspotdiff<-data_sum2
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

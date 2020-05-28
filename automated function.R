automated_endemism <- function(resolutions=resolution_vec, species_records=mydata, species="SPECIES",
                               longitude="LONGITUDE", latitude="LATITUDE", frame.raster=myraster,
                               names=names_vec, minimal_method="keep",
                               two_method="assign", three_plus_method="assign", 
                               assigned_area=0.0001, jitter_amount=0.1, difference_method="absolute", 
                               outlier_method="centroid", outlier_crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", 
                               outlier_percent=0.05, outlier_threshold, outlier_raster.res,
                               package_method, package_mltpl, package_tdi, min_occs,
                               scaling="yes", hotspot_percent=5, location_area=45599) 
#Description --
  #Calaculates Endemism Difference (ED) and Hotspot Difference (HD) for each inputted taxa and resolution
#Usage --
  #For example: automated_endemism(resolutions=resolution_vec, species_records=mydata, species="SPECIES",longitude="LONGITUDE", latitude="LATITUDE", frame.raster=myraster, names=names_vec, minimal_method="keep", two_method="assign", three_plus_method="assign", assigned_area=0.0001, difference_method="absolute", outlier_method="centroid", outlier_crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", outlier_percent=0.05, min_occs=6,  scaling="yes", hotspot_percent=5, location_area=45599)
#Arguments --
  #resolutions
    #A vector containing the resolutions to be analysed
    #Must contain resolutions that are whole interger multiples of the resolution of the frame.raster (e.g. if res(frame.raster)=0.1, resolutions could be 0.2 or 1 but not 0.45)
  #species_records
    #A data.frame with rows as individual species records, and columns that include fields for species name, longitude and latitude (see species, longitude, latitude below)
    #Or a list of data.frames
  #names
    #A vector containing the taxa names of the species_records (in the correct order)
  #species
    #What colname in the supplied species_records contains species names
  #latitude
    #What colname in the supplied species_records contains latitude values
  #longitude
    #What colname in the supplied species_records contains longitude values
  #frame.raster
    #An existing raster object the user can optionally elect to supply as the frame for calculations and mapping (with a specified res, extent and crs)
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
  #Returns a list of length 4:
    #ED_data: Data frame of ED measures
    #ED_plot: Plot of ED measures
    #HD_data: Data frame of HD measures
    #HD_plot: Plot of HD measures
#Required packages --
  #sp, raster, regeos, simba, CoordinateCleaner
#Required functions --
  #coordinate_cleaner, weighted_endemism, hotspot_difference, endemism_difference
{
  require(sp)
  require(raster)
  require(rgeos)
  require(simba)
  require(CoordinateCleaner)
  
  for (i in 1:length(resolutions)){ # checks whether the inputted resolutions are whole integer multiples of the frame raster resolution
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    check<-(resolutions[i])/(res(frame.raster)[1])
    val<-is.wholenumber(check)
    val<-as.character(val)
    if (val!="TRUE") {
      stop("Resolutions are not a whole integer multiple of the frame raster resolution")
    }
  }
  if(class(species_records)=="data.frame"){ #if a single dataframe is inputted 
    cat("Cleaning the data using the specified method ", "\n")
    species_records<-coordinate_cleaner(occurence_data=species_records, method=outlier_method, crs=outlier_crs, 
                                        percent =outlier_percent, threshold = outlier_threshold,
                                        min_occs=min_occs, raster.res=outlier_raster.res,
                                        package_method=package_method, package_mltpl=package_mltpl, 
                                        package_tdi=package_tdi)
    cat("The data has been cleaned ", "\n")
    endemdiff<-data.frame(Cell_size=numeric(0), ED=numeric(0)) #create an empty dataframe (used later for the ED function)
    hotspotdiff<-data.frame(Cell_size=numeric(0), HD=numeric(0)) #create an empty dataframe (used later for the HD function)
    if (minimal_method == "remove"){ #removes non-polygon occurring species
      cat("Removing species with minimal occurence points ", "\n")
      spec<-unique(species_records[,1])
      for (i in spec) {
        spec_data<-species_records[which(species_records$species == i),]
        spec_data<-data.frame(LONGITUDE=spec_data$longitude, LATITUDE=spec_data$latitude)
        dat<-SpatialPoints(spec_data, CRS(outlier_crs))
        species_polygon<-try(gConvexHull(dat))
        if (class(species_polygon) != "SpatialPolygons") {
          species_records<-species_records[!(species_records$species == i),]
          cat("Species ", i, "has been removed from the data set ", "\n")
        }
        if (class(species_polygon) == "SpatialPolygons") {
          species_records<-species_records
        }
      }
    }
    if (minimal_method == "keep") { #removes all non-polygon occurring species
      species_records<-species_records
    }
    for(i in 1:length(resolutions)){ 
      res_check<-as.character(isTRUE(all.equal(resolutions[i], res(frame.raster)[1])))
      if (res_check=="TRUE") {
        calc_raster<-frame.raster 
      } #creates a frame.raster for 
      if (res_check=="FALSE") {
        agg_factor<-(resolutions[i])/(res(frame.raster)[1])
        calc_raster<-aggregate(frame.raster, fact=agg_factor)
      }
      cat("Calculating weighted endemim using AOO with a raster resolution ", resolutions[i], "\n")
      aoo<-weighted_endemism(species_records=species_records, species=species,longitude=longitude,
                             latitude=latitude, frame.raster=calc_raster,
                             weight.type="cell", two_method=two_method, three_plus_method=three_plus_method, 
                             assigned_area=assigned_area, jitter_amount=jitter_amount) #calculate weighted endemism with AOO
      aooraster<-aoo$WE_raster #call upon and label the endemicm raster (map) for aoo
      if (scaling == "yes"){
        values(aooraster)<-rescale(values(aooraster), to = c(0, 1))
      }
      if (scaling == "no") {
        aooraster<-aooraster
      }
      cat("Calculating weighted endemim using EOO with a raster resolution ", resolutions[i], "\n")
      eoo<-weighted_endemism(species_records=species_records, species=species,longitude=longitude,
                             latitude=latitude, frame.raster=calc_raster,
                             weight.type="geo", two_method=two_method, three_plus_method=three_plus_method, 
                             assigned_area=assigned_area, jitter_amount=jitter_amount) #calculate weighted endemism with EOO
      eooraster<-eoo$WE_raster #call upon and label the endemicm raster (map) for eoo
      if (scaling == "yes"){
        values(eooraster)<-rescale(values(eooraster), to = c(0, 1))
      }
      if (scaling == "no") {
        eooraster<-eooraster
      }
      data<-endemism_difference(aooraster,eooraster,endemdiff, method=difference_method) #calculate the endemism difference between using AOO and EOO
      data2<-hotspot_difference(aooraster,eooraster,hotspotdiff, loc_percent=hotspot_percent, loc_area=location_area)
      endemdiff<-data #this new data will be called upon in the ED function for the next resolution (the new resolution data will be added to it)
      hotspotdiff<-data2
      } # cls for(i in ...
    taxa_summary_data1<-endemdiff # a summary of the ED values for all resolutions for that single taxa 
    Taxa<-rep(names, each=nrow(taxa_summary_data1)) 
    taxa_summary_data1<-as.data.frame(cbind(Taxa, taxa_summary_data1))
    taxa_summary_data2<-hotspotdiff # a summary of the ED values for all resolutions for that single taxa 
    taxa_summary_data2<-as.data.frame(cbind(Taxa, taxa_summary_data2))
    
    cat("Calculating ED for the taxa ", names[i], "using the specified method ", "\n")
  } # cls if(class ....
  
  if(class(species_records)=="list"){ #if a list of data frames are inputted 
    cat("Cleaning the data using the specified method ", "\n")
    total_endemdiff<-data.frame(Taxa=character(), Cell_size=numeric(0), ED=numeric(0)) #create an empty dataframe (used later for the ED function) 
    total_hotspotdiff<-data.frame(Taxa=character(), Cell_size=numeric(0), HD=numeric(0)) #create an empty dataframe (used later for the ED function) 
    for(i in 1:length(species_records)){
      occurrence_records<-species_records[[i]]
      occurrence_records<-coordinate_cleaner(occurence_data=occurrence_records, method=outlier_method, crs=outlier_crs, 
                                             percent =outlier_percent, threshold = outlier_threshold,
                                             min_occs=min_occs, raster.res=outlier_raster.res,
                                             package_method=package_method, package_mltpl=package_mltpl, 
                                             package_tdi=package_tdi)
      cat("The data for taxa ", names[i], "has been cleaned ", "\n")
      endemdiff<-data.frame(Cell_size=numeric(0), ED=numeric(0))
      hotspotdiff<-data.frame(Cell_size=numeric(0), HD=numeric(0))
      if (minimal_method == "remove"){
        spec<-unique(occurrence_records[,1])
        for (i in spec) {
          spec_data<-occurrence_records[which(occurrence_records$species == i),]
          spec_data<-data.frame(LONGITUDE=spec_data$longitude, LATITUDE=spec_data$latitude)
          dat<-SpatialPoints(spec_data, CRS(outlier_crs))
          species_polygon<-try(gConvexHull(dat))
          if (class(species_polygon) != "SpatialPolygons") {
            occurrence_records<-occurrence_records[!(occurrence_records$species == i),]
          }
        }
      }
      if (minimal_method == "keep") {
        occurrence_records<-occurrence_records
      }
      for(x in 1:length(resolutions)){
        res_check<-as.character(isTRUE(all.equal(resolutions[x], res(frame.raster)[1])))
        if (res_check=="TRUE") {
          calc_raster<-frame.raster
        }
        if (res_check=="FALSE") {
          agg_factor<-(resolutions[x])/(res(frame.raster)[1])
          calc_raster<-aggregate(frame.raster, fact=agg_factor)
        }
        cat("Calculating weighted endemim using AOO with a raster resolution ", resolutions[x], "\n")  
        aoo<-weighted_endemism(species_records=occurrence_records, species=species,longitude=longitude,
                               latitude=latitude, frame.raster=calc_raster,
                               weight.type="cell", two_method=two_method, three_plus_method=three_plus_method, 
                               assigned_area=assigned_area, jitter_amount=jitter_amount) #calculate weighted endemism with AOO
        aooraster<-aoo$WE_raster #call upon and label the endemicm raster (map) for aoo
        if (scaling == "yes"){
          values(aooraster)<-rescale(values(aooraster), to = c(0, 1))
        }
        if (scaling == "no") {
          aooraster<-aooraster
        }
        cat("Calculating weighted endemim using EOO with a raster resolution ", resolutions[x], "\n")
        eoo<-weighted_endemism(species_records=occurrence_records, species=species,longitude=longitude,
                               latitude=latitude, frame.raster=calc_raster,
                               weight.type="geo", two_method=two_method, three_plus_method=three_plus_method, 
                               assigned_area=assigned_area, jitter_amount=jitter_amount) #calculate weighted endemism with EOO
        eooraster<-eoo$WE_raster #call upon and label the endemicm raster (map) for eoo
        if (scaling == "yes"){
          values(eooraster)<-rescale(values(eooraster), to = c(0, 1))
        }
        if (scaling == "no") {
          eooraster<-eooraster
        }
        res_data1<-endemism_difference(aooraster,eooraster,endemdiff, method=difference_method) #calculate the endemism difference between using AOO and EOO
        endemdiff<-res_data1 #this new data will be called upon in the ED function for the next resolution (the new resolution data will be added to it)
        res_data2<-hotspot_difference(aooraster,eooraster,hotspotdiff, loc_percent=hotspot_percent, loc_area=location_area)
        hotspotdiff<-res_data2
      } # cls for(i in res ...
      taxadata1<-endemdiff  #a summary of the ED values for all resolutions for that sequestial taxa
      taxaname1<-rep(names[i], each=nrow(taxadata1)) #create a vector of repeated characters (of the taxa name) the same length as the number of rows in summary of EDs (essentially the number of resolutions)
      taxadata2<-hotspotdiff  #a summary of the ED values for all resolutions for that sequestial taxa
      taxaname2<-rep(names[i], each=nrow(taxadata2)) #create a vector of repeated characters (of the taxa name) the same length as the number of rows in summary of EDs (essentially the number of resolutions)
      cat("Calculating ED for taxa ", names[i], "using the specified method ", "\n")
      num_check<-as.character(isTRUE(all.equal(species_records[[i]], species_records[[1]])))
      if (num_check=="TRUE") { #if it is the first taxa in the list of data frames
        newtaxadata1<-as.data.frame(cbind(taxaname1, taxadata1)) #create a dataframe from the indiviual taxa data adding a column of repeated taxa names
        newtaxadata2<-as.data.frame(cbind(taxaname2, taxadata2))
        colnames(newtaxadata1)<-c("Taxa", "Cell_size", "ED" )
        colnames(newtaxadata2)<-c("Taxa", "Cell_size", "HD" )
        prep_taxa_data1<-rbind(total_endemdiff, newtaxadata1)
        prep_taxa_data2<-rbind(total_hotspotdiff, newtaxadata2)
        taxa_data_set1<-prep_taxa_data1
        taxa_data_set2<-prep_taxa_data2
      } #cls if
      if (num_check=="FALSE") { # if it is noth the first taxa in the list of taxa/data frames 
        addition1<-as.data.frame(cbind(taxaname1, taxadata1)) #create a dataframe from the indiviual taxa data adding a column of repeated taxa names
        addition2<-as.data.frame(cbind(taxaname2, taxadata2))
        colnames(addition1)<-c("Taxa", "Cell_size", "ED" )
        colnames(addition2)<-c("Taxa", "Cell_size", "HD" )
        prep_taxa_data1<-rbind(taxa_data_set1, addition1)
        prep_taxa_data2<-rbind(taxa_data_set2, addition2)
        taxa_data_set1<-prep_taxa_data1 #create a new dataframe adding the data from the current taxa to the data of the previous taxa
        taxa_data_set2<-prep_taxa_data2 
        } #cls else
      taxa_summary_data1<-taxa_data_set1
      taxa_summary_data2<-taxa_data_set2
    }  # cls for(i in spec ...
  } # cls if(class ....
  ED_summary_data<-taxa_summary_data1
  HD_summary_data<-taxa_summary_data2
  ED_plot<- ggplot(ED_summary_data, aes(Cell_size, ED)) + geom_line(aes(group = Taxa), colour = "grey50") + geom_point(aes(colour = Taxa)) + ylab("Endemism Difference") + xlab("Cell size (km2)") + theme_bw()
  HD_plot<- ggplot(HD_summary_data, aes(Cell_size, HD)) + geom_line(aes(group = Taxa), colour = "grey50") + geom_point(aes(colour = Taxa)) + ylab("Hotspot Difference") + xlab("Cell size (km2)") + theme_bw()
  outputs<- list(ED_data=ED_summary_data, ED_plot=ED_plot, HD_data=HD_summary_data, HD_plot=HD_plot)
  return(outputs)
  cat("Calculations complete ", "\n")
}

weighted_endemism <- function(species_records=mydata, species="SPECIES", longitude="LONGITUDE", 
                              latitude="LATITUDE", frame.raster=myraster, weight.type="cell", 
                              jitter_amount=0.1, two_method="assign", three_plus_method="assign", 
                              assigned_area=1)
#Description --
  #This is a modified version of Guerin, Ruokolainen and Lowes (2015) endemism function that can be found at https://github.com/GregGuerin/biomap/blob/master/weighted.endemism.R
  #Calculates (taxonomic/species) weighted endemism (species richness inversely weighted by species ranges) across gridded maps
#Usage --
  #For example: weighted_endemism(species_records=mydata, species="SPECIES", longitude="LONGITUDE", latitude="LATITUDE", frame.raster=myraster, weight.type="geo", shapefile.clip=aus_coastline, minimal_method="keep", one_method="assign", two_method="assign", three_plus_method="assign", assigned_area=0.0001, raster.res=res, jitter_amount=0.1)
#Arguments --
  #species_records
    #A data.frame with rows as individual species records, and columns that include fields for species name, longitude and latitude (see species, longitude, latitude below)
  #species
   #What colname in the supplied species_records contains species names?
  #latitude
    #What colname in the supplied species_records contains latitude values?
  #longitude
    #What colname in the supplied species_records contains longitude values?
  #frame.raster
    #An existing raster used for calculations and mapping (with a specified res, extent and crs)
  #weight.type: 
    # ="cell" will calculate cell-based range weights (AOO)
    # ="geo" will calculate geographic range weights (EOO)
    # ="richness" sets the weights to 1, which is equivalent to calculating simple species richness
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
#Details --
  #This implementation of weighted endemism allows alternative calculation of weights for species ranges. Weights can be calculated based on the frequency of occurrence in grid cells (AOO), or alternatively by the geographical size of the species range (EOO). 
#Value --
  #Returns a list of length 4:
    #$WE: vector of weighted endemism scores
    #$WE_Raster: raster map with endemism scores
    #$weights: a named numeric vector of weights used to calculate endemism (equivalent to range size in square kilometres if weight.type="geo")
    #$grid.matrix : a binary data.frame of species against grid cell numbers used in the function which is returned so that it can be re-used in subsequent runs to save time
#Required packages --
  #sp, raster, regeos, simba
{
  require(sp)
  require(raster)
  require(rgeos)
  require(simba)
  
  cat("Checking species occurence data (column names, NA values, outliers from the frame raster etc.)", "\n")
  if(class(species_records) != "data.frame") {
    stop("Species data must be in a data.frame")
  } #checking that the occurence data is of class data.frame 
  colnames(species_records) <- c("SPECIES", "LATITUDE", "LONGITUDE") #setting the column names of the occurence data 
  if(!("SPECIES" %in% colnames(species_records))) {stop("Cannot locate species data")} #checking if the occurence data contains a column with species names 
  if(!("LATITUDE" %in% colnames(species_records))) {stop("Cannot locate latitude data")} #checking if the occurence data contains a column with latitude values
  if(!("LONGITUDE" %in% colnames(species_records))) {stop("Cannot locate longitude data")} #checking if the occurence data contains a column with longitude values
  if(any(is.na(species_records$LONGITUDE))) {
    species_records <- species_records[-which(is.na(species_records$LONGITUDE)),] 
  } #checking for NA values in the longitude column of the occurence data and removing them if present 
  if(any(is.na(species_records$LATITUDE))) {
    species_records <- species_records[-which(is.na(species_records$LATITUDE)),] 
  } #checking for NA values in the latitude column of the occurence data and removing them if present
  coordinates(species_records) <- c("LONGITUDE", "LATITUDE") # create a matrix of longitude and latitude values
  if(!(extent(species_records)@xmin >= extent(frame.raster)@xmin & extent(species_records)@xmax <= extent(frame.raster)@xmax & extent(species_records)@ymin >= extent(frame.raster)@ymin & extent(species_records)@ymax <= extent(frame.raster)@ymax)) {
    cat("Some point locations lie outside the frame raster -- trimming these records", "\n")
    species_record_COORDS <- as.data.frame(coordinates(species_records))
    species_records <- species_records[-which(species_record_COORDS$LONGITUDE < extent(frame.raster)@xmin | species_record_COORDS$LONGITUDE > extent(frame.raster)@xmax | species_record_COORDS$LATITUDE < extent(frame.raster)@ymin | species_record_COORDS$LATITUDE > extent(frame.raster)@ymax),]
  } # checking for and removing occurence points that lie outside the region of interest (as represented by the frame raster extent)
  cat("Species occurence data has been checked ", "\n")
  
  cat("Generating the gridded occurrence matrix", "\n")
  cell_numbers <- cellFromXY(frame.raster, species_records) #create a vector of the raster cell numbers that contain each occurence point (e.g. species_records[1,] occurs in cell_numbers[1] raster cell)
  cell_occur_matrix_prep <- data.frame(cell=cell_numbers, species=species_records$SPECIES, presence=rep(1, length(cell_numbers))) # create a data frame containing only presence data of each species in each cell 
  if(any(duplicated(cell_occur_matrix_prep))) {
    cell_occur_matrix_prep <- cell_occur_matrix_prep[-which(duplicated(cell_occur_matrix_prep)),]
  } #check and remove duplicated species cell occurences (e.g. if a cell contains two occurences from the same species, one will be removed)
  if(any(is.na(cell_occur_matrix_prep$cell))) {
    cell_occur_matrix_prep <- cell_occur_matrix_prep[-which(is.na(cell_occur_matrix_prep$cell)),]
  } #check and remove cells that contain no species occurences
  cell_occur_matrix <- mama(cell_occur_matrix_prep) #transpose the prep data frame so that species represent columns and cells represent rows, adding zeros to absent cells 
  cat("Occurrence matrix generated with dimensions: ", dim(cell_occur_matrix), "\n")
  
  cat("Calculating geographic range weights", "\n")
  if(weight.type=="cell") {
    cat("Calculating cell-based range weights", "\n")
    inv_rang_cell_occur_mat <- apply(cell_occur_matrix, 2, function(x) {x/sum(x)}) # divide each presence/absence value by the total number of raster cells
    ranges_prep <- apply(cell_occur_matrix, 2, function(x) {sum(x)}) # calculate the number of cells each species occurs in (their range)
    ranges<- ranges_prep*mean(values(area(frame.raster))) #multiply by cell size to get area in km^2
  } #cls if(w.t = cell...
  if(weight.type=="richness") {
    inv_rang_cell_occur_mat <- cell_occur_matrix
    ranges <- apply(cell_occur_matrix, 2, function(x) {x = 1}) # for each species assign range values of 1 (i.e. for richness the range sizes dont matter)
  } #cls if(w.t = richness...
  if(weight.type=="geo") {	
    crdref<-crs(frame.raster)
    spp_ranges <- function(x) {
      n <- 0
      v <- rep(0, length(colnames(cell_occur_matrix))) # create a vector of zeros with the length of the total number of species
      x$SPECIES <- sub(pattern = " ", replacement = ".", x = x$SPECIES, fixed=TRUE)
      x$SPECIES <- sub(pattern = "-", replacement = ".", x = x$SPECIES, fixed=TRUE)
      for (i in 1:length(colnames(cell_occur_matrix))) {
        n <- n + 1
        temp_prep<-as.data.frame(x)
        temp_prep <- temp_prep[which(temp_prep$SPECIES == colnames(cell_occur_matrix)[i]),] # create a dataframe with only i's occurences
        temp_prep <- data.frame(LONGITUDE=temp_prep$LONGITUDE, LATITUDE=temp_prep$LATITUDE)
        temp<-SpatialPoints(temp_prep, crdref)
        species_polygon<-try(gConvexHull(temp))
        if (class(species_polygon) != "SpatialPolygons") {
          if (nrow(temp_prep) == 1) {
            v[n]<- assigned_area
            warning("Cannot compute polygon (singular occurence point), returning the assigned area ", assigned_area, " for ", colnames(cell_occur_matrix)[i])
          } # cls if (nrow(temp_prep) == 1)
          if (nrow(temp_prep) == 2) {
            if (two_method == "assign") {
              v[n]<- 2*assigned_area
              warning("Cannot compute polygon (two occurence points), returning the assigned area ", assigned_area, " for ", colnames(cell_occur_matrix)[i])
            } #cls if (two_method == "assign")
            if (two_method == "span") {
              v[n]<- (pointDistance(temp_prep[1,], temp_prep[2,], lonlat=TRUE))/1000
              warning("Cannot compute polygon (two occurence points), returning the distance between the two points (span) for ", colnames(cell_occur_matrix)[i])
            } #cls if (two_method == "span")
            if (two_method == "assign_span") {
              v[n]<- assigned_area+((pointDistance(temp_prep[1,], temp_prep[2,], lonlat=TRUE))/1000)
              warning("Cannot compute polygon (two occurence points), returning a combination of the assigned area ", assigned_area, " and the distance between the two points (span) for ", colnames(cell_occur_matrix)[i])
            } #cls if (two_method == "span_assign")
          } # cls if (nrow(temp_prep) == 2)
          if (nrow(temp_prep) >= 3) {
            if (three_plus_method == "assign") {
              mulptli<-nrow(temp_prep)
              v[n]<- mulptli*assigned_area
            } #cls if (three_method == "assign")
            if (three_plus_method == "span") {
              coord <- cbind(temp_prep$LONGITUDE, temp_prep$LATITUDE)
              dobj <- dist(coord)
              dmat <- as.matrix(dobj)
              diag(dmat) <- NA
              dmax <- max(apply(dmat,2,max,na.rm=TRUE))
              loc<-which(dmat == dmax, arr.ind = TRUE)
              loc<-c(loc[1,1], loc[2,1])
              new_dat<-SpatialPoints(rbind(coord[loc[1],], coord [loc[2],]), crdref )
              span<-pointDistance(new_dat[1,], new_dat[2,], lonlat=TRUE)
              v[n]<-span/1000
            } #cls if (three_method == "span")
            if (three_plus_method == "assign_span") {
              mulptli<-nrow(temp_prep)
              coord <- cbind(temp_prep$LONGITUDE, temp_prep$LATITUDE)
              dobj <- dist(coord)
              dmat <- as.matrix(dobj)
              diag(dmat) <- NA
              dmax <- max(apply(dmat,2,max,na.rm=TRUE))
              loc<-which(dmat == dmax, arr.ind = TRUE)
              loc<-c(loc[1,1], loc[2,1])
              new_dat<-SpatialPoints(rbind(coord[loc[1],], coord[loc[2],]), crdref )
              span<-pointDistance(new_dat[1,], new_dat[2,], lonlat=TRUE)
              v[n]<-(mulptli*assigned_area)+(span/1000)
            } #cls if (three_method == "assign_span")
            if (three_plus_method == "jitter") {
              if (length(unique(as.numeric(temp_prep$LATITUDE))) == 1){
                random_row<-sample(as.numeric(1:nrow(temp_prep)),1)
                jittered_coord<-jitter(as.numeric(temp_prep[random_row,2]), amount = jitter_amount)
                temp_prep[random_row,2]<-jittered_coord
              } #cls if (length(unique(as.numeric(temp_prep$LAT ...
              if (length(unique(as.numeric(temp_prep$LONGITUDE))) == 1){
                random_row<-sample(as.numeric(1:nrow(temp_prep)),1)
                jittered_coord<-jitter(as.numeric(temp_prep[random_row,1]), amount = jitter_amount)
                temp_prep[random_row,1]<-jittered_coord
              } #cls if (length(unique(as.numeric(temp_prep$LONG ...
              temp<-SpatialPoints(temp_prep, crdref)
              species_polygon<-gConvexHull(temp)
              polygon_area <- (area(species_polygon))/1000000 #returns size in km^2 (convert from m^2)
              v[n]<-polygon_area
            }  #cls if (three_method == "jitter")
          } # cls if (nrow(temp_prep) >= 3)
        } #cls (class(species_polygon) != "SpatialPolygons")     
        if (class(species_polygon) == "SpatialPolygons") {
          polygon_area <- (area(species_polygon))/1000000 #returns size in km^2 (convert from m^2)
          v[n]<-polygon_area
        } #cls if(class(species_polygon) == "SpatialPo...
        cat(n, "species complete:", colnames(cell_occur_matrix)[i], v[n], "\n")
      } #cls for (i in colnames...
      names(v) <- colnames(cell_occur_matrix)
      return(v)
    } #cls spp_ranges function
    ranges <- spp_ranges(species_records) # calculating species ranges 
  } #cls if weight.type == geo...
  
  inv_rang_cell_occur_mat <- cell_occur_matrix
  for(i in 1:ncol(inv_rang_cell_occur_mat)) {
    inv_rang_cell_occur_mat[,i] <- inv_rang_cell_occur_mat[,i]/ranges[which(names(ranges) == colnames(inv_rang_cell_occur_mat)[i])]
  } #cls for(i in 1:ncol...
  
  cat("Calculating weighted endemism", "\n")
  rawEndemism <- rowSums(inv_rang_cell_occur_mat) 
  WE_raster <- frame.raster
  WE_raster[] <- NA
  WE_raster[as.numeric(names(rawEndemism))] <- rawEndemism
  outputs <- list(WE = rawEndemism, WE_raster = WE_raster, weights = ranges, grid.matrix = cell_occur_matrix)
  return(outputs)
  
} #cls function

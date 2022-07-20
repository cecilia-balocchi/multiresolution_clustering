clean_crime.from_csv <- function(file_path.name, data_path = "../../data/"){
  # the names of the columns are specific for the dataset release by the Philly police
  require(sp)
  load(paste0(data_path,"shapefile/phillyblockgroup"))
  load(paste0(data_path,"shapefile/phillyblock"))
  load(paste0(data_path,"shapefile/phillytracts"))
  crime.pre <- read.csv(paste0(data_path,file_path.name),
                        sep = ",",
                        header = T,
                        stringsAsFactors = F)
  # GPS location of crime:
  temp.pointX <- crime.pre$point_x
  temp.pointY <- crime.pre$point_y
  
  ## actually they have the same NAs, but just in case
  rem.ind.points <- is.na(temp.pointX) | is.na(temp.pointY)
  ## want to make this exact same as old file. Remove crap points, keep order:
  crimes.new <- crime.pre[!rem.ind.points,c('dc_dist', 'psa', 'ucr_general',
                                            'text_general_code', 'point_x', 'point_y')]
  ## crime type oh yes - the order matters (for crimeframe)
  crimetypes <- c('murder', 'rape', 'robbery', 'assault',
                  'burglary', 'theft', 'motortheft',
                  'otherassault', 'arson', 'forgery',
                  'fraud', 'embezzlement', 'receivestolen',
                  'vandalism', 'weaponviolation',
                  'prostitution', 'sexoffense', 'drugviolation',
                  'gambling', 'familyoffense', 'dui', 'liquorlaw',
                  'publicdrunk', 'disorderly', 'vagrancy',
                  'other')
  crimeframe <- data.frame(ucr = sort(unique(crime.pre$ucr_general)),
                           crimetype = crimetypes)
  # the following reorders the crimetype for the ucs in pre.crime, so we can add crimetype to crime.new
  temp.crimetype = crimeframe$crimetype[match(crime.pre$ucr_general, crimeframe$ucr)]
  crimes.new$crimetype = temp.crimetype[!rem.ind.points]
  ## there's a timedate col, but it's in 12h... so I won't use it
  temp.date = as.Date(crime.pre$dispatch_date, tz = 'EST')
  temp.time = crime.pre$dispatch_time
  crimes.new$date = temp.date[!rem.ind.points]
  temp.timedate = as.POSIXct(paste(temp.date, temp.time), tz = 'EST')
  crimes.new$timedate = temp.timedate[!rem.ind.points]
  ## add time with any given Sunday (only so that it has a Sunday zero point
  ## for week plotting)
  temp.time <- format(crimes.new$timedate, format = "%H:%M:%S", tx = 'EST')
  given.sunday <- '2016-05-08'
  crimes.new$time <- as.POSIXct(paste(given.sunday, temp.time),  tz = "America/New_York")
  # the time zone does not work
  ## add weekday
  crimes.new$weekday <- weekdays(crimes.new$timedate)
  
  
  ##------------------------------------##
  
  ## add block, blockgroup
  ## spatial situ
  pts <- SpatialPoints(crimes.new[ ,c('point_x', 'point_y')])
  ## now over - lists
  ## you'd save time by just doing blocks and then matching blocks to
  ## blockgroups, but not much time, so copying this code was easier
  ## than writing a simple match
  
  # WARNING: here we use the order of phillyblock, of phillyblockgroup 
  # and the one of tracts. If for some reason these orders are not 'consistent'
  # there could be mistakes
  # CHANGED INTO: save also the block/blockgroup/tract identifier!
  
  # over returns a list: each element of the list corresponds to a region
  # and contains a vector with the indices of the crimes happened in that region
  # so the length of that vector is the number of crimes in that region
  
  proj4string(pts) <- phillyblock@proj4string
  crimebyblock <- sp::over(phillyblock, pts, returnList = T)
  proj4string(pts) <- phillyblockgroup@proj4string
  crimebyblockgroup <- sp::over(phillyblockgroup, pts, returnList = T)
  crimebytracts <- sp::over(tracts, pts, returnList = T)
  
  ## we want to record which block a crime happens in, in the dataframe
  crimes.new$block <- NA
  crimes.new$blockID <- NA
  crimes.new$blockgroup <- NA
  crimes.new$blockgroupID <- NA
  crimes.new$tract <- NA
  crimes.new$tractID <- NA
  
  # The identifier we'll use is GEOID10 (a factor)
  
  ## replace with:
  crimevector <- unlist(crimebyblock)
  crimeblocklength <- lapply(crimebyblock, length)
  
  if(is.factor(phillyblock@data$GEOID10)){
    levels_region <- levels(phillyblock@data$GEOID10)[phillyblock@data$GEOID10]
  } else {
    levels_region <- phillyblock@data$GEOID10
  }
  
  crimeind.temp <- c(0, cumsum(crimeblocklength))
  for(i in 1:length(crimebyblock)){
    if((crimeind.temp[i+1]-crimeind.temp[i])>0){
      temp.index <- (crimeind.temp[i] + 1): (crimeind.temp[i + 1])
      crimes.new$block[crimevector[temp.index]] = i
      crimes.new$blockID[crimevector[temp.index]] = levels_region[i]
    } 
    
    if(i %% 100 == 0) cat(i, " ")
  }
  
  
  crimevector <- unlist(crimebyblockgroup)
  crimeblockgrouplength <- lapply(crimebyblockgroup, length)
  
  levels_region <- levels(phillyblockgroup@data$GEOID10)[phillyblockgroup@data$GEOID10] 
  
  crimeind.temp <- c(0, cumsum(crimeblockgrouplength))
  for(i in 1:length(crimebyblockgroup)){
    temp.index <- (crimeind.temp[i] + 1): (crimeind.temp[i + 1])
    crimes.new$blockgroup[crimevector[temp.index]] <- rep(i, length(temp.index))
    crimes.new$blockgroupID[crimevector[temp.index]] <- levels_region[i]
    if(i %% 100 == 0) cat(i, " ")
  }
  
  crimevector <- unlist(crimebytracts)
  crimetractlength <- lapply(crimebytracts, length)
  
  levels_region <- levels(tracts@data$GEOID10)[tracts@data$GEOID10]
  
  crimeind.temp <- c(0, cumsum(crimetractlength))
  for(i in 1:length(crimebytracts)){
    temp.index <- (crimeind.temp[i] + 1): (crimeind.temp[i + 1])
    crimes.new$tract[crimevector[temp.index]] = i
    crimes.new$tractID[crimevector[temp.index]] = levels_region[i]
    if(i %% 100 == 0) cat(i, " ")
  }
  
  ##------------------------------------##
  
  ## add year
  temp.year = format(crimes.new$date, '%Y')
  crimes.new$year <- as.numeric(temp.year)
  
  ## sort the whole thing
  crimes.new <- crimes.new[order(crimes.new$timedate),]
  
  ## a small fraction (863 or 0.04%) aren't in Philly (according to our shapefile)
  ## remove them
  ## (note: each crime either has both a block and blockgroup, or neither)
  crimes.new <- crimes.new[!is.na(crimes.new$block),]
  return(crimes.new)
}



aggregate_crime.from_rdata <- function(file_path.name, ncols = 12, variable = 'blockgroup'){
  require(dplyr)
  load(file_path.name)
  violent <- crimes.new[crimes.new$ucr_general <= 400 | crimes.new$ucr_general == 800 | crimes.new$ucr_general == 1700,]
  nonviolent <- crimes.new[crimes.new$ucr_general == 600 | crimes.new$ucr_general == 500,]
  
  v.table <- table(violent[,variable],violent$year)
  nv.table <- table(nonviolent[,variable],nonviolent$year)
  
 
  # Exclude year 2018 because it's only three months
  if(variable == 'block'){
    # load("data/shapefile/phillyblock")
    list_blocks <- sort(unique(as.numeric(c(rownames(v.table), rownames(nv.table)))))
    n_v <- nrow(v.table)
    n_nv <- nrow(nv.table)
    n_tot <- length(list_blocks)
    ## let's add rows to v.table
    v.table2 <- rbind(v.table[,1:ncols], array(0, dim = c(n_tot-n_v, ncols)))
    add_rows <- setdiff(list_blocks, as.numeric(rownames(v.table)))
    rownames(v.table2)[(n_v+1):n_tot] <- as.character(add_rows)
    ind <- sort(rownames(v.table2), index.return=T)$ix
    v.table2 <- v.table2[ind,]   
    ## let's add rows to nv.table
    nv.table2 <- rbind(nv.table[,1:ncols], array(0, dim = c(n_tot-n_nv, ncols)))
    add_rows <- setdiff(list_blocks, as.numeric(rownames(nv.table)))
    rownames(nv.table2)[(n_nv+1):n_tot] <- as.character(add_rows)
    ind <- sort(rownames(nv.table2), index.return=T)$ix
    nv.table2 <- nv.table2[ind,]   
    v.table <- v.table2; nv.table <- nv.table2
  } 
  colnames(v.table) <- paste0("v_",colnames(v.table))
  colnames(nv.table) <- paste0("nv_",colnames(nv.table))
  total <- cbind(v.table[,1:ncols],nv.table[,1:ncols])
  n <- dim(v.table)[1]
  total <- cbind(as.numeric(rownames(total)),total)
  colnames(total)[1] <- variable
  total <- as.data.frame(total)
  
  varID <- paste0(variable,"ID")
  total$GEOID10 <- violent[match(total[,1], violent[,variable]),varID]
  if(variable == 'block'){
    total$GEOID10[is.na(total$GEOID10)] <- nonviolent[match(total[is.na(total$GEOID10),1], nonviolent[,variable]),varID]
    load("data/shapefile/phillyblock")
    
    GEOID_ordered <- levels(phillyblock@data$GEOID10)[phillyblock@data$GEOID10]
    
    add_geoid <- setdiff(GEOID_ordered, total$GEOID10)
    tmp_ind <- match(add_geoid, GEOID_ordered)
    add_block <- as.numeric(rownames(phillyblock@data)[tmp_ind])+1
    # add_block <- setdiff(1:length(phillyblock), total[1:n_tot,variable]) #WRONG
    total[(n_tot+1):length(phillyblock),"GEOID10"] <- add_geoid
    total[(n_tot+1):length(phillyblock),2:(1+2*ncols)] <- 0
    total[(n_tot+1):length(phillyblock),variable] <- add_block
  }
  
  c_y <-total
  # c_y <-tbl_df(total)
  t1 <- c_y[rep(1:dim(c_y)[1],ncols),]
  t2 <- arrange_(t1, .dots = variable)
  t3 <- mutate(t2,year = rep(2006:(ncols+2005),dim(c_y)[1]))
  t4 <- mutate(t3,violent = NA)
  for(y in 2006:(ncols+2005)){
    str <- paste0("v_",y)
    t4[t4$year==y,"violent"] <- c_y[,str]
    # t4[t4$year==y,"violent"] <- c_y[,str][[1]]
  }
  t5 <- mutate(t4, nonviolent = NA)
  for(y in 2006:(ncols+2005)){
    str <- paste0("nv_",y)
    t5[t5$year==y,"nonviolent"] <- c_y[,str]
    # t5[t5$year==y,"nonviolent"] <- c_y[,str][[1]]
  }
  years <- "year"
  violents <- "violent"
  nonviolents <- "nonviolent"
  geo <- "GEOID10"
  t5 <- dplyr::select_(t5,.dots = c("variable","geo","years","violents","nonviolents"))
  # t5[,variable][[1]] <- as.factor(t5[,variable][[1]])
  # t5[,geo][[1]] <- as.factor(t5[,geo][[1]])
  t5[,variable] <- as.factor(t5[,variable])
  t5[,geo] <- as.factor(t5[,geo])
  t6 <- mutate_(t5, .dots = setNames(list(variable), "X"))
  t6 <- dplyr::select(t6,X,GEOID10,year,violent,nonviolent)
  d6 <- as.data.frame(t6)
  
  
  return(d6)
}


# Inverse Hyperbolic Sine Transformation - it's a good approximation of log() but it's well defined at 0.
ihst <- function(y){
  log(y+(y^2+1)^0.5)-log(2)
}
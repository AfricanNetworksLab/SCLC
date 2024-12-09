# Calculate a network metric per unit area
# for the Spatial Conflict Life Cycle paper
# David Russell, August 2023


###### Read in libraries ######
library(tidyverse)
library(sf)
library(lubridate)


###### Read in data and clean ######
setwd("C:/Users/drussell/Dropbox (UFL)/SWAC2023-24/")


# Fishnet
fishnet <- st_read("gis/fishnet/All_Africa_fishnet_50k_proj.shp")
fishnet <- fishnet[!duplicated(fishnet$geometry),] # Get rid of duplicates
fishnet$ID <- 1:nrow(fishnet) # Make an unique ID column

# Dyadic event dataset
dyads <- read.csv("data/Networks/dyadsJune2023AFRICA.csv")

# Remove imprecisely-coded events
dyads <- dyads[dyads$geo_precision != 3,]

# Remove duplicated dyads
# This should already be done when updating actor names
# But this does it too, just in case
dyads$idAndDyadName <- paste0(dyads$dataIDs,dyads$dyadNames,dyads$coop)
dyads <- dyads[!duplicated(dyads$idAndDyadName),]

dyads$date <- dmy(dyads$event_date)

# Pare dataframe down to just essential columns
dyads <- dyads[,c("dataIDs","coop","actorsX.new","actorsY.new",
                  "latitude","longitude","date","country","year")]


###### Create sf objects and project data ######
### Prepare points data
# Indicate what the column names for longitude and latitude are
coords <- c("longitude","latitude")

# Indicate the EPSG code for the projection of the data
# (4326) is unprojected WGS84 data
crs = 4326

# Turn this into an sf object
dyads.unproj <- st_as_sf(dyads, coords = coords, crs = crs)

# Project data into Africa Lambert Conformal Conic
crsString <- "+proj=lcc +lat_1=20 +lat_2=-23 
              +lat_0=0 +lon_0=25 +x_0=0 
              +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
dyads.proj <- st_transform(dyads.unproj, crs = crsString)

# Remove unnecessary objects from memory to speed computation
rm(dyads.unproj,coords,crs,crsString)

###### Use st_join() and summarize() ######

SCLCsum <- function(dyads.proj,fishnet){
  # Point in polygon overlay (using st_join)
  dyadsInFishnet <- st_join(dyads.proj,fishnet)
  
  dyadsInFishnet <- as.data.frame(dyadsInFishnet)
  dyadsInFishnet <- subset(dyadsInFishnet,select = -c(geometry))

  # Summarize event counts, actor counts, negative ties, and positive ties
  # By fishnet cell
  dyadsInFishnet$eventCounter <- rep(1,nrow(dyadsInFishnet))
  summaryTable <- dyadsInFishnet %>% group_by(ID) %>%
    summarize(eventCount = n_distinct(dataIDs),
              actorCount = n_distinct(actorsY.new,actorsX.new),
              positiveTies = sum(eventCounter[coop == TRUE]),
              negativeTies = sum(eventCounter[coop == FALSE]))
  
  # Wherever there is a cell with only one actor,
  # recode to NA
  summaryTable$eventCount[summaryTable$actorCount == 1] <- 0
  summaryTable$positiveTies[summaryTable$actorCount == 1] <- 0
  summaryTable$negativeTies[summaryTable$actorCount == 1] <- 0
  summaryTable$actorCount[summaryTable$actorCount == 1] <- 0
  # Note how actors col has to come last...
  
  # Calculate percentage negative ties
  summaryTable$percNeg <- summaryTable$negativeTies / 
    (summaryTable$negativeTies + summaryTable$positiveTies) * 100
  
  # Make sure there is an entry for each column
  outputDF <- data.frame("ID" = fishnet$ID)
  outputDF <- merge(outputDF,summaryTable,by.x = "ID",by.y = "ID",all = T)
  
  return(outputDF)
}

SCLCoverTime <- function(points, polygons, timeColumn){
  finalDF <- data.frame("ID" = polygons$ID)
  for(time in min(timeColumn):max(timeColumn)){
    print(time)
    timePoints <- points[timeColumn == time,]
    print(nrow(timePoints))
    timeOutput <- SCLCsum(timePoints,polygons)
    for(col in 1:length(colnames(timeOutput))){
      colnames(timeOutput)[col] <- paste0(colnames(timeOutput)[col],as.character(time))
    }
    finalDF <- merge(finalDF,timeOutput, by.x = "ID", by.y = paste0("ID",time), all = T)
    print(paste0("Finished processing ",as.character(time)))
  }
  return(finalDF)
}

# Run it
SCLC <- SCLCoverTime(dyads.proj,fishnet,dyads.proj$year)

SCLC <- SCLC[!is.na(SCLC$ID),]

actorColumnsDF <- SCLC[,startsWith(colnames(SCLC),"actorCount")]
SCLC$actorSums <- rowSums(actorColumnsDF,na.rm = T)
SCLC$actorLocalAvg <- SCLC$actorSums / ncol(actorColumnsDF)

resultsDF <- data.frame("ID" = fishnet$ID)
for(i in 1:ncol(actorColumnsDF)){
  newCol <- NA
  newCol[actorColumnsDF[,i] > SCLC$actorLocalAvg] <- "High"
  newCol[actorColumnsDF[,i] <= SCLC$actorLocalAvg] <- "Low"
  resultsDF <- cbind(resultsDF,newCol)
  colnames(resultsDF)[ncol(resultsDF)] <- paste0(colnames(actorColumnsDF)[i],"rel")
}

SCLC <- merge(SCLC,resultsDF,by = "ID")

##### Write it out #####

#write.csv(SCLC,"data/Networks/SCLCoutputData/SCLCJun2023.csv",row.names = F)
write.csv(SCLC,"data/Networks/SCLCoutputData2/SCLCJun2023.csv",row.names = F)
fishnetData <- merge(fishnet,SCLC,by.x = "ID",by.y = "ID",all.x = T)
#st_write(fishnetData,"data/Networks/SCLCoutputData/SCLCJun2023.shp",append = F)
st_write(fishnetData,"data/Networks/SCLCoutputData2/SCLCJun2023.shp",append = F)



###### Classify observations by comparing to previous year ######

SCLC <- read.csv("data/Networks/SCLCoutputData/SCLCJun2023.csv")

SCLC[is.na(SCLC)] <- 0

# Define function
diffFunction <- function(inputDF,inputType){
  inputCols <- inputDF[,startsWith(colnames(inputDF),inputType)]
  outputDF <- data.frame("ID" = inputDF$ID)
  for(i in 2:ncol(inputCols)){
    diffVector <- inputCols[,i] - inputCols[,i-1]
    diffVector[inputCols[,i] == 0] <- NA
    
    col1Name <- colnames(inputCols)[i]
    col2Name <- colnames(inputCols)[i-1]
    col1Year <- substr(col1Name,nchar(col1Name) - 3, nchar(col1Name))
    col2Year <- substr(col2Name,nchar(col2Name) - 3, nchar(col2Name))
    diffVectorName <- paste0("diff_",inputType,col2Year,"_",col1Year)
    classVectorName <- paste0("class_",inputType,col2Year,"_",col1Year) 
    
    diffClass <- ifelse(diffVector > 0,
                        paste0("Increase_",inputType),
                        paste0("Decrease_",inputType))
    outputDF <- cbind(outputDF,diffVector,diffClass)
    colnames(outputDF)[(ncol(outputDF)-1):(ncol(outputDF))] <- 
      c(diffVectorName,classVectorName)
  }
  return(outputDF)
}

# Event counts
eventCountDiff <- diffFunction(SCLC,"eventCount")
#write.csv(eventCountDiff,"data/Networks/SCLCoutputData/eventCountDiff.csv")
write.csv(eventCountDiff,"data/Networks/SCLCoutputData2/eventCountDiff.csv")

# Percent negative
percNegDiff <- diffFunction(SCLC,"percNeg")
#write.csv(percNegDiff,"data/Networks/SCLCoutputData/percNegDiff.csv")
write.csv(percNegDiff,"data/Networks/SCLCoutputData2/percNegDiff.csv")


# Negative ties
negCountDiff <- diffFunction(SCLC,"negativeTies")
#write.csv(negCountDiff,"data/Networks/SCLCoutputData/negCountDiff.csv")
write.csv(negCountDiff,"data/Networks/SCLCoutputData2/negCountDiff.csv")


# Actor counts
actorCountDiff <- diffFunction(SCLC,"actorCount")
#write.csv(actorCountDiff,"data/Networks/SCLCoutputData/actorCountDiff.csv")
write.csv(actorCountDiff,"data/Networks/SCLCoutputData2/actorCountDiff.csv")


### Concatenate these
actorCountDiff <- read.csv("data/Networks/SCLCoutputData/actorCountDiff.csv")
percNegDiff <- read.csv("data/Networks/SCLCoutputData/percNegDiff.csv")


# For each pair of rows, concatenate them and give them an appropriate name
justActorDiff <- actorCountDiff[,startsWith(colnames(actorCountDiff),"class")]
justpercNegDiff <- percNegDiff[,startsWith(colnames(percNegDiff),"class")]

outputDF <- data.frame("ID" = actorCountDiff$ID)
for(i in 1:ncol(justActorDiff)){
  outputVector <- paste0(justActorDiff[,i],"_",justpercNegDiff[,i])
  outputDF <- cbind(outputDF,outputVector)
  colName <- colnames(justActorDiff)[i]
  yearName <- substr(colName,nchar(colName) - 3, nchar(colName))
  colnames(outputDF)[ncol(outputDF)] <- paste0("SCLC",yearName)
}
#write.csv(outputDF,"data/Networks/SCLCoutputData/actor_percNeg.csv",row.names = F)
write.csv(outputDF,"data/Networks/SCLCoutputData2/actor_percNeg.csv",row.names = F)
fishnetData <- merge(fishnet,outputDF,by.x = "ID",by.y = "ID",all.x = T)
#st_write(fishnetData,"data/Networks/SCLCoutputData/actor_perNeg.shp",append = F)
st_write(fishnetData,"data/Networks/SCLCoutputData2/actor_perNeg.shp",append = F)


###### Check all of this ######

# Event counts
eventCountDiff <- read.csv("data/Networks/SCLCoutputData/eventCountDiff.csv")

# Percent negative
percNegDiff <- read.csv("data/Networks/SCLCoutputData/percNegDiff.csv")

# Negative ties
negCountDiff <- read.csv("data/Networks/SCLCoutputData/negCountDiff.csv")

# Actor counts
actorCountDiff <- read.csv("data/Networks/SCLCoutputData/actorCountDiff.csv")

# Actor and Percent Neg
sclcOutput <- read.csv("data/Networks/SCLCoutputData2/actor_percNeg.csv")

nrow(sclcOutput[sclcOutput == "Increase_actorCount_NA"])
# Row 6070 SCLC2011
# Row 6061 SCLC2014
# Row 1293 SCLC2019
# Row 8016 SCLC2018
#nrow(sclcOutput[sclcOutput == "Decrease_actorCount_NA"])
#nrow(sclcOutput[sclcOutput == "NA_Increase_percNeg"])
#nrow(sclcOutput[sclcOutput == "NA_Decrease_percNeg"])

sclcOutput %>% filter_all(any_vars(. %in% c("Increase_actorCount_NA")))
sclcOutput %>% filter_all(any_vars(. %in% c("Decrease_actorCount_NA")))
sclcOutput %>% filter_all(any_vars(. %in% c("NA_Increase_percNeg")))
sclcOutput %>% filter_all(any_vars(. %in% c("NA_Decrease_percNeg")))


# I read SCLC data, dyads, and fishnet back in,
# and reproduced dyadsInFishnet dataframe (spatial join)
# using part of the SCLCsum function
View(SCLC[c(6070,6061,1293,8016),])

View(dyadsInFishnet[(dyadsInFishnet$ID  == 6070 & dyadsInFishnet$year == 2011) |
                      (dyadsInFishnet$ID  == 6061 & dyadsInFishnet$year == 2014) |
                      (dyadsInFishnet$ID  == 1293 & dyadsInFishnet$year == 2019) |
                      (dyadsInFishnet$ID  == 8016 & dyadsInFishnet$year == 2018),
                    ])
# Read in ACLED too
ACLED <- read.csv("data/ACLED/MostRecentRawData/ACLED_2000_2023_AFRICA.csv")

View(ACLED[ACLED$event_id_cnty %in% c("BFO1419","LBY283","LBY2775","ETH4603"),])
View(ACLED[is.na(ACLED$actor2) | ACLED$actor2 == "",]) #2,257 instances (of 178,986)
summary(factor(ACLED$event_type[is.na(ACLED$actor2) | ACLED$actor2 == ""])) # All explosions/rv

View(ACLED[ACLED$assoc_actor_1 == "" & ACLED$actor2 == "" & 
             ACLED$assoc_actor_2 == "",])
View(ACLED[(ACLED$assoc_actor_1 == "" & ACLED$actor2 == "" & 
             ACLED$assoc_actor_2 == "") |
             (is.na(ACLED$assoc_actor_1) & is.na(ACLED$actor2) & 
                is.na(ACLED$assoc_actor_2)),])

# Remove 4 events
sclcOutput <- read.csv("data/Networks/SCLCoutputData/actor_percNeg.csv")

sclcOutput[sclcOutput == "Increase_actorCount_NA"] <- "NA_NA"
sclcOutput[sclcOutput == "Decrease_actorCount_NA"] <- "NA_NA"
sclcOutput[sclcOutput == "NA_Decrease_percNeg"] <- "NA_NA"
sclcOutput[sclcOutput == "NA_Increase_percNeg"] <- "NA_NA"
sclcOutput[sclcOutput == "NANA"] <- "NA_NA"



write.csv(sclcOutput,"data/Networks/SCLCoutputData/actor_percNeg.csv",row.names = F)
fishnetData <- merge(fishnet,sclcOutput,by.x = "ID",by.y = "ID",all.x = T)
st_write(fishnetData,"data/Networks/SCLCoutputData/actor_perNeg.shp",append = F)

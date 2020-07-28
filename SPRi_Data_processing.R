
#################################################

'SPRI DATA ANALYSIS PIPELINE'

##################################################
"

Quick overview of the new data munging process:
  1 - Average all the spots together
  2 - Shift data by assing all numbers with the
  absolute of the most negative (minimum value of
  the data) such that the most negative one will
  become zero
  3 - unit-based normalization of the data
  4 - Take the derivative

Quick overview of the data analysis process
  1 - Determine the Rmax value and index
  2 - Find the kd value using nls
  3 - Find the ka value using nls
  4 - Calculate kD=kd/ka
"

########## Dependencies ##########
library(tidyr)
library(readr)
library(readxl)
library(dplyr)

########## Functions ##########
# Function converts 00:00:00 format to
# min.sec format
changeTime = function(hourMinSec) {
  "
  Function takes a character vector of the
  format 00:00:00
  Returns the time as a integer expressed
  in minutes
  "
  hours = as.numeric(substr(hourMinSec, 1, 2))
  minutes = as.numeric(substr(hourMinSec, 4, 5))
  seconds = as.numeric(substr(hourMinSec, 7, 8))

  hours = hours*60 # Convert hours to minutes
  seconds = seconds/60 # Convert seconds to minutes

  return(hours+minutes+seconds)
}

# Function
# Accepts file path (either relative or absolute)
# Returns list of the injection times
processedFile = function(filepath) {
  "
  Function reads in a file in one line at a time
  reading only the one line at a time into memory.
  Returns a list of times that correspond to
  positive injections
  "
  fileObj = file(filepath, "r") # File object
  i <- 1 # counter
  timeList <- list() # Main list to fill

  # Read a line from the file, on at a time
  while ( TRUE ) {
    line = readLines(fileObj, n = 1) # Read in line
    if ( length(line) == 0 ) {
      break # Stop once end of file is reached
    }
    # Capture only those lines with relevant injection times
    if ( grepl("Valve in Injection position", line)) {
      # If TRUE find match for time at the beginning of the line
      timeRegex <- regmatches(line, regexpr("?^[0-9]*:[0-9]*:[0-9]*", line))
      timeConvert <- changeTime(timeRegex) # Convert time to minutes only format
      timeList[[i]] <- timeConvert # Append converted tiem to the list
      i <- i + 1 # Increase counter
    }
  }
  return(timeList) # Return list
  close(fileObj)
}

# Function
# Accepts file path (either relative or absolute) and a list of injection times
# Returns list where keys are injection times and values are concentrations
extractConcentrations <- function(filepath, injectList) {
  "
  Function reads in a file in one line at a time
  reading only the one line at a time into memory.
  Returns a list of concentrations where each item
  corresponds to the concentration of that injection
  "
  fileObj = file(filepath, "r") # File object
  concList <- list() # List of concentrations used

  # Read a line from the file, on at a time
  while ( TRUE ) {
    line = readLines(fileObj, n = 1) # Read in line
    if ( length(line) == 0 ) {
      break # Stop once end of file is reached
    }
    # Capture only those lines with relevant injection times
    if ( grepl("Valve in Injection position", line)) {
      # If TRUE find match for time at the beginning of the line
      timeRegex <- regmatches(line, regexpr("?^[0-9]*:[0-9]*:[0-9]*", line))
      timeConvert <- changeTime(timeRegex) # Convert time to minutes only format

      if( timeConvert %in% injectList ) {
        # If the converted time is in the list of injection times extract
        # concentration information from the end of the line
        concRegx <- regmatches(line, regexpr("?\\([0-9]*[A-z]{1,2}\\s+.*\\)$", line))
        timeConvert <- as.character(timeConvert) # Convert time from numeric to character class. Otherwise list will be > 200 items long
        concList[[timeConvert]] <- concRegx # Append concentration to list with injection time as the key
      }
    }
  }
  return(concList) # Return the list of concentrations
  close(fileObj)
}

# Unit based normalization
# Need to loop through each column in the df with the exception of time
normalize <- function(x) {
  "
  performs unit-based normalization
  "
  return ((x - min(x)) / (max(x) - min(x)))
}

########## Read/format the Data ##########
# Path variables for the three data files needed
kinCommentsPath = "./Pick_Me/Kinetics_Comments_test.txt"
kinDataPath = "./Pick_Me/Kinetics_Data_test.xlsx"
spotFilePath = "./Pick_Me/Spot_File_test.xlsx"

# Named list of negatives and the spots they apply to
# Hard coded and will have to be updated for each new dataset but it's the best I can do in R
# On the plus side the loop which subtracts negatives fromtheir corresponding spots can reamin the same.
negList = list(I=c('A', 'D', 'F'), J=c('B', 'E', 'G'), K=c('C', 'H'))

# Function call -> extracts injection times and converts the format. saves as list.
kinCommentsTime = processedFile(kinCommentsPath)

# Function call -> extracts the concentration for each injection time. saves as list.
kinCommmentsConc = extractConcentrations(kinCommentsPath, kinCommentsTime)
print(kinCommmentsConc)

# Read in spot file and kinetic data file
spotFile =  read_excel(spotFilePath, col_names = TRUE)
kinFile = read_excel(kinDataPath, col_names = TRUE)

# Create column names list from the spot file
colList = list("Time") # List needs to start with 'Time' as the first column
spotNamesList = as.list(spotFile$Species) # Convert rows in 'Species' column to a list
colList = append(colList, spotNamesList) # Concatenate lists

# Label the columns in the Kinetic df object
colnames(kinFile) <- colList
kinDF <- as.data.frame(kinFile) # Convert to df

dfRowList = list() # List to hold rows in the range of the times in kinCommentsTime + 9 (i:i+9)
i <- 1 # Counter
# Iterate all rows in the dataframe
for (rows in 1:nrow(kinDF)) {
  timeRow <- kinDF[rows, "Time"] # Save the integer in the 'Time' column
  # Iterate the times kinCommentsTime list
  for (num in kinCommentsTime) {
    # Look for integer time values from the rows in the range of time:time+9
    if ( timeRow >= num & timeRow <= round(num+9.00, 2) ) {
      dfRowList[[i]] <- kinDF[rows, ] # Append entire row from df if integer in time column falls inside range
      i <- i + 1 # Increase counter
    }
  }
}

# Use rbind to form the rows in list into a df
spriDF <- do.call(rbind.data.frame, dfRowList)
rownames(spriDF) <- 1:nrow(spriDF) # reset rownames index back to 1
spriDF <- distinct(spriDF)

# Seperate different 9 minute time segments and save these as new df in a list
spriDfList <- list() # Empty list to hold df associated with each injection time
for (num in kinCommentsTime) {
  startTime <- which(spriDF$Time == num) # Index of injection start time
  endTime <- which(spriDF$Time == round(num+9.00, 2))  # Index of injection end time. Need two decimal place or weird things start to happen
  injection <- spriDF[startTime:endTime, ] # Extract injection from DF
  key <- paste('Time-', num, collapse = ',') # Create key name for DF
  spriDfList[[key]] <- injection # Save extracted DF to list with associated key name
}

########## Data Munging ##########
# Main data processing loop
# Loops over each dataframe in the spriDFList object (where each df corresponds to am injection time+9)
spriProcessed <- list() # Empty list to hold all the processed data
i <- 0 # Counter
for (df in spriDfList) {
  if (i == length(spriDfList)) {
    break # break the loop once the last df in the list is processed
  }
  # Save the time Column as a seperate vector, but do not remove it
  # Removing it before the average is taken will cause the column names to become unique, which we don't want
  timeCol <- subset(df, select = 'Time') # Extract 'Time' and save as seperate column

  # Transpose values so that all the values in the df are positive
  colMin <- which.min(apply(df, MARGIN = 2, min)) # Return the column index of the largest negative value
  rowMin <- which.min(apply(df, MARGIN = 1, min)) # Return the row index of the largest negative value
  maxMin <- df[rowMin, colMin] # location of largest negative value
  df <- df + abs(maxMin) # Increase all values by the absolute value of the largest negative integer in df

  # Take the Average of row values in like columns
  df <- sapply(split.default(df, names(df)), rowMeans)
  df <- as.data.frame(df)

  # Can now drop the time column as it has been altered by the the data munging process anyhow
  # The Time vecotor saved earlier will be added in once the derivative is taken
  drops <- c("Time")
  df <- df[ , !(names(df) %in% drops)] # Drop the 'Time' column from main df

  # Subtract the negative from the appropriate columns
  for (colName in names(df)) { # Iterate the columns in the df
    for (neg in 1:length(negList)) { # Iterate the names in the named list. The names are the names of the -ve columns
      if ( colName %in% negList[[neg]] ) { # Check if the column name under the i-th name in the column list
        negName = as.character(names(negList[neg])) # Get name of negative column
        df[[colName]] <- df[[colName]] - df[[negName]] # Subtract the negative column from the positive
      }
    }
  }

  # Remove the negative columns from the df
  drops <- c(names(negList))
  df <- df[ , !(names(df) %in% drops)]

  # Normalize the values in the df to bewteen 0 and 1
  df <- normalize(df) # Use the normalize function

  # Take the derivative of each column.
  # This will also require a new column for each
  originalColCount <- ncol(df) # number of columns
  originalColNames <- list(colnames(df)) # names of all the columns in the df
  j <- 1 # column counter
  for (cols in names(df)) {
    if (j > originalColCount) {
      break # end loop if col limit is reached
    }
    newColName <- paste(originalColNames[[1]][j], "Deriv", sep = "_") # Create column name
    newColData <- list(0) # First number in our new col needs to be zero to make the columns same length of original df
    colData <- as.list(diff(df[[cols]], lag = 1))
    df$V1 <- cbind(unlist(append(newColData, colData))) # Combine the two lists and bind to df as a new column
    names(df)[names(df) == "V1"] <- newColName # Rename the column with a descriptive name

    j <- j + 1 # Increase counter
  }

  # Add Time Column back into the df
  df$Time <- seq(from = 0, to = 540, by = 3)

  i <- i + 1 # Increases counter
  key = "key" # Generic name
  spriProcessed[[key]] <- df # Save processed df as df in list
  names(spriProcessed)[names(spriProcessed) == "key"] <- names(spriDfList)[i] # Rename df with descriptive name
}

# TODO: Work code in from kD_largest_rMax_mulitple_columns
# TODO: in the Kinetic comments, on the injection lines it lists what's being injected as well as the
#       concentration. Could extract this and put into a named list. if this was able to correspond, or could list
#       the time in the name for each entry, could use this to iterate the dfs and generate results
=======
#################################################

'SPRI DATA ANALYSIS PIPELINE'

##################################################
"
Quick overview of the new data munging process:
  1 - Average all the spots together
  2 - Shift data by assing all numbers with the
  absolute of the most negative (minimum value of
  the data) such that the most negative one will
  become zero
  3 - unit-based normalization of the data
  4 - Take the derivative

Quick overview of the data analysis process
  1 - Determine the Rmax value and index
  2 - Find the kd value using nls
  3 - Find the ka value using nls
  4 - Calculate kD=kd/ka
"

########## Dependencies ##########
library(tidyr)
library(readr)
library(readxl)
library(dplyr)

########## Functions ##########
# Function converts 00:00:00 format to
# min.sec format
changeTime = function(hourMinSec) {
  "
  Function takes a character vector of the
  format 00:00:00
  Returns the time as a integer expressed
  in minutes
  "
  hours = as.numeric(substr(hourMinSec, 1, 2))
  minutes = as.numeric(substr(hourMinSec, 4, 5))
  seconds = as.numeric(substr(hourMinSec, 7, 8))

  hours = hours*60 # Convert hours to minutes
  seconds = seconds/60 # Convert seconds to minutes

  return(hours+minutes+seconds)
}

# Function
# Accepts file path (either relative or absolute)
# Returns list of the injection times
processedFile = function(filepath) {
  "
  Function reads in a file in one line at a time
  reading only the one line at a time into memory.
  Returns a list of times that correspond to
  positive injections
  "
  fileObj = file(filepath, "r") # File object
  i <- 1 # counter
  timeList <- list() # Main list to fill

  # Read a line from the file, on at a time
  while ( TRUE ) {
    line = readLines(fileObj, n = 1) # Read in line
    if ( length(line) == 0 ) {
      break # Stop once end of file is reached
    }
    # Capture only those lines with relevant injection times
    if ( grepl("Valve in Injection position", line)) {
      # If TRUE find match for time at the beginning of the line
      timeRegex <- regmatches(line, regexpr("?^[0-9]*:[0-9]*:[0-9]*", line))
      timeConvert <- changeTime(timeRegex) # Convert time to minutes only format
      timeList[[i]] <- timeConvert # Append converted tiem to the list
      i <- i + 1 # Increase counter
    }
  }
  return(timeList) # Return list
  close(fileObj)
}

# Function
# Accepts file path (either relative or absolute) and a list of injection times
# Returns list where keys are injection times and values are concentrations
extractConcentrations <- function(filepath, injectList) {
  "
  Function reads in a file in one line at a time
  reading only the one line at a time into memory.
  Returns a list of concentrations where each item
  corresponds to the concentration of that injection
  "
  fileObj = file(filepath, "r") # File object
  concList <- list() # List of concentrations used

  # Read a line from the file, on at a time
  while ( TRUE ) {
    line = readLines(fileObj, n = 1) # Read in line
    if ( length(line) == 0 ) {
      break # Stop once end of file is reached
    }
    # Capture only those lines with relevant injection times
    if ( grepl("Valve in Injection position", line)) {
      # If TRUE find match for time at the beginning of the line
      timeRegex <- regmatches(line, regexpr("?^[0-9]*:[0-9]*:[0-9]*", line))
      timeConvert <- changeTime(timeRegex) # Convert time to minutes only format

      if( timeConvert %in% injectList ) {
        # If the converted time is in the list of injection times extract
        # concentration information from the end of the line
        concRegx <- regmatches(line, regexpr("?\\([0-9]*[A-z]{1,2}\\s+.*\\)$", line))
        timeConvert <- as.character(timeConvert) # Convert time from numeric to character class. Otherwise list will be > 200 items long
        concList[[timeConvert]] <- concRegx # Append concentration to list with injection time as the key
      }
    }
  }
  return(concList) # Return the list of concentrations
  close(fileObj)
}

# Unit based normalization
# Need to loop through each column in the df with the exception of time
normalize <- function(x) {
  "
  performs unit-based normalization
  "
  return ((x - min(x)) / (max(x) - min(x)))
}

########## Read/format the Data ##########
# Path variables for the three data files needed
kinCommentsPath = "./Pick_Me/Kinetics_Comments_test.txt"
kinDataPath = "./Pick_Me/Kinetics_Data_test.xlsx"
spotFilePath = "./Pick_Me/Spot_File_test.xlsx"

# Named list of negatives and the spots they apply to
# Hard coded and will have to be updated for each new dataset but it's the best I can do in R
# On the plus side the loop which subtracts negatives fromtheir corresponding spots can reamin the same.
negList = list(I=c('A', 'D', 'F'), J=c('B', 'E', 'G'), K=c('C', 'H'))

# Function call -> extracts injection times and converts the format. saves as list.
kinCommentsTime = processedFile(kinCommentsPath)

# Function call -> extracts the concentration for each injection time. saves as list.
kinCommmentsConc = extractConcentrations(kinCommentsPath, kinCommentsTime)
print(kinCommmentsConc)

# Read in spot file and kinetic data file
spotFile =  read_excel(spotFilePath, col_names = TRUE)
kinFile = read_excel(kinDataPath, col_names = TRUE)

# Create column names list from the spot file
colList = list("Time") # List needs to start with 'Time' as the first column
spotNamesList = as.list(spotFile$Species) # Convert rows in 'Species' column to a list
colList = append(colList, spotNamesList) # Concatenate lists

# Label the columns in the Kinetic df object
colnames(kinFile) <- colList
kinDF <- as.data.frame(kinFile) # Convert to df

dfRowList = list() # List to hold rows in the range of the times in kinCommentsTime + 9 (i:i+9)
i <- 1 # Counter
# Iterate all rows in the dataframe
for (rows in 1:nrow(kinDF)) {
  timeRow <- kinDF[rows, "Time"] # Save the integer in the 'Time' column
  # Iterate the times kinCommentsTime list
  for (num in kinCommentsTime) {
    # Look for integer time values from the rows in the range of time:time+9
    if ( timeRow >= num & timeRow <= round(num+9.00, 2) ) {
      dfRowList[[i]] <- kinDF[rows, ] # Append entire row from df if integer in time column falls inside range
      i <- i + 1 # Increase counter
    }
  }
}

# Use rbind to form the rows in list into a df
spriDF <- do.call(rbind.data.frame, dfRowList)
rownames(spriDF) <- 1:nrow(spriDF) # reset rownames index back to 1
spriDF <- distinct(spriDF)

# Seperate different 9 minute time segments and save these as new df in a list
spriDfList <- list() # Empty list to hold df associated with each injection time
for (num in kinCommentsTime) {
  startTime <- which(spriDF$Time == num) # Index of injection start time
  endTime <- which(spriDF$Time == round(num+9.00, 2))  # Index of injection end time. Need two decimal place or weird things start to happen
  injection <- spriDF[startTime:endTime, ] # Extract injection from DF
  key <- paste('Time-', num, collapse = ',') # Create key name for DF
  spriDfList[[key]] <- injection # Save extracted DF to list with associated key name
}

########## Data Munging ##########
# Main data processing loop
# Loops over each dataframe in the spriDFList object (where each df corresponds to am injection time+9)
spriProcessed <- list() # Empty list to hold all the processed data
i <- 0 # Counter
for (df in spriDfList) {
  if (i == length(spriDfList)) {
    break # break the loop once the last df in the list is processed
  }
  # Save the time Column as a seperate vector, but do not remove it
  # Removing it before the average is taken will cause the column names to become unique, which we don't want
  timeCol <- subset(df, select = 'Time') # Extract 'Time' and save as seperate column

  # Transpose values so that all the values in the df are positive
  colMin <- which.min(apply(df, MARGIN = 2, min)) # Return the column index of the largest negative value
  rowMin <- which.min(apply(df, MARGIN = 1, min)) # Return the row index of the largest negative value
  maxMin <- df[rowMin, colMin] # location of largest negative value
  df <- df + abs(maxMin) # Increase all values by the absolute value of the largest negative integer in df

  # Take the Average of row values in like columns
  df <- sapply(split.default(df, names(df)), rowMeans)
  df <- as.data.frame(df)

  # Can now drop the time column as it has been altered by the the data munging process anyhow
  # The Time vecotor saved earlier will be added in once the derivative is taken
  drops <- c("Time")
  df <- df[ , !(names(df) %in% drops)] # Drop the 'Time' column from main df

  # Subtract the negative from the appropriate columns
  for (colName in names(df)) { # Iterate the columns in the df
    for (neg in 1:length(negList)) { # Iterate the names in the named list. The names are the names of the -ve columns
      if ( colName %in% negList[[neg]] ) { # Check if the column name under the i-th name in the column list
        negName = as.character(names(negList[neg])) # Get name of negative column
        df[[colName]] <- df[[colName]] - df[[negName]] # Subtract the negative column from the positive
      }
    }
  }

  # Remove the negative columns from the df
  drops <- c(names(negList))
  df <- df[ , !(names(df) %in% drops)]

  # Normalize the values in the df to bewteen 0 and 1
  df <- normalize(df) # Use the normalize function

  # Take the derivative of each column.
  # This will also require a new column for each
  originalColCount <- ncol(df) # number of columns
  originalColNames <- list(colnames(df)) # names of all the columns in the df
  j <- 1 # column counter
  for (cols in names(df)) {
    if (j > originalColCount) {
      break # end loop if col limit is reached
    }
    newColName <- paste(originalColNames[[1]][j], "Deriv", sep = "_") # Create column name
    newColData <- list(0) # First number in our new col needs to be zero to make the columns same length of original df
    colData <- as.list(diff(df[[cols]], lag = 1))
    df$V1 <- cbind(unlist(append(newColData, colData))) # Combine the two lists and bind to df as a new column
    names(df)[names(df) == "V1"] <- newColName # Rename the column with a descriptive name

    j <- j + 1 # Increase counter
  }

  # Add Time Column back into the df
  df$Time <- seq(from = 0, to = 540, by = 3)

  i <- i + 1 # Increases counter
  key = "key" # Generic name
  spriProcessed[[key]] <- df # Save processed df as df in list
  names(spriProcessed)[names(spriProcessed) == "key"] <- names(spriDfList)[i] # Rename df with descriptive name
}

# TODO: Work code in from kD_largest_rMax_mulitple_columns
# TODO: in the Kinetic comments, on the injection lines it lists what's being injected as well as the
#       concentration. Could extract this and put into a named list. if this was able to correspond, or could list
#       the time in the name for each entry, could use this to iterate the dfs and generate results
#for (df in spriProcessed) {
  # Still need to figure out how to make this work
#}

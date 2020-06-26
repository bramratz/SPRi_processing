### kD calculations for SPRI data ###

## How the program will calculate kD given normalized values and derivative values
# DONE: Need to read in the kd and ka values from an excel workbook. Needs to sheet specific 
# DONE: Add time column as column one for for re-naming columns from workbook.
# DONE: Find the max value in each norm column. Use larget value index as rmax
# DONE: Calculate the kd value first, then calculate the ka value, save these as variables 
# DONE: Put process of finding kd and ka values into loop. Save the values in a list or df to ID
# DONE: Final kD calculation. 
# DONE: Output all values. Associate the correct values with each aptamer.
# TODO: Fix potential issues created by having columns of different length.
# TODO: robust system for determining dissociation point - old method need to be improved upon.

## Load library dependencies 
install.packages('readxl', dependencies = TRUE)
install.packages('magrittr')

library(readxl)
library(magrittr)
library(dplyr)
library(tidyr)

## Data you can change 
fileName = './SPRI_data_test_file_cot1.xlsx' # Put ./Your_File_Name.xlsx  Uses a relative path to your file
sheetName = 'Test2_125nM_LS25'
proteinConc = 0.000000125
# ***Don't alter beyond this point***

##### Data munging #####
## Read in data amd format
# Adjusted SPRI data
spriData <- read_excel(fileName, sheet = sheetName)

# Convert Data to a df to make it easier to work with 
spriDF <- as.data.frame(spriData)
head(spriDF)

# Rename columns with Norm or Deriv followed by a number
counter = 0 # keep track of number to add to column 
for (colName in 1:ncol(spriDF)) {
  counter = counter+1
  if (counter == 1) {
    names(spriDF)[colName] <- paste('Time', counter,sep="")
  }else if (counter %% 2 != 0) {
    names(spriDF)[colName] <- paste('Deriv',counter,sep="")
  }else if (counter %% 2 == 0) {
    names(spriDF)[colName] <- paste('Norm',counter,sep="")
  }
}

# Remove First row from Norm column, this will also get rid of the NA in the Deriv column
spriDF <- spriDF[-c(1),]
row.names(spriDF) <- NULL # Reset index

# Find columns with "Deriv" in column header
drops <- list() # Empty list
for (cols in names(spriDF)) { # Iterate column names
  if (grepl("Deriv", cols, fix = TRUE)) { # regular expression matching for "Deriv" in column header 
    drops <- c(drops, cols) # Append to list of columns that contain derivative values
  }else if (grepl("Time", cols, fix = TRUE)) { # Regular expression matching for "Time" in column header 
    drops <- c(drops, cols) # Append to list of columns that contain time values
  }
}

# Subset df to drop the Derivative columns 
drops <- unlist(drops)
justNormDf <- spriDF[ , !(names(spriDF) %in% drops)] # Retain only the Norm columns
justNormDf <- data.matrix(justNormDf) # Convert all values to numeric values 


# Find largest row value and its index
rMax <- 0
idx <- 0
for (rows in 1:nrow(justNormDf)) { # Iterate rows
  for (cols in 1:ncol(justNormDf)) { # Iterate columns
    if (is.na(justNormDf[rows, cols])) { # Ignore NA values in rows
      next
    }else if (justNormDf[rows, cols] > rMax) { # Only change rMax if a larger value is found
      rMax = justNormDf[rows, cols] # Update rMax
      idx = rows # Save rMax index 
    }
  }
}
lastIdx = nrow(spriDF) # Index of last value in column 
# ***Assumes that all the columns will be the same legth this is not the case... 
# ***Likely will either need to use the lowest colum last column idx or use something like 100 seconds (33 values)
#    after Rmax will be included. 

##### Calculate kd and ka values for all column pairs #####
## Calculate kd
# Extract values for kd from all columns
kdDF <- spriDF[(idx+1):lastIdx, ] # Values are from idx+1 to the end of the original df

kd = 0.001 # Default starting kd

# Need an empty list to save kd values 
# Could use a list of aptamers as well. would make it easier -> user input 
# for each column in kdDF need to do the following
  # Either save norm and deriv as variables or take straight from df with [[]]
  # Calculate fit with these values
  # save kd to list 
  # move on to the next pair of columns 

halfLen = round(ncol(kdDF)/2) # Number of Pair of columns. NOTE. This will be an odd number but will be rounded down automatically to the correct number 
startPt = 2 # 1st column in Norm/Deriv col pair
endPt = 3 # 2nd column in Norm/Deriv col pair  
kdCal = list() # List to hold values calculated
# Go through column pairs one pair at a time 
for (i in range(1:halfLen)){
  singleAptVal <- kdDF %>%
    select(all_of(startPt):all_of(endPt))
  
  kdDeriv <- (as.vector(as.numeric(singleAptVal[,2])))
  kdNorm <- (as.vector(as.numeric(singleAptVal[,1])))
  
  # Calculate fit
  kdFit = nls(kdDeriv~-kd*kdNorm, start = list(kd=kd))
  kdFitSum = summary(kdFit)
  kdCalc = kdFitSum$parameters[1] # Save calculated kd value
  kdCal = c(kdCal, kdCalc) # Add kd Value to the list 
  
  # Increase start and end point values to encompas the next pair
  startPt = startPt+2
  endPt = endPt+2
}

## Calculate ka
# Extract all values from ka from columns 
kaDF <- spriDF[1:idx, ]

ka = 100000 # Default starting ka
conc = proteinConc # Protein concentration

# Need an empty list to save ka values 
# Could use a list of aptamers as well. would make it easier -> user input 
# for each column in kdDF need to do the following:
#   Either save norm and deriv as variables or take straight from df with [[]]
#   Calculate fit with these values
#   save ka to list 
#   move on to the next pair of columns 

halfLen = round(ncol(kaDF)/2) # Number of Pair of columns. Note if this not an even number this will not work as intended
startPt = 2 # 1st column in Norm/Deriv col pair
endPt = 3 # 2nd column in Norm/Deriv col pair  
kaCal = list() # List to hold values calculated
# Go through column pairs one pair at a time 
for (i in range(1:halfLen)){
  singleAptVal <- kaDF %>%
    select(all_of(startPt):all_of(endPt)) # Grab first column pair 
  
  kaDeriv <- (as.vector(as.numeric(singleAptVal[,2])))
  kaNorm <- (as.vector(as.numeric(singleAptVal[,1])))
  
  # Calculate fit
  kaFit = nls(kaDeriv~ka*conc*rMax-(ka*conc+kdCalc)*kaNorm, start = list(ka=ka))
  kaFitSum = summary(kaFit)
  kaCalc = kaFitSum$parameters[1] # Save calculated kd value
  kaCal = c(kaCal, kaCalc) # Add kd Value to the list 
  
  # Increase start and end point values to encompas the next pair
  startPt = startPt+2
  endPt = endPt+2
}

##### Calculate kD #####
# Viola
for (i in range(1:halfLen)){
  kD = as.numeric(kdCal[i])/as.numeric(kaCal[i])
  cat("\nThis is your kd value: ", as.numeric(kdCal[i])) # Print out kd value to console
  cat("\nThis is your ka value: ", as.numeric(kaCal[i]))  # Print out ka value to console
  cat("\nThis is your kD value: ", kD) # Print out kD value to console
}


##### Interesting graphs #####

library(ggplot2)

## Graphing Norm vs. time 
x <- spriDF$Time1
y1 <- spriDF$Norm2
y2 <- spriDF$Norm4
dfGraph <- data.frame(x, y1, y2)

ggplot(dfGraph, aes(x)) +
  geom_point(aes(y=y1)) +
  geom_point(aes(y=y2)) +
  geom_smooth(aes(y=y1), colour = 'red') +
  geom_smooth(aes(y=y2), colour = 'green')

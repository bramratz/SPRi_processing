### BASIC KD CALCULATOR ### 

"
Calculates kd, ka, and kD for a SINGLE aptamer. 
Calculates the previously listed values using both the Rmax and 240 time point as the dissociation point 
"

## Load library dependencies 
install.packages('readxl', dependencies = TRUE)
install.packages('magrittr')

library(readxl)
library(magrittr)

## Data you can change 
fileName = './Pick_Me/Everything_processed.xlsx' # Put ./Your_File_Name.xlsx  Uses a relative path to your file
sheetName = 'Mark1'
proteinConc = 250e-9
# ***Don't alter beyond this point***

#####
### Data prep
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

# Find the largest value in normalized data column
rMax = max(spriDF$Norm2) # Largest value
idx = which.max(spriDF$Norm2) # Index of largest value
idx240 = which(spriDF$Time1==237) # Index of 240 time point 
lastIdx = length(spriDF$Norm2) # Index of last value in column 

##### 
### Calculating fit using Rmax
## Calculate kd
# Extract normalized values for kd. ie. those values past the index position of the max value
kdNorm <- spriDF$Norm2 %>%
  as.vector() %>%
  as.numeric()
kdNorm <- kdNorm[(idx+1):lastIdx] # Extract those values from position idx+1 to end 

# Extract derivative values for kd ie. those values past the index position of the max value
kdDeriv <-  spriDF$Deriv3 %>%
  as.vector() %>%
  as.numeric()
kdDeriv <- kdDeriv[(idx+1):lastIdx] %>% # Extract those values from position idx+1 to end 
  na.omit()

kd = 0.001 # Default starting kd

# Calculate fit
kdFit = nls(kdDeriv~-kd*kdNorm, start = list(kd=kd))
kdFitSum = summary(kdFit)
kdCalc = kdFitSum$parameters[1] # Save calculated kd value 

# Plot the original points
# first argument is the x values, second is the y values
plot(kdNorm, kdDeriv)

#This adds to the already created plot a line
# once again, first argument is x values, second is y values
lines(kdNorm, predict(kdFit))
cor(kdDeriv,predict(kdFit))

## Calculate ka
# Extract normalized values for ka. ie. those values before the index position of the max value
kaNorm <- spriDF$Norm2 %>%
  as.vector() %>%
  as.numeric()
kaNorm <- kaNorm[1:idx] # Extract those values from position 1 to idx 

# Extract derivative values for ka ie. those values before the index position of the max value
kaDeriv <-  spriDF$Deriv3 %>%
  as.vector() %>%
  as.numeric()
kaDeriv <- kaDeriv[1:idx] %>% # Extract those values from position 1 to idx 
  na.omit()

ka = 100000 # Default starting ka
conc = proteinConc # Protein concentration

# Calculate fit
kaFit = nls(kaDeriv~ka*conc*rMax-(ka*conc+kdCalc)*kaNorm, start = list(ka=ka))
kaFitSum = summary(kaFit)
kaCalc = kaFitSum$parameters[1] # Save calculated ka value 

# Plot the original points
# first argument is the x values, second is the y values
plot(kaNorm, kaDeriv)

#This adds to the already created plot a line
# once again, first argument is x values, second is y values
lines(kaDeriv, predict(kaFit))
nlcor(kaDeriv,predict(kaFit))


# Residual sum-of-squares values for Rmax models 
#kdRSS = kdFit$m$deviance()
#kaRSS = kaFit$m$deviance()

#####
### Calculating fit with 240
## Calculate kd
# Extract normalized values for kd. ie. those values past the index position of the max value
kdNorm240 <- spriDF$Norm2 %>%
  as.vector() %>%
  as.numeric()
kdNorm240 <- kdNorm240[(idx240+1):lastIdx] # Extract those values from position idx+1 to end 

# Extract derivative values for kd ie. those values past the index position of the max value
kdDeriv240 <-  spriDF$Deriv3 %>%
  as.vector() %>%
  as.numeric()
kdDeriv240 <- kdDeriv240[(idx240+1):lastIdx] %>% # Extract those values from position idx+1 to end 
  na.omit()

kd = 0.001 # Default starting kd

# Calculate fit
kdFit240 = nls(kdDeriv240~-kd*kdNorm240, start = list(kd=kd))
kdFitSum240 = summary(kdFit240)
kdCalc240 = kdFitSum240$parameters[1] # Save calculated kd value 

## Calculate ka
# Extract normalized values for ka. ie. those values before the index position of the max value
kaNorm240 <- spriDF$Norm2 %>%
  as.vector() %>%
  as.numeric()
kaNorm240 <- kaNorm240[1:idx240] # Extract those values from position 1 to idx 

# Extract derivative values for ka ie. those values before the index position of the max value
kaDeriv240 <-  spriDF$Deriv3 %>%
  as.vector() %>%
  as.numeric()
kaDeriv240 <- kaDeriv[1:idx240] %>% # Extract those values from position 1 to idx 
  na.omit()

ka = 100000 # Default starting ka
conc = proteinConc # Protein concentration

# Calculate fit
kaFit240 = nls(kaDeriv240~ka*conc*rMax-(ka*conc+kdCalc240)*kaNorm240, start = list(ka=ka))
kaFitSum240 = summary(kaFit240)
kaCalc240 = kaFitSum240$parameters[1] # Save calculated ka value 

# Residual sum-of-squares values for 240 time point models 
kdRSS240 = kdFit240$m$deviance()
kaRSS240 = kaFit240$m$deviance()

#####
### Viola
# For Rmax
kD = kdCalc/kaCalc

cat("This is your kd value using Rmax as the dissociation point: ", kdCalc) # Print out kd value to console
cat("This is your ka value using Rmax as the dissociation point: ", kaCalc)  # Print out ka value to console
cat("This is your kD value using Rmax as the dissociation point: ", kD) # Print out kD value to console

# For 240
kD = kdCalc240/kaCalc240

cat("This is your kd value using 240 seconds as the dissociation point: ", kdCalc240) # Print out kd value to console
cat("This is yourka value using 240 seconds as the dissociation point: ", kaCalc)  # Print out ka value to console
cat("This is your kD value using 240 seconds as the dissociation point: ", kD) # Print out kD value to console


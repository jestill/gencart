#########################################################################
#                                                                       #
# LOCAL MORAN'S I PERMUTATION TEST                                      #
#                                                                       #
# Author: Jamie Estill                                                  #
# Contact: jestill@uga.edu                                              #
# Started: 2/8/2004                                                     #
# Last Visisted: 2/9/2004                                               #
#                                                                       #
# A permutation test of local Moran I values written in the R           #
# statistical programming language. Takes as its input a tab delimited  #
# text file containing 3 columns with a header row. The data are in the #
# format of x,y,z where x and y are spatial coordinates of points in    #
# euclidean space while z is the value associated with each point. Note #
# that the header files must be labeled x,y, and z for the calls to the #
# array to work properly.                                               #
#                                                                       #
# Yields the following results:                                         #
#  *OResults - A vector of the local I values from the Original Data    #
#              matrix.                                                  #
#                                                                       #
#########################################################################

##########################################################
# CLEAR RUBBISH AND PRINT UPDATE                         #
##########################################################

rm(list=ls(all=TRUE))                                                    # Clear all of the current data objects
flush.console()                                                          # The flush console command allows the print data to show up right away
print ("The Local Moran I Program Has Started", quote=FALSE)

##########################################################
# LOAD NEEDED LIBRARIES                                  #
##########################################################
library(spdep)                                                           # The spatial dependency library

##########################################################
# VARIALBES                                              #
##########################################################
InputFile = "C:/RData/RiceChromOne.txt"                                  # Input file of spatial data
OutputFileSummary = "C:/RData/RiceChromOneAnalysisSummary.txt"           # Export file name for the summary of the analysis
OutputFileValues = "C:/RData/RiceChromOneAnalysisIValues.txt"            # Export file name for the entire set of permuted I values
LagMin = 0                                                               # The minimum of the lag distance
LagMax = 200000                                                          # The maximum of the lag distance
NumPerm = 99                                                              # The number of permutations to do
options(verbose=FALSE)                                                   # Set verbose to true to see if this lets me see what is going on
options(echo=TRUE)                                                       # Lets see if the echo command gives me what I want

##########################################################
# ANALYZE THE ORIGINAL DATA MATRIX                       #
##########################################################

SpatialData <- read.table(InputFile, header=T)                           # Read the data into the SpatialData vector
xycoords <- cbind(SpatialData$x, SpatialData$y)                          # Format the xy coordinates as a matrix to get the neighborhood
mynb=dnearneigh(xycoords, LagMin ,LagMax)                                # Set the neighborhood matrix
mylistw = nb2listw( mynb, glist=NULL, style="B", zero.policy = FALSE)    # Determine the neighborhood
OrigLocalIResult = localmoran(SpatialData$z, mylistw, zero.policy=FALSE, spChk=NULL)  # The values from the localmoran of the original matrix are in the OResult vector


##########################################################
# PERMUTATIONS FOR THE ENTIRE DATASET                    #
##########################################################

numRows = length(SpatialData$z)                                          # The number of rows is equal to the size of the array
numCols = NumPerm                                                        # The number of columns is equal to the number of permutations
LocalIOutput <- array(0, c(numRows,numCols))                             # Make an array filled with zeros with a number of rows as cols as defined
LocalISummary <- array(0, c(numRows,5))                                        # This will hold the max and min values


#for (i in 1:(length(SpatialData$z)))                                    # i will be used to represent each of the points in the dataset
# Replace to the above when I've worked out the conditional permutation
for (i in 1:(length(SpatialData$z)))                                     # i will be used to represent each of the points in the dataset

{
print(c("Element ", i, "of", length(SpatialData$z)), quote=FALSE)                                       # Tell me which element is currently running
    for (j in 1:NumPerm)
    {
        flush.console()                                                      # This allows me to see what is going on
        print( c("Element ", i, "Permutation ", j), quote=FALSE)             # Print the current status
        SpatialDataPerm = SpatialData                                                               # Load the spatial data to a new vector 'SptialDataPer' that we will use for a permutation
 

        #########################################################################
        # Resample the data and then load back into the SpatialDataPerm$z       #
        #########################################################################
        TempDataVector <- numeric()

        if (i > 1)                                                               # If this is not the first element in the array
        {                                                                        # then                                                                                                                                 flush.console()

          for (k in 1:i-1)                                                       # for values up until the current element load the values in the same array position
          {
              TempDataVector[k] = SpatialDataPerm$z[k]
          }
        }

        if ( i < length(SpatialData$z))                                          # If this is not the last element in the array
        {                                                                        # then
          for (k in i:(length(SpatialData$z))-1)                                 # for values after the element of interest
          {
              TempDataVector[k] =  SpatialDataPerm$z[k+1]
          }
        }


        SampTempDataVector = sample(TempDataVector)                              # Permute the TempDataVector

        if (i > 1)                                                               # If this is not the first element in the array
        {                                                                        # then
          for (k in 1:i-1)                                                       # for values up until the current element load the values in the same array position
          {
              SpatialDataPerm$z[k] = SampTempDataVector[k]
          }
        }

        if ( i < length(SpatialData$z))                                          # If this is not the last element in the array
        {                                                                        # then
          kMin = i + 1                                                           # load the values into the array
          kMax = length(SpatialData$z)
          for (k in kMin:kMax)
          {
              Offset = k-1
              SpatialDataPerm$z[k] =  SampTempDataVector[Offset]
          }
        }


        PermLocalIResult = localmoran(SpatialDataPerm$z, mylistw, zero.policy=FALSE, spChk=NULL)    # The values from the localmoran of the original matrix are in the OResult vector
        # This is the new thing that was added to get the single value that I need and load it into the 2d array
        TheIList = PermLocalIResult$Ii
        LocalIOutput[i,j] = TheIList[i]
        #LocalIOutput[i,j] = 99
    } # end of the for next loop for j


} # end of the for next loop for i which is each of the points in the dataset

##########################################################
# GET THE SUMMARY DATA                                   #
##########################################################
flush.console()
print("Calculating summary statistics")                                       # Print current status
TempVector <- numeric()
for (i in 1:(length(SpatialData$z)))
{

 for ( j in 1:NumPerm)

 {
 TempVector[j] = LocalIOutput[i,j]
 }

LocalISummary[i,1] = max(TempVector)                                          # The min value of the permuted Is
LocalISummary[i,2] = min(TempVector)                                          # The max value of the permuted Is
TheOrigIList = OrigLocalIResult$Ii                                            # Fetch the set of Is from the original list
LocalISummary[i,3] = TheOrigIList[i]                                          # The observed value of the Local I
LocalISummary[i,4] = mean(TempVector)                                         # The mean of the permuted Is
LocalISummary[i,4] = sd(TempVector)                                           # The standard deviation of the permuted Is
}

##########################################################
# EXPORT THE RELEVANT DATA TO A TEXT FILE                #
##########################################################
flush.console()
print("Exporting the summary statistics")                                     # Print Current Status
write.table(LocalISummary, OutputFileSummary)

flush.console()
print("Exporting the permuted I values")                                      # Print Current Status
write.table(LocalIOutput, OutputFileValues)

# LET ME KNOW THAT THE ANALYSIS IS COMPLETE
print ("The analysis is complete", quote=FALSE)

#########################################################################
#                                                                       #
# RANDOM REAPER aka SNPPicker                                           #
#                                                                       #
# Author: Jamie Estill                                                  #
# Contact: jestill@uga.edu                                              #
# Started: 2/11/2004                                                    #
# Last Visisted: 2/12/2004                                              #
#                                                                       #
#                                                                       #
#                                                                       #
#########################################################################

##########################################################
# CLEAR RUBBISH AND PRINT UPDATE                         #
##########################################################

rm(list=ls(all=TRUE))                                                    # Clear all of the current data objects
flush.console()                                                          # The flush console command allows the print data to show up right away
StartTime = Sys.time()
print ("The Simulation Has Started", quote=FALSE)


##########################################################
# LOAD NEEDED LIBRARIES                                  #
##########################################################
library(spdep)                                                           # The spatial dependency library

##########################################################
# VARIALBES                                              #
##########################################################
InputFile = "C:/RData/SNPChrom12.txt"                                        # Input file of spatial data
OutputFileSummary = "C:/RData/SNPChrom12_Summary.txt"                         # Export file name for the summary of the analysis
OutputFileValues = "C:/RData/SNPChrom12_Raw.txt"                              # Export file name for the entire set of permuted I values
jfind <-function(x,y) seq(along=x) [x > y]                               # Define a function that will will let me find the one I need

numPositive = 23810                                                      # The number of 'Positive' hits in the dataset
NumPerm = 3                                                              # The number of permutations to do

options(verbose=FALSE)                                                   # Set verbose to true to see if this lets me see what is going on
options(echo=TRUE)                                                       # Lets see if the echo command gives me what I want

ElementSizes <- read.table(InputFile, header=F)
numCols = NumPerm                                                        # The number of columns is equal to the number of permutations
numRows = (length(ElementSizes[[1]]))                                    # The number of rows is equal to the length of the array
Output <- array(0, c(numRows,numCols))
ElementSummary <- array(0, c(numRows,4))                                  # An array to hold the summary data

#print (numRows)

for (j in 1:NumPerm)
{
    flush.console()                                                          # The flush console command allows the print data to show up right away
    print(c("Permutation",j,"of",NumPerm))

    MyData = (ElementSizes[[1]])                                                # Load the original ElementSizes dataset into the data vector
    #print (MyData)
    for (i in 1:numPositive)
    {
    flush.console()
    cumdata = cumsum(MyData)
    RandomNumber =  runif(1, min=0, max=1)
    pick = RandomNumber * sum(MyData)

     ##########################################################
     # GETTING THE PICKED ROW BY FINDING THE VALUE BEFORE THE ROW THAT WILL BE PICKED
     ##########################################################
     #PrePick = 0
     #for (k in 1:numRows)
     #{
     #    if (cumdata[[k]]<pick) PrePick = k                         # This will stop when the current row is before the selected row, therefore the selected row is k + 1
     #}
     #PickedRow = PrePick + 1
     #print(c("Permutation: ",j," Element ",i, " of ", numPositive ), quote=FALSE)
     ##########################################################
     # END OF GET THE VALUE BEFORE THE ROW THAT WILL BE PICKED
     ##########################################################


     ##########################################################
     # A DIFFERENT APPROACH TO GET THE PICKED ROW
     ##########################################################
     Above <- jfind (cumdata,pick)                             # Find the index number of all values where the cumdata is above the pick number
     PickedRow = Above[1]                                      # Select the first element of the list where cumdata>pick, this will be the first value that meets the criteria
     ##########################################################
     # END A DIFFERENT APPROACH TO GET THE PICKED ROW
     ##########################################################

     ##########################################################
     # AND NOW FOR SOMETHING COMPLETELY DIFFERENT
     ##########################################################
     #Above <- which (cumdata>pick)                             # Find the index number of all values where the cumdata is above the pick number
     #PickedRow = Above[1]                                      # Select the first element of the list where cumdata>pick, this will be the first value that meets the criteria
     ##########################################################
     # AND NOW FOR SOMETHING COMPLETELY DIFFERENT
     ##########################################################


     MyData[PickedRow] = MyData[PickedRow] - 1                 # Subtract one from the elements that is being picked fom the MyData Vecotr and then
     Output[PickedRow,j] = Output[PickedRow,j] + 1 ;           # add one to the element in the Output vector that has been picked.

    }
}

##########################################################
# GET THE SUMMARY DATA                                   #
##########################################################
flush.console()
print("Calculating summary statistics")                                       # Print current status
TempVector <- numeric()
#numRows = (length(ElementSizes[[1]]))                                    # The number of rows is equal to the length of the array

for (i in 1:numRows)
{

for ( j in 1:NumPerm)
{
 TempVector[j] = Output[i,j]
}

ElementSummary[i,1] = max(TempVector)                                          # The min value of the permuted Is
ElementSummary[i,2] = min(TempVector)                                          # The max value of the permuted Is
ElementSummary[i,3] = mean(TempVector)                                         # The mean of the permuted Is
ElementSummary[i,4] = sd(TempVector)                                           # The standard deviation of the permuted Is

}


##########################################################
# WRITE THE RESULTS TO A TEXT FILE
##########################################################
write.table(Output, OutputFileValues)
write.table(ElementSummary, OutputFileSummary)

EndTime = Sys.time()
TotalTime = EndTime - StartTime
print("END")
print(TotalTime)


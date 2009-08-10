# This will be the network K analysis for the chromosomes
# This will initially be coded for cumulative lags since I do not know how
#  well this will work for the lag distances that I will use.

print ("The Die Has Been Cast..", quote=FALSE)                    # Let the user know that the permutation has started
flush.console()
rm(list=ls(all=TRUE))                                             # Clear out any prexisting objects

######################################
#                                    #
# IMPORTANT VARIABLES                #
#                                    #
######################################
PermModel = "all genome loci"                                     # The permutation model that is used in the analysis
Taxon = "Rice"                                                    # The name of the taxon that is being plotted, this variable is used in drawing the plots
WorkDir = "C:/RData/RiceSNP/"                                     # The working dir, this is where the input files are and where the output files will go
GenomeFileName = "RiceGenome.txt"                                 # File containing the infromation related
ObsLociName = "RiceSNP.txt"                                       # The input file that contains the data related to the observed file list
RawOutName = "RawData.txt"                                        # The name of the raw data output file
SumOutName = "SumData.txt"                                        # The name of the summary output file
NumPerm = 99                                                      # The number of permuations to do
NumBreaks = 30                                                    # The number of breaks to use in the histogram plots
# LAG DISTANCE CLASSES
LagDist <- 1000000                                                # The base lag distance to build upong
NumLags = 3                                                       # The number of lags that will need to be calculated
LagBins <- 1:NumLags                                              # Gives and array of numbers for lag numbers
LagCalc <- LagBins*LagDist                                        # Gives the numeric quantaties to use for the lag distances

# LOAD THE FUNCTIONS THAT I NEED
# MY LAPTOP PATHS
source("C:/RData/PFAM/NetK.r")                                    # The function to do the network K calculation.
# source("G:/JamieFiles/iFolder/jestill/Home/RWork/GetArgs.r")    # The function to get arguments from the R command line.

# LOAD INPUT DATA AND DETERMINE SOME SUMMARY DATA
GenomeFile = paste(WorkDir,GenomeFileName, sep="")                # Information relate to all genes/loci for the entire genome of the test organism.
Genome <- read.table(GenomeFile, header=T)                        # Load the genome data into a data table
TransGenome <- t(Genome)
NumGenomeLoci <- length(Genome$Chromosome)                        # The total number of loci in the genome
NumGenomeChrom <- max(Genome$Chromosome)                          # The total number of chromosomes in the genome file
ObsLociFile = paste(WorkDir,ObsLociName, sep="")                  # The file containing information related to the set of observed loci
ObsLoci <- read.table(ObsLociFile, header=T)                      # Load the observed loci data into a data table
NumberObsLoci <- length(ObsLoci$Chromosome)                       # The number of observed loci
NumberObsChrom <- max(ObsLoci$Chromosome)                         # The number of chromosomes in the observed loci data set

# SET VARIABLES (THESE MAY ALSO BE SPECIFIED FROM THE COMMAND LINE LATTER USING GETARGS)
RawOutputFile = paste(WorkDir,RawOutName, sep="")                 # The raw output file to hold all data
SummaryOutputFile = paste(WorkDir,SumOutName, sep="")             # The summary output file
OutRowNum = 0                                                     # The outputrow number will be used later to keep track of the number of rows
numCols = 3                                                       # The number of columns in the output file
numRows = (NumPerm * NumGenomeChrom) + 1                          # The number of rows that will be needed in the output file

# The original raw output for a single lag data frame is below
# RawOutput <- data.frame("RowNum"=0,"Chromosome"=0,"NumLoci"=0)  # The output as a data frame

# GENERATE THE DATA FRAME TO HOLD THE RAW OUTPUT

# May want to have a summary output data matrix for each chromosome to simplify things a bit


#####################################################
#                                                   #
# DEFINE DATA FRAMES                                #
# Define the summary data.                          #
#                                                   #
#####################################################
# DEFINE THE RAW OUTPUT DATA FRAME
RawOutput <- data.frame("RowNum"=0,"Chromosome"=0)                # First two cols in the data frame are RowNum and Chromosome
for (ThisLag in 1:NumLags)                                        # For each of the lags
{
 ColName = paste("Lag",ThisLag,sep="")                            # This is the name that would need to be used for this column
 RawOutput[1,ThisLag+2] = 0                                       # need to offset this by two since the first two columns will contain RowNum and Chromosome
 names(RawOutput)[ThisLag+2] <- ColName                           # Name the column after the name of the lag that this represents
}

# DEFINE THE SUMMARY OUTPUT DATAFRAM
SumOutput <- data.frame("RowNum"=0,"Chromosome"=0, "LagNum"=0, "LagMax"=0, "LagMin"=0, "LagObs"=0)              # First two cols of the summary output

#####################################################
#                                                   #
# PERMUTATION                                       #
# Permutate the values based on the selected        #
# permutation approach                              #
#                                                   #
#####################################################

for (ThisChrom in 1:NumGenomeChrom)
{
  # Reset the dataframe to hold the temp data
  # In this sampleoutput, the row number is equal to the lag number, therefore
  # I can replace the data for a specific row
  TempSampOut <- NULL



  ObsChrom <- subset(ObsLoci, ObsLoci$Chromosome==ThisChrom)
  NumObs <- length(ObsChrom$Chromosome)
  print((c("CHROMOSOME: ", ThisChrom, "NumObs", NumObs)), quote=FALSE)

  for (Perm in 1:NumPerm)
  {
   OutRowNum = OutRowNum + 1
   print (OutRowNum)
   # Get a sample number of loci per chromosome for the number of observed loci

   # THIS RANDOMIZES ACROSS THE ENTIRE GENOME
   SelectRowNums <- sample (1:NumGenomeLoci, NumberObsLoci)     # Get a random set of rows to select from the genome that are equal to the observed number of loci
   RandomGenome <- Genome[SelectRowNums, ]                      # Select those rows from the origina Genome Data frame
   SampChrom <- subset(RandomGenome, RandomGenome$Chromosome==ThisChrom)
   GenomeChrom <- subset(ObsLoci, ObsLoci$Chromosome==ThisChrom)
   FullGenomeChrom <- subset(Genome, Genome$Chromosome==ThisChrom)
   ChromLength <- max(FullGenomeChrom$MidPoint)

   RawOutput[OutRowNum,1] = OutRowNum
   RawOutput[OutRowNum,2] = ThisChrom

   # CALCULATE THE K VALUE FOR EACH OF THE LAG VALUES OF INTEREST
   for (ThisLag in 1:NumLags)                                        # For each of the lags
   {

    SampKValue <- GenKh (ThisLag*LagDist, ChromLength, SampChrom$MidPoint)
    RawOutput[OutRowNum,ThisLag+2] = SampKValue

   }

   flush.console()

  } # END OF FOR EACH PERMUTATION

} # END OF FOR EACH CHROMOSOME

# I then imported this trans table into access, and then exported it as a tab delim text file

print("The permutation is finished",quote=FALSE)

#####################################################
#                                                   #
# SUMMARY DATA                                      #
# Determine the summary data.                       #
#                                                   #
#####################################################


#####################################################
#####################################################
# OUTPUT                                            #
#####################################################
#####################################################

#####################################################
#                                                   #
# TEXT                                              #
# Raw data output and summary text file.            #
#                                                   #
#####################################################
# Raw output
write.table ( RawOutput , paste(WorkDir,RawOutName) )

# Summary output
write.table ( SumOutput , paste(WorkDir,SumOutName) )

#####################################################
#                                                   #
# PLOT: BOX PLOT OVERVIEW                           #
# One box plot for each lag.                        #
#                                                   #
#####################################################
# CAN CHANGE THIS FOR I IN 3:NUMCOLS TO GET INFORMATION FOR ALL H VALUES THAT ARE SELECTED
for (ThisLag in 1:NumLags)                                        # For each of the lags
{

 # FIRST WILL NEED TO DETERMINE HERE THE MINIMUN AND MAXIMUM Y VALUE TO SHOW
 ThisLagYMin = min(RawOutput[[ThisLag+2]])
 ThisLagYMax = max(RawOutput[[ThisLag+2]])

 for (ThisChrom in 1:NumGenomeChrom)
 {
     ObsK <- GenKh (ThisLag*LagDist, ChromLength, ObsChrom$MidPoint)
     if (ObsK < ThisLagYMin)
     { ThisLagYMin = ObsK }
     if (ObsK > ThisLagYMax)
     { ThisLagYMax = ObsK }
 }
 # DRAW THE PLOT

 MyPlot <-boxplot(RawOutput[[ThisLag+2]] ~ RawOutput$Chromosome, notch=TRUE, col="lightblue2", ylim=c(0,20000000))

 # Original below, new code above
 #MyPlot <-boxplot(RawOutput$NumLoci ~ RawOutput$Chromosome, notch=TRUE, col="lightblue2", ylim=c(0,10000000))
 title(main="Permutation from All Rice Loci", sub = paste("perm=", NumPerm," , ","lag=",ThisLag*LagDist, sep=""), xlab="CHROMOSOME", ylab="NETWORK-K")
 for (ThisChrom in 1:NumGenomeChrom)
 {
   ObsChrom <- subset(ObsLoci, ObsLoci$Chromosome==ThisChrom)
   # The following should just be calculated above and added to a separate data frame that is just looked up here, however this will work for now .. just slow
   # The current use is not so good since it assumes that the last gene position is the full space of the chromosome
   # However this assumption may be valid since the permuation is limited to the actual coding gene space and not the entire genome
   FullGenomeChrom <- subset(Genome, Genome$Chromosome==ThisChrom) # Needed to get the max value for this chromosomes
   ChromLength = ChromLength <- max(FullGenomeChrom$MidPoint)      # Needed to get the max value for this chromosome
   ObsK <- GenKh (ThisLag*LagDist, ChromLength, ObsChrom$MidPoint)

   print(paste("Chromosome",ThisChrom,"Observed K",ObsK))
   points(ThisChrom ,ObsK, pch=23, bg="red")
   # text(ThisChrom, ObsK, ObsK)                                  # The following can be added to show the value
   savePlot(paste(WorkDir,Taxon,"BoxPlot","_",NumPerm,"Perm_","Lag",ThisLag,sep=""), "wmf")
 } # End of for each chromosome
} # End of for each lag


#####################################################
#                                                   #
# PLOT: HISTOGRAMS                                  #
# One histogram for each chromosome for each lag.   #
#                                                   #
#####################################################
# To show the full histogram
for (ThisLag in 1:NumLags)                                        # For each of the lags
{

  for (ThisChrom in 1:NumGenomeChrom)
  {
    ObsChrom <- subset(ObsLoci, ObsLoci$Chromosome==ThisChrom)
    ObsK <- GenKh (ThisLag*LagDist, ChromLength, ObsChrom$MidPoint)
    SampChrom <- subset(RawOutput, RawOutput$Chromosome==ThisChrom)
    MySubSet <- subset(RawOutput, RawOutput$Chromosome==ThisChrom)
    HistPlot <-hist(RawOutput[[ThisLag+2]], breaks=NumBreaks, col="lightblue2", main=c(Taxon, ThisChrom), sub = paste("perm=", NumPerm," , ","lag=",ThisLag*LagDist, sep=""), xlab="NETWORK-K", ylab="FREQUENCY", xlim=c(0,20000000) )

    points(ObsK, 0, pch=23, bg="red", cex=2.5)
    savePlot(paste(WorkDir,Taxon,"_perm",NumPerm,"_hist_lag",ThisLag,"_chrom",ThisChrom, sep=""), "wmf")
  }
}

#####################################################
#                                                   #
# PLOT: CORRELOGRAMS                                #
# One correlogram for each chromosome               #
#                                                   #
#####################################################
for (ThisChrom in 1:NumGenomeChrom)
{
 # Get the sumsample of the summary data matrix that will be needed

 # Determine the xlim of the plot

 # Determine the ylim of the plot

 # Save the output
 savePlot(paste(WorkDir,Taxon,"_perm",NumPerm,"correl_chrom",ThisChrom, sep=""), "wmf")

}


#####################################################
#                                                   #
# WORKSPACE                                         #
# The workspace includes the objects that were      #
# generated as part of running the program          #
#                                                   #
#####################################################
save.image(paste(WorkDir,Taxon,"WorkSpace",".RData",sep=""))






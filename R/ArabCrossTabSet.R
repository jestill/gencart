library(spdep) # I need to load the spdep library if I am going to do Moran's I analysis.

# Generate Summary Data for Genome Annotations

# Variables that will need to be changed for each plot

# Chromosome 5
InputFile = "C:/RData/Arab/ArabStrandr3chr1.txt"
myTax = "Arabidopsis"
GeneData <- read.table(InputFile, header=T)


# Centromere positions, Arab Chrom 5
Ch1CentPt = 15263481
Ch2CentPt = 3635266.5
Ch3CentPt = 13875255
Ch4CentPt = 2042760
Ch5CentPt = 11814256

# Set the centromere position for the gene of interest
myCentPos = Ch1CentPt

myChrom <- GeneData$Chromosome[1]
myTitle <- paste("Gene Density", myTax, "Chromosome", myChrom, "rev 3")
numMeg <- max(GeneData$MidPoint/1000000)                                # Determine the number of megabases in a sequence
test <- hist(GeneData$MidPoint, breaks=numMeg*2, main=myTitle, xlab="Midpoint Bin (mB)", col="blue")
# The abline would be a good way to add the location of the centromere to the plot
abline(v=myCentPos,col=3,lty=3,lwd=5)

# Take the results from the bins and do a Moran's I, Show the results of the Moran's I on the graph if possible





# Create a vector of ones or zeros if I want to draw these at the bottom of the plot

# The points objects could be used to draw overgos on top of the physical map
points(OvergoData$MidPoint,GeneData$Chromosome, col="black", pch=10)


test <- hist(GeneData$MidPoint, breaks=numMeg, main=c("Gene Density", myTax,"Chromosome ",myChrom), xlab="Midpoint Bin (mB)", col="blue", ylim=numMeg)


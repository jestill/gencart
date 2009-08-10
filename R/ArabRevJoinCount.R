# This will be where I am working on using the spdep library to
# do join count tests of contiguity of +/- from genomic data sets
# For prokaryotic datasets this will need to modified to use the
# mapping onto a torus to avoid edge effects

# For the following to work the data MUST be listed in the proper order (ie. presort the list by midpoint before doing the analysis)
# LOAD THE SPDEP LIBRARY
library(spdep)

InputFile = "C:/RData/Arab/ArabStrandr1chr3.txt"
GeneData <- read.table(InputFile, header=T)
numgenes = (length(GeneData[[1]]))                                # Determine the number of genes for drawing the chromosome nb
GeneData <- GeneData[sort.list(GeneData$MidPoint),]               # Do A sort to get the data in the correct order
chrom.nb <- cell2nb( numgenes, 1, type="rook", torus=FALSE)
Result1 <-  joincount.multi(GeneData$Strand, nb2listw(chrom.nb, style="B"))                 # This may work if the GeneData$Strand is automatically used as a factor data type

InputFile = "C:/RData/Arab/ArabStrandr2chr3.txt"
GeneData <- read.table(InputFile, header=T)
numgenes = (length(GeneData[[1]]))                                # Determine the number of genes for drawing the chromosome nb
GeneData <- GeneData[sort.list(GeneData$MidPoint),]               # Do A sort to get the data in the correct order
chrom.nb <- cell2nb( numgenes, 1, type="rook", torus=FALSE)
Result2 <-  joincount.multi(GeneData$Strand, nb2listw(chrom.nb, style="B"))                 # This may work if the GeneData$Strand is automatically used as a factor data type

InputFile = "C:/RData/Arab/ArabStrandr3chr3.txt"
GeneData <- read.table(InputFile, header=T)
numgenes = (length(GeneData[[1]]))                                # Determine the number of genes for drawing the chromosome nb
GeneData <- GeneData[sort.list(GeneData$MidPoint),]               # Do A sort to get the data in the correct order
chrom.nb <- cell2nb( numgenes, 1, type="rook", torus=FALSE)
Result3 <-  joincount.multi(GeneData$Strand, nb2listw(chrom.nb, style="B"))                 # This may work if the GeneData$Strand is automatically used as a factor data type

InputFile = "C:/RData/Arab/ArabStrandr4chr3.txt"
GeneData <- read.table(InputFile, header=T)
numgenes = (length(GeneData[[1]]))                                # Determine the number of genes for drawing the chromosome nb
GeneData <- GeneData[sort.list(GeneData$MidPoint),]               # Do A sort to get the data in the correct order
chrom.nb <- cell2nb( numgenes, 1, type="rook", torus=FALSE)
Result4 <-  joincount.multi(GeneData$Strand, nb2listw(chrom.nb, style="B"))                 # This may work if the GeneData$Strand is automatically used as a factor data type

InputFile = "C:/RData/Arab/ArabStrandr5chr3.txt"
GeneData <- read.table(InputFile, header=T)
numgenes = (length(GeneData[[1]]))                                # Determine the number of genes for drawing the chromosome nb
GeneData <- GeneData[sort.list(GeneData$MidPoint),]               # Do A sort to get the data in the correct order
chrom.nb <- cell2nb( numgenes, 1, type="rook", torus=FALSE)
Result5 <-  joincount.multi(GeneData$Strand, nb2listw(chrom.nb, style="B"))                 # This may work if the GeneData$Strand is automatically used as a factor data type



# This will be where I am working on using the spdep library to
# do join count tests of contiguity of +/- from genomic data sets
# For prokaryotic datasets this will need to modified to use the
# mapping onto a torus to avoid edge effects

# For the following to work the data MUST be listed in the proper order (ie. presort the list by midpoint before doing the analysis)
# LOAD THE SPDEP LIBRARY
library(spdep)

InputFile = "C:/RData/Rice/Ricev3chr01.txt"
GeneData <- read.table(InputFile, header=T)
numgenes = (length(GeneData[[1]]))                                # Determine the number of genes for drawing the chromosome nb
GeneData <- GeneData[sort.list(GeneData$MidPoint),]               # Do A sort to get the data in the correct order
chrom.nb <- cell2nb( numgenes, 1, type="rook", torus=FALSE)
Result01 <-  joincount.multi(GeneData$Strand, nb2listw(chrom.nb, style="B"))                 # This may work if the GeneData$Strand is automatically used as a factor data type

InputFile = "C:/RData/Rice/Ricev3chr02.txt"
GeneData <- read.table(InputFile, header=T)
numgenes = (length(GeneData[[1]]))                                # Determine the number of genes for drawing the chromosome nb
GeneData <- GeneData[sort.list(GeneData$MidPoint),]               # Do A sort to get the data in the correct order
chrom.nb <- cell2nb( numgenes, 1, type="rook", torus=FALSE)
Result02 <-  joincount.multi(GeneData$Strand, nb2listw(chrom.nb, style="B"))                 # This may work if the GeneData$Strand is automatically used as a factor data type

InputFile = "C:/RData/Rice/Ricev3chr03.txt"
GeneData <- read.table(InputFile, header=T)
numgenes = (length(GeneData[[1]]))                                # Determine the number of genes for drawing the chromosome nb
GeneData <- GeneData[sort.list(GeneData$MidPoint),]               # Do A sort to get the data in the correct order
chrom.nb <- cell2nb( numgenes, 1, type="rook", torus=FALSE)
Result03 <-  joincount.multi(GeneData$Strand, nb2listw(chrom.nb, style="B"))                 # This may work if the GeneData$Strand is automatically used as a factor data type


InputFile = "C:/RData/Rice/Ricev3chr04.txt"
GeneData <- read.table(InputFile, header=T)
numgenes = (length(GeneData[[1]]))                                # Determine the number of genes for drawing the chromosome nb
GeneData <- GeneData[sort.list(GeneData$MidPoint),]               # Do A sort to get the data in the correct order
chrom.nb <- cell2nb( numgenes, 1, type="rook", torus=FALSE)
Result04 <-  joincount.multi(GeneData$Strand, nb2listw(chrom.nb, style="B"))                 # This may work if the GeneData$Strand is automatically used as a factor data type


InputFile = "C:/RData/Rice/Ricev3chr05.txt"
GeneData <- read.table(InputFile, header=T)
numgenes = (length(GeneData[[1]]))                                # Determine the number of genes for drawing the chromosome nb
GeneData <- GeneData[sort.list(GeneData$MidPoint),]               # Do A sort to get the data in the correct order
chrom.nb <- cell2nb( numgenes, 1, type="rook", torus=FALSE)
Result05 <-  joincount.multi(GeneData$Strand, nb2listw(chrom.nb, style="B"))                 # This may work if the GeneData$Strand is automatically used as a factor data type


InputFile = "C:/RData/Rice/Ricev3chr06.txt"
GeneData <- read.table(InputFile, header=T)
numgenes = (length(GeneData[[1]]))                                # Determine the number of genes for drawing the chromosome nb
GeneData <- GeneData[sort.list(GeneData$MidPoint),]               # Do A sort to get the data in the correct order
chrom.nb <- cell2nb( numgenes, 1, type="rook", torus=FALSE)
Result06 <-  joincount.multi(GeneData$Strand, nb2listw(chrom.nb, style="B"))                 # This may work if the GeneData$Strand is automatically used as a factor data type


InputFile = "C:/RData/Rice/Ricev3chr07.txt"
GeneData <- read.table(InputFile, header=T)
numgenes = (length(GeneData[[1]]))                                # Determine the number of genes for drawing the chromosome nb
GeneData <- GeneData[sort.list(GeneData$MidPoint),]               # Do A sort to get the data in the correct order
chrom.nb <- cell2nb( numgenes, 1, type="rook", torus=FALSE)
Result07 <-  joincount.multi(GeneData$Strand, nb2listw(chrom.nb, style="B"))                 # This may work if the GeneData$Strand is automatically used as a factor data type


InputFile = "C:/RData/Rice/Ricev3chr08.txt"
GeneData <- read.table(InputFile, header=T)
numgenes = (length(GeneData[[1]]))                                # Determine the number of genes for drawing the chromosome nb
GeneData <- GeneData[sort.list(GeneData$MidPoint),]               # Do A sort to get the data in the correct order
chrom.nb <- cell2nb( numgenes, 1, type="rook", torus=FALSE)
Result08 <-  joincount.multi(GeneData$Strand, nb2listw(chrom.nb, style="B"))                 # This may work if the GeneData$Strand is automatically used as a factor data type


InputFile = "C:/RData/Rice/Ricev3chr09.txt"
GeneData <- read.table(InputFile, header=T)
numgenes = (length(GeneData[[1]]))                                # Determine the number of genes for drawing the chromosome nb
GeneData <- GeneData[sort.list(GeneData$MidPoint),]               # Do A sort to get the data in the correct order
chrom.nb <- cell2nb( numgenes, 1, type="rook", torus=FALSE)
Result09 <-  joincount.multi(GeneData$Strand, nb2listw(chrom.nb, style="B"))                 # This may work if the GeneData$Strand is automatically used as a factor data type


InputFile = "C:/RData/Rice/Ricev3chr10.txt"
GeneData <- read.table(InputFile, header=T)
numgenes = (length(GeneData[[1]]))                                # Determine the number of genes for drawing the chromosome nb
GeneData <- GeneData[sort.list(GeneData$MidPoint),]               # Do A sort to get the data in the correct order
chrom.nb <- cell2nb( numgenes, 1, type="rook", torus=FALSE)
Result10 <-  joincount.multi(GeneData$Strand, nb2listw(chrom.nb, style="B"))                 # This may work if the GeneData$Strand is automatically used as a factor data type


InputFile = "C:/RData/Rice/Ricev3chr11.txt"
GeneData <- read.table(InputFile, header=T)
numgenes = (length(GeneData[[1]]))                                # Determine the number of genes for drawing the chromosome nb
GeneData <- GeneData[sort.list(GeneData$MidPoint),]               # Do A sort to get the data in the correct order
chrom.nb <- cell2nb( numgenes, 1, type="rook", torus=FALSE)
Result11 <-  joincount.multi(GeneData$Strand, nb2listw(chrom.nb, style="B"))                 # This may work if the GeneData$Strand is automatically used as a factor data type


InputFile = "C:/RData/Rice/Ricev3chr12.txt"
GeneData <- read.table(InputFile, header=T)
numgenes = (length(GeneData[[1]]))                                # Determine the number of genes for drawing the chromosome nb
GeneData <- GeneData[sort.list(GeneData$MidPoint),]               # Do A sort to get the data in the correct order
chrom.nb <- cell2nb( numgenes, 1, type="rook", torus=FALSE)
Result12<-  joincount.multi(GeneData$Strand, nb2listw(chrom.nb, style="B"))                 # This may work if the GeneData$Strand is automatically used as a factor data type


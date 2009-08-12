#-----------------------------------------------------------+
# join_count_test.r - Test join count statistics            |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill at gmail.com                         |
# STARTED: 07/2005                                          |
# UPDATED: 08/11/2009                                       |
# DESCRIPTION:                                              |
# Test using the spdep library to do join count tests of    |
# contiguity of the +/- from genomic data sets.             |
# For prokaryotic datasets, this could be modified to use   |
# the mapping onto a torus to avoid edge effects.           |
#                                                           |
#-----------------------------------------------------------+
# NOTE: For the following to work the data MUST be listed in the proper order
# (ie. presort the list by midpoint before doing the analysis)

# LOAD THE SPDEP LIBRARY
library(spdep);

#-----------------------------+
# ROOKS CASE JOIN COUNT       |
# INPUT TEXT FILE             |
#-----------------------------+
# Input is Arabidopsis chromosome 5
# This is OLD data from TIGR version 5 of the arabidopsis 
InputFile <- "../data/ArabStrandr5chr5.txt";
# Import data to the GeneData object
GeneData <- read.table(InputFile, header=T);
# Determine the number of genes for drawing the chromosome nb
numgenes = (length(GeneData[[1]]));
# Do A sort to get the data in the correct order
# Treat each locus as the midpoint
GeneData <- GeneData[sort.list(GeneData$MidPoint),];
# Generate the neighbors files
chrom.nb <- cell2nb( numgenes, 1, type="rook", torus=FALSE);
# This may work if the GeneData$Strand is automatically used as a factor data type
Result5 <- joincount.multi(GeneData$Strand, nb2listw(chrom.nb, style="B"));

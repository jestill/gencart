# To get the information for the chromosome of interest



SelectRowNums <- sample (1:NumGenomeLoci, NumberObsLoci)     # Get a random set of rows to select from the genome that are equal to the observed number of loci
RandomGenome <- Genome[SelectRowNums, ]                      # Select those rows from the origina Genome Data frame
SampChrom <- subset(RandomGenome, RandomGenome$Chromosome==ThisChrom)

GenomeChrom <- subset(ObsLoci, ObsLoci$Chromosome==ThisChrom)

FullGenomeChrom <- subset(Genome, Genome$Chromosome==ThisChrom)
ChromLength <- max(FullGenomeChrom$MidPoint)


# I WILL NEED SampChrom$MidPoint to return as column

# From the initial query
GenKh (1000000, ChromLength,SampChrom$MidPoint)
# -- result was
GenKh (1000000, ChromLength,SampChrom$MidPoint)

# USING A TEST FUNCTION THAT SHOULD RETURN THE LENGTH OF THE VECTOR

test01 <-as.vector(SampChrom$MidPoint)                    # Still kept these as rows
result01 <- WorkKh (1000000, ChromLength, test01)
# -- result was
[1] 22


test02 <- matrix(SampChrom$MidPoint, 22, 1)                           # Did something very strange
result02 <- WorkKh (1000000, ChromLength, test02)
# -- result was
[1] 22

test03 <- matrix(SampChrom$MidPoint, 22, 1)                # Looks like this may work
result03 <- WorkKh (1000000, ChromLength, test03)
# -- result was
[1] 22

result04 <- WorkKh (1000000, ChromLength, as.vector(SampChrom$MidPoint))
# -- result was
[1] 22

result05 <- GenKh (1000000, ChromLength, as.vector(SampChrom$MidPoint))
# --result was
Error in Features[i, 1] : incorrect number of dimensions

result06 <- GenKh (1000000, ChromLength, test02)
# -- result was
Error in if (EvalDistance <= CritDistance) { : 
argument is of length zero



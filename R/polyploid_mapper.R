#########################################################################
#                                                                       #
# GENETIC MAPPING IN AUTOPOLYPLOIDS                                     #
# USING MAXIMUM LIKLIHOOD                                               #
# UNDER ASYMMETRIC COUPLING                                             #
# Author: Jamie Estill                                                  #
# Contact: jestill_atthe_uga.edu                                        #
# Started: 3/10/2005                                                    #
# Last Visited: 5/3/2005                                                #
#                                                                       #
# FUNCTION NAME:                                                        #
#  ap2r: auto polyploids to recombination fraction                      #
# FUNCTION CALL:                                                        #
#  ap2r (infile, outfile, m)                                            #
#                                                                       #
# REFERENCES                                                            #
#  METHODOLOGY: Ripol et. al. 1999. 'Statistical aspects of genetic     #
#               mapping in autopolyploids'. Gene. 235: 31-41.           #
#  CHI-SQUARE VALUES: Rohlf and Sokal. 1995. 'Statistical Tables'       #
#                     W. H. Freeman and Company                         #
#                                                                       #
# USAGE:                                                                #
# ap2r >infile.txt <outfile.txt                                         #
#                                                                       #
# REQUIRED PACKAGES:                                                    #
# Bhat                                                                  #
#                                                                       #
# INPUT:                                                                #
#  - Tab delimited text file for each locus for each of the progeny     #
#    with no header that contains the following information             #
#    Locus | Dose | Progeny 1| Progeny2 | ... | Progeny N               #
#  - The loci should be scored such that                                #
#        4 = Allele is present in the progency                          #
#        1 = Allele is not present in the progeny                       #
#        0 = Allele state could not be determined                       #
#                                                                       #
# OUTPUT:                                                               #
#  - recombination values that are significantly different then 0.5     #
#    for each of the significantly paired loci in a tab delimited table #
#    use the most significant                                           #
#    Locus 1 | Locus 2 | numInfo | genconfig | recomb | pvalue          #
#      Locus 1   - The first locus from a significant pairing           #
#      Locus 2   - The second locus from a significant pairint          #
#      numInfo   - The number of informtive progeny between two loci    #
#      genconfig - The genomic configuration                            #
#      recomb    - The recombinational distance between two loci        #
#      pvalue    - The most significant p-value for the association     #
#                                                                       #
#                                                                       #
# CURRENT NEEDS                                                         #
#  - Write this as a function to allow the user to send the calculation #
#    of sets of pairwise values to nodes on a multiprocessor cluster    #
#    with input in the form of (inputFile, LowLocus, UppLocus)          #
#  - Write the function for the phenotypic probablilites to be a        #
#    variable to choose from the set of possible phenotypic             #
#    probabilities for given genotypic linkage configurations where     #
#    linkage is one of the following                                    #
#       - c11: coupling                 - r11: repulsion                #
#       - ac12: Asymmetric coupling     - r12: Repulsion 1-2            #
#       - dc: Double coupling           - cr: Coupling & Repulsion      #
#       - dr: Doubling Repulsion        - ac13: Assymetric Couping 1-3  #
#       - ac13: Assymetric Coupling 2-3                                 #
#                                                                       #
#                                                                       #
#########################################################################

# MLE
# http://www.weibull.com/LifeDataWeb/maximum_likelihood_estimation_appendix.htm   
# WITH THE N-R METHOD
# www.udayton.edu/~mathdept/EPUMD/2004/Pete.pdf 
#
# More info on MLE from
# http://www.zoo.ufl.edu/bolker/emd/notes/lect13.html
# The MLE is the same as the maximum log-likelihood or the minimum negative log-likelihood estimate.
#
# ac12 took  2.602222 hours on my machine at work for casey's data, and did not find any ac12 significant linkages

# NEED TO INCLUDE THE BHAT PROGRAM TO SOLVE THE EQUATION USING THE NEWTON-RAPHSON ALGORITHM
library(Bhat)  # The bhat package does MLE

TimeStart = Sys.time()

# GET THE INPUT DATA AND SET SOME OF THE VARIABLES
InputFile = "C:/RData/cas_map_data.txt"                               # The matrix of values associated with each locus for each of the progeny
SigOut = "C:/RData/SigOut.txt"
myM = 3                                                               # m is the ploidy = the number of homologous chromosomes

# DEFINE CRITICAL CHI-SQUARE VALUES FOR DIFFERENT ALPHA (P-VALUES) FOR ONE DEGREE OF FREEDOM
chi_05 = 3.841
chi_025 = 5.024
chi_01 = 6.635
chi_005 = 7.879              # THIS IS THE CRITICAL VALUE FOR ONE-TAILED P = 001
chi_001 = 10.828

InputValues <- read.table(InputFile, header=F, as.is=TRUE)            # Read the tab delimited mapping data set, including as.is=TRUE will allow this to correctly read the locus names instead of converting them to factors

numRowsToProcess = (nrow(InputValues))                                # The number of rows in the input data set
numColsToProcess = (nrow(InputValues))
numLocToProcess = ( ncol(InputValues) - 2 )                           # Need to subtract two here for the name of the locus and the 1 or 2
numColsMatrix = (nrow(InputValues))                                   # The number of rows in the input data set
numRowsMatrix = (nrow(InputValues))

# DEFINE THE INTERMEDIATE MATRIX TO HOLD THE DATA
# Initially this will be filled with zeros

# It may be a bad idea to do this since this will hog memory
# May just want to show the output for the significant sets
# ValMatrix <- array(0, c(numRowsMatrix,numColsMatrix,6))
# It has been check

#  ValMatrix (Row,Col, 1) = numInfo: the number of informative loci
#  ValMatrix (Row,Col, 2) = numAB: the number of progeny containing both A and B
#  ValMatrix (Row,Col, 3) = numA: the number of progeny containing only A
#  ValMatrix (Row,Col, 4) = numB: the number of progeny containing only B
#  ValMatrix (Row,Col, 5) = numNull: the number of progeny in the double null class
#  ValMatrix (Row,Col, 6) = fmin: minimum solution to the function
#  ValMatrix (Row,Col, 7) = upp: Upper 95% confidence value for the solution
#  ValMatrix (Row,Col, 8) = low: Lower 95% confidence value for the solution

# DEFINE THE MATRIX TO HOLD THE SIGNIFICANT VALUES
# INTIALLY THIS WILL BE FILLED WITH ZEROS, THIS MAY NEED TO BE SET TO A VARIABLE SIZE
SigMatrix <- array(0, c( 500, 7 ))


# function to test the sum of the four probabilities

TestSum <- function (r, m )
{

TheSum <- ProbAB(r,m) + ProbA(r,m) + ProbB(r,m) + ProbNull(r,m)
return (TheSum)

}

#############################################################
#                                                           #
# RECOMBINATION FRACTION LIKLIHOOD EQUATION                 #
#                                                           #
#############################################################
recomb <- function (x,gencon)
{
  # x is the recomb
  # gencon is the genomic configuratoin - this allows for any genomic configuration to be passed to the function

  m <- myM   # SET m to the global variable myM, this is probably redundant to do this


  # DETERMINE WHICH OF THE PHENOTYPIC PROBABILITES TO USE
  if (gencon == "ac12")     # ASYMMETRIC COUPLING 1-2 : Okay sums to 1
  {
    ProbAB <- function(r,m) { 0.5 - (0.25*((m-2)/(m-1))*r) }
    ProbA <- function(r,m) { 0.25*((m-2)/(m-1))*r }
    ProbB <- function(r,m) { (0.25*(m/(m-1))) + (0.25*((m-2)/(m-1))*r) }
    ProbNull <- function(r,m) {  (0.25*((m-2)/(m-1)))-(0.25*((m-2)/(m-1))*r) }
    #print ("ac12: Assymetric Coupling 12")
    #flush.console()
  }
  else if (gencon == "c11") # Couping 11 :  Okay sums to 1
  {
    ProbAB <- function(r,m) { 0.5 * (1-r) }
    ProbA <- function(r,m) { 0.5 * r }
    ProbB <- function(r,m) { 0.5 * r }
    ProbNull <- function(r,m) { 0.5*(1-r) }
    #print ("c1: Coupling")
    #flush.console()
  }
  else if (gencon == "r11") # Repulsion 11 :
  {
    # ERROR DOES NOT SUM TO 1
    ProbAB <- function(r,m) { (0.25*(m/(m-1))) + (0.5*(1/(m-1))*r) }
    ProbA <- function(r,m) { (0.25*(m/(m-1))) - (0.5*(1/(m-1))*r) }
    ProbB <- function(r,m) { (0.25*(m/(m-1))) - (0.5*(1/(m-1))*r) }
    ProbNull <- function(r,m) { (0.25*(m/(m-1))) + (0.5*(1/(m-1))*r) }
    #print ("r11: Repulsion")
    #flush.console()
    ProbAB <- function(r,m) { (0.25*(m/(m-1))) + (0.5*(1/(m-1))*r) }
    ProbA <- function(r,m) { (0.25*(m/(m-1))) - (0.5*(1/(m-1))*r) }
    ProbB <- function(r,m) { (0.25*(m/(m-1))) - (0.5*(1/(m-1))*r) }
    ProbNull <- function(r,m) { (0.25*(m/(m-1))) + (0.5*(1/(m-1))*r) }
  }
  else if (gencon == "r12") # Repulsion 12 : Okay sums to 1
  {
    ProbAB <- function(r,m) { ( 0.125*((3*m-4)/(m-1)) ) - ( 0.5*(1/(m-1))*r ) }
    ProbA <- function(r,m) { (0.125*((m)/(m-1))) - (0.5*(1/(m-1))*r) }
    ProbB <- function(r,m) { (0.125*((3*m)/(m-1))) + (0.5*(1/(m-1))*r)}
    ProbNull <- function(r,m) { (0.125*((m-4)/(m-1))) + (0.5*(1/(m-1))*r) }
    #print ("r1: repulsion 1-2")
    #flush.console()
  }
  else if (gencon == "dc") # Double coupling : ERROR - DOES NOT SUM TO ONE
  {
    # ERROR DOES NOT SUM TO 1
    ProbAB <- function(r,m) { (0.25*((3*m-2)/(m-1))) - (0.5*((m-2)/(m-1))*r) + (0.25*((m-2)/(m-1))*r^2) }
    ProbA <- function(r,m) { (0.5*((m-2)/(m-1))) - (0.25*((m-2)/(m-1))*r^2) }
    ProbB <- function(r,m) { (0.5*((m-2)/(m-1))) - (0.25*((m-2)/(m-1))*r^2) }
    ProbNull <- function(r,m) { (0.25*((m-2)/(m-1))) - (0.5*((m-2)/(m-1))*r) + (0.25*((m-2)/(m-1))*r^2) }
    #print ("dc: Double Coupling")
    #flush.console()
  }
  else if (gencon == "cr")
  {
    ProbAB <- function(r,m) { (0.125*((5*m-4)/(m-1))) - ( 0.125*((m-6)/(m-1))*r ) - (0.25*((1/(m-1))*r^2)) }
    ProbA <- function(r,m) { (0.125*(m/(m-1))) + (0.125*((m-6)/(m-1))*r) + (0.25*(1/(m-1))*r^2) }
    ProbB <- function(r,m) { (0.125*(m/(m-1))) + (0.125*((m-6)/(m-1))*r) + (0.25*(1/(m-1))*r^2) }
    # MAY NEED TO FIX THE FOLLOWING TO TAKE INTO ACCOUNT THE TYPO ERROR IN THE MANUSCRIPT, I DON'T KNOW IF THE FIRST MAJOR IS AN ADDITION OR A SUBTRACTION
    ProbNull <- function(r,m) { (0.125*((m-4)/(m-1))) - (0.125*((m-6)/(m-1))*r) - (0.25*(1/(m-1))*r^2)  }
    #print ("c1: Coupling & Repulsion")
    #flush.console()
  }
  else if (gencon == "dr")
  {
    ProbAB <- function(r,m) { (0.0625*((9*m^2 - 34*m + 24)/((m-1)*(m-3)))) + (0.5*((m-4)/((m-1)*(m-3)))*r) + (0.5*(1/(m-1)*(m-3))*r^2) }
    ProbA <- function(r,m) { (0.0625*((3*m^2 - 10*m)/(m-1)*(m-3))) - (0.5*((m-4)/((m-1)*(m-3)))*r) - (0.5*(1/(m-1)*(m-3))*r^2) }
    ProbB <- function(r,m) { (0.0625*((3*m^2 - 10*m)/(m-1)*(m-3))) - (0.5*((m-4)/((m-1)*(m-3)))*r) - (0.5*(1/(m-1)*(m-3))*r^2) }
    ProbNull <- function(r,m) { (0.0625*((m^2 - 10*m + 24)/((m-1)*(m-3)))) + (0.5*((m-4)/((m-1)*(m-3)))*r) + (0.5*(1/(m-1)*(m-3))*r^2) }
    #print ("dr: Coupling")
    #flush.console()
  }

  g <- (ProbAB(x,myM)^numAB)*(ProbA(x,myM)^numA)*(ProbB(x,myM)^numB)*(ProbNull(x,myM)^numNull)

  return( g )

}

#############################################################
#
# NEGATIVE LOG OF THE RECOMB FUNCTION                 #
#
#############################################################
nlrecomb <- function (x)
{

  #gencon will need to be a global variable that is set before calls to this function are made
  gencon = mygencon # set the local gencon variable to the global variable mygencon

  g <- recomb (x, gencon)
  return( -log(g))

}

#############################################################
#
# NEGATIVE LOG OF THE RECOMB FUNCTION                 #
#
#############################################################
# rhat is the estimated value
# gencon is used to determine the genotypic configuration to evaluate
lx <- function (rhat, gencon)
{

l <- ((recomb( rhat, gencon ))/(recomb( 0.5, gencon)))
return (-2*log(l))

}

#############################################################
#
# KOSAMBI'S MAPPING FUNCTION
#
#############################################################
# Converts recombination fraction to distance units in CentiMorgan
# using Kosambi's Mapping Function
kosambi <- function (r)
{

d <- (0.25 * log((1+2*r)/(1-2*r)))
return (100*d)

}


#############################################################
#
# HALDANE'S MAPPING FUNCTION
#
#############################################################
# Converts recombination fraction to distance units in CentiMorgan
# using Haldane's Mapping Function.
haldane <- function (r)
{

d <- (-0.5 * log(1-2*r))
return (100*d)

}


# INITIALIZE THE COUNTER VARIABLES
# A represent dose of first type, B of second type
numSig = 0               # The number of significant linkages that are found
numInfo = 0              # The number of informative loci
numAB = 0                # The nubmer of progeny showing both loci
numA = 0                 # The number of progeny showing only A
numB = 0                 # The number of progeny showing only B
numNull = 0              # The number of progeny showing neither

# EVALUATE EACH PAIRWISE
for (i in 1:numRowsToProcess)
{

 for (j in 1:numColsToProcess)
 {

 if (i != j) # ONLY CALCULATE THESE VALUES FOR I NOT EQUAL TO J, ALSO LIMITING TO J < I will prevent duplicate comparisons and will do this on a triangle
 {              # Begin of i not equal to j
 print (c("i: ",i,"j: ",j," Comparing:", InputValues[i,1],":",InputValues[j,1]), quote=FALSE)
 flush.console()
 # ONLY DO THE FOLLOWING FOR I NOT EQUAL TO J

 # EVALUATE THE COUPLING AND DETERMINE THE PROBABILITY MODEL TO USE

  for (k in 1:numLocToProcess)
  {

   # EVALUATE THE SITUATION AND INCREMENT THE PROPER COUNTER
   # This is pseudocode, need to write in proper R format
   # Values where i,k or j,k are equal to 0 are ignored
   if (InputValues[i,k] == 4)
   {
      if (InputValues[j,k] == 4)
      {
           numAB = numAB + 1
           numInfo = numInfo + 1
      }
      else if (InputValues[j,k] == 1)
      {
           numA = numA + 1
           numInfo = numInfo + 1
      }
   }
   else if (InputValues[i,k] == 1)
   {
       if (InputValues[j,k] == 4)
       {
           numB = numB + 1
           numInfo = numInfo + 1
       }
       else if (InputValues[j,k] == 1)
       {
           numNull = numNull + 1
           numInfo = numInfo + 1
       }
   }


  } #END OF FOR EACH K

 # SHOW THE CURRENT VALUES
 #print ( c("Info", numInfo), quote=FALSE)
 #print (c("AB", numAB), quote=FALSE)
 #print (c("A", numA), quote=FALSE)
 #print (c("B", numB), quote=FALSE)
 #print (c("Null", numNull), quote=FALSE)
 #flush.console()

 # DETERMING GENCON TO USE
 mygencon = "c11"

 # MLE ESTIMATION HERE
 # It would be nice if I could generate one that also yields a 95% confidence
 mylist <- list(label="dist", est=0.45, low=0.0, upp=0.5)
 result <- dfp( mylist, f=nlrecomb, tol=1e-01  )

 # TEST FOR SIGINIFCANCE USING THE CHI-SQUARE TEST
 test <-  lx(result$est, mygencon)

 if ( test > chi_001 )
 {
  numSig = numSig + 1
  print ("we have a significante value")
  print (c("Test Value:",test, "Chi", chi_001 ))
  print (c("Comparing:", InputValues[i,1],":",InputValues[j,1]), quote=FALSE)
  print (c("NumberInfo", numInfo), quote=FALSE)
  print (c("Model: ",mygencon,"r: ",result$est,"Hald", haldane(result$est),"Kosambi", kosambi(result$est) ), quote=FALSE)
  flush.console()
  # LOAD ANY SIGNIFICANT RESULTS TO THE SIG OUT MATRIX
  SigMatrix [numSig, 1] = InputValues[i,1]
  SigMatrix [numSig, 2] = InputValues[j,1]
  SigMatrix [numSig, 3] = numInfo
  SigMatrix [numSig, 4] = mygencon
  SigMatrix [numSig, 5] = result$est
  SigMatrix [numSig, 6] = haldane(result$est)
  SigMatrix [numSig, 7] = kosambi(result$est) )

}

if ( test < chi_001 )
 {
  print ("nothing to see here")
 }


 # LOAD RAW RESULTS TO THE VALUE MATRIX
 #ValMatrix[i,j,1] = numInfo
 #ValMatrix[i,j,2] = numAB
 #ValMatrix[i,j,3] = numA
 #ValMatrix[i,j,4] = numB
 #ValMatrix[i,j,5] = numNull

# ValMatrix[i,j,6] = rFract.fmin
# ValMatrix[i,j,7] = rFract.low
# ValMatrix[i,j,8] = rFract.upp

 # RESET THE COUNTER VARIABLES
 numInfo = 0
 numAB = 0
 numA = 0
 numB = 0
 numNull = 0

 }  # END OF FOR I NOT EQUAL TO J
 } # END OF FOR EACH J
} # END OF FOR EACH I

# OUTPUT A TEXT FILE OF THE SIGNIFICANT COMBINATIONS
write.table(SigMatrix, SigOut)
write.table(SigMatrix, "W:/common_pgml/Casey/jamie/sig_ac12_m4.txt")
# SHOW HOW MUCH TIME IT TOOK TO COMPLETE THE OPERATION
TimeEnd = Sys.time()
TimeDiff <- TimeEnd - TimeStart
print(TimeDiff)




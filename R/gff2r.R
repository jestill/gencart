##########################################################
#
# PROGRAM: GDAWGS::DataIO
# VERSION: 0.0.2
#  AUTHOR: James C. Estill
# CONTACT: jestill_#_uga.edu
# STARTED: 2/13/2005
# UPDATED: 2/14/2005
#
# DESCRIPTION:
#

# REQUIRED LIBRARIES:
#  None
#
##########################################################

##########################################################
# EXAMPLE USE                                            #
##########################################################

# Adds an isppoint column
test$ispoint <- c(TRUE,FALSE) #assignment 1

# Working out orientation
getdiff <- function (inframe)
{
  inframe$start - inframe$end
}

test <- gff2r("C:/RData/atest2.gff")
numRows = (length(test$seqname))                          # How many rows are here

#r2gff(test)                                              # Calling this should put the data.frame in standard gff format and export

##########################################################
# SUBFUNCTION                                            #
# gff2r                                                  #
#
# INPUT VARIABLES:
#   *dfName = The name of the data frame
#      Default value =
#   *inputpath = The path of the input file that will be
#    converted to a R data frame. This should be in gff format
#      Default value =
#   *header = Boolean describing if the gff frame has a header
#      Default value = FALSE
#   * May also want to designate the type of data that this will be
#       -ie. GeneticMapData
#            PhysicalMapData
#            SequenceMapData
#            RecombProximityData
#            Raw genetic mapping data
#            (It may be possible and even necessary to incorporate the ability
#             the calculate a genetic map from a given set of genetic map data)
#             This could use complied C functions from MapMaker if necessary ...
#
##########################################################

gff2r <- function (inFile, header="FALSE")
{
  gffnames = c("seqname","source","feature","start","end","score","strand","frame","group")    # Load GFF formated headers to an array
  dframe <- read.table(inFile, sep="\t", header=F, col.names=gffnames, na.strings=".")                   # Read the GFF text file into a data frame using the header array and treating decimal values as NA
 
  # Check the number of columns to make sure that there are exactly nine.
  # In the future this should occur on just a sample before the entire file is read
  numCols = length(dframe)                                                                     # See how man columns are in the file
  
  # Before continue make sure the data frame is the expected size
  if (numCols > 9)                                                                             # If there are too many columns, stop and alert the user
     stop("The input file contains too many columns")
  if (numCols < 9)                                                                             # If there are too few colums, stop and alert the user
     stop("The input file contains too few columns")

  # May need to add Integrety check for value orientation. This will make sure that every Start, is the lower value

  dframe$ispoint <- (dframe$start == dframe$end)                                               # Add an ispoint column to the data frame. This is set to true when the start and end are the same indication a point feature

  #strsplit(test3$group,"t")

  # This would be the place to parse extra data that are stored in the last column
  # Valid values in the gff3 format to fetch here include
  # ID, Name, and Parent
  # They  could be split using strsplit(dframe$group,";")

  dframe
                                                                                      # Return the data frame that was use as input
} # End of gff2r subfunction


# VERSION HISTORY:
# 2/13/2005 : 0.0.1
#  - Began coding the gff2r function
#  - Initial gff2r is to use the read.table function in R
#    assign the columns names as insert NA for . fields
# 2/14/2005 : 0.0.2
#  - Adding point, segment determinations to gff2r
#  - Only allow tabs to be used as separators, other white spaces are ignored in gff2r

##########################################################
#
# PROGRAM: gencart_core.r
#  AUTHOR: James C. Estill
# CONTACT: jestill_#_uga.edu
# STARTED: 2/13/2005
# UPDATED: 07/13/2009
#
# DESCRIPTION:
# GenCart core functions. Currently just for me to use
# but will try to keep the documentation up to date.

# REQUIRED LIBRARIES:
#  None
#
##########################################################

# The following library is used to derive substrings of tag values
# from the input GFF3 files
library(stringr);

#-----------------------------------------------------------+
# SUBFUNCTION: gff2r                                        |
#-----------------------------------------------------------+
# Load data in the GFF foramt. This easily handles data in
# the GFF2 format file.
#
# INPUT VARIABLES:
#   *inputpath = The path of the input file that will be
#    converted to a R data frame. This should be in gff format
#      Default value =
#   *header = Boolean describing if the gff frame has a header
#      Default value = FALSE
# USAGE: test <- gff2r("test.gff");

gff2r <- function (infile, header="FALSE") {

  # GFF HEADERS
  gffnames = c(
    "seqname",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attribute");

  # READ THE TABLE, using a basic read table command for now
  # This will also reseat any exsting header names to the
  # commonly used header names for GFF
  if (header == FALSE) {
    dframe <- read.table(infile, header=F, col.names=gffnames);
  } else {
    dframe <- read.table(infile, header=T, col.names=gffnames);
  }
  
  numCols = length(dframe);

  # Return the data, may want to reset the the headers to standard
  # names used by gencart
  return (dframe);
  
}


loadGFF<- function (infile, header="FALSE") {

  # As currently written the substring functions
  # require stringr
  
  # GFF HEADERS
  gffnames = c(
    "seqname",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attribute");

  # READ THE TABLE, using a basic read table command for now
  # This will also reseat any exsting header names to the
  # commonly used header names for GFF
  if (header == FALSE) {
    dframe <- read.table(infile, header=F, col.names=gffnames, 
                         stringsAsFactors=F);
  } else {
    dframe <- read.table(infile, header=T, col.names=gffnames,
                         stringsAsFactors=F);
  }


  # Length of the features
  dframe$length <- abs(dframe$end - dframe$start) + 1;

  # Some basic information about the data frame, not currently used
  numCols = length(dframe);
  numRows = length(dframe$attribute);
  
  #-----------------------------+
  # SUPPORT GFF3 ATTRIBUTE      |
  # TAGS                        |
  # CURRENTLY SUPPORTS          |
  #  - ID                       |
  #  - Name                     |
  #  - Parent                   |
  #-----------------------------+
  # Extracting ID values from GFF3
  # First look for strings terminated by semicolons starting with ID, then
  # add values terminated by end of line. The strings are then cleaned up
  # by removing the tag value and a terminating semicolon if one
  # exists
  dframe$ID <- str_extract(dframe$attribute, "ID=.*;");
  dframe$ID <- ifelse( is.na(dframe$ID),
                      str_extract(dframe$attribute, "ID=.*$") ,
                      dframe$ID);
  dframe$ID <- sub("ID=", '', dframe$ID);
  dframe$ID <- sub(";", '', dframe$ID);

  # Name values
  dframe$Name <- str_extract(dframe$attribute, "Name=.*;");
  dframe$Name <- ifelse( is.na(dframe$Name),
                      str_extract(dframe$attribute, "Name=.*$") ,
                      dframe$Name);
  dframe$Name <- sub("Name=", '', dframe$Name);
  dframe$Name <- sub(";", '', dframe$Name);

  # Parent values
  dframe$Parent <- str_extract(dframe$attribute, "Parent=.*;");
  dframe$Parent <- ifelse( is.na(dframe$Parent),
                      str_extract(dframe$attribute, "Parent=.*$") ,
                      dframe$Parent);
  dframe$Parent <- sub("Parent=", '', dframe$Parent);
  dframe$Parent <- sub(";", '', dframe$Parent);
  
  #///////////////////////////////////////////////////////////
  # TO DO :
  #///////////////////////////////////////////////////////////
  # Lookup for child parent relationships, should be able to
  # use netsted integers for these searches. But the algorithm
  # to fill this will be complicated in R
  
  # Return the GFF data frame
  return (dframe);
  
}


#///////////////////////////////////////////////////////////
# Summarize a GFF file passed to the function
# as a dataframe
#///////////////////////////////////////////////////////////
summarizeGFF <- function (dframe, genomeName = "Genome",
                          fileRoot ="genome") {

  # Get the nonredundant list of scaffolds and report length
  
  # Set some file paths (These will be PNG files)
  outExonCountImage = (paste(fileRoot, "_exon_count.png",  sep=""));
  outCDNALengthImage = (paste(fileRoot, "_gene_summary.png",  sep=""));
  outSummaryImage = (paste(fileRoot, "_gene_summary.png",  sep=""));
  outSummaryFile = (paste( fileRoot,"_genome_report.txt",  sep=""));
  
  # Get the nonredundant list of features and counts and report names
  # and values if < 100 features

  # Get the subset of just the exons to work with
  exons <- subset(dframe, dframe$feature=="exon");
  introns <- subset(dframe, dframe$feature=="intron");

  exonSummary <- summary(exons);

  # This will break the sources into individual colums, so if the
  # exons are derived from different sources, this will be relfected
  # in the output

  # Number of exons per gene
  exonsPerGene <- table(exons$Parent, exons$source);
  # The following is the list of exonsPerGene sorted
  exonCountList <- sort(unique(exonsPerGene));
  
  # Get the max vale from the unique list of the exon count
  # thus exonCountHist$breaks is the number of exons and
  #      exonCountHist$counts is the number of genes with that many exons
# Commented out the following since they are done below  
#  exonCountHist <- hist(exonsPerGene,
#                        max(unique(exonsPerGene)), col="green",
#                        main = c(genomeName, "Exons per Gene") );
  # Generate historgram of exon lengths

  # Generate histogram of gene lengths
  geneLengths <- tapply(exons$length, exons$Parent, FUN=sum);
# Commented out the following since they are done below  
#  geneLengthHist <- hist(geneLengths, 150,
#                           col="red", xlim = c(0,8000),
#                           main= c(genomeName, "cDNA Length"),
#                           xlab="cDNA Length (bp)");



  
  # Test of drawing plots ot outfile
  #  par(new=T);
  oldpar <- par(no.readonly=TRUE)
  par <- par(mfrow = c(4,1), pty = "m", bg="white");

  # Do the histograms
  exonCountHist <- hist(exonsPerGene,
                        max(unique(exonsPerGene)), col="green",
                        main = c(genomeName, "Exons per Gene"),
                        xlab="Count of exons Per Gene",
                        xlim = c(0,30) );

  geneLengthHist <- hist(geneLengths, 150,
                         col="red", xlim = c(0,8000),
                         main= c(genomeName, "cDNA Length"),
                         xlab="cDNA Length (bp)");

  exonLengthHist <-hist(exons$length, 3000,
                        xlim = c(0,2000), col="blue",
                        xlab="Exon Length (bp)",
                        main = c(genomeName, "Exon Length"));

  exonLengthHist <-hist(introns$length, 3000,
                        xlim = c(0,2000), col="brown",
                        xlab="Intron Length (bp)",
                        main = c(genomeName, "Intron Length"));
  
  # PLOT THE IMAGE
  dev.copy (png, filename=outSummaryImage, pointsize = 14,
            height=1000, width=720, bg="white" );
  dev.off();

  par(oldpar);

#  result <- list ( exonLengthSummary = summary(exons$length),
#                   exonLengthHistogram =  exonLengthHist);

  # The following just provides the summary for the values and
  # not the full histograms, althought these can be returned
  # as shown above.
  result <- list ( exonLengthSummary = summary(exons$length),
                   cDNALengthSummary = summary(geneLengths),
                   intronLengthSummary =  summary(introns$length),
                   exonCountSummary = summary(as.vector(exonsPerGene))
                  );
  
  return(result);
  
}


# Given a tag as input return the value
# assumes that input is a list that has already
# been split by ;
atrKeyVal <- function ( atrList, tag ) {

  # Set default of tag value
  tagValue <- NA;
  
  # Goal is to take a list
  # return the value for a given tag
  # ie take list of attribute values and
  # return fa

  # First split to tag value pairs
  splitList <- strsplit(atrList, "=");

  # Then fetch the value of interest
  # First make sure that grep returns a valu
  rowToFetch <- grep( tag , splitList);

  if (length(rowToFetch) != 0 ) {
    tagValue <- splitList[[rowToFetch]][2];
  }
  
  return(tagValue);
  
}




#-----------------------------------------------------------+
# DATA CATEGORIZATION                                       |
#-----------------------------------------------------------+
# This will categorize data into classes based on the
# classification method specified by the classMethod option
#
# *classVals = The values to classify
# *numClass  = The number of classes to create
# *classMethod = The classification method to use
#          - quantile
#          - bagged
#          - kmeans
#          - equal_interval
#          - pam

datcat <- function (classVals,
                    numClass = 10,
                    classMethod = "quantile",
                    numBreaks = 200,
                    plotMaxX = max(classVals),
                    plotMinX = min(classVals),
                    directToFile = FALSE,
                    copyPlot = TRUE,
                    drawPlot = TRUE,
                    flipColRamp = FALSE,
                    imageWidth = 17.15,
                    imageHeight = 6,
                    plotTitle = "Categorized Data",
                    plotXLabel = "Classified Values") {
  
  # REQUIRED LIBRARY FOR THE COLOR RAMPS
  library(colorRamps);

  #-----------------------------+
  # SET EXPORT PATHS            |
  #-----------------------------+
  outColFile = (paste("class_",classMethod, "_",
                      numClass, "_rgb.txt",  sep=""));
  outHeader = (paste("class_", classMethod, "_",
                     numClass, sep=""));
  outNAFile = (paste("class_", classMethod, "_",
                     numClass, ".NA",  sep=""));
  outCatPlotImage = (paste("class_", classMethod, "_",
                        numClass, ".tiff",  sep=""));
# Plot image as a PNG file  
#  outCopyPlotImage = (paste("class_", classMethod, "_",
#                        numClass, ".png",  sep=""));
  outCopyPlotImage = (paste("classx_", classMethod, "_",
                            numClass, ".png",  sep=""));

  #-----------------------------+
  # ESTABLISH PLOT OBJECT       |
  #-----------------------------+
  # If we are printing directly to a file object
  if (directToFile == TRUE) {
    # The following if for figures that will be sent to print
    # in PLOS genetics, supporting figures can be in any size or format
    # PNG format probably makes the most sense for these supporting
    # figures
    # PLOS GENETICS WANTS LZW COMPRESSION FOR TIFFS
    # Set the units to centimeters to facilitate drawing for PLOS Genetics
    # 1 Col figures - Width = 8.25 cm
    # 1.5 Col figures - Width = 12.06 - 12.7 cm
    # 2 Col figures - Width = 17.15 cm
    tiff(filename = outCatPlotImage, pointsize = 12,
         width = imageWidth, height = imageHeight,
         units="cm",
         bg="white", compression="lzw", res=300, antialias="n");
    
  } else {
    op <- par(bg = "white");
  }
  
  #-----------------------------+
  # COLOR PALETTE               |
  #-----------------------------+
  colRamp <- blue2red(numClass);
  #colRamp <- green2red(numClass);

  # If you want the inverse ramp
  if (flipColRamp == TRUE ) {
    altRamp <- blue2red(numClass);
    for (i in 1:numClass) {
      altRamp[i] = colRamp[(numClass+1)-i];
    }
    colRamp <- altRamp;
  }
  
  minVal = min( classVals );
  maxVal = max( classVals );


  #-----------------------------------------------------------+
  # CLASSIFY THE VALUES                                       |
  #-----------------------------------------------------------+

  #-----------------------------+
  # KMEANS CLUSTER CLASSIFY     |
  #-----------------------------+
  # Minimizes the variation of the distances from the mean within each class 
  # This is included with R and can be used when the package e1071
  # is not installed
  if (classMethod == "kmeans") {
    
    groupK <- kmeans(classVals, numClass);
    classK <- groupK$cluster;
    
    cols <- match(groupK$cluster, order(groupK$centers));
    ncl <- nrow(groupK$centers);
    
    stbrks <- matrix(unlist(tapply(classVals, factor(cols), range)),
                     ncol = 2, byrow = TRUE);
    
  }

  #-----------------------------+
  # BAGGED CLUSTERING           |
  #-----------------------------+
  if (classMethod == "bagged") {
  
    # blucst requires the e1071 package
    library(e1071);

    # Set a random number seed
    set.seed(20080624);

    newBreaks <- bclust(classVals,  iter.base=20,
                        centers=numClass, verbose=FALSE);
    cols <- match(newBreaks$cluster, order(newBreaks$centers));
    stbrks <- matrix(unlist(tapply(classVals, factor(cols), range)),
                     ncol = 2, byrow = TRUE);
  }

  #-----------------------------+
  # QUANTILES                   |
  #-----------------------------+
  if (classMethod == "quantile") {
    
    quantVal = 1/numClass;
    stbrks <- matrix (NA, ncol=2,nrow=numClass);
    quantRes <- quantile(classVals, probs = seq (0,1, quantVal),
                         names=FALSE );
    
    for (i in 1:numClass) {
      # Since quantiles represent percent below the given value we get..
      stbrks [i,1] = quantRes[i]+1; 
      stbrks [i,2] = quantRes[i+1];
    }

    breaks <- matrix (NA, ncol=1,nrow=numClass+1);
    breaks[1] = minVal;
    for (i in 1:numClass) {
      breaks[i+1] = stbrks[i,2];
    }
  
    cols <- findInterval(classVals, breaks , all.inside=TRUE);
  } # End of if classMethod is quantile

  #-----------------------------+
  # EQUAL INTERVAL              |
  #-----------------------------+
  if (classMethod == "equal_interval" ) {
    
    stbrks <- matrix (NA, ncol=2,nrow=numClass);
    stepSize = (maxVal - minVal)/numClass;
    breakTop = minVal;              # Bottom of the break
    breakBot = minVal;              # Top of the break
    
    for (i in 1:numClass) {
      breakTop = breakTop + stepSize;
      stbrks [i,1] = breakBot;
      stbrks [i,2] = breakTop;
      breakBot = breakTop + 1;
    }
    
    breaks <- matrix (NA, ncol=1,nrow=numClass+1);
    breaks[1] = minVal;
    for (i in 1:numClass) {
      breaks[i+1] = stbrks[i,2];
    }
    
    cols <- findInterval(classVals, breaks , all.inside=TRUE);
    
  }

  #-----------------------------+
  # PARTITION AROUND MEDIODS    | 
  #-----------------------------+
  # Very slow method for a large number of data points
  if (classMethod == "pam" ) {
    
    library("cluster");
    newBreaks <- pam(classVals,  numClass, diss = FALSE,
                     metric = "euclidean", cluster.only = TRUE
                     );
    
  }
  
  #-----------------------------+
  # STANDARD DEVIATIONS         |
  #-----------------------------+
  # TO DO: Get the mean and do above and below number of std deviations
  #        in this case the number of classes will need be an even number,
  #        otherwise could use odd numbers and the program will generate
  #        2x number of classes. ie, numClass = 5 will actually create
  #        5 classes. To get the "standard" deviations will probably
  #        need to know if the distribution is normal.


  
  #-----------------------------------------------------------+
  # DRAW THE PLOT                                             |
  #-----------------------------------------------------------+
  # Changing these to use the min value and max value for limits
  
  if (drawPlot == TRUE ) {

    subText =  (paste("class_", classMethod, "_", numClass,  sep=""));
  
    #-----------------------------+
    # SET UP THE PLOTTING OBJECT  |
    #-----------------------------+
    # Previous method for the plotting object 08/05/2009
    #par <- par(mfrow = c(2,1), pty = "m", bg="white", family="arial");
    par <- par(mfrow = c(2,1), pty = "m", bg="white", family="arial");

    #-----------------------------+
    # TOP HISTOGRAM PLOT          |
    #-----------------------------+
    # HISTOGRAM OBJECTS GIVES INFORMATION FOR EQUAL AREA
    # AGE CLASSES, FOR N BREAKS
    # CAN USE $breaks to delimit and use as a classifier
    # bottom, left, top, right
    # OLD MARGINS 08/05/2009 
    par(mar = c (0.5,4,1,.1) + .1);
    #par(xaxs="i");
    # NEW MARGINS
    #par(mai = c (0.1,.4,.4,.2) + .1);

    # Using xaxt="n" will supress x axis labels
    histVals <-hist(classVals, numBreaks, col="black", xlab=NULL,
                    xaxt = "n",
                    main=plotTitle,
                    xlim=c(plotMinX,plotMaxX) );
    #-----------------------------+
    # EMPIRICAL CUMULATIVE DISTN  |
    #-----------------------------+
    # Reset values for margin since this does not have a header
    # bottom, left, top, right
    # The following is the default numbers
    # OLD MARGINS 08/05/2009 
    par(mar = c (4,4,0.3,.1) + .1);
    # Commenting out yaxt, xaxt will autodraw axis
    plot.ecdf(classVals, verticals=TRUE , cex=0.1,
              main = "", xlab=plotXLabel,
              yaxt="n",
              xaxt="n",
              ylim=c(-0.2,1), xlim=c(plotMinX,plotMaxX));

    # PERCENT COMPOSITION AXIS
    #axis(1, labels=c(50,60,70,80,90), at=c(50000,60000,70000,80000,90000));
#    axis(1, labels=c(2,4,6,8,10,12,14),
#         at=c(20000,40000,60000,80000,100000,120000,140000));
    # FL LRP COUNT
#    axis(1, labels=c(1,5,10,15,20),
#         at=c(10000,50000,100000,150000,200000));
    # INSERTION DATE
    axis(1, labels=c(.25,0.5,0.75,1,1.25,1.5),
         at=c(250000,500000,750000,1000000,1250000,1500000));

    # Y-AXIS
    axis(2, labels=c(.2,.4,.6,.8,1.0), at=c(.2,.4,.6,.8,1.0));

    # DRAW RUG OF VALS, line for each datapoint
    rug(classVals);
  
    #-----------------------------+
    # ADD COLOR CLASS BOXES       |
    #-----------------------------+
    par(new=T);

    # INITIALIZE MATRICES
    xVal <- matrix ( 1, nrow = numClass, ncol = 1);
    yVal <- matrix ( -0.05, nrow = numClass, ncol = 1);
    recVals <- matrix(c(100000,0.1),
                      nrow=numClass, ncol=2, byrow=TRUE);
    
    # SET VALS
    for (i in 1:numClass) {
      xVal[i] = stbrks[i,1] + ( (stbrks[i,2]- stbrks[i,1])/2 );
      recVals[i,1] =  stbrks[i,2]- stbrks[i,1];
    }
    
    # DRAW BOXES
    symbols ( xVal, yVal, rectangles = recVals, fg=colRamp , bg=colRamp,
             axes = FALSE, xlab="",ylab="", inches=FALSE,
             ylim=c(-0.2,1), xlim=c(plotMinX,plotMaxX) );

    #-----------------------------+
    # ADD LINES TO ECDF           |
    #-----------------------------+
    xLinesLeft <- matrix (stbrks[,1], nrow = numClass, ncol=2 );
    xLinesRight <- matrix (stbrks[,2], nrow = numClass, ncol=2 );
    yLines <- matrix( c(0,1), nrow=numClass, ncol=2, byrow=TRUE);
  
    # DRAW LINES ITERATIVELY
    for (i in 1:numClass) {

      # TO DO, print the start/stop values or set these as parameters
      #        that can be returned from the function
      
      # Left side of breaks
      lines ( c (xLinesLeft[i,1],xLinesLeft[i,2]),
             c (yLines[i,1], yLines[i,2]),
             col = "gray", lty=3);
  
      # Right side of breaks
      lines ( c (xLinesRight[i,1],xLinesRight[i,2]),
             c (yLines[i,1], yLines[i,2]),
             col = "gray", lty=3);
      
    }
    
    par(new=F);
    
  } # End of if draw plot

  #-----------------------------+
  # COPY THE PLOT               |
  #-----------------------------+
  if (copyPlot == TRUE ) {
 #The following works with the previously set margins    
#    dev.copy (png, filename=outCopyPlotImage, pointsize = 12,
#                height=709, width=2025, bg="white" );
    # The following attempt to work with pointsize
#        dev.copy (png, filename=outCopyPlotImage, pointsize = 20,
#                height=709, width=2025, bg="white" );
    # Modify the margins
    #-----------------------------+
    # DRAWING TIFF                |
    #-----------------------------+
    par(cex.axis =.6, cex.lab = 0.6, cex.main = 0.6);
    dev.copy (tiff, filename="insert_date_key.tiff",
              units="cm", res=300,
              pointsize = 8,
              height=8, width=17.15, bg="white" );
    dev.off();
  }


#  #-----------------------------+
#  # CREATE OUT MATRIX           |
#  #-----------------------------+
#  # This is the classified output matrix
#  numRows <- length(classVals);
#
#  classMat <- matrix (NA, ncol=2,nrow=numRows);
#  classMat[,1] = seqData[[1]];
#  classMat[,2] = cols;
#
#  # CONVERT MATRIX TO DATA FRAME
#  classFrame <- as.data.frame(classMat);
    
  # cols is the classified data to return, classified as
  # "colors" represented by
  # May also consider returning a more compelte object
  return(cols);
  
} # End of the datcat subfunction





#-----------------------------------------------------------+
# GENOME DISTRIBUTION                                       |
#-----------------------------------------------------------+
# Currently optimized for maize LTR_retro data
# This will generate a 'heatmap' for a set of features
# for the entire genome. Takes a gff dataframe as its input
# Modifying heatmap
# * genomeGFF - Genome data in the GFF format as imported by
#               the gff2r function
# * classVals - Classification Report from dat
#
genomeHeatmap <- function (genomeGFF,
                           classVals,
                           flipColRamp = FALSE,
                           copyPlot = FALSE,
                           directToFile = TRUE,
                           alignCent = TRUE,
                           outPlotImage = "genomeHeatmap.tif",
                           mainLabel = "Genome Heatmap",
                           numCat = 10,
                           yLoc = 0,
                           plotWidth = 890,
                           imageWidth = 17.15,
                           imageHeight = 15,
                           boxHeight = 50,
                           boxWidth = 3) {


  # Append the cat data to the gff file, this will go in the
  # 10th position. Can do a last add check, just to make sure
  # and set this as the binCol
  genomeGFF$cat <- classVals;
  binCol = 10;
  
  #-----------------------------+
  # ESTABLISH PLOT              |
  #-----------------------------+
  # If we are printing directly to a file object
  if (directToFile == TRUE) {
    # The following if for figures that will be sent to print
    # in PLOS genetics, supporting figures can be in any size or format
    # PNG format probably makes the most sense for these supporting
    # figures
    # PLOS GENETICS WANTS LZW COMPRESSION FOR TIFFS
    # Set the units to centimeters to facilitate drawing for PLOS Genetics
    # 1 Col figures - Width = 8.25 cm
    # 1.5 Col figures - Width = 12.06 - 12.7 cm
    # 2 Col figures - Width = 17.15 cm
    tiff(filename = outPlotImage, pointsize = 12,
         width = imageWidth, height = imageHeight,
         units="cm",
         bg="white", compression="lzw", res=300, antialias="n");
    
  } else {
    op <- par(bg = "white");
  }

  # GENERAL PLOT PARAMATERS
  # family is the font family to use
  # mar is the margins as (bottom, left, top, right)
  # units are ---?
  par(family="arial",
      mar=c(4,0,1,0) + 0.0 );

  plot(c(0, plotWidth), c(0, 1000), type = "n", xlab="Location (MB)", ylab="",
       main=mainLabel, yaxt="n", xaxt="n", frame.plot=FALSE);

  # Draw the bottom axis, set the tickmarks and the labels
  pw <- plotWidth;
  myAxisScale <- c(1,50,100,150,200,250,300);
  magScale <- myAxisScale*boxWidth;
  magScale[1] <- 1;
  myAxisLabel <- myAxisScale;
  axis(side=1, magScale, labels=myAxisLabel);
  
  #-----------------------------+
  # SET COLOR RAMP              |
  #-----------------------------+
  colRamp <- blue2red(numCat);

  # If you want the inverse ramp
  if (flipColRamp == TRUE ) {
    altRamp <- blue2red(numCat);
    for (i in 1:numCat) {
      altRamp[i] = colRamp[(numCat+1)-i];
    }
    colRamp <- altRamp;
  }
  
  #-----------------------------+
  # DRAW FOR ALL CHROMOSOMES    |
  #-----------------------------+
  # The following is for genome with 10 chromosomes such as maize
  for (chrom in 1:10) {
    
    selChrom <- chrom;
    #yOffset <- selChrom*100;
    # The following will start at zero
    yOffset <- selChrom*100 - 100;
    cat("Processing chrom",selChrom,"\n");
    
    #-----------------------------+
    # SUBSET BY CHROMOSOME        |
    #-----------------------------+
    # Subset by genomic seqname
    distnDataSub <- subset(genomeGFF, genomeGFF$seqname==selChrom);
    numBin <- length(distnDataSub[[1]]);
    
    # INITIALIZE MATRICES
    xVal <- matrix ( boxWidth, nrow = numBin, ncol = 1);
    yVal <- matrix ( yLoc, nrow = numBin, ncol = 1);
    colVal <- matrix (1, nrow = numBin, ncol = 1);

    # GET THE COLOR VALUE
    for (i in 1:numBin) {
      startPos = i*boxWidth + 1;
      endPos = startPos + boxWidth + 1;
      xVal[i] = startPos + ( (endPos - startPos)/2 );
      
      #-----------------------------+
      # CHANGE THE FOLLOWING FOR    |
      # GINA PUTTING THESE IN       |
      # DIFFERENT COLUMNS           |
      #-----------------------------+
      colVal[i] = colRamp[distnDataSub[i,binCol]];
    }
  
  
    # DRAW THE CHROMOSOME LABEL
    text(-20, selChrom*100+25-100, selChrom);

    lastStart <- 0;

    #////////////////////////////////////////////
    # Modifiy last start to align the centromers
    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # This will need to be offset wrt chrom 1 and the bin that
    # contains the centromere multiplied by the box width
    # Data from email from J. Stein
    #
    # 1   132.9 Mb - 134.3 Mb    -- 0
    # 2   89.5 Mb - 90.9 Mb      -- 43
    # 3   94.3 Mb - 95.3 Mb      -- 39
    # 4   104.2 Mb 105.0 Mb      -- 104
    # 5   101.3 Mb -108.4Mb      -- 104
    # 6   49.8 Mb - 50.4 Mb      -- 83
    # 7   55.1 Mb - 55.5 Mb      -- 78
    # 8   45.9 Mb - 47.1 Mb      -- 87
    # 9   68.3 Mb - 69.2 Mb      -- 65
    # 10  59.3 Mb - 60.5 Mb      -- 73

    if (alignCent == TRUE) {
      if (chrom == 1) {
        # 133
        lastStart <- 0;
      }
      if (chrom == 2) {
        # 90
        lastStart <- (43 * boxWidth) - boxWidth;
      }
      if (chrom == 3) {
        # 94
        lastStart <- (39 * boxWidth) - boxWidth;
      }
      if (chrom == 4) {
        # 104
        lastStart <- (29 * boxWidth) - boxWidth;
      }
      if (chrom == 5) {
        # 104
        lastStart <- (29 * boxWidth) - boxWidth;
      }
      if (chrom == 6) {
        # 50
        lastStart <- (83 * boxWidth) - boxWidth;
      }
      if (chrom == 7) {
        # 55
        lastStart <- (78 * boxWidth) - boxWidth;
      }
      if (chrom == 8) {
        # 46
        lastStart <- (87 * boxWidth) - boxWidth;
      }
      if (chrom == 9) {
        # 68
        lastStart <- (65 * boxWidth) - boxWidth;
      }
      if (chrom == 10) {
        # 60
        lastStart <- (73 * boxWidth) - boxWidth;
      }
    }
    #////////////////////////////////////////////
    # END start modification
    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

    # DRAW COLORED BOX FOR EACH BIN IN THE CHROMOSOME
    for (i in 1:numBin) {
      xStartPos = lastStart + boxWidth;
      xEndPos = xStartPos + boxWidth;
      yStart = yOffset;
      yEnd = yOffset + boxHeight;
      rect(xStartPos, yStart, xEndPos, yEnd, col=colVal[i], border=colVal[i] );
      lastStart = xStartPos;
    }
    
    
    
  } # End of for each chrom
  
  
  if (directToFile == TRUE) {
    dev.off();
  } else {
    par(op);
  }


  #-----------------------------+
  # DRAW OUTPUT FILE            |
  #-----------------------------+
  # The following draws in landscape view, can switch for portrait view
  if (copyPlot == TRUE ) {
    dev.copy (png, filename=outPlotImage,
  #            height=3000, width=2250, bg="white" );
  # The following is for the high res            
  #            height=2250, width=3000, bg="white" );
              height=400, width=400, bg="white" );
    dev.off();
  }
  
  
} # End of genomeHeatmap function



#-----------------------------------------------------------+
# HISTORY
#-----------------------------------------------------------+
# 07/13/2009
# - Added datcat function
# - Added gff2r function
# 07/14/2009
# - Added genomeHeatmap function
#   currently best used just for drawing maize data
#
# 08/03/2009
# - Working to get image file output working for the
#   datcat in a way that will work well with TIFF files for
#   the maize genome paper.

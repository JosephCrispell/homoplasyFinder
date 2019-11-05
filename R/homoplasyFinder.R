## Tutorials
#https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/
#http://r-pkgs.had.co.nz/description.html
#https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html

## Packages to install
#install.packages("devtools")
#install.packages("digest")
#devtools::install_github("klutometis/roxygen")

## Packages to load
#library("devtools")
#library("roxygen2")

## Creating package
#packageDirectory <- "/home/josephcrispell/Desktop/Research/Homoplasy/homoplasyFinder/"
#create(packageDirectory)

## Documenting changes
#setwd(packageDirectory)
#document()

## Install
#setwd("..")
#install("homoplasyFinder")

# NOTE: Everytime the JAR file is changed you need to:
#   - Replace JAR in R package directory /inst/jar/
#   - If R open, clear workspace, terminal and close without saving
#   - Re-install R package as above

## Attaching data
# For each object of interest - save alone in RData file
# Store RData files in /data/
# Add tag into function header:
#   #' An object of class "DNAbin" containing nucleotide sequence alignment
#   #'
#   #' @name tree
#   #' @docType data
#   #' @author Joseph Crispell \email{crispelljoseph@@gmail.com}
#   #' @keywords data

## Internal functions (not available to user)
# remove @export
# include @keywords internal for internal functions

## Calling package functions
# Always use packageName::function when calling external package functions
# Add into Imports: in DESCRIPTION (separate with ", ")

## Multithreading
# added @export in for functions that are used during multithreading, these functions need to be available in the environment for this line
#       parallel::clusterExport(cl=clusterOfThreads, list("getAllelesAtPosition", "getIsolatesWithNsAtPosition", "areSetsOfIsolatesTheSame"))
#       which is necessary when multithreading is being conducted on a mac

## Usage example
# library(homoplasyFinder)
# path <- "/home/josephcrispell/Desktop/Research/Homoplasy/"
# treeFile <- paste0(path, "example-AFTER_10-08-18.tree")
# fastaFile <- paste0(path, "example_10-08-18.fasta")
# inconsistentPositions <- runHomoplasyFinderInJava(treeFile, fastaFile, path, verbose=TRUE)
# tree <- readAnnotatedTree(path)
# plotAnnotatedTree(tree, inconsistentPositions, fastaFile)

#' Run HomoplasyFinder using Java jar file. Use if you can't get rJava to install properly :-)
#'
#' This function runs HomoplasyFinder as a command line jar tool to identify positions that are potentially homoplasious. Output(s) will appear in current working directory.
#' @param treeFile The full path to the Newick formatted phylogenetic tree file
#' @param fastaFile The full path to a FASTA formatted nucleotide sequence alignment. Defaults to NULL
#' @param presenceAbsenceFile The full path to a CSV table file reporting the presence/absence of INDELs. Defaults to NULL
#' @param createFasta Flag to tell HomoplasyFinder whether to create FASTA file without any inconsistent positions. Defaults to TRUE
#' @param createAnnotatedNewickTree Flag to tell HomoplasyFinder whether to create annotated Newick formatted phylogenetic tree file. Defaults to TRUE
#' @param includeConsistentSitesInReport Flag to tell HomoplasyFinder whether to include information about the consistent sites in the report. Defaults to TRUE
#' @param verbose Flag to tell HomoplasyFinder to request detailed information. Defaults to true
#' @param multithread Flag to tell HomoplasyFinder to use multiple threads. Defaults to false
#' @keywords homoplasyFinder java
#' @export
#' @examples
#' # Find the FASTA and tree files attached to package
#' fastaFile <- system.file("extdata", "example.fasta", package = "homoplasyFinder")
#' treeFile <- system.file("extdata", "example.tree", package = "homoplasyFinder")
#' 
#' # Run the HomoplasyFinder jar tool
#' runHomosyFinderJarTool(treeFile, fastaFile)
#' 
#' # Get the current date
#' date <- format(Sys.Date(), "%d-%m-%y")
#' 
#' # Get the current working directory
#' workingDirectory <- paste0(getwd(), "/")
#' 
#' # Read in the output table
#' resultsFile <- paste0(workingDirectory, "consistencyIndexReport_", date, ".txt")
#' results <- read.table(resultsFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)
#' 
#' # Note the inconsistent positions
#' inconsistentPositions <- results[results$ConsistencyIndex < 1, "Position"]
#' 
#' # Read in the annotated tree
#' tree <- readAnnotatedTree(workingDirectory)
#' 
#' # Plot the annotated tree
#' plotAnnotatedTree(tree, inconsistentPositions, fastaFile)
runHomoplasyFinderJarTool <- function(treeFile, fastaFile=NULL, presenceAbsenceFile=NULL,  
                                      createFasta=TRUE,
                                      createAnnotatedNewickTree=TRUE,
                                      includeConsistentSitesInReport=FALSE,
                                      verbose=TRUE, multithread=FALSE){
  
  # Note the location of the HomoplasyFinder jar file
  homoplasyFinderJar <- system.file("java", "HomoplasyFinder.jar", package = "homoplasyFinder")
  
  # Build the command to be used
  command <- paste0("java -jar ", homoplasyFinderJar, " --tree ", treeFile)
  if(is.null(fastaFile) == FALSE){
    command <- paste0(command, " --fasta ", fastaFile)
  }
  if(is.null(presenceAbsenceFile) == FALSE){
    command <- paste0(command, " --presenceAbsence ", presenceAbsenceFile)
  }
  if(verbose){
    command <- paste0(command, " --verbose")
  }
  if(createFasta){
    command <- paste0(command, " --createFasta")
  }
  if(createAnnotatedNewickTree){
    command <- paste0(command, " --createAnnotatedTree")
  }
  if(includeConsistentSitesInReport){
    command <- paste0(command, " --includeConsistent")
  }
  if(multithread){
    command <- paste0(command, " --multithread")
  }
  
  # Run the command
  system(command, wait=TRUE)
}

#' Run HomoplasyFinder using Java code
#'
#' This function runs HomoplasyFinder (coded in Java) to identify positions that are potentially homoplasious
#' @param path The full path to the directory where the output files will be created
#' @param treeFile The full path to the Newick formatted phylogenetic tree file
#' @param fastaFile The full path to a FASTA formatted nucleotide sequence alignment. Defaults to "Not provided"
#' @param presenceAbsenceFile The full path to a CSV table file reporting the presence/absence of INDELs. Defaults to "Not provided"
#' @param createFasta Flag to tell HomoplasyFinder whether to create FASTA file without any inconsistent positions. Defaults to TRUE
#' @param createReport Flag to tell HomoplasyFinder whether to create report file detailing the inconsistent sites identified. Defaults to TRUE
#' @param createAnnotatedNewickTree Flag to tell HomoplasyFinder whether to create annotated Newick formatted phylogenetic tree file. Defaults to TRUE
#' @param includeConsistentSitesInReport Flag to tell HomoplasyFinder whether to include information about the consistent sites in the report. Defaults to TRUE
#' @param verbose Flag to parse to HomoplasyFinder to request detailed information. Defaults to true
#' @param multithread Flag to tell HomoplasyFinder to use multiple threads. Defaults to false
#' @keywords homoplasyFinder java
#' @export
#' @examples
#' # Find the FASTA and tree files attached to package
#' fastaFile <- system.file("extdata", "example.fasta", package = "homoplasyFinder")
#' treeFile <- system.file("extdata", "example.tree", package = "homoplasyFinder")
#' 
#' # Get the current working directory
#' workingDirectory <- paste0(getwd(), "/")
#' 
#' # Run the HomoplasyFinder java code
#' inconsistentPositions <- runHomoplasyFinderInJava(treeFile, fastaFile, workingDirectory)
#' 
#' # Get the current date
#' date <- format(Sys.Date(), "%d-%m-%y")
#'  
#' # Read in the output table
#' resultsFile <- paste0(workingDirectory, "consistencyIndexReport_", date, ".txt")
#' results <- read.table(resultsFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)
#' 
#' # Read in the annotated tree
#' tree <- readAnnotatedTree(workingDirectory)
#' 
#' # Plot the annotated tree
#' plotAnnotatedTree(tree, inconsistentPositions, fastaFile)
runHomoplasyFinderInJava <- function(treeFile, fastaFile="Not provided", presenceAbsenceFile="Not provided", path, 
                                     createFasta=TRUE,
                                     createReport=TRUE,
                                     createAnnotatedNewickTree=TRUE,
                                     includeConsistentSitesInReport=FALSE,
                                     verbose=TRUE, multithread=FALSE){

  # Add the java class path (path to the jar file in R package)
  rJava::.jaddClassPath('inst/java/HomoplasyFinder.jar')
  
  # Initialise the Java class - include package in path to class
  javaHomoplasyFinderClass <- rJava::.jnew("homoplasyFinder/HomoplasyFinder")
  
  # Run HomoplasyFinder
  result <- rJava::.jcall(javaHomoplasyFinderClass, # Loaded Java class
                          method="runHomoplasyFinderFromR", # Method in Java class to be called
                          returnSig="[I", # The return type for the method
                          # Method arguments follow
                          treeFile, fastaFile, presenceAbsenceFile, path, createFasta, createReport,
                          createAnnotatedNewickTree, includeConsistentSitesInReport,
                          verbose, multithread)
  return(result)
}

#' Read in the annotated Newick tree
#'
#' This function runs reads in the annotated Newick formatted phylogenetic tree file produced by HomoplasyFinder. Node labels of tree report where changes occured for inconsistent positions ("-" separated)
#' @param path The full path to the directory where the output files were created (must match what was used for \code{runHomoplasyFinderInJava()})
#' @param date Date \code{runHomoplasyFinderInJava()} was ran. Format: "\%d-\%m-\%y". Defaults to today.
#' @keywords tree newick annotated
#' @export
readAnnotatedTree <- function(path, date=format(Sys.Date(), "%d-%m-%y")){
  return(ape::read.tree(paste0(path, "annotatedNewickTree_", date, ".tree")))
}

#' Plotting HomoplasyFinder results
#'
#' This function plots an annotated phylogeny with the inconsistent sites identified by HomoplasyFinder. Alignment colours: A=red, C=blue, G=cyan, and T=orange. NOTE: Works best with addTextLabels() package installed.
#' @param tree An object of class "phylo" produced by \code{readAnnotatedTree()}
#' @param inconsistentPositions An integer array of the inconsistent positions identified by HomoplasyFinder. Produced by \code{runHomoplasyFinderInJava()}
#' @param fastaFile The full path to the FASTA formatted nucleotide sequence alignment
#' @param addScale A boolean value to determine whether a scale is added to plot. Defaults to FALSE
#' @param alignmentCex A multiplier value to change the size of the space alotted to plotting the homoplasy alignment. Defaults to 2
#' @param addNodeLabels A multiplier value to change the size node labels. Defaults to 1
#' @param nodeLabelCex A multiplier value to change the size of the internal node labels added to phylogeny. Defaults to 1
#' @param alignmentPositionCex A multiplier value to change the size of the position labels added above alignment. Defaults to 0.5
#' @param addSeparatorLines A boolean variable to determine whether white separator lines are plotted on alignment. Defaults to true
#' @param actualPositions An integer array of length=\code{length(inconsistentPositions)} that reports the actual positions on the genome of inconsistent positions. Useful when \code{homoplasyFinder} was run on concatenated sites. Defaults to \code{NULL} meaning it is ignored. 
#' @keywords plot tree annotated
#' @export
#' @examples
#' # Find the FASTA and tree files attached to package
#' fastaFile <- system.file("extdata", "example.fasta", package = "homoplasyFinder")
#' treeFile <- system.file("extdata", "example.tree", package = "homoplasyFinder")
#' 
#' # Get the current working directory
#' workingDirectory <- paste0(getwd(), "/")
#' 
#' # Run the HomoplasyFinder java code
#' inconsistentPositions <- runHomoplasyFinderInJava(treeFile, fastaFile, workingDirectory)
#' 
#' # Get the current date
#' date <- format(Sys.Date(), "%d-%m-%y")
#'  
#' # Read in the output table
#' resultsFile <- paste0(workingDirectory, "consistencyIndexReport_", date, ".txt")
#' results <- read.table(resultsFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)
#' 
#' # Read in the annotated tree
#' tree <- readAnnotatedTree(workingDirectory)
#' 
#' # Plot the annotated tree
#' plotAnnotatedTree(tree, inconsistentPositions, fastaFile)
plotAnnotatedTree <- function(tree, inconsistentPositions, fastaFile, addScale=FALSE, alignmentCex=2, addNodeLabels=TRUE,
                              nodeLabelCex=1, alignmentPositionCex=0.5, addSeparatorLines=TRUE, actualPositions=NULL){
  
  # Get the current plotting margins - so that we can revert to these once finished
  marginSettings <- par("mar")
  
  # Set the plotting margins
  par(mar=c(0,0,1,0.5))

  # Plot the phylogeny
  # Note plotting invisible tip labels - these provide space for plotting alignment
  ape::plot.phylo(tree, show.tip.label=TRUE, type="phylogram", align.tip.label=TRUE, 
                  tip.color=rgb(0,0,0,0), main="Homoplasies Identified", cex=1*alignmentCex)
  
  # Add a scale bar if requested
  if(addScale){
    addScaleBar()
  }
  
  # Add internal node labels indicating the inconsistent positions associated with each internal node
  if(addNodeLabels){
    addInternalNodeLabels(tree, cex=nodeLabelCex)
  }

  # Add an alignment detailing the nucleotides for each inconsistent position
  addInconsistentPositionAlignment(tree, fastaFile, inconsistentPositions, alignmentPositionCex, addSeparatorLines, actualPositions)
  
  # Revert the previous margin settings
  par("mar"= marginSettings)
}

#' Get the X and Y coordinates of each tip in the current plot.phylo plot
#'
#' Function used by \code{plotAnnotatedTree()}
#' @param tipLabels An object of class "vector" containing the labels assigned to each tip in the plotted tree. Get using (\code{tree$tip.label})
#' @keywords internal
#' @return Returns an object of class "list" containing the X and Y coordinates of each tip in the currently plotted phylogenetic tree
getTipCoordinates <- function(tipLabels){
  lastPP <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  tips <- list()
  for(i in seq_along(tipLabels)){
    tips[[as.character(tipLabels[i])]] <- c(lastPP$xx[i], lastPP$yy[i])
  }

  return(tips)
}

#' Get the X and Y coordinates of each internal node in the current plot.phylo plot
#'
#' Function used by \code{plotAnnotatedTree()}
#' @param nTips The number of tips in in the plotted tree. Calculate using (\code{length(tree$tip.label)})
#' @keywords internal
#' @return Returns an object of class "matrix" containing the X and Y coordinates of each internal node in the currently plotted phylogenetic tree
getInternalNodeCoordinates <- function(nTips){
  
  # Get all the information from the most recent plot.phylo() call
  lastPP <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  
  # Initialise a matrix to store X and Y coordinates of the internal nodes
  coordinates <- matrix(NA, nrow=length(lastPP$xx) - nTips, ncol=2)
  for(i in (nTips + 1):length(lastPP$xx)){
    coordinates[(i - nTips), 1] <- lastPP$xx[i]
    coordinates[(i - nTips), 2] <- lastPP$yy[i]
  }
  
  return(coordinates)
}

#' Adds labels to plotted tree indicating the internal nodes associated with each inconsistent position
#'
#' Function used by \code{plotAnnotatedTree()}
#' @param tree An object of class "phylo" produced by \code{readAnnotatedTree()}
#' @param cex A multiplier value to change the size of the node labels
#' @keywords internal
addInternalNodeLabels <- function(tree, cex){
  
  # Get the coordinates of the internal nodes
  internalNodeCoords <- getInternalNodeCoordinates(length(tree$tip.label))
  
  # Check if basicPlotteR package installed
  if(is.element("basicPlotteR", installed.packages()[,1])){
    
    # Get the coordinates and labels for each internal node
    xCoords <- c()
    yCoords <- c()
    labels <- c()
    for(i in 1:length(tree$node.label)){
      
      # Skip internal nodes with no label
      if(tree$node.label[i] == ""){
        next
      }
      
      # Get the coordinates of the current internal node
      coords <- internalNodeCoords[i, ]
      
      # Create the label into the inconsistent positions
      positions <- paste(strsplit(tree$node.label[i], split="-")[[1]], collapse=",")

      # Store the coordinates and label
      xCoords[length(xCoords) + 1] <- coords[1]
      yCoords[length(yCoords) + 1] <- coords[2]
      labels[length(labels) + 1] <- positions
    }
    
    # Add the labels
    basicPlotteR::addTextLabels(xCoords, yCoords, labels, col.label="white", col.line="red", 
                                col.background=rgb(0,0,0, 0.75), avoidPoints=FALSE)
    
  }else{
    # Add a label to each internal node
    for(i in 1:length(tree$node.label)){
      
      # Skip internal nodes with no label
      if(tree$node.label[i] == ""){
        next
      }
      
      # Get the coordinates of the current internal node
      coords <- internalNodeCoords[i, ]
      
      # Create the label into the inconsistent positions
      positions <- paste(strsplit(tree$node.label[i], split="-")[[1]], collapse=",")
      
      # Calculate the height and width of the label
      labelHeight <- strheight(positions, cex=cex)
      labelWidth <- strwidth(positions, cex=cex)
      
      # Add a background polygon
      polygon(x=c(coords[1] - (labelWidth * 0.52),
                  coords[1] - (labelWidth * 0.52),
                  coords[1] + (labelWidth * 0.52),
                  coords[1] + (labelWidth * 0.52)),
              y=c(coords[2] - (labelHeight * 0.6),
                  coords[2] + (labelHeight * 0.6),
                  coords[2] + (labelHeight * 0.6),
                  coords[2] - (labelHeight * 0.6)), 
              col=rgb(0,0,0, 0.75),
              border=NA)
      
      # Add each position at the current node
      text(x=coords[1], y=coords[2], labels=positions[1], col=rgb(1,1,1), cex=cex)
    }
  }
}

#' Reads in FASTA file and records which sequences is associated with each tip label
#' 
#' Function used by \code{plotAnnotatedTree()}
#' @param fastaFile The full path to the FASTA formatted nucleotide sequence alignment
#' @keywords internal
#' @return Returns an object of class "list" containing the nucleotide sequence associated with each tip in the plotted phylogenetic tree
getTipSequences <- function(fastaFile){
  
  # Read in the fasta file
  sequences <- ape::read.dna(fastaFile, format="fasta", as.character=TRUE)
  
  # Create a list recording each sequence
  output <- list()
  for(row in rownames(sequences)){
    output[[row]] <- sequences[row, ]
  }
  
  return(output)
}

#' Get the maximum X and Y coordinates used on current plot.phylo plot
#'
#' Function used by \code{plotTreeAndHomoplasySites()}
#' @param tipCoordinates An object of class "list" containing the coordinates of each tip in the current \code{plot.phylo} plot - returned by \code{getTipCoordinates()}
#' @keywords internal
#' @return Returns an object of class "vector" containing the maximum X and Y coordinates found
getMaxCoordinates <- function(tipCoordinates){
  
  max <- c(NA, NA)
  
  for(key in names(tipCoordinates)){
    
    if(is.na(max[1]) == TRUE || tipCoordinates[[key]][1] > max[1]){
      max[1] <- tipCoordinates[[key]][1]
    }
    
    if(is.na(max[2]) == TRUE || tipCoordinates[[key]][2] > max[2]){
      max[2] <- tipCoordinates[[key]][2]
    }
  }
  
  return(max)
}

#' Adds a scale bar to plotted phylogenetic tree if requested
#'
#' Function used by \code{plotAnnotatedTree()}
#' @keywords internal
addScaleBar <- function(){
  
  # Get the axis Limits
  axisLimits <- par("usr")
  
  # Add Scale bar
  xLength <- axisLimits[2] - axisLimits[1]
  yLength <- axisLimits[4] - axisLimits[3]
  points(x=c(axisLimits[1] + 0.05*xLength, axisLimits[1] + 0.15*xLength),
         y=c(axisLimits[3]+0.05*yLength, axisLimits[3]+0.05*yLength), type="l", lwd=3)
  text(x=axisLimits[1] + 0.1*xLength, y=axisLimits[3] +0.07*yLength,
       labels=round(0.1*xLength, digits=1), cex=1, xpd=TRUE)
}

#' Adds nucleotide alignment to right of plotted phylogenetic tree - nucleotides for inconsistent positions
#'
#' Function used by \code{plotAnnotatedTree()}
#' @param tree An object of class "phylo" produced by \code{readAnnotatedTree()}
#' @param fastaFile The full path to the FASTA formatted nucleotide sequence alignment
#' @param inconsistentPositions An integer array of the inconsistent positions identified by HomoplasyFinder. Produced by \code{runHomoplasyFinderInJava()}
#' @param cex A multiplier value to change the size of position labels
#' @param addSeparatorLines A boolean variable to determine whether white separator lines are plotted on alignment. Defaults to true
#' @param actualPositions An integer array of length=\code{length(inconsistentPositions)} that reports the actual positions on the genome of inconsistent positions. Useful when \code{homoplasyFinder} was run on concatenated sites. Defaults to \code{NULL} meaning it is ignored. 
#' @keywords internal
addInconsistentPositionAlignment <- function(tree, fastaFile, inconsistentPositions, cex, addSeparatorLines, actualPositions){
  
  # Get the axis Limits
  axisLimits <- par("usr")
  
  # Read in the nucleotide sequences
  tipSequences <- getTipSequences(fastaFile)
  
  # Note the locations of the tips
  tipCoordinates <- getTipCoordinates(tree$tip.label)
  maxCoords <- getMaxCoordinates(tipCoordinates)
  
  # Calculate width of space for nucleotide
  charWidth <- (axisLimits[2] - maxCoords[1]) / length(inconsistentPositions)
  
  # Set nucleotide colours
  nucleotideColours <- list("a"="red", "c"="blue", "g"="cyan", "t"="orange")
  
  # Plot FASTA alignment beside tree
  for(positionIndex in seq_along(inconsistentPositions)){
    
    # Examine each of the tips
    for(tipIndex in seq_along(tree$tip.label)){
      
      # Get XY coordinates for the current tip
      xy <- tipCoordinates[[tree$tip.label[tipIndex]]]
      
      # Get the sequence for the current tip
      sequence <- tipSequences[[tree$tip.label[tipIndex]]]
      
      # Plot a polygon for the current tip's nucleotide at the current homoplasy's position
      rect(xleft=maxCoords[1] + ((positionIndex-1) * charWidth),
           ybottom=xy[2]-0.5,
           xright=maxCoords[1] + (positionIndex * charWidth),
           ytop=xy[2] + 0.5,
           col=nucleotideColours[[sequence[inconsistentPositions[positionIndex]]]],
           border=rgb(0,0,0,0), xpd=TRUE)
    }
  }
  
  # Add separator lines
  if(addSeparatorLines){
    for(positionIndex in seq_along(inconsistentPositions)){
      points(x=c(maxCoords[1] + ((positionIndex-1) * charWidth),
                 maxCoords[1] + ((positionIndex-1) * charWidth)),
             y=c(axisLimits[3], axisLimits[4]),
             type="l", col="white", xpd=TRUE)
    }
  }
  
  # Note the positions of each homoplasy
  for(positionIndex in seq_along(inconsistentPositions)){
    
    # Check if actual position available - mapped back to genome
    position <- inconsistentPositions[positionIndex]
    if(is.null(actualPositions) == FALSE){
      position <- actualPositions[positionIndex]
    }
    
    # Add the current position
    text(x=maxCoords[1] + ((positionIndex-0.25) * charWidth),
         y=maxCoords[2],
         labels=position,
         cex=cex, srt=90, xpd=TRUE, pos=3)
  }
}


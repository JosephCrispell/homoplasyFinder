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
#packageDirectory <- file.path("C:", "Users", "Joseph Crisp", "Desktop", "UbuntuSharedFolder", "Homoplasy", "homoplasyFinder")
#create(packageDirectory)
#setwd(packageDirectory)

## Documenting changes
#setwd(packageDirectory)
#document()

## Install
#setwd("..")
#install("homoplasyFinder")

## General
# remove @export and include @keywords internal for internal functions
# save single objects into RData files of the same name into data/ folder

#' An object of class "DNAbin" containing nucleotide sequence alignment
#'
#' @name tree
#' @docType data
#' @author Joseph Crispell \email{crispelljoseph@@gmail.com}
#' @keywords data
NULL

#' An object of class "phylo" detailing a phylogenetic tree
#'
#' @name sequences
#' @docType data
#' @author Joseph Crispell \email{crispelljoseph@@gmail.com}
#' @keywords data
NULL

#' Plotting HomoplasyFinder results
#'
#' This function plots a phylogeny alongside the sites where homoplasies were identified
#' @param tree An object of class "phylo"
#' @param results An object of class "data.frame" \code{homoplasyFinder()} or \code{runHomoplasyFinderJavaTool()}#
#' @param scale A boolean value to determine whether a scale is added to plot. Defaults to TRUE
#' @keywords plot
#' @export
#' @examples
#' # Load the tree and sequence data
#' data("tree")
#' data("sequences")
#'
#' # Run homoplasy finder
#' results <- homoplasyFinder(tree, sequences)
#'
#' # Plot the results
#' plotTreeAndHomoplasySites(tree, results)
plotTreeAndHomoplasySites <- function(tree, results, scale=TRUE){

  # Load the ape library explicitely - needed to get tip coordinates
  library(ape)
  
  # Set the plotting margins
  par(mar=c(0,0,1,0.5))

  # Plot the phylogeny
  ape::plot.phylo(tree, show.tip.label=TRUE, type="phylogram", align.tip.label=TRUE, tip.color=rgb(0,0,0,0),
                  main="Homoplasies Identified")

  # Get the axisLimits
  axisLimits <- par("usr")

  # Add Scale bar
  if(scale == TRUE){
    xLength <- axisLimits[2] - axisLimits[1]
    yLength <- axisLimits[4] - axisLimits[3]
    points(x=c(axisLimits[1] + 0.1*xLength, axisLimits[1] + 0.2*xLength),
           y=c(axisLimits[3]+0.2*yLength, axisLimits[3]+0.2*yLength), type="l", lwd=3)
    text(x=axisLimits[1] + 0.15*xLength, y=axisLimits[3] +0.18*yLength,
         labels=round(0.1*xLength, digits=3), cex=1, xpd=TRUE)
  }

  # Note the locations of the tips
  tipCoordinates <- getTipCoordinates(tree$tip.label)
  maxCoords <- maxCoordinates(tipCoordinates)

  # Calculate width of space for nucleotide
  charWidth <- (axisLimits[2] - maxCoords[1]) / nrow(results)

  # Set nucleotide colours
  nucleotideColours <- list("A"="red", "C"="blue", "G"="cyan", "T"="orange")

  # Plot FASTA alignment beside tree
  for(homoplasyIndex in 1:nrow(results)){

    # Get an array of the nucleotides associated with the current position
    nucleotides = strsplit(results[homoplasyIndex, "Alleles"], split=",")[[1]]

    # Get an array of concatenated ID for each nucleotide
    concatenatedIsolates = strsplit(results[homoplasyIndex, "IsolatesForAlleles"], split=",")[[1]]

    # Examine each nucleotides
    for(nucleotideIndex in 1:length(nucleotides)){

      # Examine each isolate with the current nucleotide
      for(id in strsplit(concatenatedIsolates[nucleotideIndex], split=":")[[1]]){

        # Get XY coordinates for tip
        xy <- tipCoordinates[[id]]

        # Plot a polygon for the current tip's nucleotide at the current homoplasy's position
        polygon(x=c(maxCoords[1] + ((homoplasyIndex-1) * charWidth),
                    maxCoords[1] + ((homoplasyIndex-1) * charWidth),
                    maxCoords[1] + (homoplasyIndex * charWidth),
                    maxCoords[1] + (homoplasyIndex * charWidth)),
                y=c(xy[2], xy[2] + 1, xy[2] + 1, xy[2]),
                col=nucleotideColours[[nucleotides[nucleotideIndex]]],
                border=rgb(0,0,0,0), xpd=TRUE)
      }
    }
  }

  # Add separator lines
  for(homoplasyIndex in 1:nrow(results)){
    points(x=c(maxCoords[1] + ((homoplasyIndex-1) * charWidth),
               maxCoords[1] + ((homoplasyIndex-1) * charWidth)),
           y=c(axisLimits[3], axisLimits[4]),
           type="l", col="white")
  }

  # Note the positions of each homoplasy
  for(row in 1:nrow(results)){
    text(x=maxCoords[1] + ((row-1) * charWidth) + (0.5 * charWidth),
         y=maxCoords[2] + 0.02*yLength,
         labels=results[row, "Position"], cex=0.5, srt=90, xpd=TRUE)
  }

  # Reset plotting margins
  par(mar=c(5.1, 4.1, 4.1, 2.1))
}

#' Get the maximum X and Y coordinates used on current plot.phylo plot
#'
#' Function used by \code{plotTreeAndHomoplasySites()}
#' @param tipCoordinates An object of class "list" containing the coordinates of each tip in the current \code{plot.phylo} plot - returned by \code{getTipCoordinates()}
#' @keywords internal
#' @return Returns an object of class "vector" containing the maximum X and Y coordinates found
maxCoordinates <- function(tipCoordinates){

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

#' Get the X and Y coordinates of each tip in the current plot.phylo plot
#'
#' Function used by \code{plotTreeAndHomoplasySites()}
#' @param tipLabels An object of class "vector" containing the tip labels associated with a phylogenetic tree (\code{tree$tip.label})
#' @keywords internal
#' @return Returns an object of class "list" containing the X and Y coordinates of each tip in the currently plotted phylogenetic tree
getTipCoordinates <- function(tipLabels){
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  tips <- list()
  for(i in 1:length(tipLabels)){
    tips[[as.character(tipLabels[i])]] <- c(lastPP$xx[i], lastPP$yy[i])
  }

  return(tips)
}

#' Run HomoplasyFinder.jar file
#'
#' Function to run HomoplasyFinder.jar file (java version of HomoplasyFinder)  within R - considerably faster than R package
#' @param jarFile A character string specifying the full path to HomoplasyFinder jar file
#' @param fastaFile A character string specifying the full path to FASTA sequences file
#' @param treeFile A character string specifying the full path to HomoplasyFinder jar file
#' @param verbose Flag to pass to HomoplasyFinder.jar to request detailed progress information. Defaults to TRUE
#' @keywords HomoplasyFinder java tool jar
#' @export
#' @examples
#' # Run homoplasyFinder
#' results <- runHomoplasyFinderJavaTool(fullPathToJarFile, fullPathToFastaFile, fullPathToTreeFile)
#' 
#' # Plot the results
#' plotTreeAndHomoplasySites(read.tree(fullPathToTreeFile), results)
#' @return Returns a data.frame detailing the Positions, Alleles and Isolates associated with homoplasies identified
runHomoplasyFinderJavaTool <- function(jarFile, fastaFile, treeFile, verbose=TRUE){

  # Get the current date
  date <- format(Sys.Date(), "%d-%m-%y")

  # Check for spaces in file paths - these can cause problems. If found surround path in quotes
  jarFile <- checkForSpaces(jarFile)
  fastaFile <- checkForSpaces(fastaFile)
  treeFile <- checkForSpaces(treeFile)

  # Check whether verbose requested
  verboseFlag = 0
  if(verbose){
    verboseFlag = 1
  }

  # Run the HomoplasyFinder jar file
  system(paste("java -jar", jarFile, verboseFlag, fastaFile, treeFile, sep=" "),
         ignore.stdout=FALSE)

  # Retrieve the output from HomoplasyFinder
  file <- paste("homoplasyReport_", date, ".txt", sep="")
  results <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                        check.names=FALSE)

  return(results)
}

#' Check for spaces in file path - if found this function returns path surrounded by double quotes
#'
#' Function used by \code{runHomoplasyFinderJavaTool()}
#' @param filePath A path to a file
#' @keywords internal
#' @return Returns A character string of the file path surrounded by quotations
checkForSpaces <- function(filePath){

  if(grepl(filePath, pattern=" ") == TRUE){
    filePath <- paste("\"", filePath, "\"", sep="")
  }

  return(filePath)
}

#' HomoplasyFinder - Identify homoplasies on phylogeny using its sequence alignment
#'
#' Function to run HomoplasyFinder - R version. Note considerably slower than java version
#' @param tree An object of class "phylo"
#' @param sequencesDNABin An object containing FASTA formatted nucleotide sequences of class "DNAbin"
#' @param verbose Print detailed progress information. Defaults to TRUE
#' @keywords HomoplasyFinder tool R
#' @export
#' @examples
#' # Load the tree and sequence data
#' data("tree")
#' data("sequences")
#'
#' # Run homoplasy finder
#' results <- homoplasyFinder(tree, sequences)
#'
#' # Plot the results
#' plotTreeAndHomoplasySites(tree, results)
#'
#' @return Returns a data.frame detailing the Positions, Alleles and Isolates associated with homoplasies identified
homoplasyFinder <- function(tree, sequencesDNABin, verbose=TRUE){

  #### Parse tree and FASTA

  # Note the isolates at the tips of each node in the phylogenetic tree
  nodes <- getNodes(tree)

  # Convert the DNAbin sequence alignment to alignment class
  sequences <- ape::as.alignment(sequencesDNABin)

  # Record the alleles present in FASTA sequences
  alleles <- recordAllelesInPopulation(sequences, verbose)

  # Remove the constant sites (only one nucleotide present)
  removeConstantSites(alleles, nchar(sequences$seq[1]), verbose)

  #### Assign alleles to nodes in phylogeny

  # Assign alelles to nodes in the phylogeny - where possible
  notAssigned <- assignAllelesToNodes(nodes, alleles, sequences$nam, verbose)

  #### Report the homopalsies identified

  results <- reportHomoplasiesIdentified(notAssigned, alleles)

  return(results)
}

#' Remove constant sites
#'
#' Function used by \code{homoplasyFinder()} to remove constant sites (only one nucleotide observed across sequences)
#' @param alleles An object of class "list" recording the isolates associated with each allele produced by \code{recordAllelesInPopulation()}
#' @param nSites The number of positions in the nucleotide sequence alignment
#' @keywords internal
removeConstantSites <- function(alleles, nSites, verbose=TRUE){

  nSitesRemoved <- 0

  nucleotides <- c('a', 'c', 'g', 't')
  for(position in 1:nSites){

    # Count the number of alleles at the current position
    count <- 0
    for(nucleotide in nucleotides){
      if(is.null(alleles[[paste(position, nucleotide, sep=":")]]) == FALSE){
        count <- count + 1
      }
    }
    # If only 1 allele present - remove it as it is a constant site (same across all isolates)
    if(count == 1){
      nSitesRemoved <- nSitesRemoved + 1
      alleles[[paste(position, 'n', sep=":")]] <- NULL
      for(nucleotide in nucleotides){
        if(is.null(alleles[[paste(position, nucleotide, sep=":")]]) == FALSE){
          alleles[[paste(position, nucleotide, sep=":")]] <- NULL
        }
      }
    }
  }

  if(verbose){
    cat(paste("Removed", nSitesRemoved, "constant site(s).\n"))
  }
}

#' Summarise the results from homoplasyFinder() into results table
#'
#' Function used by \code{homoplasyFinder()} to store the details of the homoplasy sites identified in dataframe
#' @param notAssigned An object of class "vector" containing the allele IDs (Position:Nucleotide) produced by \code{assignAlellesToNodes()}
#' @param alleles An object of class "list" recording the isolates associated with each allele produced by \code{recordAllelesInPopulation()}
#' @keywords internal
#' @return Returns a data.frame detailing the Positions, Alleles and Isolates associated with homoplasies identified
reportHomoplasiesIdentified <- function(notAssigned, alleles){

  # Get the positions where homoplasies were identified
  positions <- getPositions(notAssigned)

  # Initialise a data frame to store the homoplasy information
  output <- data.frame("Position"=names(positions), "Alleles"=NA, "IsolatesForAlleles"=NA, stringsAsFactors=FALSE)

  # Report each homoplasy
  for(row in 1:length(positions)){

    # Store the position
    output[row, "Position"] <- names(positions)[row]

    # Store the nucleotides associated with the current homoplasy's position
    nucleotides <- positions[[names(positions)[row]]]
    output[row, "Alleles"] <- paste(toupper(nucleotides), collapse=",")

    # Get the isolates associated with each nucleotide observed at the current site
    isolates <- paste(alleles[[paste(names(positions)[row], nucleotides[1], sep=":")]], collapse=":")
    for(i in 2:length(nucleotides)){

      isolates <- paste(isolates, paste(alleles[[paste(names(positions)[row], nucleotides[2], sep=":")]], collapse=":"), sep=",")
    }
    output[row, "IsolatesForAlleles"] <- isolates
  }

  return(output)
}

#' Get the positions from allele IDs
#'
#' Function used by \code{reportHomoplasiesIdentified()} to extract the positions of each allele id (Position:Nucleotide)
#' @param alleles An object of class "list" recording the isolates associated with each allele produced by \code{recordAllelesInPopulation()}
#' @keywords internal
#' @return Returns a vector of the positions associated with the input allele information
getPositions <- function(alleles){

  positions <- list()
  for(allele in alleles){

    parts <- strsplit(allele, split=":")[[1]]
    position <- parts[1]
    if(is.null(positions[[position]]) == FALSE){
      positions[[position]] <- c(positions[[position]], parts[2])
    }else{
      positions[[position]] <- c(parts[2])
    }
  }

  return(positions)
}

#' Assign alleles to nodes in a phylogeny
#'
#' Function used by \code{homoplasyFinder()} to assign alleles (identified by an ID (Position:Nucleotide)) to nodes in a phylogeny
#' @param nodes An object of class "list" recording the tips (isolates) associated with each node in a phylogeny - produced by \code{getNodes()}
#' @param alleles An object of class "list" recording the isolates associated with each allele produced by \code{recordAllelesInPopulation()}
#' @param isolates An object of class "vector" containing "character" strings identifying isolates
#' @param verbose Print detailed progress information. Defaults to TRUE
#' @keywords internal
#' @return Returns an object of class "vector" with the IDs of alleles (Position:Nucleotide) not assigned to nodes on the phylogeny
assignAllelesToNodes <- function(nodes, alleles, isolates, verbose){

  # Get an array of the postions of the alleles
  positions <- getAllelePositions(alleles)

  # Initialise a hashtable to record which alleles were assigned
  nodeAlleles <- list()

  # Initialise a vector to store the assigned alleles
  assigned <- c()

  # Examine each node
  for(i in 1:length(nodes)){

    # Progress
    if(verbose){
      cat(paste("\rAssigning alleles to node", i, "of", length(nodes)))
    }

    # Get the isolates at the tips of the current node
    tips <- nodes[[names(nodes)[i]]]

    # Get the rest of the isolates
    isolatesAbove <- isolates[isolates %in% tips == FALSE]

    # Check and see if these isolates are associated with an allele
    for(position in positions){

      # Get the alleles at the current position and check that some are present
      allelesAtPosition <- getAllelesAtPosition(position, alleles, assigned)
      if(length(allelesAtPosition) == 0){
        next
      }

      # Get the isolates with an "N" at the current allele's position
      isolatesWithN <- getIsolatesWithNsAtPosition(position, alleles)

      # Remove the isolates with an "N" from those associated with the current node
      tipsWithoutNs <- tips[tips %in% isolatesWithN == FALSE]
      isolatesAboveWithoutNs <- isolatesAbove[isolatesAbove %in% isolatesWithN == FALSE]

      # Examine each of the unassigned alleles at the current position
      for(allele in allelesAtPosition){

        # Get the isolates with the current allele
        isolatesWithAllele <- alleles[[allele]]

        # Compare the two sets of alleles - if they match exactly then current allele can be assigned to node
        if(areSetsOfIsolatesTheSame(tipsWithoutNs, isolatesWithAllele) == TRUE){

          nodeAlleles[[names(nodes)[i]]] <- c(nodeAlleles[[names(nodes)[i]]], allele)
          assigned[length(assigned) + 1] <- allele

        }else if(areSetsOfIsolatesTheSame(isolatesAboveWithoutNs, isolatesWithAllele) == TRUE){
          nodeAlleles[[names(nodes)[i]]] <- c(nodeAlleles[[names(nodes)[i]]], allele)
          assigned[length(assigned) + 1] <- allele
        }
      }
    }
  }

  if(verbose){
    cat("\rFinished assigning alleles to nodes.\t\t\t\t\n")
  }

  # Note the alleles that weren't assigned
  notAssigned <- names(alleles)[names(alleles) %in% assigned == FALSE]

  # Remove the N alleles
  notAssigned <- notAssigned[grepl(notAssigned, pattern=":n") == FALSE]

  return(notAssigned)
}

#' Get a vector of the alleles found at a given position that haven't already been assigned to a node
#'
#' Function used by \code{assignAlellesToNodes()}
#' @param position A character string detailing the position of interest
#' @param alleles An object of class "list" recording the isolates associated with each allele produced by \code{recordAllelesInPopulation()}
#' @param assigned A vector of strings recording the alleles that have been assigned to nodes
#' @keywords internal
#' @return Returns a vector of unassigned alleles found at the position provided
getAllelesAtPosition <- function(position, alleles, assigned){

  # Initialise a vector to store the unassigned alleles at the current position
  output <- c()

  # Create a vector of nucleotides
  nucleotides <- c('a', 'c', 'g', 't')

  # Find the alleles present
  for(nucleotide in nucleotides){

    # Create allele
    allele <- paste(position, nucleotide, sep=":")

    # Check allele isn't assigned but is present in set of alleles (in population)
    if(allele %in% assigned == FALSE && is.null(alleles[[allele]]) == FALSE){
      output[length(output) + 1] <- allele
    }
  }

  return(output)
}

#' Get a vector of the unique positions across the alleles
#'
#' Function used by \code{assignAlellesToNodes()}
#' @param alleles An object of class "list" recording the isolates associated with each allele produced by \code{recordAllelesInPopulation()}
#' @keywords internal
#' @return Returns a vector of the unique positions across the alleles
getAllelePositions <- function(alleles){

  positions <- c()
  for(allele in names(alleles)){

    position <- strsplit(allele, split=":")[[1]][1]
    if(position %in% positions == FALSE){
      positions[length(positions) + 1] <- position
    }
  }

  return(positions)
}

#' Get a vector of the isolates (associated with sequences) that have an N at the given allele's position
#'
#' Function used by \code{assignAlellesToNodes()}
#' @param position A character string detailing the position of interest
#' @param alleles An object of class "list" recording the isolates associated with each allele produced by \code{recordAllelesInPopulation()}
#' @keywords internal
#' @return Returns a vector of the names of isolates that have an N at the position associated with the allele
getIsolatesWithNsAtPosition <- function(position, alleles){

  # Build the allele with an N
  key <- paste(position, "n", sep=":")

  # Check whether any isolates with an N were noted
  isolates <- c()
  if(is.null(alleles[[key]]) == FALSE){
    isolates <- alleles[[key]]
  }

  return(isolates)
}

#' Compare two vectors of isolate names
#'
#' Function used by \code{assignAlellesToNodes()}
#' @param a An object of class "vector" containing character strings recognising isolates
#' @param b An object of class "vector" containing character strings recognising isolates
#' @keywords internal
#' @return Returns either TRUE or FALSE depending on whether the input vectors are identical or not
areSetsOfIsolatesTheSame <- function(a, b){
  result <- FALSE
  if(length(a) == length(b) && length(intersect(a, b)) == length(a)){
    result <- TRUE
  }
  return(result)
}

#' Record the tips associated with every node in a phylogeny
#'
#' Function used by \code{homoplasyFinder()}
#' @param tree An object of class "phylo"
#' @keywords internal
#' @return Returns an object of class "list" detailing the isolates associated with every node in the input tree
getNodes <- function(tree){
  nodes <- list()

  # Nodes number 1:nTips and then onwards: two methods to calculate total number:
  # - Number of edges + 1
  # - Number of tips plus number internal nodes: length(tree$tip.label) + tree$Nnode
  for(node in 1:(length(tree$tip.label) + tree$Nnode)){
    nodes[[as.character(node)]] <- geiger::tips(tree, node)
  }

  return(nodes)
}

#' Record each allele (Position:Nucleotide) found in the set of sequences
#'
#' Function used by \code{homoplasyFinder()}
#' @param sequences An object of class "alignment" containing FASTA formatted nucleotide sequences
#' @param verbose Print detailed progress information. Defaults to TRUE
#' @keywords internal
#' @return Returns an object of class "list" detailing the isolates associated with every allele present in the sequences
recordAllelesInPopulation <- function(sequences, verbose){

  # Initialise a list to store the alleles found and the sequences they are found in
  alleles <- list()

  # Examine each of the sequences
  for(i in 1:length(sequences$nam)){

    # Progress
    if(verbose){
      cat(paste("\rReading alleles in sequence", i, "of", length(sequences$nam)))
    }

    # Split the current sequence into nucleotides
    nucleotides <- strsplit(sequences$seq[i], split="")[[1]]

    # Examine each nucleotide in the current sequence
    for(pos in 1:length(nucleotides)){

      # Build an allele ID
      id <- paste(pos, nucleotides[pos], sep=":")

      # Check if encountered this allele before
      if(is.null(alleles[[id]]) == FALSE){

        alleles[[id]] <- c(alleles[[id]], sequences$nam[i])
      }else{
        alleles[[id]] <- c(sequences$nam[i])
      }
    }
  }

  if(verbose){
    cat("\rFinished reading alleles from sequences.\t\t\t\t\n")
  }

  return(alleles)
}

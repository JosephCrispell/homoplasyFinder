<img align="left" src="HomoplasyFinder-logo.png">

# HomoplasyFinder 
## Author: Joseph Crispell
## Repository created: 25-04-18
## Licence: GPL-3
## Requires: R (>= v3.3.3) & Java (>= v10.0.1)
<br><br>An R package designed for using the Java tool, HomoplasyFinder.

## Installation
```
install.packages("devtools")
library("devtools")
install_github("JosephCrispell/homoplasyFinder")
library(homoplasyFinder)
```

## Executing
```
# Find the FASTA and tree files attached to package
fastaFile <- system.file("extdata", "example.fasta", package = "homoplasyFinder")
treeFile <- system.file("extdata", "example.tree", package = "homoplasyFinder")

# Run the HomoplasyFinder jar tool
inconsistentPositions <- runHomoplasyFinderInJava(treeFile, fastaFile, workingDirectory)
 
# Get the current date
date <- format(Sys.Date(), "%d-%m-%y")

# Get the current working directory
workingDirectory <- getwd()
 
# Read in the output table
resultsFile <- paste0(workingDirectory, "/consistencyIndexReport_" + date + ".txt")
results <- read.table(resultsFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)
 
# Read in the annotated tree
tree <- readAnnotatedTree(workingDirectory)
 
# Plot the annotated tree
plotAnnotatedTree(tree, inconsistentPositions, fastaFile)
```

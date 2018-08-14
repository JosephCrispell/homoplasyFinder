# HomoplasyFinder
## Author: Joseph Crispell
## Repository created: 25-04-18
## Licence: GPL-3
An R package for identifying and plotting homoplasies on a phylogeny

Package can be directly installed into R using:
```
install.packages("devtools")
library("devtools")
install_github("JosephCrispell/homoplasyFinder")
library(homoplasyFinder)
```

Once installed HomoplasyFinder (a Java tool) can either be run using:
```
path <- "[pathToDirectoryForOutputFilesHere]"
treeFile <- "[pathToNewickTreeFileHere]"
fastaFile <- "[pathToFASTANucleotideAlignmentFileHere]"

# Identify the positions
inconsistentPositions <- runHomoplasyFinderInJava(treeFile, fastaFile, path)

# Get the annotated phylogeny
tree <- readAnnotatedTree(path)

# Plot the annotated phylogeny
plotAnnotatedTree(tree, inconsistentPositions, fastaFile)
```
The java tool is currently available for download [with](https://github.com/JosephCrispell/Java/raw/master/ExecutableJarFiles/HomoplasyFinder-GUI.jar) and [without](https://github.com/JosephCrispell/Java/raw/master/ExecutableJarFiles/HomoplasyFinder.jar) a graphical user interface.

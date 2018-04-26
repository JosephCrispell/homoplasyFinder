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

Once installed HomoplasyFinder can either be run using only R:
```
# Load the data
data("tree")
data("sequences")

# Run HomoplasyFinder
results <- homoplasyFinder(tree, sequences)

# Plot the results
plotTreeAndHomoplasySites(tree, results)
```
OR the Java version of the tool can be run from within R (THIS IS MUCH FASTER BUT YOU'll NEED [JAVA](https://java.com/en/download/)):
```
# Run HomoplasyFinder
results <- runHomoplasyFinderJavaTool("pathToJarFile", "pathToFastaFile", "pathToTreeFile")

# Load the ape library - for reading tree file
library(ape)

# Plot the results
plotTreeAndHomoplasySites(read.tree("pathToTreeFile"), results)
```
The java tool is currently available for download [with](https://github.com/JosephCrispell/Java/raw/master/ExecutableJarFiles/HomoplasyFinder_v1.jar) and [without](https://github.com/JosephCrispell/Java/raw/master/ExecutableJarFiles/HomoplasyFinder_25-04-18.jar) a graphical user interface.

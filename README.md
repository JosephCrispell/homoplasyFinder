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

## IMPORTANT - rJava is quite a temperamental R package. If you have problems getting my package to work try:
1) Updating your Java jdk to >10.0.1 (http://www.oracle.com/technetwork/java/javase/downloads/jdk10-downloads-4416644.html)
2) Delete Java jdk from previous versions
3) Reconfigure R to use new Java version, by running this in the terminal:
  ```
  R CMD javareconf
  ```
 4) Re-install rJava from source using:
   ```
  install.packages("rJava", type = "source")
  ```
  5) Run the following code in terminal (for OSX):
  ```
  sudo ln -f -s $(/usr/libexec/java_home)/jre/lib/server/libjvm.dylib /usr/local/lib
  ```

Some useful help information:
- https://stackoverflow.com/questions/26948777/how-can-i-make-rjava-use-the-newer-version-of-java-on-osx
- https://www.r-statistics.com/2012/08/how-to-load-the-rjava-package-after-the-error-java_home-cannot-be-determined-from-the-registry/
- http://spartanideas.msu.edu/2015/06/27/the-rjava-nightmare/
- https://stackoverflow.com/questions/30738974/rjava-load-error-in-rstudio-r-after-upgrading-to-osx-yosemite/31039105#31039105

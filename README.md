# HomoplasyFinder
## Author: Joseph Crispell
## Repository created: 25-04-18
## Licence: GPL-3
## Requires: R (>= v3.3.3) & Java (>= v10.0.1)

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
## Java executables
HomoplasyFinder can be run within R, using this homoplasyFinder package, or it can be ran directly as a Java application. The java tool is currently available for download [with](https://github.com/JosephCrispell/Java/raw/master/ExecutableJarFiles/HomoplasyFinder-GUI.jar) and [without](https://github.com/JosephCrispell/homoplasyFinder/raw/master/inst/java/HomoplasyFinder.jar) a graphical user interface. Note HomoplasyFinder.jar (without GUI) is present within the homoplasyFinder R package. If you manage to get the homoplasyFinder() R package installed but you can't get rJava to work then you can run this tool directly. Just search for "HomoplasyFinder.ja" on your computer. Once you've run HomoplasyFinder.jar then you can use the functions available within the R package homoplasyFinder() to visualise the results.

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

## Java heap space error
If you running HomoplasyFinder on a really large dataset and you run into a heapscape error like this:
```
java.lang.OutOfMemoryError: Java heap space
```
You'll need to change the maximum heap size limit using the `-Xmx` argument. In the command line use:
```
java -Xmx4000m -jar HomoplasyFinder.jar
```
to change it to 4GB. If you're in R, you'll need to restart your session and add the following at the top of your script before loading any packages:
```
options(java.parameters = "-Xmx4000m")
```

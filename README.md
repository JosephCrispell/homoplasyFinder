# HomoplasyFinder
## Author: Joseph Crispell
## Repository created: 25-04-18
## Licence: GPL-3
An R package for identifying and plotting homoplasies on a phylogeny

Package can be directly installed into R using:
```
$ library("devtools")
$ install_github("JosephCrispell/homoplasyFinder")
$ library(homoplasyFinder)
```

Once installed HomoplasyFinder can either be run using only R:
```
$ # Load the data
$ data("tree")
$ data("sequences")
$
$ # Run HomoplasyFinder
$ results <- homopalsyFinder(tree, sequences)
$
$ # Plot the results
$ plotTreeAndHomoplasySites(tree, results)
```
OR the Java version of the tool can be run from within R:
```
$ # Run HomoplasyFinder
$ results <- runHomoplasyFinderJavaTool("pathToJarFile", "pathToFastaFile", "pathToTreeFile")
$
$ # Plot the results
$ plotTreeAndHomoplasySites(read.tree("pathToTreeFile"), results)
```
The java tool is currently available for download [with](https://github.com/JosephCrispell/Java/blob/master/ExecutableJarFiles/HomoplasyFinder_v1.jar) and [without](https://github.com/JosephCrispell/Java/blob/master/ExecutableJarFiles/HomoplasyFinder_25-04-18.jar) a graphical user interface.

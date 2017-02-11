library("Rcpp"); 
setwd("/home/hdong/workspace/compute/RVineAD/DVineAD")

compileAttributes(); 

setwd("/home/hdong/workspace/compute/RVineAD")

library("DVineAD")
detach("package:DVineAD", unload = TRUE)
library.dynam.unload("DVineAD", system.file(package = "DVineAD"))

system("R CMD build DVineAD")
system("R CMD INSTALL DVineAD_0.1.tar.gz")

library("DVineAD")
library.dynam()

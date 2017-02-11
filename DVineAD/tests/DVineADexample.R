rm(list = setdiff(ls(), lsf.str()))

library("DVineAD")
library("VineCopula")

d <- 4 
TruncLevel <- 2
dd <- (2*d-TruncLevel-1)*TruncLevel/2
par <- c(1, 2, 3, 0.9, 0.8)
startpar <- runif(dd)

family <- c(rep(3,dd), rep(0, d*(d-1)/2 - dd)) 
parFill <- c(par, rep(0, d*(d-1)/2 - dd))

RVM <- D2RVine(1:d, family, parFill)
simdata <- RVineSim(400, RVM) 

DVineAD_EstablishTape(simdata, par, TruncLevel)
res <- DVineMLE(startpar, trace=1)
print(res$par)

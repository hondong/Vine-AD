DVineMLE <- function(startpar, maxit=100, trace=1) { 
    numpar <- length(startpar)
    lower <- rep(1e-4, numpar) 
    upper <- rep(10, numpar)
    out <- optim(par = startpar,
          fn = DVineAD_LogLik,
          gr = DVineAD_LogLikGrad,
          method = "L-BFGS-B",
          #method = "BFGS",
          lower = lower,
          upper = upper,
          control = list(maxit = maxit, 
                         trace=trace, 
                         factr=1e+8, 
                         parscale=rep(0.05, numpar)))
    out
}

rm(list = setdiff(ls(), lsf.str()))

library("DVineAD")

library("VineCopula")
i=1
res_tab1 <- list()
for (d in c(10, 20, 30, 40, 50, 60)) {
#for (d in c(10, 20, 30)) {
#for (d in c(5)) {
    res_tab1$d[[i]] <- d
    TruncLevel <- 3 
    dd <- (2*d-TruncLevel-1)*(TruncLevel)/2
    order <- 1:d
    family <- c(rep(3,dd), rep(0, d*(d-1)/2 - dd)) 
    #par <- c(runif(dd), rep(0, d*(d-1)/2 - dd))  
    par <- c(3*runif(d-1), 2*runif(d-2), 1*runif(d-3), rep(0, d*(d-1)/2 - dd))  
    parTL <- par[1:dd]  

    startpar <- c((1+runif(dd))*par[1:dd], rep(0, d*(d-1)/2 - dd))
    startparTL <- startpar[1:dd]
    
    RVM <- D2RVine(order, family, par)
    RVMstart <- D2RVine(order, family, startpar)
    simdata <- RVineSim(200, RVM) 
    
    #print("Computing function by VineCopula")
    ptm <- proc.time()
    res_tab1$RVloglik[[i]] <- RVineLogLik(as.matrix(simdata), RVM)
    res_tab1$RVloglik_t[[i]] <- proc.time() - ptm
    
    #print("Computing gradient by VineCopula")
    ptm <- proc.time()
    res_tab1$RVgrad[[i]] <- RVineGrad(as.matrix(simdata), RVM)
    res_tab1$RVgrad_t[[i]] <- proc.time() - ptm
    
    #print("Computing hessian by VineCopula")
    #ptm <- proc.time()
    #res_tab1$RVHess[[i]] <- RVineHessian(as.matrix(simdata), RVM)
    #res_tab1$RVHess_t[[i]] <- proc.time() - ptm

    #if (d==2) #Check the value of computed hessian?
    #{
    #    RVM2 <- D2RVine(order, family, par+0.0001)
    #    grad2 <- RVineGrad(simdata,RVM2)
    #    print(c("numeric approximate hessian=",(res_tab1$RVgrad[[i]]$gradient-grad2$gradient)/0.0001))
    #    print(c("RVine computed hessian=", res_tab1$RVHess[[i]]$hessian))
    #}

    #print("Establishing Tape...")
    ptm <- proc.time()
    DVineAD_EstablishTape(simdata, parTL, TruncLevel)
    res_tab1$AD_est_t[[i]] <-proc.time() - ptm
    
    #print("Computing function by AD")
    ptm <- proc.time()
    res_tab1$ADloglik[[i]] <- DVineAD_LogLik(parTL)
    res_tab1$ADloglik_t[[i]] <- proc.time() - ptm
    
    #print("Computing gradient by AD")
    ptm <- proc.time()
    res_tab1$ADgrad[[i]] <- DVineAD_LogLikGrad(parTL)
    res_tab1$ADgrad_t[[i]] <- proc.time() - ptm

    #print("Computing hessian by AD")
    #ptm <- proc.time()
    #res_tab1$ADhess[[i]] <-  DVineAD_LogLikHess(parTL)
    #res_tab1$ADhess_t[[i]] <- proc.time() - ptm
    #if (d==2) #Check the value of computed hessian?
    #{
    #   print(c("AD computed hessian:",res_tab1$ADhess[[i]]))
    #}
    
    res_tab1$valdiff[[i]] <- abs(res_tab1$RVloglik[[i]]$loglik + res_tab1$ADloglik[[i]])/abs(res_tab1$ADloglik[[i]])
    res_tab1$graddiff[[i]] = sum(abs(res_tab1$RVgrad[[i]]$gradient + DVineVec2RVineVec(d, res_tab1$ADgrad[[i]], TruncLevel)))/sum(abs(res_tab1$ADgrad[[i]]));

    ptm <- proc.time()
    res_tab1$ADMLE[[i]] <- DVineMLE(startparTL, trace=3, maxit=100)
    res_tab1$ADMLE_t[[i]] <- proc.time() - ptm
    print(res_tab1$ADMLE_t[[i]])

    ptm <- proc.time() 
    res_tab1$RVMLE[[i]] <- RVineMLE(simdata, RVMstart, grad = TRUE, maxit=100) 
    res_tab1$RVMLE_t[[i]] <- proc.time() - ptm
    print(res_tab1$RVMLE_t[[i]])


    #print(sprintf("%3.1e %3.1e %3.1e %3.1e %3.1e %3.1e %3.1e", 
    #        res_tab1$RVloglik_t[[i]][["elapsed"]], res_tab1$RVgrad_t[[i]][["elapsed"]], 
    #        res_tab1$AD_est_t[[i]][["elapsed"]], res_tab1$ADloglik_t[[i]][["elapsed"]], res_tab1$ADgrad_t[[i]][["elapsed"]], 
    #        res_tab1$valdiff[[i]], res_tab1$graddiff[[i]]))
    print(sprintf("%3.1e %3.1e %3.1e %3.1e %3.1e %3.1e %3.1e %3.1e %3.1e", 
            res_tab1$RVloglik_t[[i]][["elapsed"]], res_tab1$RVgrad_t[[i]]["elapsed"], res_tab1$RVMLE_t[[i]]["elapsed"],
            res_tab1$AD_est_t[[i]]["elapsed"], res_tab1$ADloglik_t[[i]]["elapsed"], res_tab1$ADgrad_t[[i]]["elapsed"], 
            res_tab1$ADMLE_t[[i]]["elapsed"], res_tab1$valdiff[[i]], res_tab1$graddiff[[i]]))
    i <- i+1
}

## Code for results in table 2
i=1
res_tab2 <- list()
for (d in c(10, 20, 30, 40, 50, 60)) {
#for (d in c(10, 20, 30)) {
    res_tab2$d[[i]] <- d
    TruncLevel <- d-1 
    dd <- (2*d-TruncLevel-1)*(TruncLevel)/2
    order <- 1:d
    family <- c(rep(3,dd), rep(0, d*(d-1)/2 - dd)) 
    par <- c(runif(dd), rep(0, d*(d-1)/2 - dd))  
    parTL <- par[1:dd]  
    
    RVM <- D2RVine(order, family, par)
    simdata <- RVineSim(200, RVM) 
    
    #print("Computing function by VineCopula")
    ptm <- proc.time()
    res_tab2$RVloglik[[i]] <- RVineLogLik(as.matrix(simdata), RVM)
    res_tab2$RVloglik_t[[i]] <- proc.time() - ptm
    
    #print("Computing gradient by VineCopula")
    ptm <- proc.time()
    res_tab2$RVgrad[[i]] <- RVineGrad(as.matrix(simdata), RVM)
    res_tab2$RVgrad_t[[i]] <- proc.time() - ptm
    
    #print("Establishing Tape...")
    ptm <- proc.time()
    DVineAD_EstablishTape(simdata, parTL, TruncLevel)
    res_tab2$AD_est_t[[i]] <-proc.time() - ptm
    
    #print("Computing function by AD")
    ptm <- proc.time()
    res_tab2$ADloglik[[i]] <- DVineAD_LogLik(parTL)
    res_tab2$ADloglik_t[[i]] <- proc.time() - ptm
    
    #print("Computing gradient by AD")
    ptm <- proc.time()
    res_tab2$ADgrad[[i]] <- DVineAD_LogLikGrad(parTL)
    res_tab2$ADgrad_t[[i]] <- proc.time() - ptm

    res_tab2$valdiff[[i]] <- abs(res_tab2$RVloglik[[i]]$loglik + res_tab2$ADloglik[[i]])/abs(res_tab2$ADloglik[[i]])
    res_tab2$graddiff[[i]] = sum(abs(res_tab2$RVgrad[[i]]$gradient + DVineVec2RVineVec(d, res_tab2$ADgrad[[i]], TruncLevel)))/sum(abs(res_tab2$ADgrad[[i]]));

    print(sprintf("%3.1e %3.1e %3.1e %3.1e %3.1e %3.1e %3.1e", 
            res_tab2$RVloglik_t[[i]][["elapsed"]], res_tab2$RVgrad_t[[i]][["elapsed"]], 
            res_tab2$AD_est_t[[i]][["elapsed"]], res_tab2$ADloglik_t[[i]][["elapsed"]], res_tab2$ADgrad_t[[i]][["elapsed"]], 
            res_tab2$valdiff[[i]], res_tab2$graddiff[[i]]))

    i <- i+1
}

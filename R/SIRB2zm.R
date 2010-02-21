SIRB2zm <- function (m, Y, func, y, times, VN, VF, indBeta, aBeta, bBeta, indQ, aQ, bQ, indG, aG, bG, v, S,
tauN_sh, tauN_sc, tauF_sh, tauF_sc, indep, cred) 
    {
    rtol = 1e-06
    atol = 1e-06
    tcrit = NULL
    jacfunc = NULL
    verbose = FALSE
    hmin = 0
    hmax = Inf
    initpar = NULL
    ModelInit = NULL
    rho = environment(func)

    modY <- numeric()
    modY[seq(1,2*nrow(Y), by=2)] <- Y[,1]
    modY[seq(2,2*nrow(Y), by=2)] <- Y[,2]

    n <- nrow(Y)
    two_n <- 2*n
    Y1 <- c(Y[1:n,1],Y[1:n,2])

    s1 <- seq(4,(3*nrow(rbind(c(0,0),Y))-1),by=3)
    s2 <- seq(5,(3*nrow(rbind(c(0,0),Y))),by=3)

    index <- numeric()
    index[seq(1,2*length(s1), by=2)] <- s1
    index[seq(2,2*length(s2), by=2)] <- s2
    index2 <- c(s1,s2)

      
    S <- as.numeric(S)
    Saux <- S[2]
    S[2] <- S[3]
    S[3] <- Saux

   r <- .Call("call_sir", y, times, func, rtol, atol, rho,              
         tcrit, jacfunc, ModelInit, as.integer(verbose), hmin, 
         hmax, as.integer(m), as.integer(index),as.integer(index2),
         as.integer(two_n),as.integer(n), as.integer(indep), Y1, modY,      
         VN, VF, as.integer(indBeta), aBeta, bBeta, as.integer(indQ),  
         aQ, bQ, as.integer(indG), aG, bG, v, S, tauN_sh, tauN_sc, 
         tauF_sh, tauF_sc, PACKAGE = "B2Z")


    ur <- unlist(r[1:m])
    T <- length(times)
    logCN <- log(matrix(ur[seq(2,3*T*m,by=3)],T,m))
    logCF <- log(matrix(ur[seq(3,3*T*m,by=3)],T,m))
    Beta <- r[[(m+1)]]
    Q <- r[[(m+2)]]
    G <- r[[(m+3)]]
    tauN <- r[[(m+4)]]
    tauF <- r[[(m+5)]]
    tauNF <- r[[(m+6)]]
    l <- r[[(m+7)]]

    C <- 700 - (max(range(l)))
    weights <- l + C
    weights <- exp(weights)
    weights <- weights/sum(weights)

    draw_index <- sample(1:m, m, replace=TRUE, prob=weights)
   
    Betaout <- Beta[draw_index]
    Qout <- Q[draw_index]
    Gout <- G[draw_index]
    TauNout  <- tauN[draw_index]
    TauNFout <- tauNF[draw_index]    
    TauFout  <- tauF[draw_index]
    logCNout <- t(logCN[,draw_index])[,-1]
    logCFout <- t(logCF[,draw_index])[,-1]

    #Proportion of diferent values
    prop <- length(table(Betaout))/m

    #Computing DIC 
    Dbar <- -2*mean(l[draw_index])
    parms <- c(mean(Betaout),mean(Qout),mean(Gout), mean(TauNout), mean(TauFout), mean(TauNFout) )

    Dthetabar <- -2*(.Call("sir_likelihood", y, times, func, rtol, 
        atol, rho, tcrit, jacfunc, ModelInit, as.integer(verbose),  
        hmin, hmax, parms, as.integer(index),as.integer(index2),
        as.integer(two_n),as.integer(n), Y1, modY,      
        VN, VF,PACKAGE = "B2Z")[1])

    pD = Dbar - Dthetabar
    DIC = pD + Dbar

   #Computing ESS
    ESS = m/(1 + var(weights))

   if(indep)
     {
     r <- list(Beta=Betaout, Q=Qout, G=Gout, tauN=TauNout,
          tauF=TauFout, logCN=logCNout, logCF=logCFout,Y=Y,
          DIC=DIC, pD=pD, Dbar=Dbar, ESS=ESS, prop=prop, indep=indep,
          y0=y,times=times, weights=weights, maxw= max(weights),   
          cred=cred)
     }
   else 
     {
     r <- list(Beta=Betaout, Q=Qout, G=Gout, tauN=TauNout,
          tauNF = TauNFout, tauF=TauFout, logCN=logCNout, 
          logCF=logCFout, Y=Y, DIC=DIC, pD=pD, Dbar=Dbar, ESS=ESS, 
          prop=prop, indep=indep, y0=y,times=times, 
          weights=weights, maxw= max(weights),cred=cred)
     }

   attr(r, "class") <- "sir"
   return(r)

}



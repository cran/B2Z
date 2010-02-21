METROPB2zm <- function (NUpd, burnin, lag, initial, priors, S, v, 
                        tauN_sh, tauN_sc, tauF_sh, tauF_sc, VN, 
                        VF, funct, y0, times, Y, indep, 
                         Sigma.Cand, m, cred,indBeta, aBeta, bBeta,
                        indQ, aQ, bQ, indG, aG, bG){


    y = y0
    T = length(times)
    rtol = 1e-06
    atol = 1e-06
    tcrit = NULL
    jacfunc = NULL
    verbose = FALSE
    dllname = NULL
    hmin = 0
    hmax = Inf
    initpar = NULL
    ModelInit = NULL
    rho = environment(funct)
    Nglobal = NULL
    modY <- numeric()
    modY[seq(1, 2 * nrow(Y), by = 2)] <- Y[, 1]
    modY[seq(2, 2 * nrow(Y), by = 2)] <- Y[, 2]
    s1 <- seq(4, (3 * (nrow(Y)+1) - 1), by = 3)
    s2 <- seq(5, (3 * (nrow(Y)+1)), by = 3)
    index <- numeric()
    index[seq(1, 2 * length(s1), by = 2)] <- s1
    index[seq(2, 2 * length(s2), by = 2)] <- s2
    index2 <- c(s1, s2)
    S <- as.numeric(S)
    Saux <- S[2]
    S[2] <- S[3]
    S[3] <- Saux
    n <- nrow(Y)
    two_n <- 2 * n
    Y1 <- c(Y[1:n, 1], Y[1:n, 2])
    dim <- ifelse(indep,5,6)

    cat("|----------25%----------50%----------75%----------|\n")
    flush.console()
    cat("|")
    flush.console()

    total_it <- NUpd
    next_control_bar = control_bar = 1/50

    if(is.null(Sigma.Cand) || is.null(initial))
      {
      total_it <- total_it + m
   								
      samp.values <- matrix(nrow = m, ncol = 6)
      for(i in 1:3)
        {
        split1 <- unlist(strsplit(priors[i],  "\\(" ))
        name.dist <- split1[1]
        par.dist <- unlist(strsplit(split1,  "\\)" ))[2]
    
        x <- paste("r",name.dist,"(",m,",",par.dist,")", sep="")
        write(x,file="foo")
        samp.values[,i] <- eval(source("foo"))$value
        }    
       unlink("foo")

     if(indep)
       {
       samp.values[,4] <- 1/rgamma(m, tauN_sh, tauN_sc)
       samp.values[,5] <- 1/rgamma(m, tauF_sh, tauF_sc)
       samp.values[,6] <- 0
       }
     else
       {
       samp.values[,4:6] <- t(apply(matrix(rep(v,m),m,1), 1, riwish,     
       S=matrix(S,2,2)))
       }

     res <- .Call("call_logpost_metrop", y, times, funct, rtol, 
        atol, rho, tcrit, jacfunc, ModelInit, as.integer(verbose), 
        hmin, hmax, as.integer(m), as.integer(index), 
        as.integer(index2), as.integer(two_n), as.integer(n), 
        as.integer(indep),  Y1, modY, VN, VF, samp.values[,1], 
        samp.values[,2], samp.values[,3], samp.values[,4], 
        samp.values[,5], samp.values[,6], 
        v, S, tauN_sh, tauN_sc, tauF_sh, tauF_sc,
        as.integer(indBeta), aBeta, bBeta,
        as.integer(indQ), aQ, bQ,
        as.integer(indG), aG, bG, as.integer(total_it), 
        PACKAGE = "B2Z")


      pos <- which.max(res)
      initial <- numeric()
      initial[1] <- log(samp.values[pos,1])
      initial[2] <- log(samp.values[pos,2])
      initial[3] <- log(samp.values[pos,3])

      mat <- matrix(c(samp.values[pos,4], samp.values[pos,6], 
            samp.values[pos,6], samp.values[pos,5]),2,2)

      L <- chol(mat)

      initial[4] <- L[1]
      initial[5] <- L[4]

      if(!indep){initial[6] <- L[3]}

      sv_transformed <- matrix(0, m, dim)

      sv_transformed[,1:3] <- log(samp.values[,1:3])                
      sv_transformed[,4] <- sqrt(samp.values[,4])
      sv_transformed[,5] <- sqrt(samp.values[,5] -
                  samp.values[,6]^2/samp.values[,4]) 

       if(!indep)
         {sv_transformed[,6] <- samp.values[,6]/sqrt(samp.values[,4])}
       
        ranges <- apply(sv_transformed, 2, range)
        lower <- ranges[1,]
        upper <- ranges[2,]

        nlmobj <- nlminb(start=initial, objective=minus.logpost,     
        hessian=FALSE, VN=VN,  VF=VF, Y1=Y1, modY=modY,  funct=funct, 
        y=y, times=times, T=T,  indBeta=indBeta, aBeta=aBeta, 
        bBeta=bBeta, indQ=indQ, aQ=aQ, bQ=bQ, indG=indG, aG=aG, bG=bG, 
        S=S, v=v, tauN_sh=tauN_sh, tauN_sc=tauN_sc, tauF_sh=tauF_sh,    
        tauF_sc=tauF_sc, indep=indep, index=index, index2=index2, 
        two_n=two_n, n=n, lower = lower, upper = upper)


      if(!exists("nlmobj"))
       {
       stop("Sorry, the posterior mode was not found. Please, provide 
       the initial values for Beta, Q, G, tauN, tauNF, tauF.")
       }

      initial <- nlmobj$par
      initial[1:3] <- exp(initial[1:3])
      initial[5] <- initial[5]^2 
      initial[4] <- initial[4]^2

      if(!indep)
        {
        initial[6] <- initial[4]*initial[6]
        initial[5] <- initial[5] + initial[6]^2
        }


      if(is.null(Sigma.Cand))
       {
       sm <- hessian(func = logpost, x = nlmobj$par,
       method="Richardson", method.args=list(eps=1e-5, d=0.01,
       zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2), 
       VN=VN, VF=VF, Y1=Y1, modY=modY, 
       funct=funct, y=y, times=times, T=T, indBeta=indBeta, 
       aBeta=aBeta,  bBeta=bBeta, indQ=indQ, aQ=aQ, bQ=bQ, indG=indG, 
       aG=aG, bG=bG, S=S, v=v, tauN_sh=tauN_sh, tauN_sc=tauN_sc, 
       tauF_sh=tauF_sh, tauF_sc=tauF_sc, indep=indep,
       index=index, index2=index2, two_n=two_n, n=n)

       if(length(which(is.nan(sm)==TRUE))> 0){stop("Sorry, it was not 
          possible to estimate the estimated covariance 
          matrix for the multivariate normal proposal distribution. 
          Please, provide such matrix.")}

       try(Sigma.Cand <- -solve(sm),silent=TRUE)
       if(is.null(Sigma.Cand)){stop("Sorry, it was not possible to 
           estimate the estimated covariance matrix for the 
           multivariate normal proposal distribution. Please, provide 
           such matrix.")}

       if(length(which(diag(Sigma.Cand)<0))>0)
          {
          stop("Sorry, it was not possible to estimate the estimated 
          covariance matrix for the multivariate normal proposal 
          distribution. Please, provide such matrix.")
          }

       }

      }

    ev <- eigen(Sigma.Cand, symmetric = TRUE)
    retval <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values))%*%t(ev$vectors)
    
   
    r <- .Call("call_metrop", y, times, funct, rtol, atol, rho,              
         tcrit, jacfunc, ModelInit, as.integer(verbose), hmin, 
         hmax, as.integer(NUpd), as.integer(index), as.integer(index2),
         as.integer(two_n),as.integer(n), as.integer(indep), Y1, modY,      
         VN, VF, as.integer(indBeta), aBeta, bBeta, as.integer(indQ),  
         aQ, bQ, as.integer(indG), aG, bG, v, S, tauN_sh, tauN_sc, 
         tauF_sh, tauF_sc, as.double(retval),as.double(initial),   
         as.integer(dim), as.integer(total_it), PACKAGE = "B2Z")

 
    ur <- unlist(r[1:NUpd])
    T <- length(times)
    logCN <- matrix(ur[seq(1,length(ur),by=2)],T-1,NUpd)
    logCF <- matrix(ur[seq(2,length(ur),by=2)],T-1,NUpd)

    Beta <-  r[[NUpd+1]]
    Q <-  r[[NUpd+2]]
    G <-  r[[NUpd+3]]
    tauN <-  r[[NUpd+4]]
    tauF <-  r[[NUpd+5]]
    loglik <- r[[NUpd+7]]

    cont <- length(table(Beta))
    AR = cont/NUpd
    seq_index <- seq(burnin, NUpd, by = lag)


    if(indep)
      {
      ESS <- c(effectiveSize(Beta),effectiveSize(Q),
      effectiveSize(G),effectiveSize(tauN),effectiveSize(tauF))
      mtauNF=0
      }
    else
      {
      tauNF <- r[[NUpd+6]]   
      ESS <- c(effectiveSize(Beta),effectiveSize(Q),
      effectiveSize(G),effectiveSize(tauN),
      effectiveSize(tauNF),effectiveSize(tauF))
      mtauNF <- mean(tauNF[seq_index])
      }


    #computing DIC
    Dbar <- -2*mean(loglik[seq_index])
    parms <- c(mean(Beta[seq_index]),mean(Q[seq_index]),
    mean(G[seq_index]), mean(tauN[seq_index]), mean(tauF[seq_index]), 
    mtauNF)

    Dthetabar <- -2*(.Call("sir_likelihood", y, times, funct, rtol, 
        atol, rho, tcrit, jacfunc, ModelInit, as.integer(verbose),  
        hmin, hmax, parms, as.integer(index),as.integer(index2),
        as.integer(two_n),as.integer(n), Y1, modY,      
        VN, VF,PACKAGE = "B2Z")[1])

    pD = Dbar - Dthetabar
    DIC = pD + Dbar


    
    if(indep)
     {
     r <- list(Beta=Beta, Q=Q, G=G, tauN=tauN,
          tauNF = NULL, tauF=tauF, logCN=t(logCN), logCF=t(logCF), Y=Y,
          DIC=DIC, pD=pD, Dbar=Dbar, ESS=ESS, AR=AR, indep=indep,
          y0=y0,times=times,  
          cred=cred, initial=initial, Sigma.Cand=Sigma.Cand, NUpd=NUpd, 
          burnin=burnin, lag=lag)
     }
    else 
     {
    
     r <- list(Beta=Beta, Q=Q, G=G, tauN=tauN,
          tauNF = tauNF, tauF=tauF, logCN=t(logCN), logCF=t(logCF), 
          Y=Y,pD=pD, Dbar=Dbar,
          DIC=DIC, ESS=ESS, AR=AR, indep=indep,y0=y0,times=times, 
          cred=cred, initial=initial, Sigma.Cand=Sigma.Cand, NUpd=NUpd, 
          burnin=burnin, lag=lag)
     }


    cat("|\n")

    flush.console()

    attr(r, "class") <- "metrop"

 
    return(r)

    }

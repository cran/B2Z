IMISB2zm <-
function (N0, B, M, it.max, func, priors, S, v, 
              tauN_sh, tauN_sc, tauF_sh, tauF_sc,
              VN, VF, y0, times, Y, indep, cred)
    {
    y = y0
    T <- length(times)
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
    rho = environment(func)
    Nglobal = NULL
    options(warn = 0)
    samp.values <- matrix(nrow = N0, ncol = 3)
    for (i in 1:3) {
        split1 <- unlist(strsplit(priors[i], "\\("))
        name.dist <- split1[1]
        par.dist <- unlist(strsplit(split1, "\\)"))[2]
        x <- paste("r", name.dist, "(", N0, ",", par.dist, ")", 
            sep = "")
        write(x, file = "foo")
        samp.values[, i] <- eval(source("foo"))$value
    }
    unlink("foo")
    Beta <- samp.values[, 1]
    Q <- samp.values[, 2]
    G <- samp.values[, 3]
    if (length(which(Beta < 0)) > 0) {
        stop("Negative values for Beta were generated from its prior distribution. Try other prior distribution.")
    }
    if (length(which(Q < 0)) > 0) {
        stop("Negative values for Q were generated from its prior distribution. Try other prior distribution.")
    }
    if (length(which(G < 0)) > 0) {
        stop("Negative values for G were generated from its prior distribution. Try other prior distribution.")
    }
    modY <- numeric()
    modY[seq(1, 2 * nrow(Y), by = 2)] <- Y[, 1]
    modY[seq(2, 2 * nrow(Y), by = 2)] <- Y[, 2]
    s1 <- seq(4, (3 * nrow(rbind(c(0, 0), Y)) - 1), by = 3)
    s2 <- seq(5, (3 * nrow(rbind(c(0, 0), Y))), by = 3)
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

    total_it <- N0 + B*it.max

    fromC <- .Call("call_imis_initial", y, times, func, rtol, 
        atol, rho, tcrit, jacfunc, ModelInit, as.integer(verbose), 
        hmin, hmax, as.integer(N0), as.integer(index), as.integer(index2), 
        as.integer(two_n), as.integer(n), as.integer(indep), 
        Y1, modY, VN, VF, Beta, Q, G, v, S, tauN_sh, tauN_sc, 
        tauF_sh, tauF_sc, as.integer(total_it), PACKAGE = "B2Z")

   
    L <- fromC[[N0+4]] 

    C <- 700 - (max(range(L)))
    w <- exp(L+C)/sum(exp(L + C))

    ur <- unlist(fromC[1:N0])
    T <- length(times)
    logCN <- t(log(matrix(ur[seq(2,length(ur),by=3)],T,N0))[-1,])
    logCF <- t(log(matrix(ur[seq(3,length(ur),by=3)],T,N0))[-1,])

    tauN <- fromC[[N0+1]]
    tauF <- fromC[[N0+2]]
    tauNF <- fromC[[N0+3]]    
    pSig <- fromC[[N0+5]]
    parms <- cbind(Beta,Q,G)
    ptheta <- prior.theta(cbind(Beta,Q,G), priors)*pSig


#Importance Sampling Stage
          dimens <- ifelse(indep,5,6)

    expfrac <- numeric()
    k <- 0
    theta <- list()  
    SigmaK <- list()

    if(indep)
      {
      thetas    <- cbind(log(Beta), log(Q), log(G), sqrt(tauN), sqrt(tauF))
      }
    else
      {
      ch.sig <- t(apply(cbind(tauN,tauF,tauNF),1,chol.sigma))
      thetas    <- cbind(log(Beta), log(Q), log(G), ch.sig)
      }

    Qi <- 0

    while(Qi < (1-exp(-1)) & k <= it.max )
      {
      #a
      k <- k + 1
      Nk <- N0 + B*k
      ind.max   <- which.max(w)
      theta[[k]] <- thetas[ind.max,]
      inv.cov <- solve(cov(thetas))
      distances <- apply(thetas, 1, mahalanobis_mod, d = dimens, center = theta[[k]], invcov = inv.cov)
      o <- order(distances)[1:B]
      Mat <- thetas[o,]
      wt <- (w[o] + 1/Nk)/2
      wt <- wt/sum(wt)
      SigmaK[[k]] <- cov.wt(Mat, wt = wt)$cov

      #b
      newinput <- rmvnorm(B, theta[[k]], SigmaK[[k]])
      thetas <- rbind(thetas, newinput)
  

      #c

      newBeta <- exp(newinput[,1])
      newQ <- exp(newinput[,2])
      newG <- exp(newinput[,3])

      if(indep) 
         {
         newtauN <- newinput[,4]^2
         newtauF <- newinput[,5]^2
         newtauNF <- rep(0,B)
         }
      else  
         {
         newtauN <- newinput[,4]^2
         newtauF <- newinput[,5]^2 + newinput[,6]^2
         newtauNF <- newinput[,4]*newinput[,5]
         }


        fromC <- .Call("call_loglik_imis", y, times, func, rtol, 
        atol, rho, tcrit, jacfunc, ModelInit, as.integer(verbose), 
        hmin, hmax, as.integer(B), as.integer(index), as.integer(index2), 
        as.integer(two_n), as.integer(n), as.integer(indep), 
        Y1, modY, VN, VF, newBeta, newQ, newG, newtauN, newtauF, newtauNF, 
        v, S, tauN_sh, tauN_sc, tauF_sh, tauF_sc, as.integer(total_it), PACKAGE = "B2Z")

 
        newL <- fromC[[B+2]] 

        ur <- unlist(fromC[1:B])
        newlogCN <- t(log(matrix(ur[seq(2,length(ur),by=3)],T,B))[-1,])
        newlogCF <- t(log(matrix(ur[seq(3,length(ur),by=3)],T,B))[-1,])

        L <- c(L, newL)
        logCN <- rbind(logCN,newlogCN)
        logCF <- rbind(logCF,newlogCF)

        newpSig <- fromC[[B+1]] 
        newptheta <- prior.theta(cbind(newBeta,newQ,newG), priors)*newpSig

        ptheta <- c(ptheta, newptheta)
        if(k==1) 
          {
          invsigma <- solve(SigmaK[[1]])
          sumHs <- apply(thetas, 1, dmnorm_mod, mean=theta[[1]], d=dimens, varcov=SigmaK[[1]], invvarcov=invsigma)
          }
        else
          {
          sumHsnew <- rep(0,B)
          for(s in 1:(k-1))
             {
             invsigma <- solve(SigmaK[[s]])
             sumHsnew <- sumHsnew + apply(newinput,1, dmnorm_mod, mean=theta[[s]],d=dimens,varcov=SigmaK[[s]], invvarcov=invsigma)
             }
           
  
          sumHs <- c(sumHs, sumHsnew)
          invsigma <- solve(SigmaK[[k]])
          sumHs <- sumHs + apply(thetas, 1, dmnorm_mod, mean=theta[[k]],d=dimens, varcov=SigmaK[[k]], invvarcov=invsigma)
          }  

        qtheta <- (N0/Nk)*ptheta + (B/Nk)*sumHs
        C <- 700 - (max(range(L)))
        w <- exp(L + C)*ptheta/qtheta
        w <- w/sum(w)
        Qi <- sum(1-(1-w)^M)/Nk
        expfrac[k] <- Qi
        }


    

    V.hat <- sum((Nk*w-1)^2)/Nk
    U.hat <-    -log(prod(w^(w/log(Nk))))
    Q.hat <- sum(1-(1-w)^M)
    ESS <- 1/sum(w^2)
    maxw <- max(w)
    options(warn=-1)
    draw_index <- sample(1:Nk, M, replace=TRUE, prob=w)
  
    options(warn=0)
    if(k==(it.max+1)){
    warning("Expected fraction of unique points < (1-1/e)")}

    Betaout <- exp(thetas[draw_index,1])
    Qout <- exp(thetas[draw_index,2])
    Gout <- exp(thetas[draw_index,3])
    TauNout <- thetas[draw_index,4]^2
    logCNout <- logCN[draw_index,] 
    logCFout <- logCF[draw_index,] 
    wout <- w[draw_index]

    if(indep)
     {
     TauNFout <- 0
     TauFout <- thetas[draw_index,5]^2 
     }
    else 
     {
     TauNFout <- thetas[draw_index,4]*thetas[draw_index,5]
     TauFout <- thetas[draw_index,5]^2 + thetas[draw_index,6]^2
     }


    Dbar <- -2*mean((L[draw_index]))

   
    fromC <- .Call("call_loglik_imis", y, times, func, rtol, 
    atol, rho, tcrit, jacfunc, ModelInit, as.integer(verbose), 
    hmin, hmax, as.integer(1), as.integer(index), as.integer(index2), 
    as.integer(two_n), as.integer(n), as.integer(indep), 
    Y1, modY, VN, VF, mean(Betaout), mean(Qout), mean(Gout), 
    mean(TauNout), mean(TauFout), mean(TauNFout), v, S, tauN_sh, 
    tauN_sc, tauF_sh, tauF_sc, as.integer(total_it), PACKAGE = "B2Z")
        
    Dthetabar <- -2*fromC[[3]]

    pD = Dbar - Dthetabar
    DIC = pD + Dbar

  
   
    if(indep)
     {
     r <- list(Beta=Betaout, Q=Qout, G=Gout, tauN=TauNout,
          tauNF = NULL, tauF=TauFout, logCN=logCNout, logCF=logCFout, 
          Y=Y, DIC=DIC, pD=pD, Dbar=Dbar, ESS=ESS, Qi=Qi, 
          indep=indep,y0=y0,times=times, cred=cred, 
          expfrac=expfrac, V.hat=V.hat, U.hat=U.hat, Q.hat=Q.hat, 
          maxw=maxw, w=wout)
     }
    else 
     {
     r <- list(Beta=Betaout, Q=Qout, G=Gout, tauN=TauNout,
          tauNF = TauNFout, tauF=TauFout, logCN=logCNout, 
          logCF=logCFout, Y=Y,
          DIC=DIC, ESS=ESS, pD=pD, Dbar=Dbar, Qi=Qi, indep=indep,
          y0=y0,times=times, cred=cred, 
          expfrac=expfrac, V.hat=V.hat, U.hat=U.hat, Q.hat=Q.hat, 
          maxw=maxw, w=wout)
     }

   
    cat("=|\n")
    flush.console()

    attr(r, "class") <- "imis"
 
    return(r)
    }


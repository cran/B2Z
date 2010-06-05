BCLTB2zm <- function ( priors, S, v, 
                        tauN_sh, tauN_sc, tauF_sh, tauF_sc, VN, 
                        VF, funct, y0, times, Y, indep, 
                         m, cred,indBeta, aBeta, bBeta,
                        indQ, aQ, bQ, indG, aG, bG, size_sample){


                 
                

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

    total_it <- m + size_sample
    next_control_bar = control_bar = 1/50

   								
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

  U <- chol(mat)

      initial[4] <- U[1]
      initial[5] <- U[4]

      if(!indep){initial[6] <- U[3]}

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
        hessian=TRUE, VN=VN,  VF=VF, Y1=Y1, modY=modY,  funct=funct, 
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

      post.mode <- numeric()
      post.mode[1:3] <- exp(nlmobj$par[1:3])
      post.mode[5] <- nlmobj$par[5]^2 
      post.mode[4] <- nlmobj$par[4]^2

      if(!indep)
        {
        post.mode[6] <- nlmobj$par[4]*nlmobj$par[6]
        post.mode[5] <- nlmobj$par[5]^2 + nlmobj$par[6]^2
        }


           sm <- hessian(func = logpost_BCLT, x = post.mode,
       method="Richardson", method.args=list(eps=1e-5, d=0.01,
       zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2), 
       VN=VN, VF=VF, Y1=Y1, modY=modY, 
       funct=funct, y=y, times=times, T=T, indBeta=indBeta, 
       aBeta=aBeta,  bBeta=bBeta, indQ=indQ, aQ=aQ, bQ=bQ, indG=indG, 
       aG=aG, bG=bG, S=S, v=v, tauN_sh=tauN_sh, tauN_sc=tauN_sc, 
       tauF_sh=tauF_sh, tauF_sc=tauF_sc, indep=indep,
       index=index, index2=index2, two_n=two_n, n=n)

       if(length(which(is.nan(sm)==TRUE))> 0){stop("Sorry, it was not 
          possible to estimate the covariance 
          matrix. Try a bigger m.")}

       try(covMat <- -solve(sm),silent=TRUE)
       if(is.null(covMat)){stop("Sorry, it was not 
          possible to estimate the covariance 
          matrix. Try a bigger m.")}

       if(length(which(diag(covMat)<0))>0)
          {
          stop("Sorry, it was not 
          possible to estimate the covariance 
          matrix. Try a bigger m.")
          }


        if(!is.pos.def(covMat))
          {
         stop("Sorry, the estimated covariance 
          matrix is not numerically positive definite. Try a bigger m.")
          }

        dt <- rmvnorm(size_sample, post.mode, covMat)

        inv_pos <- which(dt[,1] < 0 | dt[,2] <0 | dt[,3] <0 | dt[,4] <0 | dt[,5] < 0)
        if(length(inv_pos)> 0)
           dt <- dt[-inv_pos, ]

dt <- cbind(dt, exp(sqrt(dt[,4])), exp(sqrt(dt[,5])))


        post_sum <- cbind(t(apply(dt,2,quantile,p=c((100-cred)/200,0.5,(100+cred)/200))),
                     matrix(apply(dt,2,mean),(dim+2),1), c(post.mode,exp(sqrt(post.mode[4])), exp(sqrt(post.mode[5]))),
                     matrix(apply(dt,2,sd),(dim+2),1))
        colnames(post_sum) <- c(paste((100-cred)/2,"%",sep=""), "Median", paste((100+cred)/2,"%",sep=""), "Mean", "Mode", "SD")
        post_sum <- round(post_sum, 4)

        B <- nrow(dt)

        Beta <- dt[,1]
        Q <- dt[,2]
        G <- dt[,3]
        TauN <- dt[,4]
        TauF <- dt[,5]

        if(indep){TauNF <- rep(0,B)}
        else{TauNF <- dt[,6]}

        fromC <- .Call("call_loglik_imis", y, times, funct, rtol, 
        atol, rho, tcrit, jacfunc, ModelInit, as.integer(verbose), 
        hmin, hmax, as.integer(B), as.integer(index), as.integer(index2), 
        as.integer(two_n), as.integer(n), as.integer(indep), 
        Y1, modY, VN, VF, Beta, Q, G, TauN, TauF, TauNF, 
        v, S, tauN_sh, tauN_sc, tauF_sh, tauF_sc, as.integer(total_it), PACKAGE = "B2Z")

 
        L <- fromC[[B+2]] 

        ur <- unlist(fromC[1:B])

        logCN <- t(log(matrix(ur[seq(2,length(ur),by=3)],T,B))[-1,])
        logCF <- t(log(matrix(ur[seq(3,length(ur),by=3)],T,B))[-1,])


        Dbar <- -2*mean(L)

        fromC <- .Call("call_loglik_imis", y, times, funct, rtol, 
        atol, rho, tcrit, jacfunc, ModelInit, as.integer(verbose), 
        hmin, hmax, as.integer(1), as.integer(index), as.integer(index2), 
        as.integer(two_n), as.integer(n), as.integer(indep), 
        Y1, modY, VN, VF, mean(Beta), mean(Q), mean(G), 
        mean(TauN), mean(TauF), mean(TauNF), v, S, tauN_sh, 
        tauN_sc, tauF_sh, tauF_sc, as.integer(total_it), PACKAGE = "B2Z")
        
        Dthetabar <- -2*fromC[[3]]

        pD = Dbar - Dthetabar
        DIC = pD + Dbar

    cat("|\n")

    flush.console()

    r <- list(post.sum = post_sum, covMat=covMat, logCN=logCN, logCF=logCF, pD=pD, Dbar=Dbar, DIC=DIC, indep=indep, y0=y0, Y=Y, Beta=Beta, Q=Q, G=G, tauN=TauN, tauF=TauF, tauNF=TauNF, times=times, cred=cred)

    attr(r, "class") <- "bclt"

    return(r)

    }

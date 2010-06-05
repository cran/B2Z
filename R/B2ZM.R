B2ZM <-
function(func = func_default, y0 = c(0,0), data = NULL, priorBeta, 
                 priorQ, priorG, 
                 v, S, tauN.sh, tauN.sc,
                 tauF.sh, tauF.sc,  V_N,
                 V_F, indep.model = FALSE, credibility = 95,
                 sampler = c("IMIS", "METROP", "GIBBS", "SIR", "BCLT"),
                 sir.control = list(m=50000),
                 bclt.control = list(m=7000, sample_size= 2000),
                 metrop.control = list(NUpd = 10000, burnin = 1000, 
                 lag = 1, initial=NULL, Sigma.Cand=NULL, m=5000),
                 imis.control = list(N0=6000, B=600, M=3000, it.max=12),
                 gibbs.control = list(NUpd = 10000, burnin = 1000, lag = 1, initial=NULL, Sigma.Cand=NULL, numstp=100),
                 figures = list(save = FALSE, type =c("ps", "eps","pdf", "png", "jpg"))) {

	if(is.null(data)) {stop("data set not found")}
	if(!is.matrix(data)) {stop("data set must be a 3-column matrix")}
 	if((is.matrix(data) & ncol(data)!=3)) {stop("data set must be a 3-column matrix")}
 	if(is.character(data)) {stop("data set must contain numeric elements")}
 	if(length(which(is.na(data)==TRUE))> 0) {stop("data set contains missing values")}
	if(length(which(data < 0))> 0) {stop("data set contains negative values")}
      
     if(!exists("priorBeta")){stop("Prior distribution for Beta was not declared")}
     if(!exists("priorQ")){stop("Prior distribution for Q was not declared")}
     if(!exists("priorG")){stop("Prior distribution for G was not declared")}



     
     if(indep.model)
       {
       S <- matrix(0,2,2)
       v <- 0

       if(!exists("tauN.sh")){stop("Shape parameter in the prior 
         distribution for tauN was not declared")}
       if(!exists("tauN.sc")){stop("Scale parameter in the prior 
         distribution for tauN was not declared")}
       if(!exists("tauF.sh")){stop("Shape parameter in the prior 
         distribution for tauF was not declared")}
       if(!exists("tauF.sc")){stop("Scale parameter in the prior 
         distribution for tauF was not declared")}

        if(!is.numeric(tauN.sh)){stop("'tauN.sh' must be numeric")}
        if(!is.numeric(tauN.sc)){stop("'tauN.sc' must be numeric")}
        if(!is.numeric(tauF.sh)){stop("'tauF.sh' must be numeric")}
        if(!is.numeric(tauF.sc)){stop("'tauF.sc' must be numeric")}

        if(tauN.sh <= 0){stop("'tauN.sh' must be greater than 0")}
        if(tauN.sc <= 0){stop("'tauN.sc' must be greater than 0")}
        if(tauF.sh <= 0){stop("'tauF.sh' must be greater than 0")}
        if(tauF.sc <= 0){stop("'tauF.sc' must be greater than 0")}

       }
     else
       {
       tauN.sh = tauN.sc = tauF.sh = tauF.sc = 0
       if(!exists("v")){stop("parameter v in the joint prior distribution for tauN, tauNF and tauF was not declared")}
       if(!exists("S")){stop("parameter S in the joint prior distribution for tauN, tauNF and tauF was not declared")}
        if(!is.matrix(S)){stop("object 'S' must be a matrix")}
        if(!is.numeric(S)){stop("object 'S' must be a numeric matrix")}
        if(nrow(S)!=ncol(S)){stop("object 'S' must be a square matrix")}
        if(nrow(S)!=2){stop("object 'S' must be a 2x2 matrix")}
        try(testingPD <- chol(S), silent=TRUE)
        if(!exists("testingPD")){stop("object 'S' must be a positive definite matrix")}

        if(v<=1){stop("'v' must be greater than 1")}

       }

      priors <- c(priorBeta, priorQ, priorG)
      name.dist <- character()
      par.dist <- character()

      for(i in 1:3)
       { 
       split1 <- unlist(strsplit(priors[i],  "\\(" ))
       name.dist[i] <- paste("d",split1[1],sep="")
       par.dist[i] <- unlist(strsplit(split1,  "\\)" ))[2]
       }

      if(!exists(name.dist[1], mode="function")){stop("invalid prior for Beta")}
      if(!exists(name.dist[2], mode="function")){stop("invalid prior for Q")}
      if(!exists(name.dist[3], mode="function")){stop("invalid prior for G")}
      if(is.na(par.dist[1])){stop("missing parameters for the prior distribution of Beta")}
      if(is.na(par.dist[2])){stop("missing parameters for the prior distribution of Q")}
      if(is.na(par.dist[3])){stop("missing parameters for the prior distribution of G")}
      
     if(!is.numeric(V_F)){stop("'V_F' must be numeric")}
     if(V_F < 0){stop("'V_F' must be greater than 0")}

     if(!is.numeric(V_N)){stop("'V_N' must be numeric")}
     if(V_N < 0){stop("'V_N' must be greater than 0")}

      indep <- indep.model  

   y0 <- y0[1:2]
   if(!is.numeric(y0)){stop("y0 must be a numeric vector")}

   if(length(which(y0 < 0))){stop("initial concentration must be greater or equal than 0")}
    
      times <- c(0,data[,1])
   
      if (min(times[-1])<=0) {stop("the observed time values must be greater than 0")}
      
      order_times <- order(times[-1])
      data[,1] <- data[order_times,1]
      data[,2] <- data[order_times,2]
      data[,3] <- data[order_times,3]
       
      times <- sort(times)
      Y <- log(data[,2:3]) 

      if (!is.function(func) && !is.character(func)) 
        stop("`func' must be a function or character vector")
     
      cred = credibility
      if(cred <=0 || cred >=100) {stop("'credibility' must be a number between 0 and 100")}      

      if(figures$save==TRUE)
         {
         if(length(which(c("ps","eps","pdf","jpg","png")==figures$type))==0)
           {
           stop("in 'figures': object 'type' must be one of the options: 'ps', 'eps', 'pdf', 'jpg' or 'png'")
           }
         }

       prval <- which.dist(priorBeta)
       indBeta <- prval[1]
       aBeta <- prval[2]
       bBeta <- prval[3]

       prval <- which.dist(priorQ)
       indQ <- prval[1]
       aQ <- prval[2]
       bQ <- prval[3]

       prval <- which.dist(priorG)
       indG <- prval[1]
       aG <- prval[2]
       bG <- prval[3]

      if(is.null(indBeta)){stop("invalid prior for Beta")}
      if(is.null(indQ)){stop("invalid prior for Q")}
      if(is.null(indG)){stop("invalid prior for G")}

      verifying_prior_Beta(indBeta,aBeta,bBeta)
	verifying_prior_Q(indQ,aQ,bQ)
	verifying_prior_G(indG,aG,bG)


      if(sampler == "SIR")
        {
        if(is.null(sir.control$m)){stop("object 'm' not found")}
        if(sir.control$m <= 0){stop("object 'm' must be greater than 0")}
        if(!is.numeric(sir.control$m)){stop("object 'm' must be numeric")}

        m <- sir.control$m        

       ans <- SIRB2zm(m, Y, func, y0, times, V_N, V_F, indBeta, aBeta, 
       bBeta, indQ, aQ, bQ, indG, aG, bG, v, as.numeric(S), tauN.sh,  
       tauN.sc, tauF.sh, tauF.sc, indep, cred) 

        }
    
     if(sampler == "METROP")
        {
        namesobj <- c("NUpd", "burnin", "lag", "initial","Sigma.Cand", "m")
        
        cont <- 0
        index <- numeric()
        indexT <- numeric()
        for(i in 1:length(names(metrop.control)))
           {
           index <- which(namesobj==names(metrop.control)[i])
           cont <- cont + length(which(namesobj==names(metrop.control)[i]))
           indexT <- c(indexT,index)
           }

        if(length(names(metrop.control))!= cont){stop("some object in list 'metrop.control' has an incorrect name")}

        chosen <- seq(1:6)[-indexT]
        for(i in chosen)
          {
          if(i==1){metrop.control$NUpd <- 10000}
          if(i==2){metrop.control$burnin <- 1000}
          if(i==3){metrop.control$lag <- 1}
          if(i==4){metrop.control$initial <- NULL}
          if(i==5){metrop.control$Sigma.Cand <- NULL}
          if(i==6){metrop.control$m <- 5000}
          }


        if(metrop.control$NUpd<=0){stop("object 'NUpd' must be greater than 0")}
        if(metrop.control$burnin<=0){stop("object 'burnin' must be greater than 0")}
        if(metrop.control$burnin >= metrop.control$NUpd){stop("object 'burnin' must be less than object 'NUpd'")}
        if(metrop.control$lag<=0){stop("object 'lag' must be greater than 0")}
        if(metrop.control$lag >= metrop.control$NUpd-metrop.control$burnin){stop("object 'lag' must be less than object 'NUpd'-'burnin'")}
        
        if(!is.null(metrop.control$initial))
           {
           initial <- metrop.control$initial
           if(indep){if(length(initial)!=5){stop("in list 'metrop.control': object 'initial' must be a vector with 5 elements")}}
           else{if(length(initial)!=6){stop("in list 'metrop.control': object 'initial' must be a vector with 6 elements")}}


           if(initial[1]<=0){stop("initial value for Beta must be greater than 0")}
           if(initial[2]<=0){stop("initial value for Q must be greater than 0")}
           if(initial[3]<=0){stop("initial value for G must be greater than 0")}
           if(initial[4]<=0){stop("initial value for Tau_N must be greater than 0")}
           if(initial[5]<=0){stop("initial value for Tau_F must be greater than 0")}

           if(indep)
            {
            initial[4] <- sqrt(initial[4])
            initial[5] <- sqrt(initial[5])
            }
           else
            {
            test1 <- matrix(c(initial[4], initial[6], initial[6], initial[5]),2,2)
            try(testingPD1 <- chol(test1), silent=TRUE)
            if(!exists("testingPD1")){stop("initial values for Tau_N, Tau_NF and TauF must define a positive definite matrix")}
            initial[4] <- testingPD1[1,1]
            initial[5] <- testingPD1[2,2]
            initial[6] <- testingPD1[1,2]
            }
          }

        if(!is.null(metrop.control$Sigma.Cand))
          {
          if(!is.matrix(metrop.control$Sigma.Cand)) {stop("object 'Sigma.Cand' must be a matrix")}
          if(nrow(metrop.control$Sigma.Cand)!= ncol(metrop.control$Sigma.Cand)) {stop("object 'Sigma.Cand' must be a square matrix")}
          if(indep){if(nrow(metrop.control$Sigma.Cand)!= 5) {stop("object 'Sigma.Cand' must be a 5x5 matrix")}}
          else{if(nrow(metrop.control$Sigma.Cand)!= 6) {stop("object 'Sigma.Cand' must be a 6x6 matrix")}}
          try(testingPD2 <- chol(metrop.control$Sigma.Cand), silent=TRUE)
          if(!exists("testingPD2")){stop("'Sigma.Cand' must define a positive definite matrix")}
          }

        if(metrop.control$m <=0){stop("'m' must be greater than 0")}

        NUpd <- metrop.control$NUpd
        burnin <- metrop.control$burnin
        lag <- metrop.control$lag
        initial <- metrop.control$initial
        Sigma.Cand <- metrop.control$Sigma.Cand
        m <- metrop.control$m


        ans <- METROPB2zm(NUpd, burnin, lag, initial, priors,  
               as.numeric(S), v, tauN.sh, tauN.sc, tauF.sh, tauF.sc, 
               V_N, V_F, func, y0, times, Y, indep, 
               Sigma.Cand, m, cred, indBeta, aBeta, bBeta,
               indQ, aQ, bQ, indG, aG, bG)

        }


if(sampler=="BCLT")
       {

       if(is.null(bclt.control$m)){stop("object 'm' not found")}
        if(bclt.control$m <= 0){stop("object 'm' must be greater than 0")}
        if(!is.numeric(bclt.control$m)){stop("object 'm' must be numeric")}

if(is.null(bclt.control$sample_size)){stop("object 'sample_size' not found")}
        if(bclt.control$sample_size <= 0){stop("object 'sample_size' must be greater than 0")}
        if(!is.numeric(bclt.control$sample_size)){stop("object 'sample_size' must be numeric")}
       

        m <- as.integer(bclt.control$m)
        sample_size <- bclt.control$sample_size
        ans <- BCLTB2zm(priors, S, v, 
                tauN.sh, tauN.sc, tauF.sh, tauF.sc, V_N, 
                  V_F, func, y0, times, Y, indep, m,         
                  cred, indBeta, aBeta, bBeta,
                  indQ, aQ, bQ, indG, aG, bG, sample_size)

        }

     if(sampler=="IMIS")
        {
        if(!is.list(imis.control)){stop("'imis.control' must be a list containg the objects: 'N0', 'B', 'M' and 'it.max'")}
        namesobj <- c("N0", "B", "M", "it.max")
        
        cont <- 0
        index <- numeric()
        indexT <- numeric()
        for(i in 1:length(names(imis.control)))
           {
           index <- which(namesobj==names(imis.control)[i])
           cont <- cont + length(which(namesobj==names(imis.control)[i]))
           indexT <- c(indexT,index)
           }

        if(length(names(imis.control))!= cont){stop("some object in list 'imis.control' has an incorrect name")}

        chosen <- seq(1:4)[-indexT]
        for(i in chosen)
          {
          if(i==1){imis.control$N0 <- 6000}
          if(i==2){imis.control$B <- 600}
          if(i==3){imis.control$M <- 3000}
          if(i==4){imis.control$it.max <- 12}
          }

        if(imis.control$N0 <=0){stop("'N0' must be greater than 0")}
        if(imis.control$B <=0){stop("'B' must be greater than 0")}
        if(imis.control$M <=0){stop("'M' must be greater than 0")}
        if(imis.control$it.max <=0){stop("'it.max' must be greater than 0")}

        N0 <- imis.control$N0
        B <- imis.control$B
        M <- imis.control$M
        it.max <- imis.control$it.max
          
       
      ans <- IMISB2zm(N0, B, M, it.max, func, priors, as.numeric(S), v, 
              tauN.sh, tauN.sc, tauF.sh, tauF.sc,
              V_N, V_F, y0, times, Y, indep, cred)


        }


     if(sampler == "GIBBS")
        {
        if(!is.list(gibbs.control)){stop("'gibbs.control' must be a list containg the objects: 'NUpd', 'burnin', 'lag', 'initial', 'Sigma.Cand', 'm'")}
        namesobj <- c("NUpd", "burnin", "lag", "initial","Sigma.Cand", "m")
        
        cont <- 0
        index <- numeric()
        indexT <- numeric()
        for(i in 1:length(names(gibbs.control)))
           {
           index <- which(namesobj==names(gibbs.control)[i])
           cont <- cont + length(which(namesobj==names(gibbs.control)[i]))
           indexT <- c(indexT,index)
           }

        if(length(names(gibbs.control))!= cont){stop("some object in list 'gibbs.control' has an incorrect name")}

        chosen <- seq(1:6)[-indexT]
        for(i in chosen)
          {
          if(i==1){gibbs.control$NUpd <- 10000}
          if(i==2){gibbs.control$burnin <- 1000}
          if(i==3){gibbs.control$lag <- 1}
          if(i==4){gibbs.control$initial <- NULL}
          if(i==5){gibbs.control$Sigma.Cand <- NULL}
          if(i==6){gibbs.control$m <- 5000}
          }


        if(gibbs.control$NUpd<=0){stop("object 'NUpd' must be greater than 0")}
        if(gibbs.control$burnin<=0){stop("burnin' must be greater than 0")}
        if(gibbs.control$burnin >= gibbs.control$NUpd){stop("object 'burnin' must be less than object 'NUpd'")}
        if(gibbs.control$lag<=0){stop("object 'lag' must be greater than 0")}
        if(gibbs.control$lag >= gibbs.control$NUpd-gibbs.control$burnin){stop("object 'lag' must be less than object 'NUpd'-'burnin'")}
        
        if(!is.null(gibbs.control$initial))
           {
           initial <- gibbs.control$initial
           if(length(initial)!=3){stop("object 'initial' must be a vector with 3 elements")}
           if(initial[1]<=0){stop("initial value for Beta must be greater than 0")}
           if(initial[2]<=0){stop("initial value for Q must be greater than 0")}
           if(initial[3]<=0){stop("initial value for G must be greater than 0")}

     
          }

        if(!is.null(gibbs.control$Sigma.Cand))
          {
          if(!is.matrix(gibbs.control$Sigma.Cand)) {stop("object 'Sigma.Cand' must be a matrix")}
          if(nrow(gibbs.control$Sigma.Cand)!= ncol(gibbs.control$Sigma.Cand)) {stop("object 'Sigma.Cand' must be a square matrix")}
          if(nrow(gibbs.control$Sigma.Cand)!= 3) {stop("object 'Sigma.Cand' must be a 3x3 matrix")}
          try(testingPD2 <- chol(gibbs.control$Sigma.Cand), silent=TRUE)
          if(!exists("testingPD2")){stop("'Sigma.Cand' must define a positive definite matrix")}
          }

        if(gibbs.control$m <=0){stop("'m' must be greater than 0")}

        NUpd <- gibbs.control$NUpd
        burnin <- gibbs.control$burnin
        lag <- gibbs.control$lag
        initial <- gibbs.control$initial
        Sigma.Cand <- gibbs.control$Sigma.Cand
        m <- gibbs.control$m

      ans <- GIBBSB2zm(NUpd, burnin, lag, initial, priors,  
               as.numeric(S), v, tauN.sh, tauN.sc, tauF.sh, tauF.sc, 
               V_N, V_F, func, y0, times, Y, indep, 
               Sigma.Cand, m, cred, indBeta, aBeta, bBeta,
               indQ, aQ, bQ, indG, aG, bG)


        }


      if(figures$save==TRUE)
           {
           if(figures$type=="ps")
               {plotps(ans)}
           else
             {
             if(figures$type=="pdf")
               {plotpdf(ans)}
             else
               {
               if(figures$type=="eps"){ploteps(ans)}
               else
                 {
                 if(figures$type=="png"){plotpng(ans)}
                 else{plotjpeg(ans)}
                 }
               }
             }
           }      

return(ans)

}


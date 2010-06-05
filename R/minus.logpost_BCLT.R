minus.logpost_BCLT <-
function(parms_full, VN, VF, Y1, modY, funct, y, 
         times, T, indBeta, aBeta, bBeta,
         indQ, aQ, bQ, indG, aG, bG, S, v, 
         tauN_sh, tauN_sc, tauF_sh, tauF_sc, 
         indep, index, index2, two_n, n)
        {
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

        Beta <- parms_full[1]
        Q <- parms_full[2]
        G <- parms_full[3]
        tauN <- parms_full[4]

        if(indep)
         {
         tauF <- parms_full[5]
         tauNF <- 0
         }
        else
         {
         tauF <- parms_full[5]
         tauNF <- parms_full[6]
         }

     
        fromC <- .Call("call_metrop_logpost", y, times, funct, rtol, 
        atol, rho, tcrit, jacfunc, ModelInit, as.integer(verbose), 
        hmin, hmax, as.integer(index), as.integer(index2), 
        as.integer(two_n), as.integer(n), as.integer(indep), 
        Y1, modY, VN, VF, Beta, Q, G, tauN, tauF, tauNF, 
        v, S, tauN_sh, tauN_sc, tauF_sh, tauF_sc,
        as.integer(indBeta), aBeta, bBeta,
        as.integer(indQ), aQ, bQ,
        as.integer(indG), aG, bG, PACKAGE = "B2Z")

        ans <- - fromC

        return(ans)
        }
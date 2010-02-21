logpost <-
function(parms_full, VN, VF, Y1, modY, funct, y, 
         times, T, indBeta, aBeta, bBeta,
         indQ, aQ, bQ, indG, aG, bG, S, v, 
         tauN_sh, tauN_sc, tauF_sh, tauF_sc, 
         indep, index, index2, two_n, n){


ans <- -minus.logpost(parms_full, VN, VF, Y1, modY, funct, y, 
         times, T, indBeta, aBeta, bBeta,
         indQ, aQ, bQ, indG, aG, bG, S, v, 
         tauN_sh, tauN_sc, tauF_sh, tauF_sc, 
         indep, index, index2, two_n, n)

return(ans)
}
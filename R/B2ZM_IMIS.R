B2ZM_IMIS <-
function(func = func_default, y0 = c(0,0), data, priorBeta, 
                 priorQ, priorG, v, S, tauN.sh, tauN.sc,
                 tauF.sh, tauF.sc,  V_N,
                 V_F, indep.model = FALSE, credibility = 95,
                 N0 = 6000, B = 600, M = 3000, it.max = 12, 
                 figures = list(save = FALSE, type =c("ps", 
                 "eps","pdf", "png", "jpg"))) {

if(!indep.model){tauN.sh = tauN.sc = tauF.sh = tauF.sc = 0}
else
{
S <- matrix(0,2,2)
v <- 0
}


ans <- B2ZM(func=func, y0 = y0, data=data, priorBeta=priorBeta, priorQ=priorQ, 
       priorG=priorG, v=v, S=S, tauN.sh=tauN.sh, 
       tauN.sc=tauN.sc, tauF.sh=tauF.sh, tauF.sc=tauF.sc, V_N= V_N,  
       V_F=V_F, indep.model = indep.model, 
       credibility = 95,  sampler = "IMIS", imis.control = list(N0=N0, 
       B=B, M=M, it.max=it.max),  figures = figures)

return(ans)

}
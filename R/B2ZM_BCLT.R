B2ZM_BCLT <-
function(func = func_default, y0, data = NULL, priorBeta, 
                 priorQ, priorG, v, S, tauN.sh, tauN.sc,
                 tauF.sh, tauF.sc,  V_N, V_F, indep.model = FALSE,  
                 credibility = 95,  m = 7000, sample_size=2000, 
                 figures = list(save = FALSE, type =c("ps", 
                 "eps","pdf", "png", "jpg"))) {

if(!indep.model){tauN.sh = tauN.sc = tauF.sh = tauF.sc = 0}
else
{
S <- matrix(0,2,2)
v <- 0
}


ans <- B2ZM(func=func, y0=y0, data=data, priorBeta=priorBeta, 
       priorQ=priorQ, priorG=priorG, v=v, S=S, tauN.sh=tauN.sh, 
       tauN.sc=tauN.sc, tauF.sh=tauF.sh, tauF.sc=tauF.sc, V_N= V_N,  
       V_F=V_F, indep.model = indep.model, 
       credibility = 95,  sampler = "BCLT", 
       bclt.control = list(m=m, sample_size=sample_size),  
       figures = figures)

return(ans)

}
B2ZM_GIBBS <-
function(func = func_default, y0, data = NULL, priorBeta, 
                 priorQ, priorG, v, S, tauN.sh, tauN.sc,
                 tauF.sh, tauF.sc,  V_N, V_F, indep.model = FALSE,  
                 credibility = 95,  NUpd = 10000, burnin = 1000, 
                 lag = 1, initial=NULL, Sigma.Cand=NULL, m = 6000, 
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
       credibility = 95,  sampler = "GIBBS", 
       gibbs.control = list(NUpd = NUpd, burnin = burnin, lag = lag,  
       initial=initial, Sigma.Cand=Sigma.Cand, m=m),  
       figures = figures)

return(ans)

}
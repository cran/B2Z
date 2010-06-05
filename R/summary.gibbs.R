summary.gibbs <-
function(object,...) {
index <- seq(object$burnin, object$NUpd, by = object$lag)
Beta <- object$Beta[index]
Q <- object$Q[index]
G <- object$G[index]
Tau_F <- object$tauF[index]
Tau_N <- object$tauN[index]
GSDF <- exp(sqrt(Tau_F))
GSDN <- exp(sqrt(Tau_N))

cred <- object$cred
{if(object$indep){
suma <- t(apply(rbind(Beta,Q,G,Tau_N,Tau_F,GSDN, GSDF),1,summary.out,cred=cred))
PostCovMat <- cov(cbind(Beta,Q,G,Tau_N,Tau_F))
 rownames(suma) <- c("Beta", "Q", "G", "Tau_N", "Tau_F", "GSD(Tau_N)",  "GSD(Tau_F)")
}
else{
Tau_NF <- object$tauNF[index]
suma <- t(apply(rbind(Beta,Q,G,Tau_N,Tau_NF,Tau_F, GSDN, GSDF),1,summary.out, cred=cred))
PostCovMat <- cov(cbind(Beta,Q,G,Tau_N,Tau_F, Tau_NF))
rownames(suma) <- c("Beta", "Q", "G", "Tau_N", "Tau_NF", "Tau_F", "GSD(Tau_N)", "GSD(Tau_F)")

}}

colnames(suma)[c(1,2,4)] <- c("Mean", "SD", "Median")

dic <- object$DIC
AR <- object$AR
Sigma.Cand <- object$Sigma.Cand
initial <- matrix(object$initial,1,3)
ESS <- matrix(object$ESS,1,3)
rownames(ESS) <- rownames(initial) <- ""
colnames(ESS) <- c("Beta", "Q", "G")
colnames(initial)<- colnames(Sigma.Cand) <- rownames(Sigma.Cand) <- c("log(Beta)", "log(Q)", "log(G)")

ans <- list(summary = suma,
            PostCovMat=PostCovMat, 
            DIC = dic, 
            ESS = ESS, 
            AcceptRate = AR,
            InitVal = initial,
            CovMat = Sigma.Cand,
            indep = object$indep)

class(ans) <- 'summary.gibbs'
ans
}


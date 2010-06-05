summary.metrop <-
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
suma <- t(apply(rbind(Beta,Q,G,Tau_N,Tau_F, GSDN, GSDF),1,summary.out,cred=cred))
PostCovMat <- cov(cbind(Beta,Q,G,Tau_N,Tau_F))
}
else{
Tau_NF <- object$tauNF[index]
suma <- t(apply(rbind(Beta,Q,G,Tau_N,Tau_NF,Tau_F, GSDN, GSDF),1,summary.out, cred=cred))
PostCovMat <- cov(cbind(Beta,Q,G,Tau_N,Tau_F, Tau_NF))
}}

colnames(suma)[c(1,2,4)] <- c("Mean", "SD", "Median")

dic <- object$DIC
pD <- object$pD
Dbar <- object$Dbar
AR <- object$AR
Sigma.Cand <- object$Sigma.Cand


if(object$indep) 
  {
  initial <- matrix(object$initial,1,5)
  ESS <- matrix(object$ESS,1,5)
  rownames(ESS) <- rownames(initial) <- ""
  colnames(initial)<-colnames(ESS) <-  c("Beta", "Q", "G", "Tau_N", "Tau_F")
  colnames(Sigma.Cand) <- rownames(Sigma.Cand) <- c("log(Beta)", "log(Q)", "log(G)", "l1", "l2")

 rownames(suma) <- c("Beta", "Q", "G", "Tau_N", "Tau_F", "GSD(Tau_N)",  "GSD(Tau_F)")

  }
else
  {
  initial <- matrix(object$initial,1,6)
  ESS <- matrix(object$ESS,1,6)
  rownames(ESS) <- rownames(initial) <- ""
  colnames(initial)<-colnames(ESS) <- c("Beta", "Q", "G", "Tau_N", "Tau_F", "Tau_NF")
  colnames(Sigma.Cand) <- rownames(Sigma.Cand) <- c("log(Beta)", "log(Q)", "log(G)", "l1", "l2", "l12")

rownames(suma) <- c("Beta", "Q", "G", "Tau_N", "Tau_NF", "Tau_F", "GSD(Tau_N)", "GSD(Tau_F)")

  }

ans <- list(summary = suma,
            PostCovMat=PostCovMat, 
            DIC = dic, 
            pD = pD,
            Dbar = Dbar,
            ESS = ESS, 
            AcceptRate = AR,
            InitVal = initial,
            CovMat = Sigma.Cand,
            indep = object$indep)

class(ans) <- 'summary.mh'
ans
}


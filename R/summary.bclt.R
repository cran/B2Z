summary.bclt <-
function(object,...) {
suma <- object$post.sum

cred <- object$cred
PostCovMat <- round(object$covMat,4)

{if(object$indep){
 rownames(suma) <- c("Beta", "Q", "G", "Tau_N", "Tau_F", "GSD(Tau_N)",  "GSD(Tau_F)")
rownames(PostCovMat) <- colnames(PostCovMat) <- c("Beta", "Q", "G", "Tau_N", "Tau_F")
}
else{
suma <- suma[c(1:4,6,5,7:8),]
rownames(suma) <- c("Beta", "Q", "G", "Tau_N", "Tau_NF", "Tau_F", "GSD(Tau_N)", "GSD(Tau_F)")

rownames(PostCovMat) <- colnames(PostCovMat) <- c("Beta", "Q", "G", "Tau_N", "Tau_F", "Tau_NF")
}}


dic <- object$DIC
pD <- object$pD
Dbar <- object$Dbar

ans <- list(summary = suma,
            PostCovMat=PostCovMat, 
            DIC = dic,
            pD = pD,
            Dbar = Dbar,
            indep = object$indep)

class(ans) <- 'summary.bclt'
ans
}


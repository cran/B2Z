plot.bclt <-
function(x,...)
  {
  options(warn=-1)
  Y <- x$Y

{   if(x$indep){parms <- cbind(x$tauN, 0, x$tauF, x$logCN, x$logCF)}
   else{parms <- cbind(x$tauN, x$tauNF, x$tauF, x$logCN, x$logCF)}}
   n <- length(x$times)- 1
   r <- apply(parms, 1, predY, n=n )


   alpha <- (100 - x$cred)/100
   q <- apply(r,1,quantile,prob=c(alpha/2,0.5,1-alpha/2))
   y0 <- x$y0
   plot(x$times, c(y0[1],Y[,1]),ylim=c(min(c(y0[1],Y[,1])), max(Y[,1])+2/7*abs(diff(range(c(y0[1],Y[,1]))))),type="l", xlab="time", ylab=expression(log(C[N])), main="Log concentrations at near field")
   lines(x$times, c(y0[1], q[1,1:n]),col="red",lty = 2)
   lines(x$times, c(y0[1], q[2,1:n]),col="blue", lwd=2, lty = 3)
   lines(x$times, c(y0[1], q[3,1:n]),col="red", lty = 2)
   legend("topright",lty=c(1,3,2), lwd= c(1,2,1), col=c("black","blue","red"),legend=c(expression(paste("data: ",log(C[N]))), "posterior median", paste((1-alpha)*100,"% posterior predictive interval",sep="")), bty="n", cex=0.85)

   devAskNewPage(ask=TRUE)  

   plot(x$times,c(y0[2],Y[,2]),ylim=c(min(c(y0[2],Y[,2])), max(Y[,2])+2/7*abs(diff(range(c(y0[2],Y[,2]))))),type="l", xlab="time", ylab=expression(log(C[F])), main="Log concentrations at far field")
   lines(x$times,c(y0[2], q[1,(n+1):(2*n)]),col="red",lty = 2)
   lines(x$times,c(y0[2], q[2,(n+1):(2*n)]),col="blue", lwd=2, lty = 3)
   lines(x$times,c(y0[2], q[3,(n+1):(2*n)]),col="red", lty = 2)
   legend("topright",lty=c(1,3,2), lwd= c(1,2,1), col=c("black","blue","red"),legend=c(expression(paste("data: ",log(C[F]))), "posterior median", paste((1-alpha)*100,"% posterior predictive interval",sep="")), bty="n", cex=0.85)


   options(warn=0)
   }


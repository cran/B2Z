prior.theta <-
function(params, priors)
 {
 n <- nrow(params)
 pBetaQG <- matrix(nrow=n,ncol=3)
 for(i in 1:3)
   {
  split1 <- unlist(strsplit(priors[i],  "\\(" ))
   name.dist <- split1[1]
   par.dist <- unlist(strsplit(split1,  "\\)" ))[2]
   xnam <- paste("d",name.dist,"(c(",paste(params[1:(n-1),i],",",sep="", collapse=""),params[n,i],"),",par.dist,")", sep="")
   write(xnam,file="foo")
   pBetaQG[,i] <- eval(source("foo"))$value
   }
 unlink("foo")

 
 return(apply(pBetaQG,1,prod))
 }
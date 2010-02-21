func_default <-
function(t,y,parms)
   {
   ydot <- matrix(0,2,1)
   ydot[1,1] <- (parms[1]/parms[2])*(y[2]-y[1]) + parms[5]/parms[2]
   ydot[2,1] <- (1/parms[3])*(parms[1]*(y[1]-y[2])-y[2]*parms[4])
   return(list(ydot))
   }


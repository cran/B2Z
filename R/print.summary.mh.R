print.summary.mh <-
function(x, digits = max(options()$digits - 4, 3),...)
  {
  cat("\nPosterior Summary Statistics: \n")
  print(x$summary, digits = digits)
 
  cat("\nPosterior Covariance Matrix: \n")
  print(x$PostCovMat, digits = digits)

  cat("\nDeviance Information Criterion: ")

  cat("\npD: ")
  dput(round(x$pD,digits))

  cat("Dbar: ")
  dput(round(x$Dbar,digits))

  cat("DIC: ")
  dput(round(x$DIC,digits))


  cat("\n\nSampler used: Metropolis")
  cat("\n--------------------------------\n")

  cat("\nEffective Sample Size (ESS): \n")
  print(x$ESS, digits=digits)

  cat("\nMCMC Acceptance Rate: ")
  dput(round(x$AcceptRate,digits))
  
  cat("\nInitial Values: \n")
  print(x$InitVal)
  
  cat("\nCovariance Matrix used in the proposal distribution (transformed scale):\n")
  print(x$CovMat)
 
  if(x$indep) {cat("\nwhere l1=sqrt(Tau_N) and l2 = sqrt(Tau_F).\n")}
  else {cat("\nwhere l1=sqrt(Tau_N), l12 = Tau_NF/sqrt(Tau_N), and l2 = sqrt(Tau_F - (Tau_NF)^2 /Tau_N).\n")}
  }


print.summary.gibbs <-
function(x, digits = max(options()$digits - 4, 3),...)
  {
  cat("\nPosterior Summaries: \n")
  print(x$summary, digits = digits)
 
  cat("\nNote: GSD is the geometric standard deviation, i.e., GSD(x) = exp(sqrt(x))\n\n")

  cat("\nPosterior Covariance Matrix: \n")
  print(x$PostCovMat, digits = digits)

  cat("\nDeviance Information Criterion (DIC): ")
  dput(x$DIC)

  cat("\n\nSampler used: Gibbs with Metropolis step")
  cat("\n--------------------------------\n")

  cat("\nEffective Sample Size (ESS): \n")
  print(x$ESS, digits=digits)

  cat("\nMCMC Acceptance Rate: ")
  dput(x$AcceptRate)
  
  cat("\nInitial Values (log scale): \n")
  print(x$InitVal)
  
  cat("\nCovariance Matrix used in the proposal distribution (transformed scale):\n")
  print(x$CovMat)
 
  }


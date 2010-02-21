chol.sigma <- function(vec)
  {
  c1 <- sqrt(vec[1])
  c2 <- vec[3]/c1
  c3 <- sqrt(vec[2] - vec[3]^2 / vec[1])
  return(c(c1,c2,c3))
  }

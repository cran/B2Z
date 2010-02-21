dmnorm_mod <- function (x, mean = rep(0, d), d, varcov, invvarcov) 
    {
    x <- matrix(x, 1, d)
    n <- 1
    X <- t(matrix(x, nrow = n, ncol = d)) - mean
    Q <- apply((invvarcov %*% X) * X, 2, sum)
    logDet <- sum(logb(abs(diag(qr(varcov)[[1]]))))
    logPDF <- as.vector(Q + d * logb(2 * pi) + logDet)/(-2)
    return(exp(logPDF))
    }
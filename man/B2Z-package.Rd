\name{B2Z-package}
\Rdversion{1.1}
\alias{B2Z-package}
\alias{B2Z}
\docType{package}
\title{
Bayesian Two-Zone Models (B2ZM)
}
\description{
The B2Z R-package fits the Bayesian two-zone model proposed by Zhang et al. (2009).  Denote \eqn{CN(t)} and \eqn{CF(t)} by the concentrations at time \eqn{t} in the near and far fields, respectively. The deterministic equations in the two-zone modeling are given by: 

\deqn{(d/dt)CN(t) = Beta/VN(CF(t)- CN(t) + G/Beta)          (1)}
\deqn{(d/dt)CF(t) = Beta/VF(CN(t) - (1 + Q/Beta)CF(t))      (2)}

where Beta is the interzonal air flow rate; Q is the supply and exhaust flow rate; and G is the contaminant emission rate. 

Denote \eqn{YN(t)} and \eqn{YF(t)} by the observed log concentrations at time \eqn{(t)} in the near and far fields, respectively. Zhang et al. (2009) assume that the joint probability distribution of \eqn{YN(t)} and \eqn{YF(t)} is a bivariate normal distribution where the mean is the solution of the differential equations (1) and (2), and the covariance matrix is composed by unknown parameters: \eqn{TauN} (variance of the measurement error at the near field), \eqn{TauNF} (covariance between the measurements at the near and far fields) and \eqn{TauF} (variance of the measurement error at the far field). If the near and far field measurements error processes are independent, then \eqn{TauNF = 0}. Both dependent (default) and independent models are provided in this R-package. 

Several options as prior distribution for \eqn{Beta}, \eqn{Q} and \eqn{G} are provided. In the independent model, the priors distributions for \eqn{TauN} and \eqn{TauF} are inverse gamma distributions; and in the dependent model the prior joint distribution of \eqn{TauN}, \eqn{TauNF} and \eqn{TauF} is the Inverse Wishart Distribution.

Four sample methods are available: Gibbs with Metropolis step, Incremental Mixture Importance Sampling (IMIS), Metropolis and Sampling Importance Resampling (SIR). In addition, 
estimation using the Bayesian central limit theorem (Laplace approximation) is also available (BCLT).  The user can fit the Bayesian two-zone model by using the function \code{B2ZM} and in this function
specifying which sampler method to use OR using one of the functions: 
\code{B2ZM_GIBBS}, \code{B2ZM_IMIS}, \code{B2ZM_METROP} and \code{B2ZM_SIR}. In any case, the output is valid as an input for the functions \code{summary} and \code{plot}.
}

\seealso{
\code{\link{B2ZM}, \link{B2ZM_BCLT}, \link{B2ZM_IMIS}, \link{B2ZM_GIBBS}, \link{B2ZM_METROP}, \link{B2ZM_SIR} }
}


\details{
\tabular{ll}{
Package: \tab B2Z\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2009-12-23\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
}
\author{
Joao Vitor Dias Monteiro, Sudipto Banerjee and Gurumurthy Ramachandran.

Maintainer: Joao Vitor Dias Monteiro <monte092@umn.edu>
}
\references{
Zhang, Y., Banerjee, S., Yang,R., Lungu,C. and Ramachandran, G. (2009).  Bayesian Modeling of Exposure and Airflow Using Two-Zone Models. \emph{The Annals of Occupational Hygiene},
\bold{53}, 409-424. \url{ http://www.biostat.umn.edu/~sudiptob/ResearchPapers/ZBYLR.pdf}   
}

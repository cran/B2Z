\name{B2ZM_BCLT}
\alias{B2ZM_BCLT}
\title{Bayesian Two-Zone Models: using Bayesian Central Limit Theorem (BCLT) - Laplace Approximation }
\description{
\code{B2ZM_BCLT} obtains random samples from the posterior distribution of the parameters and exposure concentrations for the Bayesian two-zone model proposed by Zhang et al. (2009) using the Bayesian Central Limit Theorem (Laplace Approximation). The user can choose whether the near and far field measurement error processes are dependent or not. In the independent model, 5 parameters are considered: 1) Beta: Interzonal air flow rate (m3); 2) Q: supply and exhaust flow rate (m3/min); 3) G: contaminant emission rate (mg/min); 4) Tau_N: variance of the measurement error at the near field; 5)Tau_F; variance of the measurement error at the far field. In the dependent model (default), one more parameter is considered: 6) Tau_NF: covariance between the measurements at the near and far field.  Any prior distribution for Beta, Q and G can be chosen. In the independent model, the prior distributions for Tau_N and Tau_F are inverse gamma distributions; in the dependent model, the prior joint distribution of Tau_N, Tau_NF and Tau_F is the Inverse Wishart Distribution (see the Details section for more information on the parameterization of these distributions). The output from \code{B2ZM_BCLT} is a list that belongs to the class \code{bclt}.  This output is valid as an input for the functions \code{summary} and \code{plot}.
}
\usage{

B2ZM_BCLT (data = NULL, priorBeta, priorQ, priorG, 
                      v, S, tauN.sh, tauN.sc, tauF.sh, tauF.sc,  
                      VN, VF, indep.model = FALSE, cred = 95,  
                      m = 7000, sample_size=2000, 
                      figures = list(save = FALSE, 
                      type = c("ps", "eps","pdf", "png", "jpg")))}
\arguments{
\item{data}{
	A 3-column matrix where the columns are time, concentrations at 	the near field, and concentrations at the far field, 	respectively. The time must be scaled in minutes (min), and the concentrations must be scaled in miligrams per cubic meter (mg/m3)
}
  \item{priorBeta}{
A string defining the prior distribution for the parameter Beta. To declare the prior distribution of Beta, use standard R nomenclature for probability distributions. For example, if the prior of Beta is a Uniform(0,20), declare it with "unif(0,20)"; if it is a Normal(0,1),
declare it with "norm(0,1)". DO NOT put an "d" or "r" in front the name of the distributions. The options are: "unif(a,b)", "gamma(a,b)", 
"exp(a)", "norm(a,b)", "t(a)", "weibull(a,b)", "f(a,b)", "chisq(a,b)",
"cauchy(a,b)" and "lnorm(a,b)".
}
  \item{priorQ}{
A string defining the prior distribution for Q (use the nomenclature as for \code{priorBeta}).
}
  \item{priorG}{
A string defining the prior distribution for G (use the nomenclature as for \code{priorBeta}).
}
  \item{v}{
Degrees of freedom for the Inverse Wishart distribution (prior joint distribution for Tau_N, Tau_NF and Tau_F in the dependent model).
}
  \item{S}{
A 2x2 positive definite matrix for the Inverse Wishart (prior joint distribution for Tau_N, Tau_NF and Tau_F in the dependent model).
}
  \item{tauN.sh}{
The shape parameter in the inverse gamma distribution (prior distribution for Tau_N in the independent model).
}
  \item{tauN.sc}{
The scalar parameter in the inverse gamma distribution (prior distribution for Tau_N in the independent model).
}
  \item{tauF.sh}{
The shape parameter in the inverse gamma distribution (prior distribution for Tau_F in the independent model).
}
  \item{tauF.sc}{
The scalar parameter in the inverse gamma distribution (prior distribution for Tau_F in the independent model).
}
  \item{VN}{
Volume of the near field in cubic meters (m3).
}
  \item{VF}{
Volume of the far field in cubic meters (m3).
}
  \item{indep.model}{
A logical value indicating whether the independent model should be considered. The default is FALSE.
}
  \item{cred}{
A scalar between 0 and 100 indicating the credibility level for the  posterior intervals of the parameters.
}
\item{m}{Number of sampling values from the prior distribution used to estimate a good starting value that is used in the estimation of the posterior mode and covariance matrix. See the Details section for more information.}
\item{sample_size}{Size of the sample from the posterior distribution of the parameters in model.}
 \item{figures}{
 The command \code{plot(obj)} produces several plots, where 
\code{obj} is the output from \code{B2ZM_BCLT}. Using \code{figures}, those plots are built internally and saved as eps, pdf, ps, png or jpg format. \code{figures} is a list containing the following parameters:
\describe{\item{\code{save}:}{a logical value indicating that the figures are to be saved. The default is FALSE.}
\item{\code{type}:}{a string that indicates the image file type. The default is "ps".}}
}
}
\details{


Parameterization priors:  The inverse gamma distribution with shape \emph{a} and scale \emph{b}  has mean \emph{b/(a-1)} (\emph{a}>1)  and variance \emph{(b^2)/((a-1)^2(a-2))} (\emph{a>2}). The inverse Wishart with \emph{v} degrees of freedom and scalar matrix \emph{S} has mean \emph{S/(v-p-1)}, where \emph{p} is the number of rows of \emph{S}.   


Covariance Matrix estimation: The covariance matrix is estimated as the negative inverse of the hessian matrix of the log posterior distribution at the estimated posterior mode. To estimate the posterior mode, the function \code{nlminb} is used. The values of the estimated posterior mode depends on the starting parameter values. \code{m} is the number of sampling values from the prior distributions of Beta, Q and G. The vector (among the m sampled) with largest log posterior value is used as starting parameter values. 
  
The covariance matrix is estimated using the function \code{hessian} from the package \code{numDeriv}, where the parameter vector is the estimated posterior mode. 
}
\value{
\code{B2ZM_BCLT} returns a list that belongs to the class \code{bclt}. The output from \code{B2ZM_BCLT} contains the objects:
\item{Beta}{a vector containing the sampled values from the joint posterior distribution for the parameter Beta.}
\item{Q}{a vector containing the sampled values from the joint posterior distribution for the parameter Q.}
\item{G}{a vector containing the sampled values from the joint posterior distribution for the parameter G.}
\item{tauN}{a vector containing the sampled values from the joint posterior distribution for the parameter Tau_N.}
\item{tauF}{a vector containing the sampled values from the joint posterior distribution for the parameter Tau_F.}
\item{tauNF}{a vector containing the sampled values from the joint posterior distribution for the parameter Tau_NF (if the dependent model is used).}
\item{Y}{a matrix containing the log of the observed concentrations.}
\item{DIC}{deviance information criterion.}
\item{pD}{effective number of parameters.}
\item{Dbar}{Deviance expectation.}
\item{indep}{a logical value indicating whether the independent model was used.}
\item{times}{a vector containing the times when the observed concentrations were measured.}
\item{cred}{credibility of the posterior intervals.}
 
Methods defined for B2ZM_BCLT object are \code{summary} and \code{plot}. 
}
\references{
   
Monteiro, J. V. D., Banerjee, S. and Ramachandran, G. (2011). B2Z: An R Package for Bayesian Two-Zone Models. \emph{Journal of Statistical Software} 43 (2) 1--23. \url{http://www.jstatsoft.org/v43/i02/}
  
Zhang, Y., Banerjee, S., Yang,R., Lungu,C. and Ramachandran, G. (2009).  Bayesian Modeling of Exposure and Airflow Using Two-Zone Models. \emph{The Annals of Occupational Hygiene},
\bold{53}, 409-424. \url{ http://www.biostat.umn.edu/~sudiptob/ResearchPapers/ZBYLR.pdf}
}
\author{
Joao Vitor Dias Monteiro, Sudipto Banerjee and Gurumurthy Ramachandran.
}
\seealso{
\code{\link{B2Z}, \link{B2ZM}, \link{B2ZM_IMIS},\link{B2ZM_MCMC}, \link{B2ZM_SIR} }
}
\examples{
#####################
#Independent Model #
###################

#Data 1:  100 simulated concentrations during the times between 

#0 and 4, using the parameters Beta = 5, Q = 13.8, G = 351.5, 

#VN = pi*10^-3 VF = 3.8, Tau_N = 1, Tau_NF = 0.5 

#and Tau_F = 0.64. 
\dontrun{

data(ex1)

r <- B2ZM_BCLT(data = ex1, priorBeta = "unif(0,10)",
               priorQ="unif(11,17)", priorG = "unif(281,482)", 
               S = diag(10,2), v = 4, VN = pi*10^-3, VF = 3.8, 
               m = 7000, sample_size=2000)

summary(r)
plot(r)


#####################
#Independent Model #
###################

#Data 2:  100 simulated concentrations during the times between 

#0 and 4, using the parameters Beta = 5, Q = 13.8, G = 351.5, 

#VN = pi*10^-3, VF = 3.8, Tau_N = 1, Tau_NF = 0 

#and Tau_F = 0.64. 

data(ex2)


r <- B2ZM_BCLT(data = ex2, priorBeta = "unif(0,10)",
               priorQ="unif(11,17)", priorG = "unif(281,482)",
               tauN.sh = 5 , tauN.sc = 4 , tauF.sh = 5,
               tauF.sc = 7, VN = pi*10^-3, VF = 3.8, 
               indep.model = TRUE, m = 7000, sample_size=2000)

summary(r)
plot(r)
}
}
\keyword{models}

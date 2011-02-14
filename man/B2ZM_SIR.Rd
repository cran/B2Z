\name{B2ZM_SIR}
\alias{B2ZM_SIR}
\title{Bayesian Two-Zone Models: using SIR sampler}
\description{
\code{B2ZM_SIR} obtains random samples from the posterior distribution of the parameters and exposure concentrations for the Bayesian two-zone model proposed by Zhang et al. (2009) using Sampling Importance Resampling (SIR). The user can choose whether the near and far field measurement error processes are dependent or not. In the independent model, 5 parameters are considered: 1) Beta: Interzonal air flow rate (m3); 2) Q: supply and exhaust flow rate (m3/min); 3) G: contaminant emission rate (mg/min); 4) Tau_N: variance of the measurement error at the near field; 5)Tau_F; variance of the measurement error at the far field. In the dependent model (default), one more parameter is considered: 6) Tau_NF: covariance between the measurements at the near and far field.  Any prior distribution for Beta, Q and G can be chosen. In the independent model, the prior distributions for Tau_N and Tau_F are inverse gamma distributions; in the dependent model, the prior joint distribution of Tau_N, Tau_NF and Tau_F is the Inverse Wishart Distribution (see the Details section for more information on the parameterization of these distributions). The output from \code{B2ZM_SIR} is a list that belongs to the class \code{sir}.  This output is valid as an input for the functions \code{summary} and \code{plot}.
}
\usage{

B2ZM_SIR(data, priorBeta, priorQ, priorG, 
        v, S, tauN.sh, tauN.sc,
        tauF.sh, tauF.sc,  VN,
        VF, indep.model = FALSE, 
        cred = 95,  m = 10000, 
        figures = list(save = FALSE, type =c("ps", 
        "eps","pdf", "png", "jpg"))) 
}
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
\item{m}{Number of samplings. The default is 10,000.}
 \item{figures}{
 The command \code{plot(obj)} produces several plots, where 
\code{obj} is the output from \code{B2ZM_SIR}. Using \code{figures}, those plots are built internally and saved as eps, pdf, ps, png or jpg format. \code{figures} is a list containing the following parameters:
\describe{\item{\code{save}:}{a logical value indicating that the figures are to be saved. The default is FALSE.}
\item{\code{type}:}{a string that indicates the image file type. The default is "ps".}}
}
}
\details{

Parameterization priors:  The inverse gamma distribution with shape \emph{a} and scale \emph{b}  has mean \emph{b/(a-1)} (\emph{a}>1)  and variance \emph{(b^2)/((a-1)^2(a-2))} (\emph{a>2}). The inverse Wishart with \emph{v} degrees of freedom and scalar matrix \emph{S} has mean \emph{S/(v-p-1)}, where \emph{p} is the number of rows of \emph{S}.   
}
\value{
\code{B2ZM_SIR} returns a list that belongs to the class \code{sir}. The output from \code{B2ZM_SIR} contains the objects:
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
\item{ESS}{effective sample size.}
\item{indep}{a logical value indicating whether the independent model was used.}
\item{times}{a vector containing the times when the observed concentrations were measured.}
\item{cred}{credibility of the posterior intervals.}
\item{prop}{proportion of different points in the sampled values from the joint posterior distribution.}
\item{weigths}{weights used in the SIR method.}
\item{maxw}{maximum weight used in the SIR method.}
 
Methods defined for B2ZM_SIR object are \code{summary} and \code{plot}. 
}
\references{

Zhang, Y., Banerjee, S., Yang,R., Lungu,C. and Ramachandran, G. (2009).  Bayesian Modeling of Exposure and Airflow Using Two-Zone Models. \emph{The Annals of Occupational Hygiene},
\bold{53}, 409-424. \url{ http://www.biostat.umn.edu/~sudiptob/ResearchPapers/ZBYLR.pdf}   

}
\author{
Joao Vitor Dias Monteiro, Sudipto Banerjee and Gurumurthy Ramachandran.
}
\seealso{
\code{\link{B2Z}, \link{B2ZM}, \link{B2ZM_BCLT}, \link{B2ZM_IMIS}, \link{B2ZM_MCMC}}
}
\note{
In the Examples section, \code{m} is small just for a quick demonstration. We suggest that \code{m} be greater than 15,000.
}
\examples{
##################
#Dependent Model#
################

#Data 1:  100 simulated concentrations during the times 
#between 0 and 4, using the parameters Beta = 5, Q = 13.8,
#G = 351.5, VN = pi*10^-3, VF = 3.8, Tau_N = 1, 
#Tau_NF = 0.5 and Tau_F = 0.64. 

data(ex1)

r <- B2ZM_SIR(data = ex1,  priorBeta = "unif(0,10)",   
     priorQ="unif(11,17)", priorG = "unif(281,482)", S = diag(10,2), 
     v = 4, VN = pi*10^-3, VF = 3.8, m = 100 )

#plot(r)
#summary(r)

#Saving figures with .png extension
\dontrun{
r <- B2ZM_SIR(data = ex1, priorBeta = "unif(0,10)",    
     priorQ = "unif(11,17)", priorG = "unif(281,482)", 
     S = diag(10,2), v = 4, VN = pi*10^-3, VF = 3.8, 
     m = 10000, figures = list(save = TRUE, type ="png"))
}


#####################
#Independent Model #
###################

#Data 2:  100 simulated concentrations during the times 
#between 0 and 4, using the parameters Beta = 5, Q = 13.8,
#G = 351.5, VN = pi*10^-3, VF = 3.8, Tau_N = 1, 
#Tau_NF = 0 and Tau_F = 0.64. 

\dontrun{
data(ex2)

r <- B2ZM_SIR(data = ex2, indep.model = TRUE, 
     priorBeta = "unif(0,10)", priorQ="unif(11,17)", 
     priorG = "unif(281,482)", tauN.sh = 5 , tauN.sc = 4 , 
     tauF.sh = 5, tauF.sc = 7 , VN = pi*10^-3, 
     VF = 3.8, m = 100)

plot(r)
summary(r)
}

}
\keyword{models}

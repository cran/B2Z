\name{B2ZM_BCLT}
\alias{B2ZM_BCLT}
\title{Bayesian Two-Zone Models: using Bayesian Central Limit Theorem (BCLT) - Laplace Approximation }
\description{
\code{B2ZM_BCLT} obtains random samples from the posterior distribution of the parameters and exposure concentrations for the Bayesian two-zone model proposed by Zhang et al. (2009) using the Bayesian Central Limit Theorem (Laplace Approximation). The user can choose whether the near and far field measurement error processes are independent or not. In the independent model, 5 parameters are considered: 1) Beta: Interzonal air flow rate; 2) Q: supply and exhaust flow rate; 3) G: contaminant emission rate; 4) Tau_N: variance of the measurement error at the near field; 5)Tau_F; variance of the measurement error at the far field. In the dependent model (default), one more parameter is considered: 6) Tau_NF: covariance between the measurements at the near and far field.  
Several prior distributions for Beta, Q and G are available. In the independent model, the prior distributions for Tau_N and Tau_F are inverse gamma distributions; in the dependent model, the prior joint distribution of Tau_N, Tau_NF and Tau_F is the Inverse Wishart Distribution (see the Details section for more information on the parameterization of these distributions). The output from \code{B2ZM_BCLT} is a list that belongs to the class \code{bclt}.  This output is valid as an input for the functions \code{summary} and \code{plot}.
}
\usage{

B2ZM_BCLT(func = func_default, y0, data = NULL, priorBeta, 
                 priorQ, priorG, v, S, tauN.sh, tauN.sc,
                 tauF.sh, tauF.sc,  V_N, V_F, indep.model = FALSE,  
                 credibility = 95,  m = 7000, sample_size=2000, 
                 figures = list(save = FALSE, type =c("ps", 
                 "eps","pdf", "png", "jpg"))) 
}
\arguments{
\item{func}{
A function that computes the derivatives in the dynamics of the total contaminant mass in a two-zone field at time t. The default for \code{func} contains the ordinary differential system  (ODS) presented in Zhang et al. (2009). To declare a different \code{func}, see the Details section.
}
\item{y0}{A vector of size 2, containing the initial concentrations (at time 0) in the near and far field, respectively. The default is c(0,0).
}
\item{data}{
A 3-column matrix where the columns are time, concentrations at the near field, and concentrations at the far field, respectively. 
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
  \item{V_N}{
Volume of the near field.
}
  \item{V_F}{
Volume of the far field.
}
  \item{indep.model}{
A logical value indicating whether the independent model should be considered. The default is FALSE.
}
  \item{credibility}{
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
Modifying \code{func}: The ordinary differential system component in the model is solved numerically using the R package \code{odesolve}. The declaration of \code{func} is the same as in the function \code{lsoda} of the \code{odesolve} package. The default of \code{func} is the ODS presented in Zhang et al. (2009);  it is coded in R as follows:

  \code{function(t, y, parms)}\cr
  \code{\{}\cr
    \code{ydot =  matrix(0,2,1)}\cr
    \code{ydot[1,1] = (parms[1]/parms[2])*(y[2]-y[1]) + parms[5]/parms[2]}\cr
    \code{ydot[2,1] = (1/parms[3])*(parms[1]*(y[1]-y[2])-y[2]*parms[4])
}\cr
    \code{return(list(ydot))}\cr
  \code{\}}\cr
where, Beta, V_N, V_F, Q and G are represented by \code{parms[1]}, \code{parms[2]}, \code{parms[3]},   \code{parms[4]} and \code{parms[5]}, respectively. \code{y[0]} and \code{y[1]} are the concentrations at the near and far fields at time \code{t}, respectively. \code{ydot[1,1]} and \code{ydot[2,1]} are the derivative with respect to \code{t} of the concentrations at the near and far fields, respectively.  To write different dynamics of the total contaminant mass in a two-zone field, vary the arithmetic of the terms in the function above.


Parameterization priors:  The inverse gamma distribution with shape \emph{a} and scale \emph{b}  has mean \emph{b/(a-1)} (\emph{a}>1)  and variance \emph{(b^2)/((a-1)^2(a-2))} (\emph{a>2}). The inverse Wishart with \emph{v} degrees of freedom and scalar matrix \emph{S} has mean \emph{S/(v-p-1)}, where \emph{p} is the number of rows of \emph{S}.   


Covariance Matrix estimation: The covariance matrix is estimated as the negative inverse of the hessian matrix of the log posterior distribution at the estimated posterior mode. To estimate the posterior mode, the function \code{nlminb} is used. The values of the estimated posterior mode depends on the starting parameter values. \code{m} is the number of sampling values from the prior distributions of Beta, Q and G. The vector (among the m sampled) with largest likelihood value is used as starting parameter values. 
  
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
\item{logCN}{a matrix such that the i-th row contains the log concentrations at the near field  when the i-th sampled values from the joint posterior distribution are considered.}
\item{logCF}{a matrix such that the i-th row contains the log concentrations at the far field  when the i-th sampled values from the joint posterior distribution are considered.}
\item{Y}{a matrix containing the log of the observed concentrations.}
\item{DIC}{deviance information criterion.}
\item{pD}{effective number of parameters.}
\item{Dbar}{Deviance expectation.}
\item{indep}{a logical value indicating whether the independent model was used.}
\item{y0}{a vector containing the initial concentrations for both zones.}
\item{times}{a vector containing the times when the observed concentrations were measured.}
\item{cred}{credibility of the posterior intervals.}
 
Methods defined for B2ZM_BCLT object are \code{summary} and \code{plot}. 
}
\references{
   
Zhang, Y., Banerjee, S., Yang,R., Lungu,C. and Ramachandran, G. (2009).  Bayesian Modeling of Exposure and Airflow Using Two-Zone Models. \emph{The Annals of Occupational Hygiene},
\bold{53}, 409-424. \url{ http://www.biostat.umn.edu/~sudiptob/ResearchPapers/ZBYLR.pdf}   

}
\author{
Joao Vitor Dias Monteiro, Sudipto Banerjee and Gurumurthy Ramachandran.
}
\seealso{
\code{\link{B2Z}, \link{B2ZM}, \link{B2ZM_IMIS},\link{B2ZM_GIBBS}, \link{B2ZM_METROP}, \link{B2ZM_SIR} }
}
\examples{
#####################
#Independent Model #
###################

#Data 1:  100 simulated concentrations during the times between 0 and 4, using the parameters Beta = 5, Q = 13.8, G = 351.5, VN = pi*10^-3 
#VF = 3.8, Tau_N = 1, Tau_NF = 0.5 and Tau_F = 0.64. 
\dontrun{

data(ex1)

r <- B2ZM_BCLT(data = ex1, y0 = c(0,0), priorBeta = "unif(0,10)", 	priorQ="unif(11,17)", priorG = "unif(281,482)", S = diag(10,2), 
      v = 4, V_N = pi*10^-3, V_F = 3.8, m = 7000, sample_size=2000)


summary(r)
plot(r)


#####################
#Independent Model #
###################

#Data 2:  100 simulated concentrations during the times between 0 and 4, using the parameters Beta = 5, Q = 13.8, G = 351.5, VN = pi*10^-3 
#VF = 3.8, Tau_N = 1, Tau_NF = 0 and Tau_F = 0.64. 

data(ex2)


r <- B2ZM_BCLT(data = ex2, y0 = c(0,0), priorBeta = "unif(0,10)", 	priorQ="unif(11,17)", 
               priorG = "unif(281,482)", tauN.sh = 5 , tauN.sc = 4 , tauF.sh = 5, tauF.sc = 7, V_N = pi*10^-3, 
               V_F = 3.8, indep.model = TRUE, m = 7000, sample_size=2000)

summary(r)
plot(r)
}
}
\keyword{models}

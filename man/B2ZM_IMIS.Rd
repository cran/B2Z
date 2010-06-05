\name{B2ZM_IMIS}
\alias{B2ZM_IMIS}
\title{Bayesian Two-Zone Models: using IMIS sampler}
\description{
\code{B2ZM_IMIS} obtains random samples from the posterior distribution of the parameters and exposure concentrations for the Bayesian two-zone model proposed by Zhang et al. (2009) using the Incremental Mixture Importance Sampling (IMIS). The user can choose whether the near and far field measurement error processes are independent or not. In the independent model, 5 parameters are considered: 1) Beta: Interzonal air flow rate; 2) Q: supply and exhaust flow rate; 3) G: contaminant emission rate; 4) Tau_N: variance of the measurement error at the near field; 5)Tau_F; variance of the measurement error at the far field. In the dependent model (default), one more parameter is considered: 6) Tau_NF: covariance between the measurements at the near and far field.  
Several prior distributions for Beta, Q and G are available. In the independent model, the prior distributions for Tau_N and Tau_F are inverse gamma distributions; in the dependent model, the prior joint distribution of Tau_N, Tau_NF and Tau_F is the Inverse Wishart Distribution (see the Details section for more information on the parameterization of these distributions). The output from 
\code{B2ZM_IMIS} is a list that belongs to the class \code{imis}.  This output is valid as an input for the functions \code{summary} and \code{plot}.
}
\usage{

B2ZM_IMIS(func = func_default, y0 = c(0,0), data, priorBeta, 
                 priorQ, priorG, v, S, tauN.sh, tauN.sc,
                 tauF.sh, tauF.sc,  V_N,
                 V_F, indep.model = FALSE, credibility = 95,
                 N0 = 6000, B = 600, M = 3000, it.max = 12, 
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
\item{N0}{ Initial number of inputs from the prior joint distribution of Beta, G, Q, Tau_N, Tau_F, and Tau_NF (if the dependent model is considered). The default is 6,000.}
\item{B}{ Number of inputs to be chosen in the Importance Sampling Stage. The default is 600.}
\item{M}{ Number of resamplings in the Resample Stage. The default is 3,000.}
\item{it.max}{Maximum number of iterations in the Importance Sampling Stage to be tolerated, in case the stop condition suggested by Raftery and Bao (2009) takes too many time. The default is 12.}
 \item{figures}{
 The command \code{plot(obj)} produces several plots, where 
\code{obj} is the output from \code{B2ZM_IMIS}. Using \code{figures}, those plots are built internally and saved as eps, pdf, ps, png or jpg format. \code{figures} is a list containing the following parameters:
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


Parameterization priors:  The inverse gamma and inverse Wishart distributions used in  \code{B2ZM_IMIS} are from the package \code{MCMCpack}.  The inverse gamma distribution with shape \emph{a} and scale \emph{b}  has mean \emph{b/(a-1)} (\emph{a}>1)  and variance \emph{(b^2)/((a-1)^2(a-2))} (\emph{a>2}). The inverse Wishart with \emph{v} degrees of freedom and scalar matrix \emph{S} has mean \emph{S/(v-p-1)}, where \emph{p} is the number of rows of \emph{S}.   

}
\value{
\code{B2ZM_IMIS} returns a list that belongs to the class \code{imis}. The output from \code{B2ZM_IMIS} contains the objects:

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
\item{ESS}{effective sample size.}
\item{indep}{a logical value indicating whether the independent model was used.}
\item{y0}{a vector containing the initial concentrations for both zones.}
\item{times}{a vector containing the times when the observed concentrations were measured.}
\item{cred}{credibility of the posterior intervals.}
\item{expfrac}{expected fraction of unique points.}
\item{V.hat}{variance of the rescaled importance weights.}
\item{U.hat}{entropy of importance weights relative to uniformity.}
\item{Q.hat}{expected number of unique points after re-sampling.}
\item{maxw}{maximum importance weight.}
\item{w}{importance weights.}

Methods defined for B2ZM_IMIS object are \code{summary} and \code{plot}. 
}
\references{
   
Zhang, Y., Banerjee, S., Yang,R., Lungu,C. and Ramachandran, G. (2009).  Bayesian Modeling of Exposure and Airflow Using Two-Zone Models. \emph{The Annals of Occupational Hygiene},
\bold{53}, 409-424.   \url{ http://www.biostat.umn.edu/~sudiptob/ResearchPapers/ZBYLR.pdf}

Raftery, A. E. and Bao, L. (2009).   Estimating and Projecting Trends in HIV/AIDS Generalized Epidemics Using Incremental Mixture Importance Sampling.  \url{http://www.stat.washington.edu/research/reports/2009/tr560.pdf}.
}
\note{
In the Examples section, \code{N0}, \code{B}, \code{M} and \code{it.max} are small just for a quick demonstration. We suggest that \code{N0} = 6,000; \code{B} = 600; \code{M}= 3,000  and \code{it.max} greater than 10.
}
\author{
Joao Vitor Dias Monteiro, Sudipto Banerjee and Gurumurthy Ramachandran.
}
\seealso{
\code{\link{B2Z}, \link{B2ZM},\link{B2ZM_BCLT},\link{B2ZM_GIBBS},\link{B2ZM_METROP},\link{B2ZM_SIR}} 
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

r <- B2ZM_IMIS(data = ex1, y0 = c(0,0), priorBeta = "unif(0,10)",         
     priorQ="unif(11,17)", priorG = "unif(281,482)", S = diag(10,2), 
     v = 4, V_N = pi*10^-3, V_F = 3.8, N0 = 500, B = 50, M = 250, 
     it.max = 4)

plot(r)
summary(r)

#Saving figures with .pdf extension
\dontrun{
r <- B2ZM_IMIS(data = ex1, y0 = c(0,0), priorBeta = "unif(0,10)", 
     priorQ = "unif(11,17)", priorG = "unif(281,482)", S = diag(10,2), 
     v = 4, V_N = pi*10^-3, V_F = 3.8, N0 = 6000, B = 600, M = 3000, 
     it.max = 12, figures = list(save = TRUE, type ="pdf"))
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

r <- B2ZM_IMIS(data = ex2, y0 = c(0,0), indep.model = TRUE, 
     priorBeta = "unif(0,10)", priorQ="unif(11,17)", 
     priorG = "unif(281,482)", tauN.sh = 5 , tauN.sc = 4 , 
     tauF.sh = 5, tauF.sc = 7 , V_N = pi*10^-3, V_F = 3.8, 
     N0 = 5000, B = 500, M =  3000,  it.max = 12)

plot(r)
summary(r)}


}
\keyword{models}

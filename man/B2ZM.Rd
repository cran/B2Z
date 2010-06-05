\name{B2ZM}
\alias{B2ZM}
\title{
  Bayesian Two-Zone Models
}
\description{
\code{B2ZM} obtains random samples from the posterior distribution of parameters and exposure concentrations for the Bayesian two-zone model proposed by Zhang et al. (2009). Four sample methods are available: Metropolis, Gibbs with Metropolis step, Incremental Mixture Importance Sampling (IMIS), and Sampling Importance Resampling (SIR). In addition, 
estimation using the Bayesian central limit theorem (Laplace approximation) is also available. The user can choose whether the near and far field measurement error processes are dependent or not. In the independent model, 5 parameters are considered: 1) Beta: Interzonal air flow rate; 2) Q: supply and exhaust flow rate; 3) G: contaminant emission rate; 4) Tau_N: variance of the measurement error at the near field; 5)Tau_F; variance of the measurement error at the far field. In the dependent model (default), one more parameter is considered: 6) Tau_NF: covariance between the measurements at the near and far field.  Any prior distribution for Beta, Q and G can be chosen. In the independent model, the prior distributions for Tau_N and Tau_F are inverse gamma distributions; in the dependent model, the prior joint distribution of Tau_N, Tau_NF and Tau_F is the Inverse Wishart Distribution (see the Details section for more information on the parameterization of these distributions). The output from \code{B2ZM} is a list that belongs to one of the following classes: \code{metrop}, \code{gibbs}, \code{imis} and \code{sir}.  This output is valid as an input for the functions \code{summary} and \code{plot}. 
}
\usage{
B2ZM(func = func_default, y0 = c(0,0), data = NULL, 
     priorBeta, priorQ, priorG, v, S, tauN.sh, tauN.sc,
     tauF.sh, tauF.sc,  V_N, V_F, indep.model = FALSE, 
     credibility = 95, sampler = c("IMIS", "METROP", "GIBBS", "SIR", "BCLT"),
     sir.control = list(m=50000),
     bclt.control = list(m = 7000, sample_size = 2000),
     metrop.control = list(NUpd = 10000, burnin = 1000, lag = 1,   
     initial=NULL, Sigma.Cand=NULL, m=5000),
     imis.control = list(N0=6000, B=600, M=3000, it.max=12),
     gibbs.control = list(NUpd = 10000, burnin = 1000, lag = 1, 
     initial=NULL, Sigma.Cand=NULL, numstp=100),
     figures = list(save = FALSE, type =c("ps", "eps","pdf", "png", 
     "jpg")))
}
\arguments{
   \item{func}{
	A function that computes the derivatives in the dynamics of the 	total contaminant mass in a two-zone field at time t. The default 	for \code{func} contains the ordinary differential system  (ODS) 	presented in Zhang et al. (2009). To declare a different \code{	func}, see the Details section.
	}
   \item{y0}{
      A vector of size 2, containing the initial concentrations (at 	time 0) in the near and far field, respectively. The default is c	(0,0).
	}
   \item{data}{
	A 3-column matrix where the columns are time, concentrations at 	the near field, and concentrations at the far field, 	respectively. 
	}
  \item{priorBeta}{
	A string defining the prior distribution for the parameter Beta. 	To declare the prior distribution of Beta, use standard R 	nomenclature for probability distributions. For example, if the 	prior of Beta is a Uniform(0,20), declare it with "unif(0,20)"; 	if it is a Normal(0,1), declare it with "norm(0,1)". DO NOT put 	an "d" or "r" in front the name of the distributions. The options 	are: "unif(a,b)", "gamma(a,b)", "exp(a)", "norm(a,b)", "t(a)", 	"weibull(a,b)", "f(a,b)", "chisq(a,b)", "cauchy(a,b)" and 
	"lnorm(a,b)".
	}
  \item{priorQ}{
	A string defining the prior distribution for Q (use the 	nomenclature as for \code{priorBeta}).
	}
  \item{priorG}{
	A string defining the prior distribution for G (use the 	nomenclature as for \code{priorBeta}).
	}
  \item{v}{
	Degrees of freedom for the Inverse Wishart distribution (prior 	joint distribution for Tau_N, Tau_NF and Tau_F in the dependent 	model).
	}
  \item{S}{
	A 2x2 positive definite matrix for the Inverse Wishart (prior 	joint distribution for Tau_N, Tau_NF and Tau_F in the dependent 	model).
	}
  \item{tauN.sh}{
	The shape parameter in the inverse gamma distribution (prior 	distribution for Tau_N in the independent model).
	}
  \item{tauN.sc}{
	The scalar parameter in the inverse gamma distribution (prior 	distribution for Tau_N in the independent model).
	}
  \item{tauF.sh}{
	The shape parameter in the inverse gamma distribution (prior 	distribution for Tau_F in the independent model).
	}
  \item{tauF.sc}{
	The scalar parameter in the inverse gamma distribution (prior 	distribution for Tau_F in the independent model).
	}
  \item{V_N}{
	Volume of the near field.
	}
  \item{V_F}{
	Volume of the far field.
	}
  \item{indep.model}{
	A logical value indicating whether the independent model should 	be considered. The default is FALSE.
	}
  \item{credibility}{
	A scalar between 0 and 100 indicating the credibility level for 	the  posterior intervals of the parameters.
	}
  \item{sampler}{
	A string indicating which sampler method should be used. The 	options are: "IMIS", "METROP", "GIBBS"  and "SIR". The default is 	"IMIS".
	}
 \item{sir.control}{
   	A list containing the SIR sampler specifications:
 	\describe{ 
	   \item{\code{m}:}{
		Number of samplings. The default is 50,000.
            }
         }
      }
 \item{bclt.control}{
   	A list containing the BCLT specifications:
 	\describe{ 
	   \item{\code{m}:}{
		Number of sampling values from the prior distribution used to estimate a good starting value that is used in the estimation of the posterior mode and covariance matrix. See the Details section for more information.
            }
	   \item{\code{sample_size}:}{
Size of the sample from the posterior distribution of the parameters in model.            }

         }
      }
  \item{metrop.control}{
   	A list containing the Metropolis sampler specifications:
   	\describe{
   	   \item{\code{NUpd}:}{Number of updates. The default is 10,000.}
         \item{\code{burnin}:}{Period of burn-in. The default is 	    1,000.}
         \item{\code{lag}:}{Thin interval. The default is 1.}
         \item{\code{initial}:}{A vector containing the initial values 	   for the parameters: Beta, Q, G, Tau_N and Tau_F (exactly in 	   this order) if the independet model is used.  If the dependent 	   model is considered, \code{initial} is a vector containing the 	   initial values for Beta, Q, G, Tau_N, Tau_F and Tau_NF. The 	   default is NULL, and in both cases, the initial points are the 	   estimated posterior mode.}
         \item{\code{Sigma.Cand}:}{If the independent model is            
         considered, \code{Sigma.Cand} is a 5x5 covariance matrix for 	   the normal multivariate proposal distribution of log(Beta),          
         log(Q), log(G), sqrt(Tau_N), sqrt(Tau_F);  if the dependent 
         model is considered,\code{Sigma.Cand} is a 6x6 covariance 
         matrix for the normal multivariate proposal distribution of 
         log(Beta), log(Q), log(G), sqrt(Tau_N), 
         sqrt(Tau_F - (Tau_NF)^2 / Tau_N) and Tau_NF/sqrt(Tau_N). In        
         both cases, \code{Sigma.Cand} must be a positive definite 
         matrix.  The default is NULL, and in this case, the covariance 
         matrix used is the negative inverse of the hessian matrix of 
         the log posterior distribution at the estimated mode   
         posterior.}
	  \item{\code{m}:}{Number of sampling values from the prior 
        distribution used to estimate a good starting value that is 
        used in the estimation of \code{Sigma.Cand}. It is used only if 
        \code{Sigma.Cand} is not declared. See the Details section for 
        more information.}
        }
      }
 \item{imis.control}{
  A list containing the IMIS sampler specifications:
   \describe{\item{\code{N0}:}{Initial number of inputs from the prior joint distribution of Beta, G, Q, Tau_N, Tau_F, and Tau_NF (if the dependent model is considered). The default is 6,000.}
   \item{\code{B}:}{ Number of inputs to be chosen in the Importance Sampling Stage. The default is 600.}
   \item{\code{M}:}{ Number of resamplings in the Resample Stage. The default is 3,000.}
   \item{\code{it.max}:}{Maximum number of iterations in the Importance Sampling Stage to be tolerated, in case the stop condition suggested by Raftery and Bao (2009) takes too many time. The default is 12.}}
For details see Raftery and Bao (2009).}
  \item{gibbs.control}{
  A list containing the Gibbs Sampling with Metropolis step sampler specifications:
   \describe{
   \item{\code{NUpd}:}{Number of updates. The default is 10,000.}
   \item{\code{burnin}:}{Period of burn-in. The default is 1,000.}
   \item{\code{lag}:}{Thin interval. The default is 1.}
   \item{\code{initial}:}{A vector containing the initial values for the parameters: Beta, Q and G (exactly in this order).}
\item{\code{Sigma.Cand}:}{A 3x3 covariance matrix for the normal multivariate proposal distribution of log(Beta), log(Q) and log(G). The default is NULL, and in this case, the covariance matrix used is the negative inverse of the hessian matrix of the log posterior distribution at the estimated posterior mode.}
\item{\code{m}:}{Number of sampling values from the prior distribution used to estimate a good starting value that is used in the estimation of \code{Sigma.Cand}. It is used only if \code{Sigma.Cand} is not declared. See the Details section for more information.}}
}
 \item{figures}{
 The command \code{plot(obj)} produces several plots, where 
\code{obj} is the output from \code{B2ZM}. Using \code{figures}, those plots are built internally and saved as eps, pdf, ps, png or jpg format. \code{figures} is a list containing the following parameters:
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


Parameterization priors:  The inverse gamma and inverse Wishart distributions used in  \code{B2ZM} are from the package \code{MCMCpack}.  The inverse gamma distribution with shape \emph{a} and scale \emph{b}  has mean \emph{b/(a-1)} (\emph{a}>1)  and variance \emph{(b^2)/((a-1)^2(a-2))} (\emph{a>2}). The inverse Wishart with \emph{v} degrees of freedom and scalar matrix \emph{S} has mean \emph{S/(v-p-1)}, where \emph{p} is the number of rows of \emph{S}.   


Covariance Matrix for the proposal distribution: If the covariance matrix for the multivariate normal poposal distribution is NULL for the Metropolis sampler the covariance matrix used is the negative inverse of the hessian matrix of the log posterior distribution at the estimated posterior mode. To estimate the posterior mode, the function \code{nlm} is used. The values of the estimated posterior mode depends on the starting parameter values. \code{m} is the number of sampling values from the prior distributions of Beta, Q and G. The vector (among the m sampled) with largest likelihood value is used as starting parameter values. 
  
The covariance matrix is estimated using the function \code{hessian} from the package \code{numDeriv}, where the parameter vector is 
the estimated posterior mode. 

Covariance Matrix estimation using BCLT: The covariance matrix is estimated as the negative inverse of the hessian matrix of the log posterior distribution at the estimated posterior mode. To estimate the posterior mode, the function \code{nlminb} is used. The values of the estimated posterior mode depends on the starting parameter values. \code{m} is the number of sampling values from the prior distributions of Beta, Q and G. The vector (among the m sampled) with largest likelihood value is used as starting parameter values. 

}
\value{
\code{B2ZM} returns a list that belongs to one of the four classes: \code{metrop}, \code{gibbs}, \code{imis} and \code{sir}. The output from \code{B2ZM} contains the objects:
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
\item{ESS}{effective sample size.}
\item{indep}{a logical value indicating whether the independent model was used.}
\item{y0}{a vector containing the initial concentrations for both zones.}
\item{times}{a vector containing the times when the observed concentrations were measured.}
\item{cred}{credibility of the posterior intervals.}
\item{prop}{proportion of different points in the sampled values from the joint posterior distribution (only for 'sir' class).}
\item{weigths}{weights used in the SIR method (only for 'sir' class).}
\item{maxw}{maximum weight used in the SIR method (only for 'sir' class).}
\item{expfrac}{expected fraction of unique points (only for 'imis' class)}
\item{V.hat}{variance of the rescaled importance weights (only for 'imis' class)}
\item{U.hat}{entropy of importance weights relative to uniformity (only for 'imis' class)}
\item{Q.hat}{expected number of unique points after re-sampling (only for 'imis' class)}
\item{maxw}{maximum importance weight (only for 'imis' class)}
\item{w}{importance weights (only for 'imis' class)}
\item{AR}{acceptance rate in the metropolis step (only for 'metrop' or 'gibbs' classes).}
\item{Sigma.Cand}{covariance matrix used in the proposal distribution (only for 'metrop' or 'gibbs' classes).}
\item{initial}{a vector containing the initial points (only for 'metrop' or 'gibbs' classes).}
\item{NUpd}{number of updates (only for 'metrop' or 'gibbs' classes).}
\item{burnin}{burn-in period (only for 'metrop' or 'gibbs' classes).}
\item{lag}{thin interval (only for 'metrop' or 'gibbs' classes).}
 
Methods defined for B2ZM objects are \code{summary} and \code{plot}. 
}
\references{
   
Zhang, Y., Banerjee, S., Yang,R., Lungu,C. and Ramachandran, G. (2009).  Bayesian Modeling of Exposure and Airflow Using Two-Zone Models. \emph{The Annals of Occupational Hygiene},
\bold{53}, 409-424. \url{ http://www.biostat.umn.edu/~sudiptob/ResearchPapers/ZBYLR.pdf}    

Raftery, A. E. and Bao, L. (2009).   Estimating and Projecting Trends in HIV/AIDS Generalized Epidemics Using Incremental Mixture Importance Sampling.  \url{http://www.stat.washington.edu/research/reports/2009/tr560.pdf}.
}
\author{
Joao Vitor Dias Monteiro, Sudipto Banerjee and Gurumurthy Ramachandran.
}
\seealso{
\code{\link{B2Z}, \link{B2ZM_GIBBS}, \link{B2ZM_BCLT}, \link{B2ZM_IMIS}, \link{B2ZM_METROP}, \link{B2ZM_SIR} }
}
\examples{

##################
#Dependent Model#
################

#Data 1:  100 simulated concentrations during the times 
#between 0 and 4, using the parameters Beta = 5, Q = 13.8, 
#G = 351.5, VN = pi*10^-3, VF = 3.8, Tau_N = 1, 
#Tau_NF = 0.5 and Tau_F = 0.64. 
\dontrun{
data(ex1)

#######
#BCLT#
#####

r <- B2ZM(data = ex1, y0 = c(0,0), priorBeta = "unif(0,10)", 	
         priorQ="unif(11,17)", priorG = "unif(281,482)", S = diag(10,2),  
         v = 4, V_N = pi*10^-3, V_F = 3.8, sampler = "BCLT", 
         bclt.control=list( m = 7000, sample_size=2000))


summary(r)
plot(r)


########
#GIBBS#
######

r <- B2ZM(data = ex1, y0 = c(0,0), priorBeta = "unif(0,10)",   
     priorQ="unif(11,17)", priorG = "unif(281,482)", S = diag(10,2), 
     v = 4, V_N = pi*10^-3, V_F = 3.8, sampler = "GIBBS",     
     gibbs.control = list(NUpd = 10000, burnin = 1000, 
     lag = 1, m = 5000) )

summary(r)
plot(r)


#######
#IMIS#
#####

r <- B2ZM(data = ex1, y0 = c(0,0), priorBeta = "unif(0,10)",  
     priorQ="unif(11,17)", priorG = "unif(281,482)", S = diag(10,2), 
     v = 4, V_N = pi*10^-3, V_F = 3.8, sampler="IMIS", 
     imis.control = list( N0 = 6000, B = 600, M = 3000, it.max = 12))

summary(r)
plot(r)

#########
#METROP#
#######

r <- B2ZM(data = ex1, y0 = c(0,0), priorBeta = "unif(0,10)",  
     priorQ="unif(11,17)", priorG = "unif(281,482)", S = diag(10,2), 
     v = 4, V_N = pi*10^-3, V_F = 3.8, sampler = "METROP", 
     metrop.control = list(NUpd = 10000, burnin = 1000, 
     lag = 1, m = 5000) )

summary(r)
plot(r)


######
#SIR#
####

r <- B2ZM(data = ex1, y0 = c(0,0), priorBeta = "unif(0,10)",  
     priorQ="unif(11,17)", priorG = "unif(281,482)", S = diag(10,2), 
     v = 4, V_N = pi*10^-3, V_F = 3.8, sampler="SIR", 
     sir.control = list(m = 10000) )

plot(r)
summary(r)

#######################################################################

#####################
#Independent Model #
###################

#Data 2:  100 simulated concentrations during the times 
#between 0 and 4, using the parameters Beta = 5, Q = 13.8, 
#G = 351.5, VN = pi*10^-3, VF = 3.8, Tau_N = 1, 
#Tau_NF = 0 and Tau_F = 0.64. 

data(ex2)

#######
#BCLT#
#####

r <- B2ZM(data = ex2, y0 = c(0,0), priorBeta = "unif(0,10)", 	
     priorQ="unif(11,17)", priorG = "unif(281,482)", 
     tauN.sh = 5 , tauN.sc = 4 , tauF.sh = 5, tauF.sc = 7, 
     V_N = pi*10^-3,  V_F = 3.8, sampler = "BCLT", 
     indep.model = TRUE, bclt.control=list(m = 7000,  
     sample_size=2000))

summary(r)
plot(r)


########
#GIBBS#
######

r <- B2ZM(data = ex2, y0 = c(0,0), indep.model = TRUE, 
     priorBeta = "unif(0,10)", priorQ="unif(11,17)", 
     priorG = "unif(281,482)", tauN.sh = 5 , tauN.sc = 4 , tauF.sh = 5, 
     tauF.sc = 7 , V_N = pi*10^-3, V_F = 3.8, sampler = "GIBBS",
     gibbs.control = list(NUpd = 10000,  burnin = 1000, lag = 1, 
     m = 10000))

summary(r)
plot(r)


#######
#IMIS#
#####

r <- B2ZM(data = ex2, y0 = c(0,0), indep.model = TRUE, 
     priorBeta = "unif(0,10)", priorQ="unif(11,17)", 
     priorG = "unif(281,482)", tauN.sh = 5 , tauN.sc = 4 , tauF.sh = 5, 
     tauF.sc = 7 , V_N = pi*10^-3, V_F = 3.8, sampler = "IMIS",
     imis.control = list(N0 = 5000, B = 500, M =  3000,  it.max = 12))

summary(r)
plot(r)

#########
#METROP#
#######

r <- B2ZM(data = ex2, y0 = c(0,0), indep.model = TRUE, 
     priorBeta = "unif(0,10)", priorQ="unif(11,17)", 
     priorG = "unif(281,482)", tauN.sh = 5 , tauN.sc = 4 , tauF.sh = 5, 
     tauF.sc = 7 , V_N = pi*10^-3, V_F = 3.8, sampler = "METROP",
     metrop.control = list(NUpd = 10000,  burnin = 1000, lag = 1, 
     m = 10000))

summary(r)
plot(r)


######
#SIR#
####

r <- B2ZM(data = ex2, y0 = c(0,0), indep.model = TRUE, 
     priorBeta = "unif(0,10)", priorQ="unif(11,17)", 
     priorG = "unif(281,482)", tauN.sh = 5 , tauN.sc = 4 , tauF.sh = 5, 
     tauF.sc = 7 , V_N = pi*10^-3, V_F = 3.8, sampler = "SIR",
     sir.control = list(m = 10000))

plot(r)
summary(r)
}



}
\keyword{models}

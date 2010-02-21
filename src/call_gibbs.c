#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts

#include <Rdefines.h>
#include <Rinternals.h>

#include <Rmath.h>
#include "pacoteR_utils.h"



void F77_NAME(loglikelihood)(double *, double *, double *, double *,
                             double *, int *, int *, double *);

SEXP call_gibbs(SEXP y, SEXP times, SEXP func, SEXP rtol,
		SEXP atol, SEXP rho, SEXP tcrit, SEXP jacfunc, SEXP initfunc,
		SEXP verbose, SEXP hmin, SEXP hmax,    
        SEXP NUpd_R, SEXP index_R,SEXP index2_R, SEXP two_n_R, SEXP n_R,
        SEXP indep_R, SEXP Y_R, SEXP modY_R, SEXP VN_R, SEXP VF_R, 
        SEXP indBeta_R, SEXP aBeta_R, SEXP bBeta_R,
        SEXP indQ_R, SEXP aQ_R, SEXP bQ_R,
        SEXP indG_R, SEXP aG_R, SEXP bG_R,
        SEXP v_R, SEXP S_R, SEXP tauN_sh_R, SEXP tauN_sc_R, SEXP tauF_sh_R, 
        SEXP tauF_sc_R, SEXP retval_R, SEXP initial_R, SEXP total_it_R){
     
     SEXP r, ans, parms, logCNCF, Betaout, Qout, Gout, tauNout, tauFout, tauNFout, loglikout;
     int NUpd, indBeta, indQ, indG ;
     int i, j, k, indep, n, two_n, *index, *index2, total_it;
     double aBeta, bBeta, aQ, bQ, aG, bG;
     double *Beta, *Q, *G, *tauN, *tauF, *tauNF, *loglik,  loglik1;
     double *Y, *Ytild, *canYtild, *modY, *mu, *sigma, VN, VF;
     double v, *S, *S1, tauN_sh, tauN_sc, tauF_sh, tauF_sc, taus[4];
     double *Mu_Cand, *retval, logprior, current, candidate, ratio, u;
     double tauN_sh_S, tauN_sc_S, tauF_sh_S,tauF_sc_S, v1, sum1, sum2, sum3;
     double cBeta, cQ, cG;

     NUpd = INTEGER(NUpd_R)[0];       
     two_n = INTEGER(two_n_R)[0];
     n = INTEGER(n_R)[0];
     total_it = INTEGER(total_it_R)[0];

     double cand[3];
     S = (double *) R_alloc(4, sizeof(double));
     S1 = (double *) R_alloc(4, sizeof(double));
     sigma = (double *) R_alloc(4, sizeof(double));
     mu = (double *) R_alloc(two_n, sizeof(double));
     Y = (double *) R_alloc(two_n, sizeof(double));
     Ytild = (double *) R_alloc(two_n, sizeof(double));
     canYtild = (double *) R_alloc(two_n, sizeof(double));
     modY = (double *) R_alloc(two_n, sizeof(double));
     Mu_Cand = (double *) R_alloc(3, sizeof(double));
     retval = (double *) R_alloc((int)9, sizeof(double));
     
     index = (int *) R_alloc(two_n, sizeof(int));
     index2 = (int *) R_alloc(two_n, sizeof(int));

     indep = INTEGER(indep_R)[0];
     VN = REAL(VN_R)[0];
     VF = REAL(VF_R)[0];

     v = REAL(v_R)[0];
     for(i=0;i<4;i++){S[i] = REAL(S_R)[i];}
     tauN_sh = REAL(tauN_sh_R)[0];  tauN_sc = REAL(tauN_sc_R)[0];
     tauF_sh = REAL(tauF_sh_R)[0];  tauF_sc = REAL(tauF_sc_R)[0];
     indBeta = INTEGER(indBeta_R)[0]; aBeta = REAL(aBeta_R)[0]; 
     bBeta = REAL(bBeta_R)[0];
     indQ = INTEGER(indQ_R)[0]; aQ = REAL(aQ_R)[0]; bQ = REAL(bQ_R)[0];
     indG = INTEGER(indG_R)[0]; aG = REAL(aG_R)[0]; bG = REAL(bG_R)[0];
          
     for(i = 0; i < two_n; i++ )
       {
       Y[i] = REAL(Y_R)[i];
       modY[i] = REAL(modY_R)[i];
       index[i] = INTEGER(index_R)[i];
       index2[i] = INTEGER(index2_R)[i];
       }
 
     for(i =0; i < 9; i++)
       retval[i] = REAL(retval_R)[i];

     PROTECT(ans = allocVector(VECSXP, (NUpd + 7) )); 
     PROTECT(Betaout = NEW_NUMERIC(NUpd)); PROTECT(Qout = NEW_NUMERIC(NUpd)); 
     PROTECT(Gout = NEW_NUMERIC(NUpd)); PROTECT(tauNout = NEW_NUMERIC(NUpd)); 
     PROTECT(tauFout = NEW_NUMERIC(NUpd));PROTECT(tauNFout = NEW_NUMERIC(NUpd)); 
     PROTECT(loglikout = NEW_NUMERIC(NUpd)); 

     Beta = NUMERIC_POINTER(Betaout); Q = NUMERIC_POINTER(Qout); 
     G = NUMERIC_POINTER(Gout); tauN = NUMERIC_POINTER(tauNout); 
     tauF = NUMERIC_POINTER(tauFout); tauNF = NUMERIC_POINTER(tauNFout);   
     loglik = NUMERIC_POINTER(loglikout);
    

//INITIAL VALUES

     Beta[0] = cBeta = REAL(initial_R)[0];
     Q[0] = cQ =REAL(initial_R)[1];
     G[0] = cG =REAL(initial_R)[2];

     parms = NEW_NUMERIC(5);
     REAL(parms)[0] = Beta[0];
     REAL(parms)[1] = VN;
     REAL(parms)[2] = VF;
     REAL(parms)[3] = Q[0];
     REAL(parms)[4] = G[0];

     r = call_lsoda(y, times, func, parms, rtol, 
           atol, rho, tcrit, jacfunc, initfunc,
		   verbose, hmin, hmax);

     logCNCF = NEW_NUMERIC(two_n);
     for(j=0; j < two_n; j++)
        {
        REAL(logCNCF)[j] = mu[j] = log(REAL(r)[index[j]]); 
        Ytild[j] = log(REAL(r)[index2[j]]);
        }

     if(indep==1)
       {
       sum1  = sum2 = 0;
       for(k=0; k < n; k++)
         {
         sum1 = sum1 + (Y[k] - Ytild[k])*(Y[k] - Ytild[k]);
         sum2 = sum2 + (Y[n+k] - Ytild[n+k])*(Y[n+k] - Ytild[n+k]);
         }
                 
       tauN_sh_S = tauN_sh + n/(double)2;
       tauN_sc_S = tauN_sc + 0.5*sum1;

       tauF_sh_S = tauF_sh +  n/(double)2;
       tauF_sc_S = tauF_sc + 0.5*sum2;

       tauN[0] = 1.0/rpriors(2,tauN_sh_S,tauN_sc_S);
       tauF[0] = 1.0/rpriors(2,tauF_sh_S,tauF_sc_S);

       sigma[1] = sigma[2] = 0;
       sigma[0] = tauN[0];
       sigma[3] = tauF[0];
       }
     else
       {
       v1 = v + n;
       sum1  = sum2 = sum3 = 0;
       for(k=0; k < n; k++)
         {
         sum1 = sum1 + (Y[k] - Ytild[k])*(Y[k] - Ytild[k]);
         sum2 = sum2 + (Y[k] - Ytild[k])*(Y[n+k] - Ytild[n+k]);
         sum3 = sum3 + (Y[n+k] - Ytild[n+k])*(Y[n+k] - Ytild[n+k]);
         }
       
       S1[0] = S[0] + sum1;
       S1[1] = S1[2] = S[1] + sum2;
       S1[3] = S[3] + sum3;
        
       rinvwish(v1,S1, taus);
       tauN[0]  =  taus[0];
       tauF[0]  =  taus[3];
       tauNF[0] =  taus[1];
          
       sigma[1] = sigma[2] = tauNF[0];
       sigma[0] = tauN[0];
       sigma[3] = tauF[0];
       }

     logprior = dpriors(Beta[0], indBeta, aBeta, bBeta, (int)1);
     logprior = logprior + dpriors(Q[0], indQ, aQ, bQ, (int)1) + dpriors(G[0], indG, aG, bG, (int)1);
   
     loglik1 = 0;
     F77_CALL(loglikelihood) (mu, sigma, Y, Ytild, modY, &n, &two_n, &loglik1);
     loglik[0] = loglik1;  
     
     current = logprior + loglik1;
     SET_VECTOR_ELT(ans, 0, logCNCF);
     
      
//Gibbs with Metropolis step updates
    double control_bar, next_control_bar;
     next_control_bar = control_bar = 1/(double)50;
     
     for(i = 1; i < NUpd; i++)
       {
       logprior = 0;

       for(j = 0; j < 3; j++)
         cand[j]=0;
       
       Mu_Cand[0] = log(cBeta);
       Mu_Cand[1] = log(cQ);
       Mu_Cand[2] = log(cG);
        
       rmtvnorm(3, Mu_Cand, retval, cand);
      
       parms = NEW_NUMERIC(5);
       REAL(parms)[0] = exp(cand[0]);
       REAL(parms)[1] = VN;
       REAL(parms)[2] = VF;
       REAL(parms)[3] = exp(cand[1]);
       REAL(parms)[4] = exp(cand[2]);

       r = call_lsoda(y, times, func, parms, rtol, 
           atol, rho, tcrit, jacfunc, initfunc,
		   verbose, hmin, hmax);

       for(j=0; j < two_n; j++)
          {
          mu[j] = log(REAL(r)[index[j]]); 
          canYtild[j] = log(REAL(r)[index2[j]]);
          }
  
       logprior =  dpriors(exp(cand[0]), indBeta, aBeta, bBeta, (int)1);
       logprior = logprior + dpriors(exp(cand[1]), indQ, aQ, bQ, (int)1) + dpriors(exp(cand[2]), indG, aG, bG, (int)1);

   
       loglik1 = 0;
       F77_CALL(loglikelihood) (mu, sigma, Y, canYtild, modY, &n, &two_n, &loglik1);
  
       candidate = logprior + loglik1;
       ratio = exp((candidate - current));  
          
      if(ratio >= 1)
       {
       Beta[i] = exp(cand[0]);
       Q[i] = exp(cand[1]);
       G[i] = exp(cand[2]);
       
       logCNCF = NEW_NUMERIC(two_n);
       for(j=0; j < two_n; j++)
         {REAL(logCNCF)[j] = mu[j]; 
         Ytild[j]=canYtild[j];}
       
       loglik[i] = loglik1;
       current = candidate;
       }
      else
       {
       GetRNGstate();
       u = runif(0,1);
       PutRNGstate();

       if(u < ratio)
         {
         Beta[i] = exp(cand[0]);
         Q[i] = exp(cand[1]);
         G[i] = exp(cand[2]);
       
         logCNCF = NEW_NUMERIC(two_n);
         for(j=0; j < two_n; j++)
           {REAL(logCNCF)[j] = mu[j]; 
            Ytild[j]=canYtild[j];}
       
         loglik[i] = loglik1;
         current = candidate;
         }
        else
         {
         Beta[i] = REAL(Betaout)[i-1];
         Q[i] = REAL(Qout)[i-1];
         G[i] = REAL(Gout)[i-1];
         loglik[i] = REAL(loglikout)[i-1];
         } 
       }


     if(indep==1)
       {
       sum1  = sum2 = 0;
       for(k=0; k < n; k++)
         {
         sum1 = sum1 + (Y[k] - Ytild[k])*(Y[k] - Ytild[k]);
         sum2 = sum2 + (Y[n+k] - Ytild[n+k])*(Y[n+k] - Ytild[n+k]);
         }
                 
       tauN_sh_S = tauN_sh + n/(double)2;
       tauN_sc_S = tauN_sc + 0.5*sum1;

       tauF_sh_S = tauF_sh +  n/(double)2;
       tauF_sc_S = tauF_sc + 0.5*sum2;

       tauN[i] = 1.0/rpriors(2,tauN_sh_S,tauN_sc_S);
       tauF[i] = 1.0/rpriors(2,tauF_sh_S,tauF_sc_S);

       sigma[1] = sigma[2] = 0;
       sigma[0] = tauN[i];
       sigma[3] = tauF[i];
       }
     else
       {
       v1 = v + n;
       sum1  = sum2 = sum3 = 0;
       for(k=0; k < n; k++)
         {
         sum1 = sum1 + (Y[k] - Ytild[k])*(Y[k] - Ytild[k]);
         sum2 = sum2 + (Y[k] - Ytild[k])*(Y[n+k] - Ytild[n+k]);
         sum3 = sum3 + (Y[n+k] - Ytild[n+k])*(Y[n+k] - Ytild[n+k]);
         }
       
        S1[0] = S[0] + sum1;
        S1[1] = S1[2] = S[1] + sum2;
        S1[3] = S[3] + sum3;
        
        rinvwish(v1,S1,taus);
        tauN[i]  =  taus[0];
        tauF[i]  =  taus[3];
        tauNF[i] =  taus[1];
          
        sigma[1] = sigma[2] = tauNF[i];
        sigma[0] = tauN[i];
        sigma[3] = tauF[i];
       }


       SET_VECTOR_ELT(ans, i, logCNCF);
       
       if(i/(double)total_it >= next_control_bar)
         {
         Rprintf("=");
         R_FlushConsole();
         next_control_bar = next_control_bar + control_bar;
         }
       }

     SET_VECTOR_ELT(ans, NUpd, Betaout);
     SET_VECTOR_ELT(ans, (NUpd + 1), Qout);
     SET_VECTOR_ELT(ans, (NUpd + 2), Gout);
     SET_VECTOR_ELT(ans, (NUpd + 3), tauNout);
     SET_VECTOR_ELT(ans, (NUpd + 4), tauFout);
     SET_VECTOR_ELT(ans, (NUpd + 5), tauNFout);
     SET_VECTOR_ELT(ans, (NUpd + 6), loglikout);

     UNPROTECT(8);
     return (ans);
     }



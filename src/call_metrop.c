#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts

#include <Rdefines.h>
#include <Rinternals.h>

#include <Rmath.h>
#include "pacoteR_utils.h"


void F77_NAME(loglikelihood)(double *, double *, double *, double *,
                             double *, int *, int *, double *);

SEXP call_metrop(SEXP y, SEXP times, SEXP func, SEXP rtol,
		SEXP atol, SEXP rho, SEXP tcrit, SEXP jacfunc, SEXP initfunc,
		SEXP verbose, SEXP hmin, SEXP hmax,    
        SEXP NUpd_R, SEXP index_R,SEXP index2_R, SEXP two_n_R, SEXP n_R,
        SEXP indep_R, SEXP Y_R, SEXP modY_R, SEXP VN_R, SEXP VF_R, 
        SEXP indBeta_R, SEXP aBeta_R, SEXP bBeta_R,
        SEXP indQ_R, SEXP aQ_R, SEXP bQ_R,
        SEXP indG_R, SEXP aG_R, SEXP bG_R,
        SEXP v_R, SEXP S_R, SEXP tauN_sh_R, SEXP tauN_sc_R, SEXP tauF_sh_R, 
        SEXP tauF_sc_R, SEXP retval_R, SEXP initial_R, SEXP d_R, SEXP total_it_R){
     
     SEXP r, ans, parms, logCNCF, Betaout, Qout, Gout, tauNout, tauFout, tauNFout, loglikout;
     int NUpd, indBeta, indQ, indG ;
     int i, j, indep, n, two_n, d, *index, *index2, total_it;
     double aBeta, bBeta, aQ, bQ, aG, bG;
     double *Beta, *Q, *G, *tauN, *tauF, *tauNF, *loglik,  loglik1;
     double *Y, *Ytild, *modY, *mu, *sigma, VN, VF;
     double v, *S, tauN_sh, tauN_sc, tauF_sh, tauF_sc;
     double *Mu_Cand, *retval, logprior, current, candidate, ratio, u;

     NUpd = INTEGER(NUpd_R)[0];       
     two_n = INTEGER(two_n_R)[0];
     n = INTEGER(n_R)[0];
     d = INTEGER(d_R)[0];
     total_it = INTEGER(total_it_R)[0];


     double cand[d];
     S = (double *) R_alloc(4, sizeof(double));
     sigma = (double *) R_alloc(4, sizeof(double));
     mu = (double *) R_alloc(two_n, sizeof(double));
     Y = (double *) R_alloc(two_n, sizeof(double));
     Ytild = (double *) R_alloc(two_n, sizeof(double));
     modY = (double *) R_alloc(two_n, sizeof(double));
     Mu_Cand = (double *) R_alloc(d, sizeof(double));
     retval = (double *) R_alloc((int)d*d, sizeof(double));
     
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
 
     for(i =0; i < d*d; i++)
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
    
    double cBeta, cQ, cG, ctauN, ctauF, ctauNF = 0.0;
//INITIAL VALUES

     Beta[0] = cBeta = REAL(initial_R)[0];
     Q[0] = cQ =REAL(initial_R)[1];
     G[0] = cG =REAL(initial_R)[2];
     tauN[0] = ctauN = REAL(initial_R)[3];
     tauF[0] = ctauF = REAL(initial_R)[4];

     sigma[0] = tauN[0];
     sigma[3] = tauF[0];

     if(indep==1)
       {
       sigma[1] = sigma[2] = 0;
       logprior =  tauN_sh * log(tauN_sc) - lgammafn(tauN_sh) - (tauN_sh + 1) * log(sigma[0]) - (tauN_sc/sigma[0]);
       logprior += tauF_sh * log(tauF_sc) - lgammafn(tauF_sh) - (tauF_sh + 1) * log(sigma[3]) - (tauF_sc/sigma[3]);
       }
     else
       {
       tauNF[0] = ctauNF = REAL(initial_R)[5];
       sigma[1] = sigma[2] = tauNF[0];
       logprior = log(dinvwish(sigma, v, S));
       }

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

       logprior = logprior + dpriors(Beta[0], indBeta, aBeta, bBeta, (int)1);
       logprior = logprior + dpriors(Q[0], indQ, aQ, bQ, (int)1) + dpriors(G[0], indG, aG, bG, (int)1);
   
     loglik1 = 0;
     F77_CALL(loglikelihood) (mu, sigma, Y, Ytild, modY, &n, &two_n, &loglik1);
     loglik[0]=loglik1;  
     
     current = logprior + loglik1;
     SET_VECTOR_ELT(ans, 0, logCNCF);
 
 

//Metropolis Updates         
     double control_bar, next_control_bar;
     next_control_bar = control_bar = 1/(double)50;

     for(i = 1; i < NUpd; i++)
       {
       logprior = 0;

       for(j = 0; j <d; j++)
         cand[j]=0;
       
       if(indep==1)
         {
         Mu_Cand[0] = log(cBeta);
         Mu_Cand[1] = log(cQ);
         Mu_Cand[2] = log(cG);
         Mu_Cand[3] = sqrt(ctauN);
         Mu_Cand[4] = sqrt(ctauF); 
         
         rmtvnorm(5, Mu_Cand, retval, cand);
        
         sigma[0] = cand[3]*cand[3];
         sigma[1] = sigma[2] = 0;
         sigma[3] = cand[4]*cand[4];
         
         logprior =  tauN_sh * log(tauN_sc) - lgammafn(tauN_sh) - (tauN_sh + 1) * log(sigma[0]) - (tauN_sc/sigma[0]);
         logprior = logprior + tauF_sh * log(tauF_sc) - lgammafn(tauF_sh) - (tauF_sh + 1) * log(sigma[3]) - (tauF_sc/sigma[3]);
         }
       else
         {
         Mu_Cand[0] = log(cBeta);
         Mu_Cand[1] = log(cQ);
         Mu_Cand[2] = log(cG);
         Mu_Cand[3] = sqrt(ctauN);
         Mu_Cand[4] = sqrt((ctauF - ctauNF*ctauNF / ctauN)); 
         Mu_Cand[5] = ctauNF/sqrt(ctauN);

         rmtvnorm(6, Mu_Cand, retval, cand);
   
         sigma[0] = cand[3]*cand[3];
         sigma[1] = sigma[2] = cand[3]*cand[5];
         sigma[3] = cand[5]*cand[5] + cand[4]*cand[4];
         
         logprior = log(dinvwish(sigma, v, S));
         }

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
          Ytild[j] = log(REAL(r)[index2[j]]);
          }
  
      logprior = logprior + dpriors(exp(cand[0]), indBeta, aBeta, bBeta, (int)1);
      logprior = logprior + dpriors(exp(cand[1]), indQ, aQ, bQ, (int)1) + dpriors(exp(cand[2]), indG, aG, bG, (int)1);
   
      loglik1 = 0;
      F77_CALL(loglikelihood) (mu, sigma, Y, Ytild, modY, &n, &two_n, &loglik1);
  
      candidate = logprior + loglik1;
      ratio = exp((candidate - current));  

          
     if(ratio >= 1)
       {
       Beta[i] = exp(cand[0]);
       Q[i] = exp(cand[1]);
       G[i] = exp(cand[2]);
       tauN[i] = cand[3]*cand[3];
              
       if(indep==1)
         {tauF[i] = cand[4]*cand[4];}
       else
         {
         tauF[i] = cand[5]*cand[5] + cand[4]*cand[4];
         tauNF[i] = cand[3]*cand[5];
         }
       
       logCNCF = NEW_NUMERIC(two_n);
       for(j=0; j < two_n; j++)
         REAL(logCNCF)[j] = mu[j]; 
       
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
          tauN[i] = cand[3]*cand[3];
              
          if(indep==1)
           {tauF[i] = cand[4]*cand[4];}
          else
           {
           tauF[i] = cand[5]*cand[5] + cand[4]*cand[4];
           tauNF[i] = cand[3]*cand[5];
           }
       
          logCNCF = NEW_NUMERIC(two_n);
          for(j=0; j < two_n; j++)
            REAL(logCNCF)[j] = mu[j]; 
       
          loglik[i] = loglik1;
          current = candidate;
          }
        else
          {
          Beta[i] = REAL(Betaout)[i-1];
          Q[i] = REAL(Qout)[i-1];
          G[i] = REAL(Gout)[i-1];
          tauN[i] = REAL(tauNout)[i-1];
          tauF[i] = REAL(tauFout)[i-1];
          tauNF[i] = REAL(tauNFout)[i-1];
          loglik[i] = REAL(loglikout)[i-1];
          } 
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

SEXP call_metrop_logpost(SEXP y, SEXP times, SEXP func, SEXP rtol,
		SEXP atol, SEXP rho, SEXP tcrit, SEXP jacfunc, SEXP initfunc,
		SEXP verbose, SEXP hmin, SEXP hmax,    
        SEXP index_R,SEXP index2_R, SEXP two_n_R, SEXP n_R,
        SEXP indep_R, SEXP Y_R, SEXP modY_R, SEXP VN_R, SEXP VF_R, 
        SEXP Beta_R, SEXP Q_R, SEXP G_R, 
        SEXP tauN_R, SEXP tauF_R, SEXP tauNF_R,
        SEXP v_R, SEXP S_R, SEXP tauN_sh_R, SEXP tauN_sc_R, SEXP tauF_sh_R, SEXP tauF_sc_R,
        SEXP indBeta_R, SEXP aBeta_R, SEXP bBeta_R,
        SEXP indQ_R, SEXP aQ_R, SEXP bQ_R,
        SEXP indG_R, SEXP aG_R, SEXP bG_R){
     
     SEXP r, ans, parms;
     int i,j, indep, n, two_n, indBeta, indQ, indG, *index, *index2;
     double aBeta, bBeta, aQ, bQ, aG, bG, loglik1;
     double *Y, *Ytild, *modY, *mu, *sigma, VN, VF;
     double v, *S, tauN_sh, tauN_sc, tauF_sh, tauF_sc, logprior;

          
     two_n = INTEGER(two_n_R)[0];
     n = INTEGER(n_R)[0];

     S = (double *) R_alloc(4, sizeof(double));
     sigma = (double *) R_alloc(4, sizeof(double));
     mu = (double *) R_alloc(two_n, sizeof(double));
     Y = (double *) R_alloc(two_n, sizeof(double));
     Ytild = (double *) R_alloc(two_n, sizeof(double));
     modY = (double *) R_alloc(two_n, sizeof(double));
     index = (int *) R_alloc(two_n, sizeof(int));
     index2 = (int *) R_alloc(two_n, sizeof(int));

     indep = INTEGER(indep_R)[0];
     VN = REAL(VN_R)[0];
     VF = REAL(VF_R)[0];

     v = REAL(v_R)[0];
     for(i=0;i<4;i++){S[i]=REAL(S_R)[i];}
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

     PROTECT(ans = allocVector(REALSXP, 1)); 

     parms = NEW_NUMERIC(5);
     REAL(parms)[0] = REAL(Beta_R)[0];
     REAL(parms)[1] = VN;
     REAL(parms)[2] = VF;
     REAL(parms)[3] = REAL(Q_R)[0];
     REAL(parms)[4] = REAL(G_R)[0];

     r = call_lsoda(y, times, func, parms, rtol, 
           atol, rho, tcrit, jacfunc, initfunc,
		   verbose, hmin, hmax);

     for(j=0; j < two_n; j++)
        {
        mu[j] = log(REAL(r)[index[j]]); 
        Ytild[j] = log(REAL(r)[index2[j]]);
        }
  
    sigma[0] = REAL(tauN_R)[0];
    sigma[1] = sigma[2] = REAL(tauNF_R)[0];
    sigma[3] = REAL(tauF_R)[0];
 
   logprior = 0; 
    if(indep==1)
       {
       logprior =  tauN_sh * log(tauN_sc) - lgammafn(tauN_sh) - (tauN_sh + 1) * log(sigma[0]) - (tauN_sc/sigma[0]);
       logprior = logprior + tauF_sh * log(tauF_sc) - lgammafn(tauF_sh) - (tauF_sh + 1) * log(sigma[3]) - (tauF_sc/sigma[3]);
       }
    else
       {
       logprior = log(dinvwish(sigma, v, S));
       }

   logprior = logprior + log(dpriors(REAL(parms)[0], indBeta, aBeta, bBeta, (int)0) + 1E-200);
   logprior = logprior + log(dpriors(REAL(parms)[3], indQ, aQ, bQ, (int)0) + 1E-200);
   logprior = logprior + log(dpriors(REAL(parms)[4], indG, aG, bG, (int)0) + 1E-200);

   loglik1 = 0;
  
   F77_CALL(loglikelihood) (mu, sigma, Y, Ytild, modY, &n, &two_n, &loglik1);
   REAL(ans)[0] = loglik1 + logprior;
  
   UNPROTECT(1);
   return (ans);
}


SEXP call_logpost_metrop(SEXP y, SEXP times, SEXP func, SEXP rtol,
		SEXP atol, SEXP rho, SEXP tcrit, SEXP jacfunc, SEXP initfunc,
		SEXP verbose, SEXP hmin, SEXP hmax,    
        SEXP size_R, SEXP index_R,SEXP index2_R, SEXP two_n_R, SEXP n_R,
        SEXP indep_R, SEXP Y_R, SEXP modY_R, SEXP VN_R, SEXP VF_R, 
        SEXP Beta_R, SEXP Q_R, SEXP G_R, 
        SEXP tauN_R, SEXP tauF_R, SEXP tauNF_R,
        SEXP v_R, SEXP S_R, SEXP tauN_sh_R, SEXP tauN_sc_R, SEXP tauF_sh_R, SEXP tauF_sc_R,
        SEXP indBeta_R, SEXP aBeta_R, SEXP bBeta_R,
        SEXP indQ_R, SEXP aQ_R, SEXP bQ_R,
        SEXP indG_R, SEXP aG_R, SEXP bG_R,
        SEXP total_it_R){
     
     SEXP r, ans, parms;
     int i,j, size, indep, n, two_n, *index, *index2, total_it;
     int indBeta, indQ, indG;
     double loglik1, logprior;
     double *Y, *Ytild, *modY, *mu, *sigma, VN, VF;
     double v, *S, tauN_sh, tauN_sc, tauF_sh, tauF_sc;
     double aBeta, bBeta, aQ, bQ, aG, bG;

  
     two_n = INTEGER(two_n_R)[0];
     n = INTEGER(n_R)[0];
     size = INTEGER(size_R)[0];
     total_it = INTEGER(total_it_R)[0];

     S = (double *) R_alloc(4, sizeof(double));
     sigma = (double *) R_alloc(4, sizeof(double));
     mu = (double *) R_alloc(two_n, sizeof(double));
     Y = (double *) R_alloc(two_n, sizeof(double));
     Ytild = (double *) R_alloc(two_n, sizeof(double));
     modY = (double *) R_alloc(two_n, sizeof(double));
     index = (int *) R_alloc(two_n, sizeof(int));
     index2 = (int *) R_alloc(two_n, sizeof(int));

     indep = INTEGER(indep_R)[0];
     VN = REAL(VN_R)[0];
     VF = REAL(VF_R)[0];

     v = REAL(v_R)[0];
     for(i=0;i<4;i++){S[i]=REAL(S_R)[i];}
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

    PROTECT(ans = allocVector(REALSXP, size)); 

    double control_bar, next_control_bar;
    next_control_bar = control_bar = 1/(double)50;

     for(i = 0; i < size; i++)
       {
       parms = NEW_NUMERIC(5);
       REAL(parms)[0] = REAL(Beta_R)[i];
       REAL(parms)[1] = VN;
       REAL(parms)[2] = VF;
       REAL(parms)[3] = REAL(Q_R)[i];
       REAL(parms)[4] = REAL(G_R)[i];

       r = call_lsoda(y, times, func, parms, rtol, 
           atol, rho, tcrit, jacfunc, initfunc,
		   verbose, hmin, hmax);

       for(j=0; j < two_n; j++)
          {
          mu[j] = log(REAL(r)[index[j]]); 
          Ytild[j] = log(REAL(r)[index2[j]]);
          }
 
       sigma[0] = REAL(tauN_R)[0];
       sigma[1] = sigma[2] = REAL(tauNF_R)[0];
       sigma[3] = REAL(tauF_R)[0];
 
      logprior = 0; 
      if(indep==1)
       {
       logprior =  tauN_sh * log(tauN_sc) - lgammafn(tauN_sh) - (tauN_sh + 1) * log(sigma[0]) - (tauN_sc/sigma[0]);
       logprior = logprior + tauF_sh * log(tauF_sc) - lgammafn(tauF_sh) - (tauF_sh + 1) * log(sigma[3]) - (tauF_sc/sigma[3]);
       }
      else
       {
       logprior = log(dinvwish(sigma, v, S));
       }

      logprior = logprior + log(dpriors(REAL(parms)[0], indBeta, aBeta, bBeta, (int)0) + 1E-500);
      logprior = logprior + log(dpriors(REAL(parms)[3], indQ, aQ, bQ, (int)0) + 1E-500);
      logprior = logprior + log(dpriors(REAL(parms)[4], indG, aG, bG, (int)0) + 1E-500);

      loglik1 = 0;
      F77_CALL(loglikelihood) (mu, sigma, Y, Ytild, modY, &n, &two_n, &loglik1);
      REAL(ans)[i] = loglik1 + logprior;

      if((i+1)/(double)total_it >= next_control_bar)
         {
         Rprintf("=");
         R_FlushConsole();
         next_control_bar = next_control_bar + control_bar;
         }

       }
     UNPROTECT(1);
     return (ans);
     }

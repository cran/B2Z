#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts

#include <Rdefines.h>
#include <Rinternals.h>

#include <Rmath.h>
#include "pacoteR_utils.h"


void F77_NAME(loglikelihood)(double *, double *, double *, double *, double *, int *, int *, double *);

SEXP call_sir(SEXP y, SEXP times, SEXP func, SEXP rtol,
		SEXP atol, SEXP rho, SEXP tcrit, SEXP jacfunc, SEXP initfunc,
		SEXP verbose, SEXP hmin, SEXP hmax,    
        SEXP m_R, SEXP index_R,SEXP index2_R, SEXP two_n_R, SEXP n_R,
        SEXP indep_R, SEXP Y_R, SEXP modY_R, SEXP VN_R, SEXP VF_R, 
        SEXP indBeta_R, SEXP aBeta_R, SEXP bBeta_R,
        SEXP indQ_R, SEXP aQ_R, SEXP bQ_R,
        SEXP indG_R, SEXP aG_R, SEXP bG_R,
        SEXP v_R, SEXP S_R, SEXP tauN_sh_R, SEXP tauN_sc_R, SEXP tauF_sh_R, SEXP tauF_sc_R){
     
     SEXP r, ans, parms, Beta, Q, G, tauN, tauF, tauNF, loglik;
     int  indBeta, aBeta, bBeta, indQ, aQ, bQ;
     int indG, aG, bG, i,j, m, indep, n, two_n, *index, *index2;
     double *pBeta, *pQ, *pG, *ptauN, *ptauF, *ptauNF, *ploglik,  loglik1;
     double *Y, *Ytild, *modY, *mu, *sigma, taus[4], VN, VF;
     double v, *S, tauN_sh, tauN_sc, tauF_sh, tauF_sc;
          
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
     m = INTEGER(m_R)[0];
     VN = REAL(VN_R)[0];
     VF = REAL(VF_R)[0];

     v = REAL(v_R)[0];
     for(i=0;i<4;i++){S[i]=REAL(S_R)[i];}
     tauN_sh = REAL(tauN_sh_R)[0];  tauN_sc = REAL(tauN_sc_R)[0];
     tauF_sh = REAL(tauF_sh_R)[0];  tauF_sc = REAL(tauF_sc_R)[0];
     indBeta = INTEGER(indBeta_R)[0]; aBeta = REAL(aBeta_R)[0]; bBeta = REAL(bBeta_R)[0];
     indQ = INTEGER(indQ_R)[0]; aQ = REAL(aQ_R)[0]; bQ = REAL(bQ_R)[0];
     indG = INTEGER(indG_R)[0]; aG = REAL(aG_R)[0]; bG = REAL(bG_R)[0];
          
     for(i = 0; i < two_n; i++ )
       {
       Y[i] = REAL(Y_R)[i];
       modY[i] = REAL(modY_R)[i];
       index[i] = INTEGER(index_R)[i];
       index2[i] = INTEGER(index2_R)[i];
       }

     PROTECT(ans = allocVector(VECSXP, (m+7) )); 
     PROTECT(Beta = NEW_NUMERIC(m)); PROTECT(Q = NEW_NUMERIC(m)); PROTECT(G = NEW_NUMERIC(m)); PROTECT(tauN = NEW_NUMERIC(m)); 
     PROTECT(tauF = NEW_NUMERIC(m)); PROTECT(tauNF = NEW_NUMERIC(m)); PROTECT(loglik = NEW_NUMERIC(m));

     pBeta = NUMERIC_POINTER(Beta); pQ = NUMERIC_POINTER(Q); pG = NUMERIC_POINTER(G); ptauN = NUMERIC_POINTER(tauN); 
     ptauF = NUMERIC_POINTER(tauF); ptauNF = NUMERIC_POINTER(tauNF);   ploglik = NUMERIC_POINTER(loglik);

    if (m < 50)
       {
       Rprintf("|-");
       for(i=0;i<m-1;i++)
         Rprintf("-");
       
       Rprintf("--|\n");
       }
    else
       {
       Rprintf("|----------25%%----------50%%----------75%%----------|\n");
       }
       
    R_FlushConsole();
    Rprintf("|");
    R_FlushConsole();
    
    double control_bar, next_control_bar;
    next_control_bar = control_bar = 1/(double)50;
    
     for(i = 0; i < m; i++)
       {
       parms = NEW_NUMERIC(5);
       REAL(parms)[0] = pBeta[i] = rpriors(indBeta, aBeta, bBeta);
       REAL(parms)[1] = VN;
       REAL(parms)[2] = VF;
       REAL(parms)[3] = pQ[i] = rpriors(indQ, aQ, bQ);
       REAL(parms)[4] = pG[i] = rpriors(indG, aG, bG);

       r = call_lsoda(y, times, func, parms, rtol, 
           atol, rho, tcrit, jacfunc, initfunc,
		   verbose, hmin, hmax);

       for(j=0; j < two_n; j++)
          {
          mu[j] = log(REAL(r)[index[j]]); 
          Ytild[j] = log(REAL(r)[index2[j]]);
          }
  
       if(indep == 1)
         {
         ptauN[i] = 1.0/rpriors(2,tauN_sh,tauN_sc);  
         ptauF[i] = 1.0/rpriors(2,tauF_sh,tauF_sc);
         ptauNF[i]= 0;
         }
       else
         {
         rinvwish(v,S,taus);
         ptauN[i] = taus[0];
         ptauF[i] = taus[3];
         ptauNF[i] = taus[1];
         }
   
  sigma[0] = ptauN[i];
  sigma[1] = sigma[2] = ptauNF[i];
  sigma[3] = ptauF[i];

  
  loglik1 = 0;
  
  F77_CALL(loglikelihood) (mu, sigma, Y, Ytild, modY, &n, &two_n, &loglik1);
  ploglik[i] = loglik1;

  SET_VECTOR_ELT(ans, i, r);
  
  if((i+1)/(double)m >= next_control_bar)
    {
    Rprintf("=");
    R_FlushConsole();
    next_control_bar = next_control_bar + control_bar;
    }
    
  }

  SET_VECTOR_ELT(ans, m, Beta);
  SET_VECTOR_ELT(ans, (m+1), Q);
  SET_VECTOR_ELT(ans, (m+2), G);
  SET_VECTOR_ELT(ans, (m+3), tauN);
  SET_VECTOR_ELT(ans, (m+4), tauF);
  SET_VECTOR_ELT(ans, (m+5), tauNF);
  SET_VECTOR_ELT(ans, (m+6), loglik);

  Rprintf("|\n");
  R_FlushConsole();

UNPROTECT(8);
return (ans);
}


SEXP sir_likelihood(SEXP y, SEXP times, SEXP func, SEXP rtol,
		SEXP atol, SEXP rho, SEXP tcrit, SEXP jacfunc, SEXP initfunc,
		SEXP verbose, SEXP hmin, SEXP hmax,    
        SEXP parms_R, SEXP index_R,SEXP index2_R, SEXP two_n_R, SEXP n_R,
        SEXP Y_R, SEXP modY_R, SEXP VN_R, SEXP VF_R){
     
     SEXP r, ans, parms;
     int i, j, n, two_n, *index, *index2;
     double *Y, *Ytild, *modY, *mu, *sigma, VN, VF, loglik1;
          
     two_n = INTEGER(two_n_R)[0];
     n = INTEGER(n_R)[0];

     sigma = (double *) R_alloc(4, sizeof(double));
     mu = (double *) R_alloc(two_n, sizeof(double));
     Y = (double *) R_alloc(two_n, sizeof(double));
     Ytild = (double *) R_alloc(two_n, sizeof(double));
     modY = (double *) R_alloc(two_n, sizeof(double));
     index = (int *) R_alloc(two_n, sizeof(int));
     index2 = (int *) R_alloc(two_n, sizeof(int));

     VN = REAL(VN_R)[0];
     VF = REAL(VF_R)[0];
                    
     for(i = 0; i < two_n; i++ )
       {
       Y[i] = REAL(Y_R)[i];
       modY[i] = REAL(modY_R)[i];
       index[i] = INTEGER(index_R)[i];
       index2[i] = INTEGER(index2_R)[i];
       }

     PROTECT(ans = allocVector(REALSXP,2));
     REAL(ans)[0] = REAL(ans)[1] = 0;
     
     parms = NEW_NUMERIC(5);
     REAL(parms)[0] = REAL(parms_R)[0];
     REAL(parms)[1] = VN;
     REAL(parms)[2] = VF;
     REAL(parms)[3] = REAL(parms_R)[1];
     REAL(parms)[4] = REAL(parms_R)[2];

     r = call_lsoda(y, times, func, parms, rtol, 
           atol, rho, tcrit, jacfunc, initfunc,
		   verbose, hmin, hmax);

     for(j=0; j < two_n; j++)
        {
        mu[j] = log(REAL(r)[index[j]]); 
        Ytild[j] = log(REAL(r)[index2[j]]);
        }
  
     sigma[0] = REAL(parms_R)[3];
     sigma[1] = sigma[2] = REAL(parms_R)[5];
     sigma[3] = REAL(parms_R)[4];
   
     loglik1 = 0;
  
     F77_CALL(loglikelihood) (mu, sigma, Y, Ytild, modY, &n, &two_n, &loglik1);
     REAL(ans)[0] = loglik1;
    
    
     UNPROTECT(1);
  
     return (ans);
     }

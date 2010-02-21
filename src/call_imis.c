#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts

#include <Rdefines.h>
#include <Rinternals.h>

#include <Rmath.h>
#include "pacoteR_utils.h"



void F77_NAME(loglikelihood)(double *, double *, double *, double *, double *, int *, int *, double *);

SEXP call_imis_initial(SEXP y, SEXP times, SEXP func, SEXP rtol,
		SEXP atol, SEXP rho, SEXP tcrit, SEXP jacfunc, SEXP initfunc,
		SEXP verbose, SEXP hmin, SEXP hmax,    
        SEXP size_R, SEXP index_R,SEXP index2_R, SEXP two_n_R, SEXP n_R,
        SEXP indep_R, SEXP Y_R, SEXP modY_R, SEXP VN_R, SEXP VF_R, 
        SEXP Beta_R, SEXP Q_R, SEXP G_R, 
        SEXP v_R, SEXP S_R, SEXP tauN_sh_R, SEXP tauN_sc_R, SEXP tauF_sh_R, SEXP tauF_sc_R,
        SEXP total_it_R){
     
     SEXP r, ans, parms, tauN, tauF, tauNF, loglik, priorsigma;
     int i,j, size, indep, n, two_n, *index, *index2, total_it;
     double *ptauN, *ptauF, *ptauNF, *ploglik, *ppriorsigma, loglik1;
     double *Y, *Ytild, *modY, *mu, *sigma, taus[4], VN, VF;
     double v, *S, tauN_sh, tauN_sc, tauF_sh, tauF_sc, lp1, lp2;
          
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
          
     for(i = 0; i < two_n; i++ )
       {
       Y[i] = REAL(Y_R)[i];
       modY[i] = REAL(modY_R)[i];
       index[i] = INTEGER(index_R)[i];
       index2[i] = INTEGER(index2_R)[i];
       }

     PROTECT(ans = allocVector(VECSXP, (size + 5) )); 
     PROTECT(tauN = NEW_NUMERIC(size)); PROTECT(tauF = NEW_NUMERIC(size)); PROTECT(tauNF = NEW_NUMERIC(size)); 
     PROTECT(loglik = NEW_NUMERIC(size));  PROTECT(priorsigma = NEW_NUMERIC(size));


     ptauN = NUMERIC_POINTER(tauN); ptauF = NUMERIC_POINTER(tauF); ptauNF = NUMERIC_POINTER(tauNF);   
     ploglik = NUMERIC_POINTER(loglik);  ppriorsigma = NUMERIC_POINTER(priorsigma);


     Rprintf("|----------25%%----------50%%----------75%%----------|\n");
       
     R_FlushConsole();
     Rprintf("|");
     R_FlushConsole();

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
  
     if(indep == 1)
         {
         ptauN[i] = 1.0/rpriors(2,tauN_sh,tauN_sc);  
         ptauF[i] = 1.0/rpriors(2,tauF_sh,tauF_sc);
         ptauNF[i]= 0;

         sigma[0] = ptauN[i];
         sigma[1] = sigma[2] = 0;
         sigma[3] = ptauF[i];
  
         lp1 =  tauN_sh * log(tauN_sc) - lgammafn(tauN_sh) - (tauN_sh + 1) * log(sigma[0]) - (tauN_sc/sigma[0]);
         lp2 =  tauF_sh * log(tauF_sc) - lgammafn(tauF_sh) - (tauF_sh + 1) * log(sigma[3]) - (tauF_sc/sigma[3]);

         ppriorsigma[i] = exp(lp1 + lp2);
         }

       else
         {
         rinvwish(v,S,taus);
         ptauN[i]  =  taus[0];
         ptauF[i]  =  taus[3];
         ptauNF[i] =  taus[1];
         sigma[0] = ptauN[i];
         sigma[1] = sigma[2] = ptauNF[i];
         sigma[3] = ptauF[i];
         ppriorsigma[i] = dinvwish(sigma, v, S);
         }
   
 
  loglik1 = 0;
  
  F77_CALL(loglikelihood) (mu, sigma, Y, Ytild, modY, &n, &two_n, &loglik1);
  ploglik[i] = loglik1;

  SET_VECTOR_ELT(ans, i, r);
  
   if((i+1)/(double)total_it >= next_control_bar)
    {
    Rprintf("=");
    R_FlushConsole();
    next_control_bar = next_control_bar + control_bar;
    }
    
  }

  SET_VECTOR_ELT(ans, size, tauN);
  SET_VECTOR_ELT(ans, (size+1), tauF);
  SET_VECTOR_ELT(ans, (size+2), tauNF);
  SET_VECTOR_ELT(ans, (size+3), loglik);
  SET_VECTOR_ELT(ans, (size+4), priorsigma);
   
UNPROTECT(6);
return (ans);
}



SEXP call_loglik_imis(SEXP y, SEXP times, SEXP func, SEXP rtol,
		SEXP atol, SEXP rho, SEXP tcrit, SEXP jacfunc, SEXP initfunc,
		SEXP verbose, SEXP hmin, SEXP hmax,    
        SEXP size_R, SEXP index_R,SEXP index2_R, SEXP two_n_R, SEXP n_R,
        SEXP indep_R, SEXP Y_R, SEXP modY_R, SEXP VN_R, SEXP VF_R, 
        SEXP Beta_R, SEXP Q_R, SEXP G_R, 
        SEXP tauN_R, SEXP tauF_R, SEXP tauNF_R,
        SEXP v_R, SEXP S_R, SEXP tauN_sh_R, SEXP tauN_sc_R, SEXP tauF_sh_R, SEXP tauF_sc_R,
        SEXP total_it_R){
     
     SEXP r, ans, parms, loglik, priorsigma;
     int i,j, size, indep, n, two_n, *index, *index2, total_it;
     double *ploglik, *ppriorsigma, loglik1;
     double *Y, *Ytild, *modY, *mu, *sigma, VN, VF;
     double v, *S, tauN_sh, tauN_sc, tauF_sh, tauF_sc, lp1, lp2;
          
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
          
     for(i = 0; i < two_n; i++ )
       {
       Y[i] = REAL(Y_R)[i];
       modY[i] = REAL(modY_R)[i];
       index[i] = INTEGER(index_R)[i];
       index2[i] = INTEGER(index2_R)[i];
       }

     PROTECT(ans = allocVector(VECSXP, (size + 2) )); 
     PROTECT(loglik = NEW_NUMERIC(size));  PROTECT(priorsigma = NEW_NUMERIC(size));

     ploglik = NUMERIC_POINTER(loglik);  ppriorsigma = NUMERIC_POINTER(priorsigma);

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
 
   
       sigma[0] = REAL(tauN_R)[i];
       sigma[1] = sigma[2] = REAL(tauNF_R)[i];
       sigma[3] = REAL(tauF_R)[i];

       if(indep == 1)
         {  
         lp1 =  tauN_sh * log(tauN_sc) - lgammafn(tauN_sh) - (tauN_sh + 1) * log(sigma[0]) - (tauN_sc/sigma[0]);
         lp2 =  tauF_sh * log(tauF_sc) - lgammafn(tauF_sh) - (tauF_sh + 1) * log(sigma[3]) - (tauF_sc/sigma[3]);

         ppriorsigma[i] = exp(lp1 + lp2);
         }
       else
         {
         ppriorsigma[i] = dinvwish(sigma, v, S);
         }

       loglik1 = 0;
  
       F77_CALL(loglikelihood) (mu, sigma, Y, Ytild, modY, &n, &two_n, &loglik1);
       ploglik[i] = loglik1;

       SET_VECTOR_ELT(ans, i, r);

       if((i+1)/(double)total_it >= next_control_bar)
         {
         Rprintf("=");
         R_FlushConsole();
         next_control_bar = next_control_bar + control_bar;
         }

       }

       SET_VECTOR_ELT(ans, size, priorsigma);
       SET_VECTOR_ELT(ans, (size + 1), loglik);



UNPROTECT(3);
return (ans);
}


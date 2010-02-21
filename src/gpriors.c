#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts

#include <Rdefines.h>
#include <Rinternals.h>

#include <Rmath.h>


void rinvwish(double v, double *S, double *answer){
double det, Sinv[4], CC[4], Z[4], ZCC[4], rwish[4];
det = S[0]*S[3] - S[1]*S[2];

Sinv[0] = (1/det)*S[3];
Sinv[1] = -(1/det)*S[1];
Sinv[2] = -(1/det)*S[2];
Sinv[3] = (1/det)*S[0];

CC[0] = sqrt(Sinv[0]);
CC[1] = Sinv[1]/CC[0];
CC[2] = 0;
CC[3] = sqrt(Sinv[3]-(Sinv[2]*Sinv[2])/Sinv[0]);
	
GetRNGstate();
Z[0] = sqrt(rchisq(v));
Z[1] = rnorm(0,1);
Z[2] = 0;
Z[3] =sqrt(rchisq(v - 1.0));
PutRNGstate();

ZCC[0]= Z[0]*CC[0] + Z[1]*CC[2];
ZCC[1]= Z[0]*CC[1] + Z[1]*CC[3];
ZCC[2]= Z[2]*CC[0] + Z[3]*CC[2];
ZCC[3]= Z[2]*CC[1] + Z[3]*CC[3];

rwish[0] = ZCC[0]*ZCC[0] + ZCC[2]*ZCC[2];
rwish[1] = rwish[2] = ZCC[0]*ZCC[1] + ZCC[2]*ZCC[3];
rwish[3] = ZCC[1]*ZCC[1] + ZCC[3]*ZCC[3];

det = rwish[0]*rwish[3] - rwish[1]*rwish[2];
answer[0] = (1/det)*rwish[3];
answer[1] = -(1/det)*rwish[1];
answer[2] = -(1/det)*rwish[2];
answer[3] = (1/det)*rwish[0];

}

void rmtvnorm(int d, double *mean, double *retval, double *answer){
double mat[d];

for(int i=0;i< d;i++)
   {
   GetRNGstate();
   mat[i] = rnorm(0,1);
   PutRNGstate();

   for(int j=0; j<d ; j++)
      answer[j] = answer[j] + retval[i + j*d]*mat[i] + mean[j]/(double)d;
   }
}


double dinvwish(double *W, double v, double *S)
   {
   double detS, detW, gammapart;       
   double num, denom, invW[4], tracehold;    
   int k =2;
   gammapart = 1;
   
   for (int i=1;i<=2;i++) {gammapart = gammapart * gammafn((double)0.5*(v + 1 - i));}
   
   denom = gammapart * pow(2.0,(double)(v * 0.5*k)) * pow(M_PI,(double)(k * 0.25*(k - 1)));
   detS = S[0]*S[3] - S[1]*S[2];
   detW = W[0]*W[3] - W[1]*W[2];
  
   invW[0] = (1/detW)*W[3];
   invW[1] = -(1/detW)*W[1];
   invW[2] = -(1/detW)*W[2];
   invW[3] = (1/detW)*W[0];
   
   tracehold = S[0]*invW[0] + S[1]*invW[2] + S[2]*invW[1] + S[3]*invW[3];
   
   num = pow(detS,(double)0.5*v) * pow(detW,(double)(-0.5*(v + k + 1))) * exp(-0.5 * tracehold);
   return(num/denom);
   }


double rpriors(int ind, double a, double b)
  {
  double val = 0.0;
  switch(ind)
  {
  case 1:
  GetRNGstate();
  val=runif(a,b);
  PutRNGstate();
  break;
   
  case 2:
  GetRNGstate();
  val=rgamma(a,1.0/b);
  PutRNGstate(); 
  break;

  case 3:
  GetRNGstate();
  val=rexp(1/(a));
  PutRNGstate(); 
  break;

 
  case 4:
  PutRNGstate(); 
  val = rnorm(a,b);
  PutRNGstate(); 
  break;
  
  case 5:
  GetRNGstate();
  val=rt(a);
  PutRNGstate();
  break;
   
  case 6:
  GetRNGstate();
  val=rweibull(a,b);
  PutRNGstate(); 
  break;
  
  case 7:
  GetRNGstate();
  val = rf(a,b);
  PutRNGstate(); 
  break;
  
  case 8:
  PutRNGstate(); 
  val = rchisq(a);
  PutRNGstate(); 
  break;

  case 9:
  PutRNGstate(); 
  val = rcauchy(a,b);
  PutRNGstate(); 
  break;

  case 10:
  GetRNGstate();
  val = rlnorm(a,b);
  PutRNGstate(); 
  break;

  }
  return(val);
 }

double dpriors(double x, int ind, double a, double b, int give_log)
  {
  double val = 0.0;
  switch(ind)
  {
  case 1:
  val=dunif(x, a, b, give_log);
  break;
   
  case 2:
  val=dgamma(x, a, b, give_log);
  break;

  case 3:
  val=dexp(x, 1/(a), give_log);
  break;

 
  case 4:
  val = dnorm(x, a, b, give_log);
  break;
  
  case 5:
  val=dt(x, a, give_log);
  break;
   
  case 6:
  val=dweibull(x, a, b, give_log);
  break;
  
  case 7:
  val = df(x, a, b, give_log);
  break;
  
  case 8:
  val = dchisq(x, a, give_log);
  break;

  case 9:
  val = dcauchy(x, a, b, give_log);
  break;

  case 10:
  val = dlnorm(x, a, b, give_log);
  break;
  }
  return(val);
 }


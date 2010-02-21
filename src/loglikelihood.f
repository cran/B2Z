      subroutine loglikelihood (mu,sigma,Y,Ytild,modY,n,two_n,loglik)
      integer n, two_n, i, j
	double precision mu, sigma, Y, Ytild, modY, loglik
	double precision x, det, inv_sigma, trace_sigma
	double precision eigen1, eigen2, Smat, Svet, distval, logdet
	dimension mu(two_n),sigma(2,2),Y(n,2),Ytild(n,2),modY(two_n)
	dimension x(two_n), inv_sigma(2,2), Smat(n,2), Svet(two_n)
	x = modY - mu
	det = sigma(1,1)*sigma(2,2) - sigma(2,1)*sigma(1,2)
      inv_sigma(1,1) =(1/det)*sigma(2,2) 
      inv_sigma(1,2) =-(1/det)*sigma(1,2)
      inv_sigma(2,1) =-(1/det)*sigma(2,1)
      inv_sigma(2,2) =(1/det)*sigma(1,1)
      trace_sigma = sigma(1,1) + sigma(2,2)
	eigen1 = 0.5*(trace_sigma + sqrt(trace_sigma**2-4*det))
      eigen2 = 0.5*(trace_sigma - sqrt(trace_sigma**2-4*det))
      Smat = matmul((Y-Ytild),inv_sigma)
   	j=1
      do i = 1, n
      Svet(j) = Smat(i,1)
      j = j+1
      Svet(j) = Smat(i,2)
      j = j+1
      end do
	distval = dot_product(Svet,x) 
      logdet = n*(log(eigen1) + log(eigen2))
      loglik = -(two_n * 1.837877 + logdet + distval)/2.0 
      return
      end
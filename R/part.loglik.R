part.loglik <- function (par, VN, VF, funct, y0, times, Y1, 
               modY,  tauN, tauF, tauNF, index, index2, n, 
               two_n, total_it, next_control_bar){

    y = y0
    T = length(times)
    rtol = 1e-06
    atol = 1e-06
    tcrit = NULL
    jacfunc = NULL
    verbose = FALSE
    dllname = NULL
    hmin = 0
    hmax = Inf
    initpar = NULL
    ModelInit = NULL
    rho = environment(funct)
    Nglobal = NULL
				
    ans <-  .Call("call_metrop_partialloglik", y, times, funct, rtol,
        atol, rho, tcrit, jacfunc, ModelInit, as.integer(verbose), 
        hmin, hmax, as.integer(index), as.integer(index2), 
        as.integer(two_n), as.integer(n), Y1, modY, VN, VF, 
        par[1], par[2], par[3],tauN, tauF, tauNF, PACKAGE="B2Z")


    if(par[4]/total_it >= par[5])
       {
       cat("=")
       flush.console()
       }

return(ans)
}
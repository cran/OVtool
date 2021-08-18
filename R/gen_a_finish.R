#### gen_a_finish fn ####
gen_a_finish <- function(a_res, my_estimand, wts){

  if(my_estimand == "ATE"){
    U1 <- stats::rnorm(a_res$n1)*sqrt(a_res$ve1) + a_res$b1*a_res$Rstar_1 +
      a_res$es*(1-a_res$pi)
    U0 <- stats::rnorm(a_res$n0)*sqrt(a_res$ve0) + a_res$b0*a_res$Rstar_0 -
      a_res$es*a_res$pi

    U <- rep(NA, a_res$n)
    ind <- a_res$ind
    U[ind] <- U1
    U[-ind] <- U0

    a <- rep(NA, a_res$n)

    a[ind] <- 1 / wts[ind] + (1 - (1 / wts[ind])) *
      (1 / exp(-a_res$es^2 * (1 - 2 * a_res$pi) / 2 + a_res$es * U[ind]))

    a[-ind] <- 1 / wts[-ind] + (1 - (1 / wts[-ind])) *
      (1 / exp(a_res$es^2 * (1 - 2 * a_res$pi) / 2 - a_res$es * U[-ind]))

  } else if(my_estimand == "ATT"){
    U0 <- stats::rnorm(a_res$n0)*sqrt(a_res$ve0) + a_res$b0*a_res$Rstar_0 -
      a_res$es*a_res$pi

    U <- rep(NA, a_res$n)
    ind <- a_res$ind
    U[ind] <- 1
    U[-ind] <- U0

    a <- rep(NA, a_res$n)

    a[ind] = 1
    a[-ind] = exp(-a_res$es^2 * (1 - 2 * a_res$pi) / 2 + a_res$es * U[-ind])

  }  else{
    stop("OVtool currently only works under the ATE or ATT setting.")
  }

  return(a)
}



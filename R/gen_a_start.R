#### gen_a_start fn ####
gen_a_start <- function(y, tx, residuals, es, rho, my_estimand){
  ind <- which(tx == 1)

  cdf1 = EnvStats::ecdfPlot(residuals[ind], discrete = F, plot.it = F)
  cdf0 = EnvStats::ecdfPlot(residuals[-ind], discrete = F, plot.it = F)
  Rstar1 = stats::qnorm(cdf1$Cumulative.Probabilities[rank(residuals[ind], ties.method = 'random')])
  Rstar0 = stats::qnorm(cdf0$Cumulative.Probabilities[rank(residuals[-ind], ties.method = 'random')])

  n1 <- sum(tx)
  n0 <- sum(1-tx)
  n <- n1 + n0
  c1 <- stats::cov(Rstar1,y[ind])
  c0 <- stats::cov(Rstar0,y[-ind])

  pi <- mean(tx)

  v1 <- stats::var(y[ind])
  v0 <- stats::var(y[-ind])

  sd1 <- sqrt(v1)
  sd0 <- sqrt(v0)

  vU <- 1 + es^2*pi*(1-pi)

  Y <- y - mean(y)

  vY <- stats::var(y)

  b1lim <- sd1/c1
  b0lim <- sd0/c0

  A <- rho * sqrt(vU * vY)

  Q <- es*(1-pi)*mean(y[ind])*pi - es*pi*mean(y[-ind])*(1-pi)

  # Take b1
  alpha = (A - Q)/((1-pi)*c0)
  beta = (-c1*pi)/((1-pi)*c0)

  if(beta > 0){
    b1low = max(-b1lim, ((-b0lim - alpha) / beta))
    b1high = min(b1lim, ((b0lim - alpha) / beta))
  }
  if(beta < 0){
    b1low = max(-b1lim, ((b0lim - alpha) / beta))
    b1high = min(b1lim, ((-b0lim - alpha) / beta))
  }

  # Set b0 == b1
  b1 = (A - Q) / ( ((1-pi)*c0) * (1 + ((c1*pi) / ((1-pi)*c0))) )
  b0 <- (A-b1*c1*pi-Q)/((1-pi)*c0)
  ve1 <- 1 - b1^2 * stats::var(Rstar1)
  ve0 <- 1 - b0^2 * stats::var(Rstar0)
  if(!(abs(b0) <= b0lim)) stop("b0 is too large in absolute value. Try reducing the size of the grid.")

  return(a_res = list(n1 = n1, ve1 = ve1, b1 = b1, es=es, pi = pi,
                      n0 = n0, ve0 = ve0, b0 = b0, n = n, ind = ind,
                      y = y, rho = rho, Rstar_1 = Rstar1, Rstar_0 = Rstar0,
                      residuals = residuals))
}


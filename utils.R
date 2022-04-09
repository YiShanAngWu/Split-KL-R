## Helper functions
square_diff <- function(a,b){
  return((a-b)^2)
}

KLGauss <- function(mean1, mean2, sigma2){
  return(sum(square_diff(mean1,mean2))/(2*sigma2))
}

bisection <- function(fn, a, b, tol = 1e-8) {
  # If the signs of the function at the evaluated points, a and b, stop the function and return message.
  if (fn(a) * fn(b) > 0)
    stop('signs of f(a) and f(b) are the same')
  p <-  (a + b)/2;
  err = abs(fn(p));
  while (err > tol){
    ifelse(fn(a)*fn(p)<0, b <- p, a <- p)
    p <- (a + b)/2; 
    err <- abs(fn(p));
  }
  return(p)
}

kl_inv_sup <- function(Emp, RHS){
  if (Emp==0)
    # will result in numerical issues
    Emp <- 1e-8
  fn <- function(p) Emp*log(Emp/p) + (1-Emp)*log((1-Emp)/(1-p)) - RHS;
  if (fn(1-1e-8) <= 0)
    return(1-1e-8)
  p <- bisection(fn, Emp, 1-1e-8)
  return(p)
}

kl_inv_inf <- function(Emp, RHS){
  if (Emp==0)
    # will result in numerical issues
    Emp <- 1e-8
  fn <- function(p) Emp*log(Emp/p) + (1-Emp)*log((1-Emp)/(1-p)) - RHS;
  if (fn(1e-8) <= 0)
    return(1e-8)
  p <- bisection(fn, 1e-8, Emp)
  return(p)
}

bin_inv_sup <- function(n, k, delta){
  #print("bin_inv_sup")
  fn <- function(p) pbinom(k, n, p) - delta;
  if (fn(1)>=0)
    return(1)
  p <- bisection(fn, k/n, 1)
  return(p)
}

get_sample <- function(type, mean, variance2, n_samples=1){
  return(t(mvrnorm(n = n_samples, mu = mean, Sigma = variance2*diag(d))))
}
## End of helper functions
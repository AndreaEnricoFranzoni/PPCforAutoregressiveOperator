#function to evaluate mean squared error MSE on the predictions
mse <- function(x1,x2)
{
  n = length(x1)
  err = numeric(n)
  
  for (i in 1:n) {
    err[i] = (x1[i] - x2[i])^2
  }
  
  return((sum(err)/n))
}


#function to evaluate average distance AD on the predictions
ad <- function(x1,x2)
{
  n = length(x1)
  err = numeric(n)
  
  for (i in 1:n) {
    err[i] = abs(x1[i] - x2[i])
  }
  
  return((sum(err)/n))
}


#En as described in Kokoszka
En <- function(x1,x2,t.grid)
{
  err.vals    <- (x1 - x2)^2
  err.fun     <- approxfun(x=t.grid, y=err.vals)
  en          <- sqrt(integrate(err.fun, lower=range(t.grid)[1], upper=range(t.grid)[2], rel.tol=.Machine$double.eps^.05)$value)
  
  return(en)
}


#Rn as described in Kokoszka
Rn <- function(x1,x2,t.grid)
{
  err.vals    <- abs(x1 - x2)
  err.fun     <- approxfun(x=t.grid, y=err.vals)
  rn          <- integrate(err.fun, lower=range(t.grid)[1], upper=range(t.grid)[2], rel.tol=.Machine$double.eps^.05)$value
  
  return(rn)
}
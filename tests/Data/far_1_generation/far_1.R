#innovation
gaussian_k <- function(y, t.grid, a)
{
  f_int_t <- approxfun(x=t.grid, y=y*exp(-(t.grid^2)/2))
  ft_1 <- a*exp(-(t.grid^2)/2)*(integrate(f_int_t, lower=range(t.grid)[1], upper=range(t.grid)[2], rel.tol=.Machine$double.eps^.1)$value)
  
  return(ft_1)
}


identity_k <- function(y, t.grid, a)
{
  f_int_t = approxfun(x=t.grid, y=y)
  ft_1 = a*(integrate(f_int_t, lower=range(t.grid)[1], upper=range(t.grid)[2], rel.tol=.Machine$double.eps^.1)$value)
  
  return(ft_1)
}


sloping_plane_t_k <- function(y, t.grid, a)
{
  f_int_t = approxfun(x=t.grid, y=y)
  ft_1 = a*t.grid*(integrate(f_int_t, lower=range(t.grid)[1], upper=range(t.grid)[2], rel.tol=.Machine$double.eps^.1)$value)
  
  return(ft_1) 
}


sloping_plane_s_k <- function(y, t.grid, a)
{
  f_int_t <- approxfun(x=t.grid, y=y*t.grid)
  ft_1 <- a*(integrate(f_int_t, lower=range(t.grid)[1], upper=range(t.grid)[2], rel.tol=.Machine$double.eps^.1)$value)
  
  return(ft_1)
}




innovation <- function(kernel_id)
{
  if(kernel_id=="gaussian") return(gaussian_k)          #gaussian kernel
  if(kernel_id=="identity") return(identity_k)          #identity kernel
  if(kernel_id=="sp_t") return(sloping_plane_t_k)       #sloping plane t
  if(kernel_id=="sp_s") return(sloping_plane_s_k)       #sloping plane s
}




#noise
eps_1 <- function(t.grid)
{
  w <- Wiener(n=1, pts=t.grid)
  w <- as.vector(w)
  e <- w - t.grid*tail(w, n=1)
  
  return(e)
}


eps_2 <- function(t.grid)
{
  set.seed(23032000)
  eps1 <- rnorm(1)
  eps2 <- rnorm(2)
  lambda <- 0.5
  
  e <- eps1*sqrt(2)*sin(2*pi*t.grid) + eps2*sqrt(2*lambda)*cos(2*pi*t.grid)
  
  return(e)
}


eps_3 <- function(t.grid, a=1)
{
  e <- eps_2(t.grid) + a*eps_1(t.grid)
  
  return(e)
}


noise <- function(noise_id)
{
  if(noise_id=="1") return(eps_1)                   #brownian brridges
  if(noise_id=="2") return(eps_2)                   #combination of sin and cos
  if(noise_id=="3") return(eps_3)                   #combination of the first two errors
}




far_1 <- function(kernel_id, noise_id, n, t.grid, a, burnin)
{
  
  f.sample <- matrix(data=0, nrow=length(t.grid), ncol=n)
  
  #innovation function
  innovation_f <- innovation(kernel_id)
  #noise function
  noise_f <- noise(noise_id)
  
  #f at instant 0 is equivalent to the initial error
  e0 <- noise_f(t.grid)
  yt <- e0
  
  for(i in 1:(n+burnin))
  {
    # Noise at t+1
    et_1 <- noise_f(t.grid)
    
    # Innovations at t+1
    ft_1.vals <- innovation_f(yt, t.grid, a)
    
    #new value of the function
    yt <- ft_1.vals + et_1
    if(i>burnin) f.sample[, i-burnin] <- yt
  }
  
  return(f.sample)
}



feat_far_1_process <- function(kernel_id, norm)
{
  if(kernel_id=="gaussian")
  { 
    if(norm==0.5)
    {
      a    <- (1/2)*(1/0.7468)
      name <- "Gaussian, norm=0.5"
    }
    if(norm==0.8)
    {
      a    <- (4/5)*(1/0.7468)
      name <- "Gaussian, norm=0.8"
    }
  }
  if(kernel_id=="identity")
  {
    if(norm==0.5)
    {
      a    <- 1/2
      name <- "Identity, norm=0.5"
    }
    if(norm==0.8)
    {
      a    <- 4/5
      name <- "Identity, norm=0.8"
    }
  }
  if(kernel_id=="sp_t")
  {
    if(norm==0.5)
    {
      a    <- sqrt(3)/2
      name <- "Sloping plane t, norm=0.5"
    }
    if(norm==0.8)
    {
      a    <- sqrt(3)*4/5
      name <- "Sloping plane t, norm=0.8"
    }
  }
  if(kernel_id=="sp_s")
  {
    if(norm==0.5)
    {
      a    <- sqrt(3)/2
      name <- "Sloping plane s, norm=0.5"
    }
    if(norm==0.8)
    {
      a    <- sqrt(3)*4/5
      name <- "Sloping plane s, norm=0.8"
    }
  }
  
  #kernel = kernel_id
  #constant = a
  #kern_name = name
  l=list(kernel_id,a,name)
  names(l) = c("kernel","constant","name")
  #l["kernel_id"] <- kernel_id
  #l["constant"]  <- a
  #l["name"]      <- name
  
  return(l)
}

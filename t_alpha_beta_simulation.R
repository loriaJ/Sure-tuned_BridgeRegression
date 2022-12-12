log_a_x <- function(x,alpha){
  (1-alpha)^(-1) * alpha * log(sin(alpha*x)) + log(sin((1-alpha)*x)) - (1-alpha)^(-1) *log(sin(x))
}

log_B_x <- function(x,alpha){
  log(sin(x)) - alpha * log(sin(alpha*x)) - (1-alpha)*log(sin((1-alpha)*x))
}

B_x <- function(x,alpha){
  (sin(x))/(sin(alpha*x)^alpha * sin((1-alpha)*x)^(1-alpha))
}

log_2pi <- log(2*base::pi)

## verify correct simulations using moments of rt_alpha_beta
rt_alpha_beta <- function(n=1,beta0,alpha){
  
  sigma0 <- 1/sqrt(beta0*(1-alpha))
  b_0 <- alpha^(-alpha) * (1-alpha)^(-(1-alpha))
  beta_alpha <- beta0/alpha
  if(sigma0 >= sqrt(2*pi)){
    u0 <- runif(1,0,base::pi)
    v0 <- runif(1,0,1)
    x <- u0
    w <- B_x(x,alpha = alpha)
    while(! v0<= (w/b_0)^(beta_alpha)){
      u0 <- runif(1,0,base::pi)
      v0 <- runif(1,0,1)
      x <- u0
      w <- B_x(x,alpha = alpha)
    }
  } else { 
    N0 <- rnorm(1,0,1)
    v0 <- runif(1,0,1)
    x <- sigma0*abs(N0)
    w <- B_x(x,alpha)
    while(!((x <= base::pi) & (exp(-N0^2/2)*v0 <=
                               (w/b_0)^beta_alpha))){
      N0 <- rnorm(1)
      v0 <- runif(1,0,1)
      x <- sigma0*abs(N0)
      w <- B_x(x,alpha)
    }
  }
  
  g <- rgamma(1,shape = 1+beta_alpha * (1-alpha),scale = 1)
  w^(-1/alpha)*g^(-(1-alpha)/(alpha))
}


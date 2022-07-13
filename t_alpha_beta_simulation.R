

log_a_x <- function(x,alpha){
  (1-alpha)^(-1) * alpha * log(sin(alpha*x)) + log(sin((1-alpha)*x)) - (1-alpha)^(-1) *log(sin(x))
}


log_B_x <- function(x,alpha){
  log(sin(x)) - alpha * log(sin(alpha*x)) - (1-alpha)*log(sin((1-alpha)*x))
}

B_x <- function(x,alpha){
  (sin(x))/(sin(alpha*x)^alpha * sin((1-alpha)*x)^(1-alpha))
}


# B_x_b0 <- function(x,alpha){
#   (rwfec::sinc(x))/(rwfec::sinc(alpha*x)^alpha * rwfec::sinc((1-alpha)*x)^(1-alpha))
# }


log_2pi <- log(2*base::pi)

# rzolotarev <- function(n=1,alpha,b){
#   # sigma0 <- 1/(b*alpha*(1-alpha))^(1/2)
#   log_sigma0 <- -log(b) -log(alpha) -log(1-alpha)
#   
#   if(log_sigma0 >= log_2pi){
#     u0 <- runif(1,0,base::pi)
#     v0 <- runif(1,0,1)
#     
#   } else { 
#     
#     
#   }
#     
# }


# rt_alpha_beta <- function(n=1,beta0,alpha){
#   
#   log_sigma0 <- -log(beta0) -log(1-alpha)
#   log_b_0 <- -alpha*log(alpha) - (1-alpha)*log((1-alpha))
#   # inv_b_0 <- alpha^(alpha) * (1-alpha)^((1-alpha))
#   beta_alpha <- beta0/alpha
#   if(log_sigma0 >= log_2pi){
#     u0 <- runif(1,0,base::pi)
#     v0 <- runif(1,0,1)
#     x <- u0
#     log_w <- log_B_x(x,alpha = alpha)
#     while(! log(v0)<= beta_alpha * (log_w - log_b_0)){
#       u0 <- runif(1,0,base::pi)
#       v0 <- runif(1,0,1)
#       x <- u0
#       log_w <- log_2B_x(x,alpha)
#     }
#   } else { 
#     N0 <- rnorm(1)
#     v0 <- runif(1,0,1)
#     x <- abs(N0)
#     w <- B_x(x,alpha)
#     while(!((x <= base::pi) & (-N0^2/2 + log(v0) <=
#                                beta_alpha*(log_w - log_b_0)))){
#     while(!((x <= base::pi) & (-N0^2/2 + log(v0) <=
#                                beta_alpha*(log_w - log_b_0)))){
#       N0 <- rnorm(1)
#       v0 <- runif(1,0,1)
#       x <- exp(log_sigma)*abs(N0)
#       log_w <- log_B_x(x,alpha)
#     }
#   }
#   
#   g <- rgamma(1,shape = 1+beta_alpha * (1-alpha),scale = 1)
#   exp(-1/(alpha) * log_w - 1/(alpha*(1-alpha)) * log(g))
# }


## verify correct simulations using moments of rt_alpha_beta
rt_alpha_beta <- function(n=1,beta0,alpha){
  
  sigma0 <- 1/sqrt(beta0*(1-alpha))
  # log_b_0 <- -alpha*log(alpha) - (1-alpha)*log((1-alpha))
  b_0 <- alpha^(-alpha) * (1-alpha)^(-(1-alpha))
  beta_alpha <- beta0/alpha
  if(sigma0 >= sqrt(2*pi)){
    u0 <- runif(1,0,base::pi)
    v0 <- runif(1,0,1)
    x <- u0
    # log_w <- log_B_x(x,alpha = alpha)
    w <- B_x(x,alpha = alpha)
    # while(! log(v0)<= beta_alpha * (log_w - log_b_0)){
    while(! v0<= (w/b_0)^(beta_alpha)){
      u0 <- runif(1,0,base::pi)
      v0 <- runif(1,0,1)
      x <- u0
      w <- B_x(x,alpha = alpha)
    }
    #print(w)    
    
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
    #print(w)
  }
  
  g <- rgamma(1,shape = 1+beta_alpha * (1-alpha),scale = 1)
  # exp(-1/(alpha) * log_w - 1/(alpha*(1-alpha)) * log(g))
  
  # 1/(w*g^(1-alpha))^(1/alpha)
  w^(-1/alpha)*g^(-(1-alpha)/(alpha))
}



# for(i in 1:19){
#   alp <- i/20
#   beta0 <- 0.5
#   r00 <- replicate(1000,mean(1/replicate(1000,rt_alpha_beta(1,beta0,alp))))
#   ## verify correct simulations using moments of rt_alpha_beta
#   # hist(r00)
#   true_mean <- gamma(1+beta0)*gamma(1+(1+beta0)/alp)/(gamma(1+beta0/alp)*gamma(1+1+beta0))
#   print(
#   abs(mean(r00) - true_mean)/true_mean
#   )
# }

# alp <- 0.1
# beta0 <- 0.5
# r00 <- replicate(1000,mean(1/replicate(1000,rt_alpha_beta(1,beta0,alp))))
# ## verify correct simulations using moments of rt_alpha_beta
# hist(r00)
# mean(r00)
# gamma(1+beta0)*gamma(1+(1+beta0)/alp)/(gamma(1+beta0/alp)*gamma(1+1+beta0))
# 
# hist(r00)


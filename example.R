source('multivariate_sure_n_p.R')
set.seed(1030)
n <- 10
p <- 30
beta <- c(10,0,-10,rep(0,p-3))+ rnorm(p,0,0.5)

xx0 <- matrix(rnorm(n*p),ncol = p,nrow = n)
y <- xx0 %*% beta + rnorm(n,0,1)

results <- t_multivariate_sure_schains_n_p(alpha = 0.3,x = xx0,
					   y = y,n_sim = 3000,sigma2=1,
					   nu_vect = c(10^c(-10:3)),
					   nchains = 1)
plot(results$beta_estim$beta_mean,x=beta,
     ylab = 'Estimated coefficients',
     xlab = 'True coefficients')

									   

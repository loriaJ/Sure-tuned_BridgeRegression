source('multivariate_sure_no_svd.R')
set.seed(1030)
n <- 10
p <- 30
beta <- c(10,0,-10,rep(0,p-3))+ rnorm(p,0,0.5)

xx0 <- matrix(rnorm(n*p),ncol = p,nrow = n)
y <- xx0 %*% beta + rnorm(n,0,1)

results <- t_multivariate_sure_schains(alpha = 0.3,x = xx0,
									   y = y,n_sim = 3000,sigma2=1,
									   nu_vect = c(10^c(-10:3)),
									   nchains = 4)
plot(results$beta_estim,x=beta,
     ylab = 'Estimated coefficients',
     xlab = 'True coefficients')

									   
source('t_alpha_beta_simulation.R')


## for the case of n < p:
t_moments_p_n <- function(y,nu0,xlambd_xt,sigma2){
  n_samp <- dim(xlambd_xt)[1]
  n <- length(y)
  expectations <- matrix(NA,nrow = n,ncol = n_samp)
  variances <- array(NA,dim = c(n,n,n_samp))
  
  det_1 <- rep(NA,n_samp)
  p_i <- rep(NA,n_samp)
  
  expectations <- array(dim = c(n,n_samp))
  variances <- array(dim = c(n,n,n_samp))
  
  for(i in 1:n_samp){
    A <- nu0 * xlambd_xt[i,,]
    premp <- A
    diag(premp) <- diag(premp) + sigma2
    premp <- chol(premp)
    # premp <- (chol(A + sigma_mat))
    det_1[i] <- sum(log(diag(premp)))
    
    premp <- chol2inv(premp)
    prempy <- premp %*% y
    expectations[,i] <- A %*% prempy 
    variances[,,i] <- 
      sigma2 * premp %*% A
    p_i[i] <- 0.5*t(y)%*% prempy
    
    variances[,,i] <- variances[,,i] + expectations[,i] %*% t(expectations[,i])
  }
  
  s0 <- - p_i - det_1
  s1 <- sum(exp(s0-max(s0,na.rm = T)))
  # s2 <- sum(exp(-s0+max(s0,na.rm = T)))
  # print(dim(expectations))
  xbeta_estim <- expectations %*% exp(s0-max(s0,na.rm = T))/s1
  variances <- apply(variances,MARGIN = c(1,2),`%*%`,exp(s0-max(s0,na.rm = T)))/s1
  # print(dim(xbeta_estim))
  # variances0 <- variances - xbeta_estim %*% t(xbeta_estim)
  # print(c(sum(diag(variances0)),
  #         sum(diag(variances)) - sum(xbeta_estim^2) ))
  #print()
  #variances2 <- sum(xbeta_estim^2)
  
  return(list(xbeta_estim = xbeta_estim,
              variance = variances,
              var2 = sum(xbeta_estim^2),
              post_prob = s1))
}


## for the case of n < p, and measuring variances
t_moments_p_n_mod <- function(y,nu0,sigma_mat,x,n_sim=1000,alpha=0.3){
  n <- length(y)
  p <- dim(x)[2]
  xlambd_xt <- array(dim= c(n_sim,n,n))
  lambda_tx <- array(dim= c(n_sim,p,n))
  lambda_t <- array(dim= c(n_sim,p,p))
  
  for(j in 1:n_sim){
    samp <- replicate(rt_alpha_beta(n = 1,
                                    alpha = alpha/2,
                                    beta0 = 0.5),n = p)
    # samp <- diag(1/samp)
    lambda_t[j,,] <- samp
    lambda_tx[j,,] <- (1/samp) * t(x)
    xlambd_xt[j,,] <- x %*% lambda_tx[j,,]
    
  }
  n_samp <- dim(xlambd_xt)[1]
  expectations <- matrix(NA,nrow = p,ncol = n_samp)
  variances <- array(NA,dim = c(p,p,n_samp))
  
  det_1 <- rep(NA,n_samp)
  p_i <- rep(NA,n_samp)
  
  for(i in 1:n_samp){
    A <- nu0 * xlambd_xt[i,,]
    premp <- (chol(A + sigma_mat))
    det_1[i] <- sum(log(diag(premp)))
    
    premp <- chol2inv(premp)
    prempy <- premp %*% y

    expectations[,i] <- lambda_tx[i,,] %*% prempy 
    variances[,,i] <- 
      nu0 * lambda_t[i,,] - lambda_tx[i,,] %*% premp %*% t(lambda_tx[i,,]) 
    p_i[i] <- 0.5*t(y)%*% prempy
    
    variances[,,i] <- variances[,,i] + expectations[,i] %*% t(expectations[,i])
  }
  
  s0 <- - p_i - det_1
  post_prob <- exp(s0-max(s0,na.rm = T))
  s1 <- sum(post_prob)
  x_expectations <- expectations %*% post_prob / s1
  x_variances <- apply(variances,MARGIN = c(1,2),`%*%`,post_prob) /s1
  
  return(list(x_expectations = x_expectations,
              x_variances = x_variances,
              expectations=expectations,
              variances=variances,
              post_prob = s1*exp(max(s0,na.rm = T))))
}


compute_betas <- function(xlambd_xt,lambda_tx,sigma2,nu0,y){
  dims <- dim(lambda_tx)
  nchains <- dims[1]
  nsim <- dims[2]
  p <- dims[3]
  n <- dims[4]
  # sigma_inv <- solve(sigma_mat)
  p_i <- rep(NA,nsim)
  det_1 <- rep(NA,nsim)
  exp_1 <- matrix(NA,ncol = nsim,nrow = p)
  weights <- rep(NA,nchains)
  beta_estim <- matrix(NA,nrow = nchains,ncol = p)
  for(ch in 1:nchains){
    for(i in 1:nsim){
      v_t_nu <- xlambd_xt[ch,i,,]
      diag(v_t_nu) <- diag(v_t_nu)  + 1/nu0 * sigma2
      cxtx <- chol(v_t_nu)
      inv_cxtx <- chol2inv(cxtx)
      inv_cxtxy <- inv_cxtx %*% y
      det_1[i] <- sum(log(diag(cxtx)))
      
      p_i[i] <- 0.5*t(y)%*% inv_cxtxy/nu0
      
      exp_1[,i] <- lambda_tx[ch,i,,] %*% inv_cxtxy
       
    }
    s0 <- -p_i - det_1
    weights[ch] <- sum(exp(s0-max(s0,na.rm = T)))
    beta_estim[ch,] <- exp_1 %*% exp(s0-max(s0,na.rm = T))/weights[ch]
  }
  beta_matrix <- beta_estim
  
  beta_mean <- colMeans(beta_estim)
  return(list(beta_estim = beta_estim,
              beta_mean = beta_mean,
              weights = weights))
}



t_multivariate_sure_schains_n_p <- function(alpha,x,y,n_sim=10000,
                                            sigma2,nu_vect=10^(-20:3),
                                            nchains = 1,seed = 10){
  
  n <- dim(x)[1]
  p <- dim(x)[2]
  xlambd_xt <- array(dim= c(nchains,n_sim,n,n))
  lambda_tx <- array(dim= c(nchains,n_sim,p,n))
  for(ch in 1:nchains){
    for(j in 1:n_sim){
      samp <- replicate(rt_alpha_beta(n = 1,
                                      alpha = alpha/2,
                                      beta0 = 0.5),n = p)
      # samp <- diag(1/samp)
      lambda_tx[ch,j,,] <- (1/samp) * t(x)
      xlambd_xt[ch,j,,] <- x %*% lambda_tx[ch,j,,]
    }
  }
  
  
  # sigma_mat <- diag(sigma2,n)
  sure_f <- function(nu0,just_sure = TRUE){
    xbeta_mat <- matrix(data = NA,nrow = nchains,ncol = n)
    xbeta_var <- array(data = NA,dim = c(nchains,n,n))
    var2 <- rep(NA,nchains)
    m_i <- list()
    post_probs <- rep(NA,nchains)
    for(i in 1:nchains){
      m_i[[i]] <- t_moments_p_n(y = y,nu0 = nu0,
                                xlambd_xt = xlambd_xt[i,,,],
                                sigma2 = sigma2)
      xbeta_mat[i,] <- m_i[[i]]$xbeta_estim
      xbeta_var[i,,] <- m_i[[i]]$variance
      var2[i] <- m_i[[i]]$var2
      post_probs[i] <- m_i[[i]]$post_prob 
    }
    xbeta_gr <- colMeans(xbeta_mat)
    xbeta_var_gr <- apply(xbeta_var,MARGIN = c(2,3),FUN= mean)
    if(nchains > 1){
      var_xbeta <- sum(diag(var(xbeta_mat)))
    } else{
      var_xbeta <- 0
    }
    sure <- sum((y - xbeta_gr)^2) + 2*(sum(diag(xbeta_var_gr)) - mean(var2)) + 2*var_xbeta
    if(just_sure){
      print(c(nu0,sure ))
      return(sure)
    } else {
      var_xbeta <- array(data = NA,dim = c(nchains,n,n))
      for(i in 1:nchains){
        var_xbeta[i,,] <- (m_i[[i]]$variance)  
      }
      dim(xbeta_var_gr)
      var_exp_beta <- var(xbeta_mat) 
      dim(var_exp_beta)
      var_exp_beta <- var_exp_beta + xbeta_var_gr
      list(sure = sure,
           xbeta_estim = xbeta_gr,
           var_beta = var_exp_beta,
           sure_decomp = c(sse = sum((y -xbeta_gr)^2),
                           var_1 = sum(diag(xbeta_var_gr)),
                           var_2 = sum(diag(var(xbeta_mat)))),
           nu0 = nu0,
           mi = m_i)
    }
  }
  sure_vect <- sapply(nu_vect,sure_f)
  n0 <- which.min(sure_vect)
  n0_inf <- max(1,n0-1)
  n0_sup <- min(length(nu_vect),n0+1)
  opt0 <- optimise(function(log_nu)sure_f(exp(log_nu)),
                   interval = log(nu_vect[c(n0_inf,n0_sup)]))
  
  # optim_sol <- sure_f(nu0=exp(opt0$minimum),
  #                     just_sure = FALSE)
  beta_estim <- compute_betas(xlambd_xt = xlambd_xt,
                              lambda_tx = lambda_tx,
                              sigma2 = sigma2,
                              nu0 = exp(opt0$minimum),y = y)
  
  estims <- list(
    nu = exp(opt0$minimum),
    sure = (opt0$objective),
    beta_estim = beta_estim,
    # sure_decomp = optim_sol$sure_decomp,
    # mse_opt = mean((y-optim_sol$xbeta_estim)^2),
    mse_beta = mean((y- x %*% beta_estim$beta_mean)^2)
  )
  return(estims)
}


## method for obtaining the coefficients using cross validation
compute_betas_cv <- function(x,y,nu,sigma2=1,sampl_T){
  n_sim <- dim(sampl_T)[1]
  p <- dim(sampl_T)[2]
  n <- dim(x)[1]
  exp_1 <- matrix(NA,ncol = n_sim,nrow = p)
  sigma_mat <- diag(sigma2,n)
  p_i <- rep(NA,n_sim)
  det_1 <- rep(NA,n_sim)
  for(i in 1:n_sim){
    lambda_t <- diag(1/sampl_T[i,],p)
    A <- x %*% lambda_t %*% t(x)
    v_t_nu <- A + 1/nu * sigma_mat
    cxtx <- chol(v_t_nu)
    inv_cxtx <- chol2inv(cxtx)
    inv_cxtxy <- inv_cxtx %*% y
    det_1[i] <- sum(log(diag(cxtx)))
    
    p_i[i] <- 0.5*t(y)%*% inv_cxtxy/nu
    exp_1[,i] <- c(lambda_t %*% t(x) %*% inv_cxtxy)
  }
  s0 <- -p_i - det_1
  weights <- sum(exp(s0-max(s0,na.rm = T)))
  beta_estim <- exp_1 %*% exp(s0-max(s0,na.rm = T))/weights
  
  return(beta_estim)
}


cross_validation_betas <- function(alpha,x,y,n_sim=1000,
                                   sigma2,nu_vect=10^(-20:3),
                                   n_folds = 10){
  n <- dim(x)[1]
  p <- dim(x)[2]
  rand_samp <- matrix(replicate(
    n_sim*p,
    rt_alpha_beta(n = 1,
                  alpha = alpha/2,
                  beta0 = 0.5)),
    ncol = p)

  cv_partition <- suppressWarnings(split(sample(n),1:n_folds))
  cv_error_beta <- function(nu0){
    errors_i <- rep(NA,n_folds)
    for(i in 1:n_folds){
      x0 <- x[-(cv_partition[[i]]),]
      y0 <- y[-(cv_partition[[i]])]
      x1 <- x[(cv_partition[[i]]),]
      y1 <- y[(cv_partition[[i]])]
      errors_i[i] <- mean((y1 - x1%*% compute_betas_cv(x=x0,y=y0,nu=nu0,sigma2=1,
                                                       sampl_T = rand_samp))^2)
    }
    print(c(nu0,sum(errors_i)))
    if(is.na(sum(errors_i))){
      print(errors_i)
    } 
    return(sum(errors_i))
  }
  
  err_vect <- sapply(nu_vect,cv_error_beta)
  n0 <- which.min(err_vect)
  n0_inf <- max(1,n0-1)
  n0_sup <- min(length(nu_vect),n0+1)
  opt0 <- optimise(function(log_nu)cv_error_beta(exp(log_nu)),
                   interval = log(nu_vect[c(n0_inf,n0_sup)]))
  return(list(betas = compute_betas_cv(x = x,y = y,nu = exp(opt0$minimum),sigma2 = 1,sampl_T = rand_samp),
              nu = exp(opt0$minimum)))
}

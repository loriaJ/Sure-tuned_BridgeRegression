source('t_alpha_beta_simulation.R')



t_ab_post_moments_nu <- function(x,y,nu0,rand_samp,alpha=0.3,sigma_mat){
  n <- dim(x)[1]
  p <- dim(x)[2]
  n_samp <- dim(rand_samp)[1]
  p_i <- rep(NA,n_samp)
  det_1 <- rep(NA,n_samp)
  sqrti <- rep(NA,n_samp)
  exp_1 <- matrix(NA,ncol = n_samp,nrow = p)
  exp_bb <- array(NA,dim=c(p,p,n_samp))
  exp_bb2 <- array(NA,dim=c(p,p,n_samp))
  var_b_t <- array(NA,dim=c(p,p,n_samp))
  diag_n <- diag(1,n)
  tr_0 <- rep(NA,n_samp)
  nu_inv <- nu0^(-1)
  nu_x_sigma_inv_x <- nu0 * t(x)%*%solve(sigma_mat,x)

  for(i in 1:n_samp){
    samp <- rand_samp[i,]
    lambda_t <- diag(1/(samp))
    lambda_t_inv <- diag(samp)
    
    cl2 <- chol(lambda_t_inv + nu_x_sigma_inv_x)
    bs2 <- chol2inv(cl2)
    
    
    v_t_nu <- x%*% lambda_t %*% t(x) + nu_inv * sigma_mat
    
    { 
      cxtx <- chol(v_t_nu)
      inv_cxtx <- chol2inv(cxtx)
      inv_cxtxy <- inv_cxtx %*% y
      det_1[i] <- sum(log(diag(cxtx)))
    }
    
    p_i[i] <- 0.5*t(y)%*% inv_cxtxy*nu_inv
    
    exp_1[,i] <- lambda_t %*% t(x) %*% inv_cxtxy # * nu_inv
    
    exp_bb2[,,i] <- nu0 * bs2 + exp_1[,i] %*% t(exp_1[,i])
    
    exp_1[,i] <- exp_1[,i] 
  }
  
  s0 <- - p_i - det_1
  s1 <- sum(exp(s0-max(s0,na.rm = T)))
  s2 <- sum(exp(-s0+max(s0,na.rm = T)))
  beta_estim <- exp_1 %*% exp(s0-max(s0,na.rm = T))/s1
  exp_bb2 <- apply(exp_bb2,MARGIN = c(1,2),`%*%`,exp(s0-max(s0,na.rm = T)))/s1
  var_b <- exp_bb2 - beta_estim %*% t(beta_estim)
  
  xtx <- t(x) %*% x
 
  return(list(s0 = s0,
              mean = beta_estim, 
              exp_bb2 = exp_bb2,
              var2 = var_b,
              post_prob = s1,
              post_prob_inv = s2))
}



t_multivariate_seq_sure <- function(alpha,x,y,n_sim=10000,sigma2,nu_vect=10^(-20:3),y_scale =1){
  
  n <- dim(x)[1]
  p <- dim(x)[2]
  rand_samp <- matrix(replicate(
    n_sim*p,
    rt_alpha_beta(n = 1,
                  alpha = alpha/2,
                  beta0 = 0.5)),
    ncol = p)
  
  sigma_mat <- diag(sigma2,n)
  sure_f <- function(nu0){
    # print(nu0)
    m_i <- t_ab_post_moments_nu(x = x,y = y,nu0 = nu0,
                                rand_samp = rand_samp,sigma_mat = sigma_mat)
    
    var_b <- m_i$var2
    d0 <- diag(x%*%var_b%*% t(x))
    print(c(nu0,sum(d0)))
    log(sum((y-x%*%m_i$mean)^2)*y_scale^2 +  2 *(sum(d0))*y_scale^2)
  }
  
  sure_vect <- sapply(nu_vect,sure_f)
  n0 <- which.min(sure_vect)
  n0_inf <- max(1,n0-1)
  n0_sup <- min(length(nu_vect),n0+1)
  opt0 <- optimise(sure_f,
                   interval = nu_vect[c(n0_inf,n0_sup)])
  # opt0 <- optimise(sure_f,interval = range(nu_vect))
  
  optim_sol <- t_ab_post_moments_nu(x=x,y=y,nu0=opt0$minimum,
                                    rand_samp = rand_samp,sigma_mat = sigma_mat)
  
  estims <- list(
    nu = opt0$minimum,
    sure = exp(opt0$objective),
    beta_estim = optim_sol$mean,
    var_beta = optim_sol$var_sure,
    post_probs = optim_sol$post_prob,
    mse = sum((y-x%*%optim_sol$mean)^2)# ,
    # penalty = #2 *(sum(diag(xtx %*% optim_sol$var2 %*% xtx)) + tr_pen)
  )
  return(estims)
}


t_multivariate_seq_sure_optim <- function(alpha,x,y,n_sim=10000,sigma2,nu_vect=10^(-20:3)){
  
  n <- dim(x)[1]
  p <- dim(x)[2]
  rand_samp <- matrix(replicate(
    n_sim*p,
    rt_alpha_beta(n = 1,
                  alpha = alpha/2,
                  beta0 = 0.5)),
    ncol = p)
  
  xtx <- t(x)%*% x
  sigma_mat <- diag(sigma2,n)
  sure_f <- function(nu0){
    
    m_i <- t_ab_post_moments_nu(x = x,y = y,nu0 = nu0,
                             rand_samp = rand_samp,sigma_mat = sigma_mat)
    
    pred_y <- x%*%m_i$mean
    var_b <- m_i$var2
    d0 <- diag( x %*% var_b %*% t(x))
    print(c(nu0,sum(d0)))
    log(sum((y-pred_y)^2) +  2 *(sum(d0)))
  }
  
  sure_vect <- sapply(nu_vect,sure_f)
  n0 <- which.min(sure_vect)
  n0_inf <- max(1,n0-1)
  n0_sup <- min(length(nu_vect),n0+1)
  opt0 <- optim(sure_f,par = sure_vect[n0],method = 'L-BFGS-B',
                lower = nu_vect[n0_inf],upper = nu_vect[n0_sup],control = list(parscale=nu_vect[n0]))
  # opt0 <- optimise(sure_f,interval = range(nu_vect))
  
  optim_sol <- t_ab_post_moments_nu(x=x,y=y,nu0=opt0$par,
                                 rand_samp = rand_samp,sigma_mat = sigma_mat)
  
  estims <- list(
    nu = opt0$par,
    sure = exp(opt0$value),
    beta_estim = optim_sol$mean,
    var_beta = optim_sol$var_sure,
    post_probs = optim_sol$post_prob,
    mse = sum((y-x%*%optim_sol$mean)^2)
  )
  return(estims)
}




t_multivariate_sure_schains <- function(alpha,x,y,n_sim=10000,sigma2,nu_vect=10^(-20:3),nchains = 2){
  
  if(nchains <= 1){
    return(t_multivariate_sure(alpha,x,y,n_sim,sigma2,nu_vect))
  }
  
  n <- dim(x)[1]
  p <- dim(x)[2]
  rand_samp <- array(replicate(
    n_sim*p*nchains,
    rt_alpha_beta(n = 1,
                  alpha = alpha/2,
                  beta0 = 0.5)),
    dim = c(n_sim,p,nchains))
  
  sigma_mat <- diag(sigma2,n)
  sure_f <- function(nu0,just_sure = TRUE){
    beta_mat <- matrix(data = NA,nrow = nchains,ncol = p)
    xbeta_var <- array(data = NA,dim = c(nchains,n,n))
    m_i <- list()
    post_probs <- rep(NA,nchains)
    for(i in 1:nchains){
      m_i[[i]] <- t_ab_post_moments_nu(x = x,y = y,nu0 = nu0,
                                       rand_samp = rand_samp[,,i],
                                       sigma_mat = sigma_mat)
      beta_mat[i,] <- m_i[[i]]$mean
      xbeta_var[i,,] <- x%*%m_i[[i]]$var2%*%t(x)
      post_probs[i] <- m_i[[i]]$post_prob 
    }
    beta_gr <- colMeans(beta_mat)
    xbeta_var_gr <- apply(xbeta_var,MARGIN = c(2,3),FUN= mean)
    var_xbeta <- sum(diag(var(beta_mat%*%t(x))))
    sure <- sum((y -x%*%beta_gr)^2) + 2*sum(diag(xbeta_var_gr)) + 2*var_xbeta
    if(just_sure){
      return(sure)
    } else {
      var_beta <- array(data = NA,dim = c(nchains,p,p))
      for(i in 1:nchains){
        var_beta[i,,] <- (m_i[[i]]$var2)  
      }
      var_beta <- apply(var_beta,MARGIN = c(2,3),FUN = mean)
      var_exp_beta <- ((var(beta_mat)))
      list(sure = sure,
           beta_estim = beta_gr,
           var_beta = var_beta + var_exp_beta,
           nu0 = nu0,
           mi = m_i)
    }
  }
  sure_vect <- sapply(nu_vect,sure_f)
  n0 <- which.min(sure_vect)
  n0_inf <- max(1,n0-1)
  n0_sup <- min(length(nu_vect),n0+1)
  opt0 <- optimise(sure_f,
                   interval = nu_vect[c(n0_inf,n0_sup)])

  optim_sol <- sure_f(nu0=opt0$minimum,
                      just_sure = FALSE)
  
  estims <- list(
    nu = opt0$minimum,
    sure = (opt0$objective),
    beta_estim = optim_sol$beta_estim,
    var_beta = optim_sol$var_beta,
    post_probs = optim_sol$post_prob,
    mi = optim_sol$mi,
    mse = sum((y-x%*%optim_sol$beta_estim)^2)# ,
    # penalty = #2 *(sum(diag(xtx %*% optim_sol$var2 %*% xtx)) + tr_pen)
  )
  return(estims)
}





t_multivariate_sure_schains_surface <- function(alpha,x,y,n_sim=10000,sigma2,nu_vect=10^(-20:3),nchains = 2){
  
  if(nchains <= 1){
    return(t_multivariate_sure(alpha,x,y,n_sim,sigma2,nu_vect))
  }
  
  n <- dim(x)[1]
  p <- dim(x)[2]
  rand_samp <- array(replicate(
    n_sim*p*nchains,
    rt_alpha_beta(n = 1,
                  alpha = alpha/2,
                  beta0 = 0.5)),
    dim = c(n_sim,p,nchains))
  
  sigma_mat <- diag(sigma2,n)
  sure_f <- function(nu0,just_sure = TRUE){
    # print(nu0)
    beta_mat <- matrix(data = NA,nrow = nchains,ncol = p)
    xbeta_var <- array(data = NA,dim = c(nchains,n,n))
    post_probs <- rep(NA,nchains)
    m_i <- list()
    for(i in 1:nchains){
      m_i[[i]] <- t_ab_post_moments_nu(x = x,y = y,nu0 = nu0,
                                       rand_samp = rand_samp[,,i],
                                       sigma_mat = sigma_mat)
      beta_mat[i,] <- m_i[[i]]$mean
      xbeta_var[i,,] <- x%*%m_i[[i]]$var2%*%t(x)
      post_probs[i] <- m_i[[i]]$post_prob
    }
    
    ## v1 
    # beta_gr <- t((post_probs %*% beta_mat)/sum(post_probs))
    # # xbeta_var_gr <- apply(xbeta_var,
    # #                       MARGIN = c(2,3),
    # #                       FUN = function(s)(post_probs^2 %*% s)/sum(post_probs)^2)
    # # var_xbeta <- (colMeans(((t(beta_mat)-c(beta_gr))*post_probs))/sum(post_probs))
    # 
    # ## v2:
    # xbeta_var_gr <- apply(xbeta_var,
    #                       MARGIN = c(2,3),
    #                       FUN = function(s)(post_probs %*% s)/sum(post_probs))
    # beta_difs <- (t(beta_mat)-c(beta_gr))
    # var_xbeta <- sum(diag(x%*%beta_difs %*% diag(post_probs/sum(post_probs)) %*% t(beta_difs) %*% t(x)))
    # 
    # 
    # x%*%(t(beta_mat) - c(beta_gr))^2*post_probs^2/sum(post_probs)^2
    # v3:

    beta_gr <- colMeans(beta_mat)
    xbeta_var_gr <- apply(xbeta_var,MARGIN = c(2,3),FUN= mean)
    var_xbeta <- sum(diag(var(beta_mat%*%t(x))))

    sure <- sum((y -x%*%beta_gr)^2) + 2*sum(diag(xbeta_var_gr)) + 2*var_xbeta
    if(just_sure){
      return(sure)
    } else {
      var_beta <- array(data = NA,dim = c(nchains,p,p))
      for(i in 1:nchains){
        var_beta[i,,] <- (m_i[[i]]$var2)  
      }
      var_beta <- apply(var_beta,MARGIN = c(2,3),FUN = mean)
      var_exp_beta <- ((var(beta_mat)))
      list(sure = sure,
           beta_estim = beta_gr,
           var_beta = var_beta + var_exp_beta,
           nu0 = nu0,
           mi = m_i)
    }
  }
  
  
  
  sure_vect <- lapply(nu_vect,sure_f,just_sure=FALSE)
  
  estims <- list(
    nu_vect,
    sure_vect
  )
  return(estims)
}









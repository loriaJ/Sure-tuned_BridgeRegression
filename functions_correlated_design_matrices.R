# modular prediction error functions

### functions to run the prediction schemes:

### MSE error comparison
# source('multivariate_sure_minimization_copy.R')
source('multivariate_sure_no_svd.R')

set.seed(2106)
library(BayesBridge)
library(tidyverse)
pred_error <- function(beta_v,alpha_tests,sigma2=1,n0,p,x0,type,descr_list,n_sims=100){

  n <- dim(x0)[1]
  
  y_obs <- x0 %*% beta_v + rnorm(0,n = n,sd=1)
  
  time_bayesbridge_stb <- vector(mode='list',length = length(alpha_tests))
  time_bayesbridge_tri <- vector(mode='list',length = length(alpha_tests))
  time_stb_exp_var <- vector(mode='list',length = length(alpha_tests))
  time_stb_exp_var_chains <- vector(mode='list',length = length(alpha_tests))

  # to save the results and compare outputs
  bayesBridge_stb <- vector(mode='list',length = length(alpha_tests))
  estims_post_exp_var <- vector(mode='list',length = length(alpha_tests))
  estims_post_exp_var_chains <- vector(mode='list',length = length(alpha_tests))
  bayesBridge_stb_tau <- rep(NA,length(alpha_tests))
  y <- y_obs 
  x <- x0 
  
  
  for(i in seq_along(alpha_tests)){
    time_bayesbridge_stb[[i]] <- system.time(
      bb_stb <- bridge.reg.stb(y = y,X = x,alpha = alpha_tests[i],
                               sig2.true = sigma2,nsamp = 10000,
                               ortho = FALSE) )
    
    bayesBridge_stb[[i]] <- colMeans(bb_stb$beta)
    bayesBridge_stb_tau[i] <- mean(bb_stb$tau)
    
    time_stb_exp_var[[i]]  <- system.time(
      estims_post_exp_var[[i]] <- t_multivariate_seq_sure(alpha = alpha_tests[i],
                                                          x = x,y = y,
                                                          n_sim = 10000,
                                                          sigma2 = sigma2,
                                                          nu_vect = 10^c(-10:3)))
    
    time_stb_exp_var_chains[[i]]  <- system.time(
      estims_post_exp_var_chains[[i]] <- t_multivariate_sure_schains(alpha = alpha_tests[i],
                                                                     x = x,y = y,
                                                                     n_sim = 1000,
                                                                     sigma2 = sigma2,nchains = 20,
                                                                     nu_vect = 10^c(-10:3)))
  }
  
  timing_df <- bind_rows(
    time_bayesbridge_stb %>% 
      map(`as.vector`) %>% 
      map(`[`,1:3) %>% 
      map(t) %>% 
      map_dfr(as.data.frame) %>% 
      mutate(alpha_col=alpha_tests,
             type='bayesbridge_stb'),
    time_stb_exp_var %>% 
      map(`as.vector`) %>% 
      map(`[`,1:3) %>% 
      map(t) %>% 
      map_dfr(as.data.frame) %>% 
      mutate(alpha_col=alpha_tests,
             type='stb_proposed')#,
    ) %>% 
    rename(user.self = V1,
           sys.self = V2,
           elapsed = V3)
  x1_list <- generate_x1s(n0,p = p,n_sims = n_sims,
                          type = type,descr_list = descr_list)

  post_mean_bbstb <- matrix(NA,nrow = length(alpha_tests),
                            ncol = n_sims)
  post_mean_stb <- matrix(NA,nrow = length(alpha_tests),
                          ncol = n_sims)
  post_mean_stb_chains <- matrix(NA,nrow = length(alpha_tests),
                          ncol = n_sims)
  sd_x0 <- sqrt(diag(var(x0)))
  for(j in 1:n_sims){
    x1 <- x1_list[[j]]
    x1_s <- x1 
    y_new <- x1 %*% beta_v + rnorm(n,sd = 1)
    
    y_new2 <- y_new
    for(i in seq_along(alpha_tests)){
      post_mean_bbstb[i,j] <- sum((y_new2 - x1_s%*% bayesBridge_stb[[i]])^2)
      post_mean_stb[i,j] <- sum((y_new2 - x1_s %*% estims_post_exp_var[[i]]$beta_estim)^2)
      post_mean_stb_chains[i,j] <- sum((y_new2 - x1_s %*% estims_post_exp_var_chains[[i]]$beta_estim)^2)
    }
    
  }
  
  shrinkage_param <- data.frame('nu_stb_proposed' = map_dbl(estims_post_exp_var,'nu'),
                                'tau_bayesbridge_stb' = bayesBridge_stb_tau,
                                'alpha' = alpha_tests
  )
  
  shrinkage_param$tau_stb_proposed <- shrinkage_param$nu_stb_proposed ^ (1/shrinkage_param$alpha)
  
  return(list(estims_post_exp_var = estims_post_exp_var,
              estims_post_exp_var_chains = estims_post_exp_var_chains,
              Pred_MSE = 
                data.frame(alpha = alpha_tests,
                           BB_stb = post_mean_bbstb,
                           Prop_stb = post_mean_stb,
                           Prop_stb_chains = post_mean_stb_chains
                           ),
              timing = timing_df,
              Shrinkage_parameter = shrinkage_param))
}



### Pred error:
generate_x1s <- function(n0,p,n_sims,type='mvnormal',descr_list = list(rho=corr1)){
  # types: 
  ## mvnormal
  ## wishart post
  ## factor
  x1_list <- list()
  if(type == 'mvnormal'){
    corr1 <- descr_list$rho
    for(i in 1:n_sims){
      
      x1_list[[i]] <- mvtnorm::rmvnorm(n0+1,
                                       sigma = diag(1-corr1,p) + corr1)
    }
    
  } else if (type == 'wishart'){
    sigma_w <- descr_list$sigma_w
    for(i in 1:n_sims){
      x1_list[[i]] <- mvtnorm::rmvnorm(n = n0+1,
                                       sigma = sigma_w)
    }
    
  } else if (type == 'factors'){
    k <- descr_list$k
    B <- descr_list$B
    sd_1 <- descr_list$sd_1
    sd_2 <- descr_list$sd_2
    
    
    for(i in 1:n_sims){
      F_matrix <- matrix(rnorm((n0+1) * k,mean = 0,sd=sd_1),
                         nrow = k,ncol = n0+1)
      
      x1_list[[i]] <- t(B %*% F_matrix) + matrix(rnorm(sd = sd_2,n = p*(n0+1)),nrow = n0+1)
    }
    
  }
  
  return(x1_list)
  
}

replicate_pred_error3 <- function(n,p,sigma2,alpha_tests = seq(0.3,1.9,by=0.1),
                                  n_rep = 100, b1 = 10,type='mvnormal',
                                  descr_list = list(rho = 0.9),signal_perc = 0.1){
  
  noise_perc <- 1-signal_perc 
  beta_v <- rnorm(mean=c(rep(0,noise_perc*p),rep(b1, signal_perc*p)),
                  n = p,sd = 0.5)
  beta_v <- sample(beta_v)
  if(type == 'wishart'){
    descr_list$sigma_w <- rWishart(1,p,diag(1,p))[,,1]
  }
  if(type =='factors'){
    # descr_list$k <- 8
    descr_list$B <- matrix(1,nrow = p,ncol = descr_list$k)
    
  }
  print(names(descr_list))
  x0 <- generate_x1s(n0 = n-1,p = p,n_sims = 1,type = type,
                     descr_list = descr_list)[[1]]
    print(dim(x0))

  list_results <- pred_error(alpha_tests = alpha_tests, sigma2 = sigma2, p =p, n0 = n-1,x0=x0,
                             type=type,descr_list=descr_list,
                             beta_v=beta_v,n_sims = n_rep)
  
  results <- list(list_results,
                  beta = beta_v)
  
  return(results)
  
}


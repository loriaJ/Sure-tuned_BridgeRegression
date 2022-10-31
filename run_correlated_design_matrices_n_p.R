# run functions design matrices:

## script to run the prediction schemes:
source('functions_correlated_design_matrices_n_p_cv.R') 
repsx1 <- list()


n <- 100 ## c(8,16,32,64,128,256,512)
p <- 50 ## c(3,10,100,200,500)
n_rep <- 1#00 # 

b1 <- 10 
cor10 <- 0.9
type <-  1
alpha <- 0.3 # seq(3,19,2)/10
nsim <-1
invisible(eval(parse(text=commandArgs(TRUE)))) 
set.seed(25+nsim) 
if(type == 1){
  type <- 'mvnormal'
  descr_list = list(rho = cor10) 
}
if(type == 2){
  type <- 'factors'
  descr_list = list(sd_1 = 0.5,sd_2=0.1,k=8)
}
if(type == 3 ){
  type <- 'wishart'
  descr_list = list()
}
repsx1 <- replicate_pred_error3(n,p,sigma2=1,n_rep = n_rep,b1=b1,
                                type = type,
                                descr_list =  descr_list,
                                alpha_tests = alpha)


if( type == 'mvnormal'){
  img_name <- paste('simulations_2022_09_02_n0',n,'p',p,
                    'b1',b1,'type_',type,'cor',cor10,'alpha',
                    alpha,'n_sim',nsim,'.RData',sep = '_')
} else{
  img_name <- paste('simulations_2022_09_02_n0',n,'p',p,
                    'b1',b1,'type_',type,'.RData',sep = '_')
}

save.image(img_name)




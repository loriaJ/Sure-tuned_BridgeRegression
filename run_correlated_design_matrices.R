# run functions design matrices:

## script to run the prediction schemes:
source('functions_correlated_design_matrices.R') 
repsx1 <- list()
# for(cor0 in 1:10){
set.seed(25)
n <- 100 ## c(8,16,32,64,128,256,512)
p <- 50 ## c(3,10,100,200,500)
n_rep <- 100 # 
#corss_list <- #0.6 # corr1 # 0:9/10
b1 <- 10 # c(0,0.5)
cor10 <- 0.9
type <-  1
alpha <- seq(3,19,2)/10
invisible(eval(parse(text=commandArgs(TRUE)))) 
# cor0 <- cor10
# if( type %in% c('mvnormal','factors','wishart'))
cor0 <- cor10
if(type == 1){
  type <- 'mvnormal'
  descr_list = list(rho = cor0) 
}
if(type == 2){
  type <- 'factors'
  descr_list = list(sd_1 = 0.5,sd_2=0.1,k=8)
  # descr_list$k <- 8
}
if(type == 3 ){
  type <- 'wishart'
  descr_list = list()
}
repsx1 <- replicate_pred_error3(n,p,sigma2=1,n_rep = n_rep,b1=b1,type = type,
                                descr_list =  descr_list,
                                # corr1 = cor0,
                                # alpha_tests = c(0.1,0.3),
                                alpha_tests = alpha)


repsx1[[1]]$Pred_MSE %>%
  gather(-alpha,key = 'Method',value = 'SSE') %>%
  separate(Method,into= c('Method','n_sim'),sep='[.]') %>%
  filter(Method != 'BBtri') %>%
  group_by(alpha,Method) %>%
  summarise(Mean_SSE = mean(SSE),
            SD_SSE = sd(SSE)) %>%
  ggplot(aes(x=alpha,y = Mean_SSE,color = Method)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(x=alpha,ymax = Mean_SSE + SD_SSE,ymin = Mean_SSE - SD_SSE))


if( type == 'mvnormal'){
  img_name <- paste('simulations_2022_02_15_n0',n,'p',p,
                    'b1',b1,'type_',type,'cor',cor10,'alpha',alpha,'.RData',sep = '_')
} else{
  img_name <- paste('simulations_2022_02_15_n0',n,'p',p,
                    'b1',b1,'type_',type,'.RData',sep = '_')
}

save.image(img_name)








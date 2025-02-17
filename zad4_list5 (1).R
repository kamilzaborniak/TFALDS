# EM algorithm

iter_task4 <- function(n,m, eps){
  Alternative <- rbinom(n,1,eps)
  lambdas<- 1/(Alternative*2.5 +3)
  sapply(lambdas, function(i){rexp(m,i)})-> samples
  
  T_test<- colSums(samples)
  #hist(T_test)
  
  # EM algorithm
  k<- 1
  eps_est<-numeric(500)
  eps_est[1]<-0.005
  mu_est<-numeric(500)
  mu_est[1]<-mean(T_test/m)
  d<-1
  
  while(d>10^(-5) & k<500){
    #print(k)
    pi<- dgamma(T_test, m, scale=mu_est[k])*eps_est[k]/
      (dgamma(T_test, m, scale=mu_est[k])*eps_est[k] +
         dgamma(T_test, m, rate=1/3)*(1-eps_est[k]))
    
    mu_est[k+1]<-sum(pi*(T_test/m))/sum(pi)
    #mu_est[1:5]
    eps_est[k+1]<- mean(pi)
    #eps_est[1:5]
    d<- (mu_est[k+1]-mu_est[k])^2+(eps_est[k+1]-eps_est[k])^2
    k<- k+1
  }
  
  data.frame(m=m, 
             n = n,
             eps = eps,
             eps_est = eps_est[k],
             mu_est = mu_est[k]) %>% 
    return()
}

sim4<-function(n,m,eps){
  plyr::rlply(1000,
              iter_task4(n,m,eps)) %>% 
    do.call("rbind", .) %>% as.data.frame() %>%
    return()
}

sapply(N, function(n){
  sapply(M, function(m){
    sapply(E, function(e){
      sim4(n,m,e)
    }, simplify = F) %>% 
      do.call("rbind",.)
  }, simplify = F) %>% 
    do.call("rbind",.)
}, simplify = F) %>% 
  do.call("rbind",.) -> zad4

ggplot2::ggplot(zad4, aes(x=mu_est))+
  geom_histogram()+
  facet_grid(vars(eps), vars(n,m), labeller ="label_both", scales = "free")+
  geom_vline(xintercept=5.5, colour="red")+
  xlab(TeX("$\\hat{\\mu}$"))->plt1
ggplot2::ggplot(zad4,aes(x=eps_est))+
  geom_histogram()+
  facet_grid(vars(eps), vars(n,m), labeller ="label_both", scales = "free")+
  geom_vline(aes(xintercept=eps), colour="red")+
  xlab(TeX("$\\hat{\\epsilon}$"))->plt2


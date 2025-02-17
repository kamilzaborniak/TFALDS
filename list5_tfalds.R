library(mutoss)
library(ggplot2)
library(plyr)
library(dplyr)

{
 # n<-200 #i
 # m<-20  #i
 # e<-eps<- 0.05 #ii
 # q<- 0.1*sqrt(200/m) #i

BFDR <- function(c_bfdr, eps, size_sample, alpha){
    ((1 - eps) * (1 - pgamma(c_bfdr, size_sample, 1/3))/
      (1 - ((1 - eps) * pgamma(c_bfdr, size_sample, 1/3) + 
              eps * pgamma(c_bfdr, size_sample, 1/5.5))) - alpha) %>%
      return()
  }
p_val<- function(sample){
  m<-length(sample)
  t<- sum(sample)
  p<- 1-pgamma(test_statistics, m, 1/3)
  return(p)
}
#c_BFDR <- uniroot(BFDR,c(qgamma(0.95, m, 1/3)*0.1,qgamma(0.95,m, 1/3)*2), e = eps, m =m, q=q)$root
# b_treshold<- BFDR_threshold<- uniroot(BFDR, c(qgamma(0.95, m, 1/3)*0.1,
#                                    qgamma(0.95, m, 1/3)*2), 
#                          eps=e, size_sample=m, alpha=q)$root

iter<-function(n,m,eps,q, b_treshold){
  Alternative<- rbinom(n,1,eps)
  index_A<- which(Alternative>0)
  index_0<- which(Alternative==0)
  lambdas<- 1/(Alternative*2.5 +3)
  sapply(lambdas, function(i){rexp(m,i)})-> samples
  test_stat<- sapply(1:n, function(i) sum(samples[,i]))
  p_vals<- sapply(1:n, function(i)(1-pgamma(test_stat[i], m, 1/3))) %>% unlist()
  
  bonf_test<- bonferroni(p_vals, q, silent = TRUE)$rejected
  R_bonf<- sum(bonf_test)
  TD_bonf<- sum(bonf_test[index_A])
  FD_bonf<- sum(bonf_test[index_0])
  FnoD_bonf<- sum(bonf_test[index_A]*-1+1)
  FDP_bonf<- ifelse(sum(R_bonf)>0, FD_bonf/R_bonf, 0)
  TDP_bonf<- ifelse(sum(Alternative)>0, TD_bonf/ sum(Alternative), 0)
  
  bh_test<- BH(p_vals, q, silent = TRUE)$rejected
  R_bh<- sum(bh_test)
  TD_bh<- sum(bh_test[index_A])
  FD_bh<- sum(bh_test[index_0])
  FnoD_bh<- sum(bh_test[index_A]*-1+1)
  FDP_bh<- ifelse(sum(R_bh)>0, FD_bh/R_bh, 0)
  TDP_bh<- ifelse(sum(Alternative)>0, TD_bh/ sum(Alternative), 0)
  
  BFDR_test<- test_stat > b_treshold
  R_BFDR_test<- sum(BFDR_test)
  TD_BFDR_test<- sum(BFDR_test[index_A])
  FD_BFDR_test<- sum(BFDR_test & !Alternative)
  FnoD_BFDR_test<- sum(BFDR_test[index_A]*-1+1)
  FDP_BFDR_test<- FD_BFDR_test/R_BFDR_test
  FDP_BFDR_test<- ifelse(is.nan(FDP_BFDR_test),0,FDP_BFDR_test)
  TDP_BFDR_test<- ifelse(sum(Alternative)>0, TD_BFDR_test/ sum(Alternative), 0)
  
  df<- data.frame(FD_bonf, FDP_bonf, FnoD_bonf, TDP_bonf,
                  FD_bh, FDP_bh, FnoD_bh, TDP_bh,
                  FD_BFDR_test, FDP_BFDR_test, FnoD_BFDR_test, TDP_BFDR_test)
  return(df)
}

sim<- function(no_test=200,size_smpl=20,epsilon=0.1,q=0.01){
  #df<- data.frame(n,m,eps,q)
  treshold<- stats::uniroot(BFDR, c(qgamma(0.95, size_smpl, 1/3)*0.1,
                             qgamma(0.95, size_smpl, 1/3)*2), 
                     eps=epsilon,
                     size_sample=size_smpl, alpha=q)$root
  set.seed(410)
  rlply(1000,
        iter(n=no_test, m=size_smpl,eps=epsilon,q=q, b_treshold=treshold)) %>% 
    do.call("rbind", .) %>% colMeans()%>% t() %>% as.data.frame()  -> result

  meanFD_bonf=result$FD_bonf
  meanFnoD_bonf=result$FnoD_bonf
  cost1_bonf= meanFD_bonf+ meanFnoD_bonf
  cost2_bonf= 2*meanFD_bonf+ meanFnoD_bonf
  cost3_bonf= meanFD_bonf+ 2*meanFnoD_bonf
  
  meanFD_bh=result$FD_bh
  meanFnoD_bh=result$FnoD_bh
  cost1_bh= meanFD_bh+ meanFnoD_bh
  cost2_bh= 2*meanFD_bh+ meanFnoD_bh
  cost3_bh= meanFD_bh+ 2*meanFnoD_bh
  
  meanFD_BFDR_test=result$FD_BFDR_test
  meanFnoD_BFDR_test=result$FnoD_BFDR_test
  cost1_BFDR_test= meanFD_BFDR_test+ meanFnoD_BFDR_test
  cost2_BFDR_test= 2*meanFD_BFDR_test+ meanFnoD_BFDR_test
  cost3_BFDR_test= meanFD_BFDR_test+ 2*meanFnoD_BFDR_test
  
  data.frame(n=rep(no_test,3), m=rep(size_smpl,3), eps=rep(epsilon,3), q=rep(q,3),
             test= c("Bonferroni", "BH", "BFDR test"),
             FDR=c(result$FDP_bonf, result$FDP_bh, result$FDP_BFDR_test), 
             Power=c(result$TDP_bonf, result$TDP_bh, result$TDP_BFDR_test),
             cost1=c(cost1_bonf,cost1_bh,cost1_BFDR_test), 
             cost2=c(cost2_bonf, cost2_bh,cost2_BFDR_test),
             cost3=c(cost3_bonf, cost3_bh, cost3_BFDR_test),
             row.names = NULL) -> df 
    return(df)
}

N<-c(200,1000)
M<-c(20,100)
E<-c(0.01, 0.05, 0.1, 0.2)
#Q<-c(0.1, 0.1*sqrt(200/m))
}
#sim(200,20,0.01,0.1)
############################################## ZAD 2
{
sapply(N, function(n){
  sapply(M, function(m){
    sapply(E, function(e){
      Q<-c(0.1, 0.1*sqrt(200/m))
      sapply(Q, function(q){
        sim(no_test=n, size_smpl=m, epsilon=e, q=q)},
             simplify = F) %>% 
        do.call("rbind",.)},
      simplify = F) %>% 
      do.call("rbind",.)},
    simplify = F) %>%
    do.call("rbind",.)
}, simplify = F) %>% do.call("rbind",.) -> zad2
}

ggplot(zad2, aes(x=as.factor(eps), y=FDR, colour=test))+
  geom_point()+
  facet_grid(vars(n), vars(m,q=round(q,3)), labeller = "label_both")+
  geom_hline(aes(yintercept = q))+
  xlab("e")

ggplot(zad2, aes(x=as.factor(eps), y=Power, colour=test))+
  geom_point()+
  facet_grid(vars(n), vars(m,q=round(q,3)), labeller = "label_both")+
  xlab("e")

####################################### ZAD 3
{
  tau<- function(eps, c0, cA){
    t<- 33/5 * log(11/6 * c0/cA *(1-eps)/eps)
    return(t)
  }
  power_exact <- function(m,tau_val) {
    1 - pgamma(tau_val,  m, 1 / 5.5) %>%
      return()
  }
  BFDR_exact <- function(m, eps, tau_val){
    (1 - pgamma(tau_val, m, 1/3)) * (1 - eps) / 
      ((1 - pgamma(tau_val, m, 1/3) * (1 - eps)) + 
         ((1 - pgamma(tau_val, m, 1 / 5.5)) * eps)) %>%
      return()
  }
  cost_bfdr_exact<- function(m,eps,c0,cA,tau_val){
    c0 * (1 - eps) * (1 - pgamma(tau_val, m, 1/3)) + 
      cA * eps * (1 - pgamma(tau_val, m, 1/5.5)) %>%
    return()
  }
  sim_t3<- function(n,size_sample,e,c0,cA){
    tauuuu<- tau(eps=e,c0,cA)
    cost<-cost_bfdr_exact(m=size_sample,
                          eps=e,
                          c0,
                          cA,
                          tau_val = tauuuu)
    data.frame(n=n,
      m=c(size_sample), eps=c(e), c0=c(c0), 
               cA=c(cA), #test=c("Bayesian classifier"), 
               tau=c(tauuuu), 
               BFDR=c(BFDR_exact(m=size_sample,
                                 eps=e,tau_val = tauuuu)), 
               Power=c(power_exact(m=size_sample,tau_val=tauuuu)), 
               cost=c(cost),
               cost_total=n*cost)%>%
      return()
  }
  
  }


####### task 4

iter4<-function(n,m,eps){
  samples<- replicate(n, rexp(m,1/3))
  alternatives<-(rbinom(n, 1, eps))
  id_alter<- which(alternatives>0) # tam gdzie 5.5 mamy zapisane w alternatives
  sapply(id_alter, function(i){samples[,i]<-rexp(m,1/5.5)})
  
  T_test <- sapply(1:n, function(i)sum(samples[,i])) # test statistics pochodza z rozkladu gamma(m, 1/3) dla H_0
  
  k<-0
  epsi0<-0
  epsi1<-0.005
  mu0<-0
  mu1<-3
  d<- (mu0-mu1)^2+(epsi0-epsi1)^2
  
  while(d>10^(-5) & k<5000){
    ####expectation step
    epsi0=epsi1
    mu0=mu1
    pi=dgamma(T_test, m, 1/mu1)*epsi1/
      (dgamma(T_test, m, 1/mu1)*epsi1 + dgamma(T_test, m, 1/3)*(1-epsi1))
    
    mu0= mu1
    mu1<-sum(pi*T_test)/(sum(pi)*m)
    epsi1<- mean(pi)
    
    d<- (mu0-mu1)^2+(epsi0-epsi1)^2
    k<- k+1
    #print(d)
  }
  
  data.frame(m=m, 
             n = n,
             eps = eps,
             eps_est = epsi1,
             mu_est = mu1) %>% 
  return()
  
}
sim4<-function(n,m,eps){
  rlply(100,
        iter4(n,m,eps)) %>% 
    do.call("rbind", .) %>% as.data.frame() %>%
    return()
}

sapply(N, function(n){
  sapply(M, function(m){
    sapply(E, function(eps){
            sim4(n,m,eps)
    }, simplify = F)%>% 
      do.call("rbind",.)
  }, simplify = F)%>% 
    do.call("rbind",.)
}, simplify = F)%>% 
  do.call("rbind",.) -> zad4

ggplot(zad4, aes(x=mu_est))+
  geom_histogram()+
  facet_grid(vars(n,m), vars(eps), labeller ="label_both", scales = "free_x")->plt1
ggplot(zad4,aes(x=eps_est))+
  geom_histogram()+
  facet_grid(vars(n,m), vars(eps), labeller ="label_both", scales = "free_x")->plt2


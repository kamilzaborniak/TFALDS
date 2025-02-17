library(Sim.DiffProc)
library(mutoss)

p_val<- function(sample){
  return(c(1-pnorm(sample)))
}

zad1<- function(p_vec, alpha=0.05, no.Tr.H1){
  n<- length(p_vec)
  
  df_result<- data.frame(test=c(), FWER=c(), FDR=c(), Power=c())
  bnf <- bonferroni(p_vec, alpha, silent=TRUE)$rejected
  hlm <- holm(p_vec, alpha, silent=TRUE)$rejected
  hcb <- hochberg(p_vec, alpha, silent=TRUE)$rejected
  sdk <- sidak(p_vec, alpha, silent=TRUE)$rejected
  bhq <- mutoss::BH(p_vec, alpha, silent=TRUE)$rejected
  
  df_result<- data.frame(V_bnf=bnf[(no.Tr.H1+1):n]%>% sum() >=1,
                         FDP_bnf= ifelse(sum(bnf)>0, sum(bnf[(no.Tr.H1+1):n])/sum(bnf), 0),
                         P_bnf= bnf[1:no.Tr.H1]%>% mean(),
                         
                         V_hlm=hlm[(no.Tr.H1+1):n]%>% sum() >=1,
                         FDP_hlm= ifelse(sum(hlm)>0, sum(hlm[(no.Tr.H1+1):n])/sum(hlm), 0),
                         P_hlm= hlm[1:no.Tr.H1]%>% mean(),
                         
                         V_hcb=hcb[(no.Tr.H1+1):n]%>% sum() >=1,
                         FDP_hcb= ifelse(sum(hcb)>0, sum(hcb[(no.Tr.H1+1):n])/sum(hcb), 0),
                         P_hcb= hcb[1:no.Tr.H1]%>% mean(),
                         
                         V_sdk=sdk[(no.Tr.H1+1):n]%>% sum() >=1,
                         FDP_sdk=ifelse(sum(sdk)>0, sum(sdk[(no.Tr.H1+1):n])/sum(sdk), 0),
                         P_sdk=sdk[1:no.Tr.H1]%>% mean(),
                         
                         V_bhq= bhq[(no.Tr.H1+1):n]%>% sum() >=1,
                         FDP_bhq= ifelse(sum(bhq)>0, sum(bhq[(no.Tr.H1+1):n])/sum(bhq), 0),
                         P_bhq= bhq[1:no.Tr.H1]%>% mean()
  )
  return(df_result)
}

zad3<- function(sample, alpha=0.05, no.Tr.H1){
  # strong sense
  n<- length(sample)
  p_vec<- sample %>% p_val()
  # 2 step Fisher's
      # 1. step - global null test
  if (min(p_vec)<=alpha/n) { Bonf_2st_test<- p_vec<=alpha # 2. step if glob. H_0 rejected, test each hypothesis at lvl. alpha
  } else { Bonf_2st_test<- p_vec*0 }# if glob. null is true, all H_0i are true
   
  # Chi^2 test
  if (sum(sample^2)>=qchisq(1-alpha, n)) { chi2_2st_test<- p_vec<=alpha 
  } else {chi2_2st_test<- p_vec*0}
  
  #weak sense
  smpl<- rnorm(n)
  p_v<- smpl %>% p_val()
  
  if (min(p_v)<=alpha/n) { Bonf_2st_test_H0<- p_v<=alpha # 2. step if glob. H_0 rejected, test each hypothesis at lvl. alpha
  } else { Bonf_2st_test_H0<- p_v*0 }# if glob. null is true, all H_0i are true
  
  if (sum(smpl^2)>=qchisq(1-alpha, n)) { chi2_2st_test_H0<- p_v<=alpha 
  } else {chi2_2st_test_H0<- p_v*0}
  
  result<- data.frame(V_bnf  = Bonf_2st_test[(no.Tr.H1+1):n]%>% sum() >=1,
                      V_bnf_weak= Bonf_2st_test_H0%>% sum() >=1,
                      FDP_bnf= ifelse(sum(Bonf_2st_test)>0,
                                      sum(Bonf_2st_test[(no.Tr.H1+1):n])/sum(Bonf_2st_test),
                                      0),
                      P_bnf  = Bonf_2st_test[1:no.Tr.H1]%>% mean(),
                      
                      V_chi2 = chi2_2st_test[(no.Tr.H1+1):n]%>% sum() >=1,
                      V_chi2_weak= chi2_2st_test_H0 %>% sum() >=1,
                      FDP_chi2= ifelse(sum(chi2_2st_test)>0,
                                      sum(chi2_2st_test[(no.Tr.H1+1):n])/sum(chi2_2st_test),
                                      0),
                      P_chi2  = chi2_2st_test[1:no.Tr.H1]%>% mean()
                      )
  return(result)
}

rep_z3<- function(n, means, alpha, no.Tr.H1){
  rlply(1000, {rnorm(n, means) %>% 
      zad3(alpha , no.Tr.H1)}) %>% 
    do.call("rbind", .) %>%
    colMeans() %>%
    return()
}

#zad1
{n<-20
means_1a<- c(1.2*sqrt(2*log(n)), rep(0,n-1))
rlply(1000, {rnorm(n, means_1a) %>% 
    p_val() %>% 
    zad1(alpha = 0.05, 1)}) %>% 
  do.call("rbind", .) -> zad1a
colMeans(zad1a)

means_1b<- c(rep(1.02*sqrt(2*log(n/10)),5), rep(0,n-5))
rlply(1000, {rnorm(n, means_1b) %>% 
    p_val() %>% 
    zad1(alpha = 0.05, 5)}) %>% 
  do.call("rbind", .) -> zad1b
colMeans(zad1b)

means_1c<- c(sqrt(2*log(20/(1:10))), rep(0,n-10))
rlply(1000, {rnorm(n, means_1c) %>% 
    p_val() %>% 
    zad1(alpha = 0.05, 10)}) %>% 
  do.call("rbind", .) -> zad1c
colMeans(zad1c)
}
# zad2
{
n<- 5000

means_2a<- c(1.2*sqrt(2*log(n)), rep(0,n-1))
rlply(1000, {rnorm(n, means_2a) %>% 
    p_val() %>% 
    zad1(alpha = 0.05, 1)}) %>% 
  do.call("rbind", .) -> zad2a
colMeans(zad2a)

means_2b<- c(rep(1.02*sqrt(2*log(n/200)),100), rep(0,n-100))
rlply(1000, {rnorm(n, means_2b) %>% 
    p_val() %>% 
    zad1(alpha = 0.05, 100)}) %>% 
  do.call("rbind", .) -> zad2b
colMeans(zad2b)

means_2c<- c(rep(sqrt(2*log(n/200)),100), rep(0,n-100))
rlply(1000, {rnorm(n, means_2c) %>% 
    p_val() %>% 
    zad1(alpha = 0.05, 100)}) %>% 
  do.call("rbind", .) -> zad2c
colMeans(zad2c)

means_2d<- c(rep(1.002*sqrt(2*log(n/2000)),1000), rep(0,n-1000))
rlply(1000, {rnorm(n, means_2d) %>% 
    p_val() %>% 
    zad1(alpha = 0.05, 1000)}) %>% 
  do.call("rbind", .) -> zad2d
colMeans(zad2d)
}
#zad3
{
  n<-20
  #a)
  means3a_20<- c(1.2*sqrt(2*log(n)), rep(0,n-1))
  rep_z3(n,means3a_20,0.05,1)
  #b)
  means3b_20<- c(rep(1.02*sqrt(2*log(n/10)),5), rep(0,n-5))
  rep_z3(n,means3b_20,0.05,5)
  #c)
  means3c_20<- c(sqrt(2*log(20/(1:10))), rep(0,n-10))
  rep_z3(n,means3c_20,0.05,10)
  #d)
  means3d_20<- c(rep(1.002*sqrt(2*log(n/8)),4), rep(0,n-4))
  rep_z3(n,means3d_20,0.05,4)

  n<-5000
  #a)
  means3a_5000<- c(1.2*sqrt(2*log(n)), rep(0,n-1))
  rep_z3(n,means3a_5000,0.05,1)
  #b)
  means3b_5000<- c(rep(1.02*sqrt(2*log(n/10)),5), rep(0,n-5))
  rep_z3(n,means3b_5000,0.05,5)
  #c)
  means3c_5000<- c(sqrt(2*log(20/(1:10))), rep(0,n-10))
  rep_z3(n,means3c_5000,0.05,10)
  #d)
  means3d_5000<- c(rep(1.002*sqrt(2*log(n/2000)),1000), rep(0,n-1000))
  rep_z3(n,means3d_5000,0.05,1000)
  
  
  
}

# task 4
{
  n<- 5000
  no_trajectories<- 1000

  U<-function(p_vec){
      n<- length(p_vec)
      t<- seq(0,1, by=1/n)
      CDF<- ecdf(p_vec)
      
      sapply(t, function(i){
                sqrt(n)* (CDF(i)-i)
          }) %>%
        data.frame(t=t,
                   Process=rep("U",n+1),
                   Trajectory=.) %>%
      return()
    }
  B<-function(n=5000){
      BB(n,t0=0,T=1, x0=0,y=0) %>%
      c() %>%
      data.frame(t=seq(0,1,by=1/n),
                 Process=rep("B",n+1),
                 Trajectory=.) %>%
      return()
  }
  iter<- function(n=5000){
    rnorm(n) %>% 
      p_val() %>%
      U() -> u
    B(n) -> b
    return(rbind(u,b))
  }
  sapply(1:5, 
         function(j) data.frame(i_0=rep(j,n+1)) %>%
           cbind( iter(n=5000)),
         simplify = F) %>% 
    t() %>%
    do.call("rbind", .) -> zad4
  
  
  ggplot(zad4, aes(x=t, y=Trajectory, group=i_0, col=Process))+
    geom_path()+
    ylab("Value")
  
  replicate(1000, {rnorm(5000) %>%
                    p_val() %>%
                    U()-> df; df$Trajectory %>%
                    abs() %>% max()}) %>% 
    quantile(probs = c(0.8, 0.9, 0.95)) -> U_quant
  replicate(1000, {BB(n,t0=0,T=1, x0=0,y=0) %>% 
                    c() %>% abs() %>% 
                    max()}) %>% 
    quantile(probs = c(0.8, 0.9, 0.95)) -> BB_quant
  
}












p_vec<- rep(0.05, 20)

holm(p_vec, alpha, silent=TRUE)$rejected
hochberg(p_vec, alpha, silent=TRUE)$rejected

### działa tak jak ma działać i jak było mówione na wykładzie







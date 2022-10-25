#Change of notation:
#notations in the paper:     rho1 rho2   alpha1 alpha2 
#notations used in the code: rho0 alpha0 rho1   alpha1 

#Sample size calculation
sample_size<-function(eff,t1,t2,n1,n2,rho0,rho1,alpha0,alpha1,sigma1){
  alpha=0.05
  beta=0.2
  kappa=sqrt(alpha1*rho1)
  sigma2=sigma1*rho1/alpha1
  vb<-sigma1/n1/t1*((1+(n1-1)*rho0)+n1*(t1-1)*rho1)+sigma2/n2/t2*((1+(n2-1)*alpha0)+n2*(t2-1)*alpha1)-2*sqrt(sigma1*sigma2)*sqrt(rho1*alpha1)
  site <- ceiling((qnorm(1-alpha/2) + qnorm(1-beta))^2*vb/(eff^2))
  return(site)
}

#power calculation
power<-function(site,eff,t1,t2,n1,n2,rho0,rho1,alpha0,alpha1,sigma1){
  alpha=0.05
  kappa=sqrt(alpha1*rho1)
  sigma2=sigma1*rho1/alpha1
  vb<-sigma1/n1/t1/site*((1+(n1-1)*rho0)+n1*(t1-1)*rho1)+sigma2/n2/t2/site*((1+(n2-1)*alpha0)+n2*(t2-1)*alpha1)-2*sqrt(sigma1*sigma2)*sqrt(rho1*alpha1)/site
  pw<-pnorm(eff/sqrt(vb)-qnorm(1-alpha/2))
  return(pw)
}

#############################################################################
#simulation for predicted power

#data generation
data_gen<-function(site=5,eff=0.2,t1=2,t2=3,n1=10,n2=13,rho0=0.05,rho1=0.03,alpha0=0.05,alpha1=0.02,sigma1=1){
  #sigma2 is determined by sigma1, alpha1 and rho1
  sigma2=sigma1*rho1/alpha1
  intercept=0.5
  N=site*(t1*n1+t2*n2)
  patid <- 1:N
  trt<-rep(c(rep(1,t1*n1),rep(0,t2*n2)),site)
  sigma_lambda<-sigma1*rho1
  sigma_theta1<-rho0*sigma1-sigma_lambda
  sigma_theta2<-alpha0*sigma2-sigma_lambda
  
  sigma_e1<-sigma1-sigma_lambda-sigma_theta1
  sigma_e2<-sigma2-sigma_lambda-sigma_theta2
  
  strata<-rep(1:site,each=(t1*n1+t2*n2))
  tid_s<-c(rep(1:t1,each=n1),rep((t1+1):(t1+t2),each=n2))
  tid<-c()
  for(i in 1:site){
    tid<-c(tid,((i-1)*(t1+t2)+tid_s))
  }
  
  site_effect<-rep(rnorm(site,0,sqrt(sigma_lambda)),each=(t1*n1+t2*n2))
  
  therapist_effect1<-rnorm(site*(t1+t2),0,sqrt(sigma_theta1))
  therapist_effect2<-rnorm(site*(t1+t2),0,sqrt(sigma_theta2))
  
  therapist_effect<-rep(therapist_effect1,as.numeric(table(tid)))*trt+rep(therapist_effect2,as.numeric(table(tid)))*(1-trt)
  
  epsilon<-c()
  for (i in 1:site){
    epsilon_s <-c(rnorm(t1*n1,0,sqrt(sigma_e1)),rnorm(t2*n2,0,sqrt(sigma_e2)))
    epsilon<-c(epsilon,epsilon_s)
  }
  
  y=intercept-eff*trt+site_effect+therapist_effect+epsilon
  y0=intercept+site_effect+therapist_effect+epsilon
  
  df<-data.frame(y=y,y0=y0,trt=trt,tid=tid,sid=strata)
  return(df)
}


#simulation function
sim_func<-function(eff=0.2,t1=2,t2=2,n1=10,n2=10,rho0=0.05,rho1=0.03,alpha0=0.05,alpha1=0.02,sigma1=1,nsim=1000,seed=1112){
  require(nlme)
  #model fitting
  set.seed(seed)
  pval<-rep(NA,nsim)
  pval0<-rep(NA,nsim)
  site<-ceiling(sample_size(eff=eff,t1=t1,t2=t2,n1=n1,n2=n2,rho0=rho0,rho1=rho1,alpha0=alpha0,alpha1=alpha1,sigma1=sigma1))
  pw<-power(site=site,eff=eff,t1=t1,t2=t2,n1=n1,n2=n2,rho0=rho0,rho1=rho1,alpha0=alpha0,alpha1=alpha1,sigma1=sigma1)
  
  #count how many singular fitting cases
  count.singular <- 0
  
  for (s in 1:nsim){
    
    exit <- FALSE
    while (exit==FALSE){
      df<-data_gen(eff=eff,site=site,t1=t1,t2=t2,n1=n1,n2=n2,rho0=rho0,rho1=rho1,alpha0=alpha0,alpha1=alpha1,sigma1=sigma1)
      fit0 = try(lme(y0 ~ trt, random = list(sid = ~ 1, tid = ~trt-1),weights = varIdent(form= ~ 1 | trt),
                     data=df),silent=T)
      fit1 = try(lme(y ~ trt, random = list(sid = ~ 1, tid = ~trt-1),weights = varIdent(form= ~ 1 | trt),
                     data=df),silent=T)
      if(class(fit0)!="try-error"&class(fit1)!="try-error"){
        exit <- TRUE
        
        count.singular <- count.singular + 1
      }
    }
    
    #z-test
    testan <- coef(summary(fit0))[2,1]/coef(summary(fit0))[2,2]
    pval0[s] <- (min((1-pnorm(testan)),pnorm(testan)))*2
  
    testaa <- coef(summary(fit1))[2,1]/coef(summary(fit1))[2,2]
    pval[s] <- (min((1-pnorm(testaa)),pnorm(testaa)))*2
    
  }
  
  count.singular <- count.singular - nsim
  tp1<-round(mean(pval0<0.05,na.rm=T) , 3)
  empower<-round(mean(pval<0.05,na.rm=T) , 3)
  pw <- round(pw,3)
  
  return(c(rho0,alpha0,rho1,alpha1,site,tp1,pw,empower,count.singular))
}

rho0 <- c(0.02,0.04,0.08,0.20)
T1 <- c(2,4,6,10)
T2 <- 1

#table 3 column 1
fidata <- c()
for (i in 1:length(rho0)){
    for (u in 1:length(T1)){
      fidata <- rbind(fidata,sim_func(eff=0.2,t1=T1[u],t2=1,n1=60/T1[u],n2=60,
                                      rho0=rho0[i],rho1=0.25*rho0[i],alpha0=0.5*rho0[i],alpha1=0.5*rho0[i],
                                      sigma1=1,nsim=5000,seed=920784642))
      
  }
}
#table 3 column 2
fidata <- c()
for (i in 1:length(rho0)){
  for (u in 1:length(T1)){
    fidata <- rbind(fidata,sim_func(eff=0.2,t1=T1[u],t2=1,n1=60/T1[u],n2=60,
                                    rho0=rho0[i],rho1=0.25*rho0[i],alpha0=0.75*rho0[i],alpha1=0.75*rho0[i],
                                    sigma1=1,nsim=5000,seed=920784642))
    
  }
}






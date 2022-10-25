#Sample size calculation
sample_size<-function(eff,t1,t2,n1,n2,rho0,alpha0,sigma1,sigma2){

  alpha=0.05
  beta=0.1
  rho1 <- alpha0*sigma2/sigma1
  alpha1=alpha0
  vb<-sigma1/n1/t1*((1+(n1-1)*rho0)+n1*(t1-1)*rho1)+sigma2/n2/t2*((1+(n2-1)*alpha0)+n2*(t2-1)*alpha1)-2*sqrt(sigma1*sigma2)*sqrt(rho1*alpha1)
  site <- ceiling((qnorm(1-alpha/2) + qnorm(1-beta))^2*vb/(eff^2))
  #return(n0)
  #return(site)
  percentage <- rho1/rho0
  return(c(rho0,alpha0,rho1,alpha1,t1,n1,t2,n2,percentage,site))
}

alpha0<- c(rep(0.01,4),rep(0.02,4),rep(0.04,3),rep(0.08,2))
rho0 <- c(0.02,0.05,0.1,0.2,0.04,0.05,0.10,0.20,0.08,0.10,0.20,0.16,0.20)
T1 <- 48
T2 <- 1
N1 <- 16.16667
N2 <- 776
sigma1 <- 0.0169
sigma2 <- c(0.10^2,0.11^2,0.12^2,0.13^2,0.14^2)
eff <- 0.008

fidata <- c()
for (k in 1:length(sigma2)){
    for (i in 1:length(rho0)){
      fidata <- rbind(fidata,sample_size(eff=eff,t1=T1,t2=T2,n1=N1,n2=N2,
                                         rho0=rho0[i],alpha0=alpha0[i],sigma1=sigma1,sigma2=sigma2[k]))
    }
}
fidata

colnames(fidata) <- c("rho0","alpha0","rho1","alpha1","T1","N1","T2","N2","Percentage","N.sites")
fidata

setwd("/Users/deckard/Desktop/Fan Li Project/Project 3/simulation result")
write.csv(fidata,"CODA_NEW.csv")


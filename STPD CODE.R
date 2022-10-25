#Sample size calculation
sample_size<-function(eff,t1,t2,n1,n2,rho0,alpha0,rho1,sigma1,sigma2){
  
  alpha=0.05
  beta=0.1

  alpha1=rho1*sigma1/sigma2
  vb<-sigma1/n1/t1*((1+(n1-1)*rho0)+n1*(t1-1)*rho1)+sigma2/n2/t2*((1+(n2-1)*alpha0)+n2*(t2-1)*alpha1)-2*sqrt(sigma1*sigma2)*sqrt(rho1*alpha1)
  site <- ceiling((qnorm(1-alpha/2) + qnorm(1-beta))^2*vb/(eff^2))
  #return(n0)
  #return(site)
  percentage <- rho1/rho0
  return(c(rho0,alpha0,rho1,alpha1,t1,n1,t2,n2,percentage,site))
}

#############################################################################

#objective loss function
obf <- function(J,t1,t2,n1,n2,rho0,alpha0,rho1,sigma1,sigma2){

  alpha1=rho1*sigma1/sigma2
  f <- sigma1*(1+(n1-1)*rho0+n1*(t1-1)*rho1)/t1/n1+sigma2*(1+(n2-1)*alpha0+n2*(J-t1-1)*alpha1)/(J-t1)/n2
  return(f)
}

#optimal number of therapists
opthera <- function(eff,J,n1,n2,rho0,alpha0,rho1,sigma1,sigma2){

  alpha1=rho1*sigma1/sigma2
  T1op <- J*sqrt(sigma1)*sqrt(n2*(1-rho0+n1*(rho0-rho1)))/
    ((sqrt(sigma1)*sqrt(n2*(1-rho0+n1*(rho0-rho1))))+(sqrt(sigma2*n1*(1-alpha0+n2*(alpha0-alpha1)))))
  T1op1 <- ceiling(T1op)
  T1op2 <- T1op1 - 1
  obf1 <- obf(J,T1op1,(J-T1op1),n1,n2,rho0,alpha0,rho1,sigma1,sigma2)
  obf2 <- obf(J,T1op2,(J-T1op2),n1,n2,rho0,alpha0,rho1,sigma1,sigma2)
  if (obf1 >= obf2){
    T1op <- T1op2
  } else {T1op <- T1op1}
  T2op <- J-T1op
  opsite <- sample_size(eff=eff,t1=T1op,t2=T2op,n1=n1,n2=n2,
                        rho0,alpha0,rho1,sigma1,sigma2)[10]
  return(c(T1op,T2op,opsite))
  
}

rho1 <- 0.007
rho0<- c(rep(0.01,4),rep(0.05,4),rep(0.10,4),rep(0.20,4))
alpha0 <- rep(c(0.01,0.05,0.10,0.20),4)
T1 <- 66
T2 <- 79
N1 <- 2.2
N2 <- 1.7
sigma1 <- 169
sigma2 <- c(11^2,13^2,14^2)
eff <- 1.3

fidata <- c()
for (k in 1:length(sigma2)){
  for (i in 1:length(rho0)){
    fidata <- rbind(fidata,sample_size(eff=eff,t1=T1,t2=T2,n1=N1,n2=N2,
                                       rho0=rho0[i],alpha0=alpha0[i],rho1=rho1,sigma1=sigma1,sigma2=sigma2[k]))
  }
}
fidata

sample_size(eff=eff,t1=77,t2=68,n1=N1,n2=N2,
            rho0=rho0[1],alpha0=alpha0[3],rho1=rho1,sigma1=sigma1,sigma2=sigma2[1])

J <- T1+T2
opdata <- c()
for (k in 1:length(sigma2)){
  for (i in 1:length(rho0)){
      opdata <- rbind(opdata,opthera(eff=eff,J=J,n1=N1,n2=N2,
                                     rho0=rho0[i],alpha0=alpha0[i],rho1=rho1,sigma1=sigma1,sigma2=sigma2[k]))
  }
}

colnames(opdata) <- c("optimalT1","optimalT2","optimal n.site")
opdata
colnames(fidata) <- c("rho0","alpha0","rho1","alpha1","T1","N1","T2","N2","Percentage","N.sites")
fidata

setwd("/Users/deckard/Desktop/Fan Li Project/Project 3/simulation result")
write.csv(fidata,"CODA_NEW.csv")
write.csv(opdata,"STPD optimal therapists_NEW.csv")

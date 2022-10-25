#function calculating the ratio
ratio<-function(t1,t2,n1,n2,rho0,rho1,alpha0,alpha1,sigma1){
  kappa=sqrt(alpha1*rho1)
  sigma2=sigma1*rho1/alpha1 
  vb1<-sigma1/n1/t1*((1+(n1-1)*rho0)+n1*(t1-1)*rho1)+sigma2/n2/t2*((1+(n2-1)*alpha0)+n2*(t2-1)*alpha1)-2*sqrt(sigma1*sigma2)*sqrt(rho1*alpha1)
  vb2<-sigma1/n1/t1*((1+(n1-1)*rho0))+sigma2/n2/t2*((1+(n2-1)*alpha0))
  
  if(alpha1 >0){
    ratio <- vb1/vb2
  } else {ratio <- 1}

  return(ratio)
}

library(latex2exp)

t1 <- 4
n1 <- 15
t2 <- 6
n2 <- 10
sigma1 <- 1

x_seq<-seq(0.01,1,length=100)
y_seq<-seq(0.01,1,length=100)
par(mar=c(5, 5, 5, 5),mfrow=c(2,2))
mat_power<-outer(x_seq,y_seq,Vectorize(function(x,y)ratio(t1=4,t2=6,n1=15,n2=10,rho0=0.01,rho1=(x*0.01),
                                                          alpha0=0.01,alpha1=(y*0.01),sigma1=1)))
#panel A
contour(x_seq,y_seq,mat_power,labcex = 1,cex.lab = 1.5,cex.axis = 1.5,col="blue",lwd=2)
mtext(TeX("$\\alpha_1/\\rho_1$"),side=1,line=3,cex=1.5)
mtext(expression(frac(alpha[2],rho[2])),side=2,line=3,cex=1.5,las=1)
mtext(TeX("$(A):\\rho_1=0.01,\\rho_2=0.01$"),side=3,line=1,cex=1.5)

#panel B
mat_power<-outer(x_seq,y_seq,Vectorize(function(x,y)ratio(t1,t2,n1,n2,rho0=0.2,rho1=(x*0.2),alpha0=0.2,alpha1=(y*0.2),sigma1)))
contour(x_seq,y_seq,mat_power,labcex = 1,cex.lab = 1.5,cex.axis = 1.5,col="blue",lwd=2)
mtext(TeX("$\\alpha_1/\\rho_1$"),side=1,line=3,cex=1.5)
mtext(expression(frac(alpha[2],rho[2])),side=2,line=3,cex=1.5,las=1)
mtext(TeX("$(B):\\rho_1=0.2,\\rho_2=0.2$"),side=3,line=1,cex=1.5)

#panel C
mat_power<-outer(x_seq,y_seq,Vectorize(function(x,y)ratio(t1,t2,n1,n2,rho0=0.01,rho1=x*0.01,alpha0=0.2,alpha1=y*0.2,sigma1)))
contour(x_seq,y_seq,mat_power,labcex = 1,cex.lab = 1.5,cex.axis = 1.5,col="blue",lwd=2)
mtext(TeX("$\\alpha_1/\\rho_1$"),side=1,line=3,cex=1.5)
mtext(expression(frac(alpha[2],rho[2])),side=2,line=3,cex=1.5,las=1)
mtext(TeX("$(C):\\rho_1=0.01,\\rho_2=0.2$"),side=3,line=1,cex=1.5)

#panel D
mat_power<-outer(x_seq,y_seq,Vectorize(function(x,y)ratio(t1,t2,n1,n2,rho0=0.2,rho1=x*0.2,alpha0=0.01,alpha1=y*0.01,sigma1)))
contour(x_seq,y_seq,mat_power,labcex = 1,cex.lab = 1.5,cex.axis = 1.5,col="blue",lwd=2)
mtext(TeX("$\\alpha_1/\\rho_1$"),side=1,line=3,cex=1.5)
mtext(expression(frac(alpha[2],rho[2])),side=2,line=3,cex=1.5,las=1)
mtext(TeX("$(D):\\rho_1=0.2,\\rho_2=0.01$"),side=3,line=1,cex=1.5)
## Ex 4
## data

setwd(paste0(getwd(), "/SDS383D-Schwob/data"))
data <- read.csv("mathtest.csv",header=TRUE)

library(dplyr)
library(invgamma)
library(mrs)

y_school <- data %>%
  group_by(school) %>%
  summarize(bar_y = mean(mathscore), n_data = n())
## plot
jpeg("bary_n.jpeg",width=600,height=600)
plot(y_school$n_data,y_school$bar_y, xlab="bar_y",ylab="n")
dev.off()

P <- nrow(y_school)
Ni <- y_school$n_data
ybar <- y_school$bar_y
M <- 10000
# initialize
theta <- kappa <- matrix(NA,ncol=P,nrow=M)
mu <- sigma2 <- tau2 <- rep(NA,M)
theta[1,] <- rep(1,P)
mu[1] <- sigma2[1] <- tau2[1] <- 1
kappa[1,] <- 1/(tau2[1]*Ni + 1)

pb <- mrs.pb("Testing Yunshan's Code", M)

for (m in 2:M){

    pb$tick()
  tilde_sig <- 1/(Ni/sigma2[m-1] + 1/(tau2[m-1]*sigma2[m-1]))
  tilde_mu <- tilde_sig*(Ni*ybar/sigma2[m-1] + mu[m-1]/(tau2[m-1]*sigma2[m-1]))
  for (i in 1:P){
    theta[m,i] <- rnorm(1,mean=tilde_mu[i], sd=sqrt(tilde_sig[i]))
  }
  
  mu[m] <- rnorm(1,mean=mean(theta[m,]),sd=sqrt(tau2[m-1]*sigma2[m-1]/P))
  
  b_sig <- 0
  for (i in 1:P){
    b_sig <- b_sig + sum((data[which(data$school==i),2] - theta[m,i])^2)/2 + (theta[m,i] - mu[m])^2/tau2[m-1]/2
  }
  a_sig <- (sum(Ni) +P)/2
  sigma2[m] <- rinvgamma(1,a_sig,b_sig)
  
  b_tau <- 1/2 + sum((theta[m,] - mu[m])^2)/sigma2[m]/2
  a_tau <- (P+1)/2
  tau2[m] <- rinvgamma(1,a_tau,b_tau)
  
  kappa[m,] <- 1/(tau2[m]*Ni + 1)
}

burnin <- 1000
iters <- burnin:M
jpeg("mu_tra.jpeg",width=600,height=600)
plot(iters,mu[iters],xlab="iterations",ylab=expression(mu),bty="l",type="l")
dev.off()
jpeg("sig_tra.jpeg",width=600,height=600)
plot(iters,sigma2[iters],xlab="iterations",ylab=expression(sigma^2),bty="l",type="l")
dev.off()
jpeg("tau_tra.jpeg",width=600,height=600)
plot(iters,tau2[iters],xlab="iterations",ylab=expression(tau^2),bty="l",type="l")
dev.off()

kappa_post <- apply(kappa[iters,],2,mean)
jpeg("kappa.jpeg",width=600,height=600)
plot(Ni,kappa_post,xlab="school sample size",ylab=expression(kappa_i))
dev.off()

var_school <- rep(NA,100)
for (i in 1:P){
  var_school[i] <- var(data[which(data$school==i),2])
}

jpeg("var_school.jpeg",width=600,height=600)
hist(var_school, breaks = 20)
dev.off()
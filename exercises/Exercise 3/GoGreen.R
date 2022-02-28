## title: Analyis of Green Buildings
## author: Michael R. Schwob


###
### Problem Set-up
###

set.seed(702)
setwd(paste0(getwd(), "/SDS383D-Schwob/data"))
library(mvtnorm)

data <- read.csv("greenbuildings.csv")
n <- dim(data)[1]

## Response  =  revenue/sq.ft.
y <- data$Rent/data$leasing_rate*100
y[which(data$leasing_rate == 0)] <- 0 # should we set to 0 or remove entirely ??

## design matrix
X <- cbind(rep(1:n), data$green_rating, data$City_Market_Rent, data$age, data$class_a,  data$class_b)
p <- dim(X)[2] # number of parameters

## Initialize Variables
Lambda <- diag(n)
kappa <- c(.1, 1)
K <- diag(p)
for(i in 1:n){
    if(i <= n/2){
        K[i, i] <- kappa[1]
    }else{
        K[i, i] <- kappa[2]
    }
}

## Gibbs Initialization
M <- 10000
beta.save <- matrix(0, p, M) # of parameters by # of Iterations
d <- eta <- 1
nu.star <- n + d
Lambda.star <- t(X)%*%Lambda%*%X + K
eta.star <- eta + t(y)%*%Lambda%*%y - t(y)%*%Lambda%*%X%*%solve(t(X)%*%Lambda%*%X+K)%*%t(t(y)%*%Lambda%*%X)
mu.star <- solve(Lambda.star)%*%(t(X)%*%Lambda%*%y)

###
### Gibbs Sampler (to get beta)
###

for(i in 1:M){

    ###
    ### Sample beta
    ###

    beta.save[, i] <- rmvt(1, sigma = Lambda.star*nu.star/eta.star, df = nu.star, delta = mu.star)

}
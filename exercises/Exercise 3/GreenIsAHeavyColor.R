## title: A Heavy-tailed Error Model for Green Buildings
## author: Michael R. Schwob

###
### Set-up
###

set.seed(702)
library(ggplot2)
library(mvtnorm)
M <- 10000 # number of iterations

###
### Get Data
###

setwd(paste0(getwd(), "/SDS383D-Schwob/data"))
data <- read.csv("greenbuildings.csv")
n <- dim(data)[1]

## Response  =  revenue*ft^2/100
y <- data$Rent*data$leasing_rate/100

## design matrix
X <- cbind(rep(1:n), data$green_rating, data$City_Market_Rent, data$age, data$class_a,  data$class_b)
p <- dim(X)[2] # number of parameters

###
### Initializations
###

## Priors
K <- rep(0.001, p)
K <- diag(K)
m <- rep(1, p)
d <- 1
eta <- 1
h <- 2

## Initial Values
lambda.vec <- rgamma(n, h/2, h/2)
Lambda <- diag(lambda.vec)
omega <- rgamma(d/2, eta/2)
beta <- rmvnorm(1, m, solve(omega*K))

## Save Matrices
lambda.save <- matrix(0, n, M)
lambda.save[, 1] <- lambda.vec
omega.save <- rep(0, M)
omega.save[1] <- omega
beta.save <- matrix(0, p, M)
beta.save[, 1] <- beta


## Computations
nu.star <- n + d
Lambda.star <- t(X)%*%Lambda%*%X + K
mu.star <- solve(Lambda.star)%*%(t(X)%*%Lambda%*%y + t(K)%*%m)
eta.star <- eta + t(y)%*%Lambda%*%y + t(m)%*%K%*%m - (t(y)%*%Lambda%*%X + t(m)%*%K)%*%solve(t(X)%*%Lambda%*%X + K)%*%t(t(y)%*%Lambda%*%X + t(m)%*%K)
Sigma.star <- drop(eta.star/nu.star)*solve(Lambda.star)

###
### Gibbs Sampler
###

for(i in 2:M){

    ###
    ### Sample beta
    ###

    tmp.Sig <- solve(omega[i-1]*t(X)%*%Lambda%*%X + omega[i-1]%*%K) ## remember to update Lambda after lambda_i's
    tmp.mn <- tmp.Sig%*%(omega[i-1]*t(X)%*%Lambda%*%y + omega[i-1]*K%*%m)
    beta.save[, i] <- rmvnorm(1, tmp.mn, tmp.Sig)

    ###
    ### Sample omega
    ###

    

    ###
    ### Sample lambda_i
    ###

}
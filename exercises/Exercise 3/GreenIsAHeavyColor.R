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
omega <- rgamma(1, d/2, eta/2)
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
eta.star <- eta + t(y)%*%Lambda%*%y + t(m)%*%K%*%m - (t(y)%*%Lambda%*%X + t(m)%*%K)%*%solve(t(X)%*%Lambda%*%X + K)%*%t(t(y)%*%Lambda%*%X + t(m)%*%K)

## work on the whole multiply by a diagonal matrix thing

## Progress Bar
pb <- progress_bar$new(format = " downloading [:bar] :percent eta: :eta", total = M, clear = FALSE)

###
### Gibbs Sampler
###

for(i in 2:M){
    
    pb.tick()

    ###
    ### Sample beta
    ###

    tmp.Sig <- solve(omega.save[i-1]*t(X)%*%Lambda%*%X + omega.save[i-1]*K)
    tmp.mn <- tmp.Sig%*%(omega.save[i-1]*t(X)%*%Lambda%*%y + omega.save[i-1]*K%*%m)
    beta.save[, i] <- rmvnorm(1, tmp.mn, tmp.Sig)

    ###
    ### Sample omega
    ###

    omega.save[i] <- rgamma(1, nu.star/2, rate = eta.star/2)

    ###
    ### Sample lambda_i
    ###

    for(j in 1:n){
        lambda.save[j, i] <- rgamma(1, (h+1)/2, rate = (omega.save[i]*(y[i] - X[i, ]%*%beta.save[, i])^2 + h)/2) ## may need tranposes
    }

    ## New computations
    Lambda <- diag(lambda.save[, i])
    eta.star <- eta + t(y)%*%Lambda%*%y + t(m)%*%K%*%m - (t(y)%*%Lambda%*%X + t(m)%*%K)%*%solve(t(X)%*%Lambda%*%X + K)%*%t(t(y)%*%Lambda%*%X + t(m)%*%K)

}
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
y <- data$Rent*data$leasing_rate/100

## design matrix
X <- cbind(rep(1:n), data$green_rating, data$City_Market_Rent, data$age, data$class_a,  data$class_b)
p <- dim(X)[2] # number of parameters

## Priors
K <- rep(0.001, p)
K <- diag(K)
m <- rep(1, p)
d <- 1
eta <- 1
Lambda <- diag(n)

## Initialize Variables
nu.star <- n + d
Lambda.star <- t(X)%*%Lambda%*%X + K
mu.star <- solve(Lambda.star)%*%(t(X)%*%Lambda%*%y + t(K)%*%m)
eta.star <- eta + t(y)%*%Lambda%*%y + t(m)%*%K%*%m - (t(y)%*%Lambda + t(m)%*%K)%*%X%*%solve(t(X)%*%Lambda%*%X+K)%*%t(t(y)%*%Lambda%*%X + t(m)%*%K)

Sigma.star <- nu.star/eta.star*Lambda.star




###
### Beta update
###

beta.save[, i] <- rmvt(1, sigma = Lambda.star*nu.star/eta.star, df = nu.star, delta = mu.star)
## title: Analyis of Green Buildings
## author: Michael R. Schwob

###
### Set-up
###

set.seed(702)
setwd("/home/mikel/Desktop/Code/SDS383D-Schwob/data")
library(mvtnorm) # sample from multivariate t distribution
library(bayestestR) # calculate Bayesian credible intervals

###
### Structure Data
###

data <- read.csv("greenbuildings.csv")
n <- dim(data)[1]

## Response  =  revenue*ft^2/100
y <- data$Rent*data$leasing_rate/100

## design matrix
X <- cbind(rep(1, n), data$green_rating, data$City_Market_Rent, data$age, data$class_a,  data$class_b) # fixed this
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
Lambda <- diag(n)

## Computations Using (c)
nu.star <- n + d
Lambda.star <- t(X)%*%Lambda%*%X + K
mu.star <- solve(Lambda.star)%*%(t(X)%*%Lambda%*%y + t(K)%*%m)
eta.star <- eta + t(y)%*%Lambda%*%y + t(m)%*%K%*%m - (t(y)%*%Lambda%*%X + t(m)%*%K)%*%solve(t(X)%*%Lambda%*%X + K)%*%t(t(y)%*%Lambda%*%X + t(m)%*%K)
Sigma.star <- drop(eta.star/nu.star)*solve(Lambda.star)

###
### Beta updates
###

betas <- rmvt(n = 1000, sigma = Sigma.star, df = nu.star, delta = mu.star)

###
### Obtain 95% Intervals
###

## Using our method
ci(betas[, 2], method = "HDI")

## Using lm() method
fit <- lm(y ~ 0 + X)
confint(fit)[2, ]

###
### Residual Analysis
###

res <- y - X%*%apply(betas, 2, mean)
png("hist.png")
hist(res, main = "Histogram of Model Residuals", breaks = 50, col = "black", border = "white")
dev.off()
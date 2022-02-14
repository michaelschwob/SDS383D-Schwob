###
### Fitting GLMs
###

set.seed(702)
setwd(paste0(getwd(), "/SDS383D-Schwob/data"))

## Read and Structure Data
data <- read.csv("wdbc.csv",header=FALSE)

y <- as.numeric(data[,2]=="M")
X <- as.matrix(data[,3:12])
X <- apply(X,2,scale)
X <- cbind(rep(1,nrow(data)),X)

y_train <- y[1:400]
X_train <- X[1:400,]

y_test <- y[401:nrow(data)]
X_test <- X[401:nrow(data),]

p <- 11

## Write log-likelihood function
loglike <- function(y,theta){
  ll <-  sum(y*theta - log(1+exp(theta)))
  return(ll)
}

###
### Line Search Implementation
###

## Function to find optimal step (via line search)
f_gamma <- function(gamma){
  theta <- X%*%(beta+gamma*g)
  ll <-  sum(y*theta - log(1+exp(theta)))
  return(-ll)
}

## Initializations
beta_iter <- beta <- matrix(rep(0.1,11),ncol=1)
ll_iter <- c()
tol <- 1e-6
conv <- FALSE
niter <- 0

## Gradient Descent Loop
while (!conv){
  g <- 0
  for (i in 1:nrow(X_train)){
    g <- g + (y_train[i] - exp(sum(X_train[i,]*beta))/(1+exp(sum(X_train[i,]*beta))))*X_train[i,]
  }
  g <- matrix(g,ncol=1)
  gamma <- optimize(f_gamma,c(1e-6,0.3))$minimum
  beta <- beta + gamma*g
  beta_iter <- cbind(beta_iter,beta)
  theta <- X_train%*%beta
  ll <- loglike(y_train,theta)
  ll_iter <- c(ll_iter,ll)
  if (abs(sum(g^2)) < tol ){
    conv <- TRUE
  }
  niter <- niter + 1
}

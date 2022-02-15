###
### Fitting GLMs
###

set.seed(702)
setwd(paste0(getwd(), "/SDS383D-Schwob/data"))

## Read and Structure Data
data <- read.csv("wdbc.csv",header=FALSE)

y <- as.numeric(data[,2] == "M")
X <- as.matrix(data[, 3:12])
X <- apply(X, 2, scale)
X <- cbind(rep(1,nrow(data)), X)

p <- 11

## Write log-likelihood function
loglike <- function(y, theta){
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
tol <- 1e-5
conv <- FALSE
ll.old <- -100000

## Gradient Descent Loop
while (!conv){
  g <- 0
  for (i in 1:nrow(X)){
    g <- g + (y[i] - exp(sum(X[i,]*beta))/(1+exp(sum(X[i,]*beta))))*X[i,]
  }
  g <- matrix(g,ncol=1)
  gamma <- optimize(f_gamma, c(1e-6, 0.3))$minimum
  beta <- beta + gamma*g
  beta_iter <- cbind(beta_iter,beta)
  theta <- X%*%beta
  ll <- loglike(y,theta)
  ll_iter <- c(ll_iter,ll)
  if (abs(ll.old - ll) < tol){
    conv <- TRUE
  }

  ll.old <- ll

}

###
### Compare with glm()
###

comparison <- glm(y ~ X, family = binomial())


###
### Plot Log-likelihood values
###

png("gradientdescent.png")
plot(ll_iter, type = "l", main = "Convergence of Log-likelihood", ylab = "Log-likelihood")
abline(h = logLik(comparison), col = "red")
legend("bottomright", legend = c("Gradient Descent Inference", "glm()"), lty = 1, col = c("black", "red"))
dev.off()

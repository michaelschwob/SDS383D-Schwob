## title: Newton Method
## author: Michael R. Schwob

###
### Set Libraries and Data
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

## Initializations
beta.iter <- beta <- matrix(rep(0.1,11),ncol=1)
theta <- X%*%beta
ll.iter <- c()
tol <- 1e-5
conv <- FALSE
ll.old <- -100000

###
### Newton-Raphson Method
###

while (!conv){

  ## Obtain gradient
  g <- 0
  for (i in 1:nrow(X)){
    g <- g + (y[i] - exp(sum(X[i,]*beta))/(1+exp(sum(X[i,]*beta))))*X[i,]
  }
  g <- matrix(g, ncol=1)

  ## Obtain Hessian
  W.vec <- exp(theta)/(1+exp(theta))^2
  W <- diag(as.vector(W.vec)) # calculate W
  H <- -t(X)%*%W%*%X # simplify computation since W is diagonal matrix
  
  ## Compute new beta
  beta <- beta - solve(H)%*%g
  beta.iter <- cbind(beta.iter, beta)

  ## Update theta
  theta <- X%*%beta

  ## Compute log-likelihood
  ll <- loglike(y, theta)
  ll.iter <- c(ll.iter, ll)

  ## Test for convergence
  if (abs(ll.old - ll)/abs(ll.old) < tol){
    conv <- TRUE
  }

  ll.old <- ll

}

###
### Compare with glm()
###

comparison <- glm(y ~ 0 + X, family = binomial()) # fit model without intercept, since X includes the intercept

beta - comparison$coefficients

###
### Plot Log-likelihood values
###

png("newtonGraphic.png")
plot(ll.iter, type = "l", main = "Convergence of Log-likelihood", ylab = "Log-likelihood")
abline(h = logLik(comparison), col = "red")
legend("bottomright", legend = c("Newton's Method Inference", "glm()"), lty = 1, col = c("black", "red"))
dev.off()

###
### Obtain Standard Errors via Hessian Matrix
###

inv.Hess <- -solve(H)
sd.errors <- rep(0, p)
for(i in 1:p){
    sd.errors[i] <- sqrt(inv.Hess[i, i])
}

sd.errors - summary(comparison)$coefficients[, 2]
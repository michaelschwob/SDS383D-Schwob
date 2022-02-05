## title: Solving for Beta
## author: Michael R. Schwob

library(matlib) # for inv() function
library(matrixcalc) # for LU decomposition
library(microbenchmark) # for comparing the two methods

###
### Inversion Method
###

## Solves for beta by inverting X'WX and post-multiplying Wy
invertFun <- function(X, y, W){
  beta <- inv(t(X)%*%W%*%X)%*%t(X)%*%W%*%y
  return(beta)
}

###
### Using LU Decomposition
###

solveFun <- function(X, y, W){
  
  ## Pseudo-code step 1 ; factor X'WX as LU
  decomp <- lu.decomposition(t(X)%*%W%*%X)
  L <- decomp$L
  U <- decomp$U
  
  ## Pseudo-code step 2 ; solve for z in Lz = X'Wy
  z <- forwardsolve(L, t(X)%*%W%*%y)
  
  ## Pseudo-code step 3 ; solve for beta.hat in U*beta.hat = z
  beta <- backsolve(U, z)

  return(beta)
}

###
### Simulate Data and Implement
###

N <- c(10, 100, 500, 800)
P <- c(2, 50, 100, 200)

for(i in 1:length(N)){

  ## Initialize variables
  n <- N[i]
  p <- P[i]
  W <- diag(n) # W = I_n
  X <- matrix(rnorm(n*p), n, p) # design matrix filled with random values
  y <- rnorm(n, 0.3*X[,1]+0.5*X[,2], 1) # instilled some kind of dependence on X
  
  ## Implementation
  assign(paste0("benchmark", i), microbenchmark(invertFun(X, y, W), solveFun(X, y, W), times = 10)) # save benchmark information
}

###
### Print Results
###

for(i in 1:length(N)){
  print(paste0("Benchmark when N=", N[i]," and P=", P[i]))
  print(get(paste0("benchmark", i)))
}
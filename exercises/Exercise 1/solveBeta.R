## title: Solving for Beta
## author: Michael R. Schwob

library(matlib) # for inv() function
library(matrixcalc) # for LU decomposition
library(microbenchmark) # for comparing the two methods

###
### Inversion Method
###

invertFun <- function(X, y, W){
  beta <- inv(t(X)%*%W%*%X)%*%t(X)%*%W%*%y
  return(beta)
}

###
### Using LU Decomposition
###

solveFun <- function(X, y, W){
  
  ## Pseudo-code step 1
  decomp <- lu.decomposition(t(X)%*%W%*%X)
  L <- decomp$L
  U <- decomp$U
  
  ## Pseudo-code step 2
  y <- forwardsolve(L, t(X)%*%W%*%y)
  
  ## Pseudo-code step 3
  beta <- backsolve(U, y)
  return(beta)
}

###
### Simulate Data and Implement
###

N <- c(800)
P <- c(200)

for(i in 1:length(N)){

  ## Initialize variables
  n <- N[i]
  p <- P[i]
  W <- diag(n) # identity matrix for W for
  X <- matrix(rnorm(n*p), n, p) # design matrix
  y <- rnorm(n, 0.3*X[,1]+0.5*X[,2], 1)
  
  ## Implementation
  assign(paste0("benchmark",i),microbenchmark(invertFun(X, y, W), solveFun(X, y, W), times=10)) # save benchmark information
}

###
### Print Results
###

for(i in 1:length(N)){
  print(paste0("Benchmark when N=",N[i]," and P=",P[i]))
  print(get(paste0("benchmark",i)))
}

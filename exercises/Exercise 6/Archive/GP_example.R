library(lhs)
library(plgp)
eps <- sqrt(.Machine$double.eps)

X <- randomLHS(40, 2)
X[, 1] <- (X[, 1] - 0.5)*6 + 1
X[, 2] <- (X[, 2] - 0.5)*6 + 1
y <- X[, 1]*exp(-X[, 1]^2 - X[, 2]^2)

xx <- seq(-2, 4, length = 40) # vector of length 40
XX <- expand.grid(xx, xx) # 1600 x 2

D <- distance(X) # 40 x 40
Sigma <- exp(-D) # 40 x 40

DXX <- distance(XX) # 1600 x 1600
SXX <- exp(-DXX) + diag(eps, ncol(DXX)) # 1600 x 1600

DX <- distance(XX, X) # 1600 x 40
SX <- exp(-DX) # 1600 x 40

Si <- solve(Sigma) # 40 x 40
mup <- SX %*% Si %*% y # 1600 x 40 x 40 x 40 x 40 x n
Sigmap <- SXX - SX %*% Si %*% t(SX) # 1600 x 1600 - 1600 x 40 x 40 x 40 x 40 x 1600

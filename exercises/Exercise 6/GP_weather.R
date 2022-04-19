###
### Gaussian Processes - An Application to Weather
###

###
### Set-up
###

## Load the Goodies
library(mrs)
library(plgp) # for distance function
mrs.seed()
mrs.load()

## Get Data
setwd(paste0(getwd(), "/SDS383D-Schwob/data"))
data <- read.csv("weather.csv")
y1 <- data$pressure
y2 <- data$temperature
X <- data[, c(3, 4)]
n <- length(y1)

y <- y2 # change to y2 to get temperature visuals

## Set Directory
setwd("/home/mikel/Desktop/Code/SDS383D-Schwob/exercises/Exercise 6")

## Hyperparameters
t2 <- 10^(-6) # suggested by James
s2 <- 1

###
### Functions
###

## Function of Squared exponential covariance function
Exp2Sigma <- function(X, b, tau1sq, tau2sq){
  eucDist = as.matrix(dist(X,diag=T,upper=T))
  kron.delta = diag(nrow(X))
  tau1sq*exp(-.5*(eucDist/b)^2) + tau2sq*kron.delta
}

## Log-marginal likelihood function
margin.likelihood <- function(b, t1){
    C <- Exp2Sigma(X, b, t1, t2)
    Sig <- C + s2*diag(n)
    result <- -1/2*t(y)%*%solve(Sig)%*%y - 1/2*log(det(C + s2*diag(n))) - n/2*log(2*pi)
    return(result)
}

###
### Obtain Optimal Hyperparameters
###

## Set-up
M <- 50 # number of points in each grid
t1.grid <- seq(0.01, 200, length.out = M)
b.grid <- seq(0.01, 5, length.out = M)
grid <- expand.grid(b.grid, t1.grid)
y <- as.matrix(y)

## Iterate through grid of points
log.vals <- rep(0, dim(grid)[1])
pb <- mrs.pb("Evaluating log-marginal likelihood: ", dim(grid)[1])
for(i in 1:dim(grid)[1]){
    pb$tick()
    log.vals[i] <- margin.likelihood(grid[i, 1], grid[i, 2]) # b, tau
}

## Save optimal hyperparameters
b.hat <- grid[which(log.vals == max(log.vals, na.rm = TRUE)), 1]
t1.hat <- grid[which(log.vals == max(log.vals, na.rm = TRUE)), 2]

###
### Obtain Posterior Means and Standard Deviations
###

## Obtain regular grid to obtain posterior means and standard deviations
M <- 100
x.lat <- seq(min(X[, 2]), max(X[, 2]), length.out = M)
x.long <- seq(min(X[, 1]), max(X[, 1]), length.out = M)
X.new <- expand.grid(x.long, x.lat) 

t1 <- t1.hat
b <- b.hat

## Function to calculate C(x,x^*)
C_tilde <- function(X, Xstar, b, tau1sq, tau2sq){
  Dist1 = (X[,1] - Xstar[,1])
  Dist2 = (X[,2] - Xstar[,2])
  tau1sq*exp(-.5*((Dist1/as.numeric(b))^2 + (Dist2/as.numeric(b))^2)) #+ tau2sq*colSums(X == Xstar)
}

C <- Exp2Sigma(X, b, t1, t2)
Cinv <- solve(C + s2*diag(n))
fhat <- var <- rep(0, M)
for (i in 1:nrow(X.new)){
  Ct <- C_tilde(X, X.new[i,], b, t1, t2)
  Cstar <- C_tilde(X.new[i,], X.new[i,], b, t1, t2)
  fhat[i] <- Ct%*%Cinv%*%as.matrix(y,ncol=1)
  var[i] <- Cstar - Ct%*%Cinv%*%matrix(Ct,ncol=1)
}

###
### Plotting Results
###

## 2-d Heat Plot
pdf("y2_new.pdf")
par(mfrow = c(1, 2))
cols <- heat.colors(100)
image(x.lat, x.long, matrix(fhat, ncol = M), xlab = "Latitude", ylab = "Longitude", col = cols)
points(X[, 2], X[, 1])
image(x.lat, x.long, matrix(sqrt(var), ncol = M), xlab = "Latitude", ylab = "Longitude", col = cols)
points(X[, 2], X[, 1])
dev.off()

## Perspective Plot
pdf("persp2_new.pdf")
persp(x.lat, x.long, matrix(fhat, ncol = M), zlim = c(min(fhat)-1, max(fhat)+1), xlab = "Latitude", ylab = "Longitude", zlab = "Posterior Preditive Mean", theta = -30, phi = 30)
dev.off()

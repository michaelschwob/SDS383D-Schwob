###
### Gaussian Processes - An Application to Weather
###

###
### Set-up
###

## Load the Goodies
library(mrs)
library(plgp) # for distance function
library(image) # for plotting
mrs.seed()
mrs.load()

## Get Data
setwd(paste0(getwd(), "/SDS383D-Schwob/data"))
data <- read.csv("weather.csv")
y1 <- data$pressure
y2 <- data$temperature
X <- data[, c(4, 3)]
n <- length(y1)

y <- y1 # change to y2 to get temperature visuals

## Set Directory
setwd("/home/mikel/Desktop/Code/SDS383D-Schwob/exercises/Exercise 6")

## Hyperparameters
t2 <- 10^(-6) # suggested by James
s2 <- 1

###
### Covariance Functions
###

## Create covariance matrix with Matern function (squared exponential)
matern.func <- function(X, b, tau1sq, tau2sq){
    eucDist <- distance(X)
	kron.delta = diag(nrow=dim(X)[1])
	tau1sq*exp(-.5*(eucDist/b)^2) + tau2sq
}

## Create covariance matrix with Matern function (squared exponential) given distance matrix
matern.dist <- function(X, b, tau1sq, tau2sq) {
	kron.delta = diag(nrow=dim(X)[1])
	tau1sq*exp(-.5*(X/b)^2) + tau2sq
}

###
### Obtain Optimal Hyperparameters
###

## Set-up
M <- 250 # number of points in each grid
t1.grid <- seq(0.001, 100, length.out = M)
b.grid <- seq(0.0001, 100, length.out = M)
grid <- expand.grid(b.grid, t1.grid)
y <- as.matrix(y)

## Log-marginal likelihood function
margin.likelihood <- function(b, t1){
    C <- matern.func(X, b, t1, t2) # get C
    Sig <- C + s2*diag(n)
    result <- -1/2*t(y)%*%solve(Sig)%*%y - 1/2*log(det(Sig)) - n/2*log(2*pi)
    return(result)
}

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
M <- 50
x.lat <- seq(min(X[, 1]), max(X[, 1]), length.out = M)
x.long <- seq(min(X[, 2]), max(X[, 2]), length.out = M)
X.new <- expand.grid(x.lat, x.long) 

## Get optimal Sigma under given data
Sigma <- matern.func(X, b.hat, t1.hat, t2)
S.inv <- solve(Sigma + s2*diag(ncol(Sigma)))

## Get optimal Sigma for new data
DX.new <- distance(X.new) # distance matrix for all new points ; 2500 x 2500
SX.new <- matern.dist(DX.new, b.hat, t1.hat, t2) # C for new points ; 2500 x 2500

## Get optimal Sigma across new data and given data
DX <- distance(X.new, X) # 2500 x 157
SX <- matern.dist(DX, b.hat, t1.hat, t2) # 2500 x 157

## Get posterior means and variances ; from GP (B)
mu.post <- SX%*%S.inv%*%y
var.post <- SX.new - SX%*%S.inv%*%t(SX)
sd.post <- sqrt(diag(var.post))
sd.post <- ifelse(is.na(sd.post), 8, sd.post) # to smooth across NA values

###
### Plotting Results
###

## 2-d Heat Plot
pdf("y1.pdf")
par(mfrow = c(1, 2))
cols <- heat.colors(128)
image(x.lat, x.long, matrix(mu.post, ncol = M), xlab = "Latitude", ylab = "Longitude", col = cols)
points(X[, 1], X[, 2])
image(x.lat, x.long, matrix(sd.post, ncol = M), xlab = "Latitude", ylab = "Longitude", col = cols)
points(X[, 1], X[, 2])
dev.off()

## Perspective Plot
pdf("persp1.pdf")
persp(x.lat, x.long, matrix(mu.post, ncol = M), xlab = "Latitude", ylab = "Longitude", zlab = "Posterior Preditive Mean", theta = -30, phi = 30)
dev.off()
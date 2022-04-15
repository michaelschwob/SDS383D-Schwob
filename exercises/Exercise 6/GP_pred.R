###
### Gaussian Processes - An Application to Utilities
###

###
### Set-up
###

## Load the Goodies
library(mrs)
library(ggpubr)
#library(GauPro)
mrs.seed()
mrs.load()

## Get Data
setwd(paste0(getwd(), "/SDS383D-Schwob/data"))
data <- read.csv("utilities.csv")
y <- data$gasbill/data$billingdays
x <- data$temp

n <- length(x)

## Set Directory
setwd("/home/mikel/Desktop/Code/SDS383D-Schwob/exercises/Exercise 6")

## Hyperparameters
t2 <- 10^(-10)
#t2 <- 0 # suggested by James
b <- seq(10, 100, length.out = 10)
t1 <- seq(0.1, 10, length.out = 10)

b <- 5
t1 <- 0.001

###
### Covariance Functions
###

## Create covariance matrix with Matern function (squared exponential)
matern.func <- function(x, b, tau1sq, tau2sq) {
	eucDist = as.matrix(dist(x,diag=T,upper=T))
	kron.delta = diag(nrow=length(x))
	tau1sq*exp(-.5*(eucDist/b)^2) + tau2sq*kron.delta
}

## Create covariance vector with Matern function (squared exponential) # !! new
matern.vec <- function(x, x.star, b, tau1sq, tau2sq) {
    eucDist <- abs(x - x.star)
    kron.delta <- 0 # this is off the main diagonal
	tau1sq*exp(-.5*(eucDist/b)^2) + tau2sq*kron.delta
}

## Matern function with parameter 5/2 (matrix)
matern2.func <- function(x, b, tau1sq, tau2sq) {
	eucDist = as.matrix(dist(x,diag=T,upper=T))
	kron.delta = diag(nrow=length(x))

    tau1sq*(1 + sqrt(5)*eucDist/b + 5*eucDist^2/(3*b^2))*exp(-sqrt(5)*eucDist/b) + tau2sq*kron.delta
}

## Matern function with parameter 5/2 (vector)
matern2.vec <- function(x, x.star, b, tau1sq, tau2sq) {
	eucDist <- abs(x - x.star)
	kron.delta = 0

    tau1sq*(1 + sqrt(5)*eucDist/b + 5*eucDist^2/(3*b^2))*exp(-sqrt(5)*eucDist/b) + tau2sq*kron.delta
}

###
### Implementation
###

counter <- 0
pb <- mrs.pb("Gauss me", length(t1)*length(b))

#gp <- GauPro(x, y)

for(k in 1:length(t1)){
    for(j in 1:length(b)){

        pb$tick()

        post.mean <- lb <- ub <- rep(0, n) # initialize posterior vectors

        for(i in 1:n){ # for each data point ; c = current, o = other

            Sig.cc <- matern.func(x[i], b[j], t1[k], t2) # should just be 1 + t2
            Sig.cc <- 1
            Sig.oo <- matern.func(x[-i], b[j], t1[k], t2)
            Sig.oc <- as.matrix(matern.vec(x[-i], x[i], b[j], t1[k], t2))

            # we don't have means, right ??
            mu.c <- 0
            mu.o <- rep(0, n-1)

            post.mean[i] <- mu.c + t(Sig.oc)%*%solve(Sig.oo)%*%(y[-i] - mu.o)

            Sig.bar <- Sig.cc - t(Sig.oc)%*%solve(Sig.oo)%*%Sig.oc

            #diff <- 1.96*sqrt(Sig.bar)*(n-1) # is this right??
            diff <- 1.96*sqrt(Sig.bar) # is this right??

            lb[i] <- post.mean[i] - diff
            ub[i] <- post.mean[i] + diff

        }

        plot.df <- data.frame(mean = post.mean, y = y, x = x, lb = lb, ub = ub)

        plot <- ggplot(plot.df, aes(x = x, y = y)) + geom_point(color = "red", alpha = 0.8) + geom_ribbon(mapping = aes(ymin = lb, ymax = ub), alpha = 0.5) + geom_line(mapping = aes(x = x, y = post.mean), color = "#BF5700") + theme_classic() + ggtitle("Gaussian Process for Utilities")
        assign(paste0("p", counter), plot)

        counter <- counter + 1

    }
}




###
### Gaussian Processes - An Application to Utilities
###

###
### Set-up
###

## Load the Goodies
library(mrs)
library(ggpubr)
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
t2 <- 0 # suggested by James
b <- c(10, 50, 100)
t1 <- c(1, 10, 20)

s2 <- 2

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

###
### Implementation
###

counter <- 1
pb <- mrs.pb("Gauss me", length(t1)*length(b))

for(k in 1:length(t1)){
    for(j in 1:length(b)){

        pb$tick()

        post.mean <- lb <- ub <- rep(0, n) # initialize posterior vectors

        for(i in 1:n){ # for each observation

            other.data <- x[-i]

            rho.ii <- matern.func(x[i], b[j], t1[k], t2)    # constant
            rho.col <- as.matrix(matern.vec(other.data, x[i], b[j], t1[k], t2))  # column: (n-1)x1
            rho.row <- t(rho.col)  # row: 1x(n-1) 
            rho.Sig <- matern.func(other.data, b[j], t1[k], t2) # matrix: (n-1)x(n-1)
            rho.Sig <- rho.Sig + diag(s2, nrow = n-1) # add jitter

            post.mean[i] <- rho.row%*%solve(rho.Sig)%*%y[-i]

            post.cov <- rho.ii - rho.row%*%solve(rho.Sig)%*%rho.col

            lb[i] <- post.mean[i] - 1.96*sqrt(post.cov)
            ub[i] <- post.mean[i] + 1.96*sqrt(post.cov)

        }

        plot.df <- data.frame(mean = post.mean, y = y, x = x, lb = lb, ub = ub)

        plot <- ggplot(plot.df, aes(x = x, y = y)) + geom_point(color = "red", alpha = 0.8) + geom_ribbon(mapping = aes(ymin = lb, ymax = ub), alpha = 0.5) + geom_line(mapping = aes(x = x, y = post.mean), color = "#BF5700") + theme_classic() + ggtitle(paste0("Gaussian Process: b = ", b[j], ", t1 = ", t1[k]))
        assign(paste0("p", counter), plot)

        counter <- counter + 1

    }
}

plot <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3, ncol = 3)
ggsave("post_conf.png", plot)

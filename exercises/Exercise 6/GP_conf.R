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
t2 <- 10^(-6) # suggested by James
b <- c(3, 10, 15)
t1 <- c(1, 5, 10)

s2 <- 0.61 # the MAP

###
### Covariance Functions
###

## Create covariance matrix with Matern function (squared exponential)
matern.func <- function(x, b, tau1sq, tau2sq) {
	eucDist = as.matrix(dist(x,diag=T,upper=T))
	kron.delta = diag(nrow=length(x))
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
        Sig <- matern.func(x, b[j], t1[k], t2)

        post.mean <- solve(diag(n) + s2*solve(Sig))%*%y
        diff.Sig <- solve(diag(n)/s2 + solve(Sig)) 
        lb <- post.mean - 1.96*sqrt(diag(diff.Sig))
        ub <- post.mean + 1.96*sqrt(diag(diff.Sig))

        plot.df <- data.frame(mean = post.mean, y = y, x = x, lb = lb, ub = ub)

        plot <- ggplot(plot.df, aes(x = x, y = y)) + geom_point(color = "red", alpha = 0.8) + geom_ribbon(mapping = aes(ymin = lb, ymax = ub), alpha = 0.5) + geom_line(mapping = aes(x = x, y = post.mean), color = "#BF5700") + theme_classic() + ggtitle(paste0("Gaussian Process: b = ", b[j], ", t1 = ", t1[k]))
        assign(paste0("p", counter), plot)

        counter <- counter + 1

    }
}

plot <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3, ncol = 3)
#ggsave("post_conf.png", plot)
plot

###
### (E) : Marginal Likelihood to Fing Optimal Hyperparameters
###

###
### Obtain Optimal Hyperparameters
###

M <- 500 # number of points in each grid
t1.grid <- seq(0.001, 100, length.out = M)
b.grid <- seq(0.0001, 100, length.out = M)
grid <- expand.grid(b.grid, t1.grid)
y <- as.matrix(y)

margin.likelihood <- function(b, t1){
    C <- matern.func(x, b, t1, t2) # get C
    Sig <- C + s2*diag(n)
    result <- -1/2*t(y)%*%solve(Sig)%*%y - 1/2*log(det(Sig)) - n/2*log(2*pi)
    return(result)
}

log.vals <- rep(0, dim(grid)[1])
pb <- mrs.pb("Evaluating log-marginal likelihood: ", dim(grid)[1])

for(i in 1:dim(grid)[1]){
    pb$tick()
    log.vals[i] <- margin.likelihood(grid[i, 1], grid[i, 2]) # b, tau
}

b.hat <- grid[which(log.vals == max(log.vals, na.rm = TRUE)), 1]
t1.hat <- grid[which(log.vals == max(log.vals, na.rm = TRUE)), 2]

b.hat <- 61.52308
t1.hat <- 39.67996

###
### Compute Posterior Mean for f|y and Plot
###

Sig <- matern.func(x, b.hat, t1.hat, t2)

post.mean <- solve(diag(n) + s2*solve(Sig))%*%y
diff.Sig <- solve(diag(n)/s2 + solve(Sig))

plot.df <- data.frame(mean = post.mean, y = y, x = x)

plot <- ggplot(plot.df, aes(x = x, y = y)) + geom_point(color = "red", alpha = 0.8) + geom_line(mapping = aes(x = x, y = post.mean), color = "#BF5700") + theme_classic() + ggtitle(paste0("Gaussian Process: b = ", b.hat, ", t1 = ", t1.hat))
plot
ggsave("optimal_hyperparameters.png", plot)
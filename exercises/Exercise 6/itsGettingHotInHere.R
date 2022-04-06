###
### Fitting Local Linear Estimators with a Gaussian Kernel
###

###
### Structure Data
###

setwd(paste0(getwd(), "/SDS383D-Schwob/data"))
data <- read.csv("utilities.csv")

Y <- data$gasbill/data$billingdays
X <- data$temp
H.mat <- X%*%solve(t(X)%*%X)%*%t(X)
n <- length(X)

###
### Set-up
###

library(mrs)
mrs.load()
mrs.seed()
setwd("/home/mikel/Desktop/Code/SDS383D-Schwob/exercises/Exercise 6")

###
### Functions
###

## Weight Function
weight.func <- function(x.new, x.old, h){ # indexed by i is the actual "old" data
    int0 <- (x.new - x.old)/h
    s1 <- sj.func(x.new, h, 1)
    s2 <- sj.func(x.new, h, 2)
    weight <- gauss.kernel(int0)*(s2 - (x.old - x.new)*s1)
    return(weight)
}

## s function
sj.func <- function(x.new, h, j){
    sum <- 0
    for(i in 1:n){
        int <- (x.new - X[i])/h
        sum <- sum + gauss.kernel(int)*(X[i] - x.new)^j
    }
    return(sum)
}

## Kernel Function
gauss.kernel <- function(x){
    sol <- exp(-(x)^2/2)/sqrt(2*pi)
    return(sol)
}

##
## (E) Fit Local Linear Estimators
##

## Problem Set-up
H <- c(0.1, 0.2, 0.5, 1) # list of h values
H <- seq(0.1, 1, length.out = 20)
LOOCV <- rep(0, length(H)) # initialize average squared error in prediction (LOOCV) matrix

M <- 1000 # number of "new" points to fit kernel-estimator
x.grid <- seq(min(X), max(X), length.out = M) # "new points"; must be within range of X data
Y.matrix <- matrix(0, M, length(H)) # save estimated y for plotting

## Initialize Progress Bar
pb <- mrs.pb("bills, bills, bills...", length(H))

## For each value of h
for(k in 1:length(H)){

    pb$tick() # update progress bar

    h <- H[k] # pick h
    smooth.y <- rep(0, M) # initialize vector for smoothed y

    ## Fit kernel-regression estimator to training data
    for(j in 1:M){ # for each value in x.grid
        num <- den <- 0 # initialize sums
        for(i in 1:n){ # for each observation
            num <- num + weight.func(x.grid[j], X[i], h)*Y[i]
            den <- den + weight.func(x.grid[j], X[i], h)
        }
        smooth.y[j] <- num/den
    }

    ## Approximate Function
    approx.func <- approxfun(x.grid, smooth.y, rule = 2) # rule = 2 is needed to extrapolate (get rid of NAs)

    ## Obtain LOOCV from test data
    LOOCV.tmp <- 0
    for(i in 1:n){ # for each observation
        num.tmp <- Y[i] - approx.func(X[i]) 
        den.tmp <- 1 - H.mat[i, i]
        LOOCV.tmp <- LOOCV.tmp + (num.tmp/den.tmp)^2
    }
    LOOCV[k] <- LOOCV.tmp/n

    Y.matrix[, k] <- smooth.y
}

## Prepare Data for Plotting
name.vec <- 0
for(i in 1:length(H)){
    tmp <- paste0("h = ", H[i])
    name.vec <- c(name.vec, tmp)
}
name.vec <- name.vec[-1]
colnames(Y.matrix) <- name.vec
plot.df <- data.frame(cbind(x.grid, Y.matrix))
plot.df <- plot.df %>% gather(key = "h", value = "value", -1)
point.df <- data.frame(test.y = Y, test.x = X)

## Plot
plot <- ggplot(plot.df, aes(x = x.grid, y = value)) + geom_line(aes(color = h)) + theme_classic() + ggtitle("Local Linear Estimation") + xlab("x") + ylab("Estimated Value") + geom_point(point.df, mapping = aes(x = test.x, y = test.y), alpha = 0.4)
ggsave("imFairlyLocal.png", plot)

## Get LOOCV
LOOCV

## Fit data using new function
optim.h <- which(LOOCV == min(LOOCV)) # get optimal h
approx.func <- approxfun(x.grid, Y.matrix[, optim.h], rule = 2)
y.fit <- rep(0, n)
for(i in 1:n){
    y.fit[i] <- approx.func(X[i])
}

## Plot Fit vs Actual
plot.df <- data.frame(fit = y.fit, x = X, actual = Y)
plot2 <- ggplot(plot.df, aes(x = X)) + geom_line(aes(y = fit)) + geom_point(aes(y = actual), color = "red") + ggtitle("Fit Y vs. Actual Y") + xlab("X") + ylab("Y") + theme_classic()
ggsave("fitvsactual.png", plot2)

###
### (F) Inspect Residuals
###

resid <- Y- y.fit
resid.df <- data.frame(x = X, resid = resid)
plot.resid <- ggplot(resid.df, aes(x = x, y = resid)) + geom_point(color = "red") + ggtitle("Residual Plot") + xlab("X") + ylab("Y") + theme_classic() + geom_hline(yintercept = 0, linetype = "dashed")
ggsave("residuals.png", plot.resid)

###
### (G) Construct an Approximate Point-wise 95% CI
###

## Obtain confidence intervals
ci.low <- ci.high <- rep(0, n) # initialize confidence intervals
ci.high <- rep(10, n)
s2.hat <- sum((Y-y.fit)^2)/(n - 2*sum(diag(H.mat)) + sum(diag(t(H.mat)%*%H.mat)))
fact <- 1.96*sqrt(s2.hat)
for(i in 1:n){
    ci.low[i] <- y.fit[i] - fact
    ci.high[i] <- y.fit[i] + fact
}

## Confidence Interval Data Frame
ci.df <- data.frame(low = ci.low, high = ci.high, x = X)

## Plot
plot.df <- data.frame(fit = y.fit, x = X, actual = Y)
plot2 <- ggplot(plot.df, aes(x = X)) + geom_line(aes(y = fit)) + geom_point(aes(y = actual), color = "red") + ggtitle("Fit Y vs. Actual Y") + xlab("X") + ylab("Y") + theme_classic() + geom_segment(ci.df, mapping = aes(x = x, xend = x, y = low, yend = high), color = "#BF5700", size = 1.5)
ggsave("ci_plot.png", plot2)

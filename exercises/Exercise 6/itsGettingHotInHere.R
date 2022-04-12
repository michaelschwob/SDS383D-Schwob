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
H <- seq(0.1, 10, length.out = 20)
LOOCV <- rep(0, length(H)) # initialize average squared error in prediction (LOOCV) matrix

x.grid <- X
M <- length(x.grid)
Y.matrix <- matrix(0, M, length(H)) # save estimated y for plotting

## Initialize Progress Bar
pb <- mrs.pb("bills, bills, bills...", length(H))

H.mat <- matrix(0, n, M) # !!
weightz <- rep(0, n)

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
            weightz[i] <- weight.func(x.grid[j], X[i], h)
        }
        weightz <- weightz/sum(weightz)

        for(i in 1:n){
            H.mat[i, j] <- weightz[i]
        }

        smooth.y[j] <- num/den
    }

    ## Obtain LOOCV from test data
    LOOCV.tmp <- 0
    for(i in 1:n){ # for each observation
        num.tmp <- Y[i] - smooth.y[i]
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
optim.h <- which(LOOCV == min(LOOCV, na.rm = TRUE)) # get optimal h
y.fit <- Y.matrix[, optim.h]

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

## Get weight matrix
Weight.mat <- matrix(0, n, n)
for(i in 1:n){ # for each observation
    for(j in 1:n){ # compare to other observations
        Weight.mat[i, j] <- weight.func(X[j], X[i], optim.h) # x.new, x.old
    }
    Weight.mat[i, ] <- Weight.mat[i, ]/sum(Weight.mat[i, ])
}

## Normalize Weight Matrix
s2.hat <- t(Y - Weight.mat%*%Y)%*%(Y - Weight.mat%*%Y)/(n - 2*sum(diag(Weight.mat)) + sum(diag(t(Weight.mat)%*%Weight.mat)))

## Obtain confidence intervals
ci.low <- ci.high <- rep(0, n) # initialize confidence intervals
ci.high <- rep(10, n)

for(i in 1:n){
    fact <- 1.96*sqrt(s2.hat*norm(Weight.mat[i, ], type = "2"))
    ci.low[i] <- y.fit[i] - fact
    ci.high[i] <- y.fit[i] + fact
}

## Confidence Interval Data Frame
ci.df <- data.frame(low = ci.low, high = ci.high, x = X)

## Plot
plot.df <- data.frame(fit = y.fit, x = X, actual = Y)
plot2 <- ggplot(plot.df, aes(x = X)) + geom_line(aes(y = fit)) + geom_point(aes(y = actual), color = "red") + ggtitle("Fit Y vs. Actual Y") + xlab("X") + ylab("Y") + theme_classic() + geom_segment(ci.df, mapping = aes(x = x, xend = x, y = low, yend = high), color = "#BF5700", size = 1)
ggsave("ci_plot.png", plot2)

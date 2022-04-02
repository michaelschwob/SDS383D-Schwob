###
### A Study of Cross Validation with Smoothers
###

library(mrs)
mrs.load()
mrs.seed()

setwd(paste0(getwd(), "/SDS383D-Schwob/exercises/Exercise 6"))

##
## (A)
##

source("smooth_like_Jagger.R")

X.1 <- rnorm(700, 0, 40) # training data
X.2 <- rnorm(300, 0, 20) # testing data; must be within range of X.1
X <- c(X.1, X.2) # combine training and testing data
Y <- f.func(X)
data <- cbind(scale(Y), scale(X))
train.data <- data[1:700, ]
test.data <- data[701:1000, ]

##############################################################

H <- c(1, 2, 5, 1)
M <- 1000
x.grid <- seq(-10, 10, length.out = M)
MSPE <- rep(0, length(H))

pb <- mrs.pb("I need more validation. ", length(H))

for(k in 1:length(H)){

    pb$tick()

    h <- H[k]
    smooth.y <- test.y <- rep(0, M)

    ## Fit kernel-regression estimator to training data
    for(j in 1:length(x.grid)){
        tmp.sum <- 0
        for(i in 1:dim(train.data)[1]){
            tmp.sum <- tmp.sum + weight.ind(train.data[i, 2], x.grid[j], h)*train.data[i, 1]
        }
        smooth.y[j] <- tmp.sum
    }

    ## Approximate Function
    approx.func <- approxfun(x.grid, smooth.y)

    ## Obtain MSPE from test data
    MSPE.tmp <- 0
    for(i in 1:dim(test.data)[1]){
        MSPE.tmp <- MSPE.tmp + (test.data[i, 1] - approx.func(test.data[i, 2]))^2
    }
    MSPE[k] <- MSPE.tmp/dim(test.data)[1]
    
}
plot(x.grid,smooth.y)
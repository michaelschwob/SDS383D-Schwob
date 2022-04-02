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
#X.1 <- scale(X.1)
#X.2 <- rnorm(300, 0, 20) # testing data; must be within range of X.1
X.2 <- sample(X.1, 300, replace = TRUE) # hear me out !!
X <- c(X.1, X.2) # combine training and testing data
noise <- rnorm(1000, 0, 3)
Y <- f.func(X) + noise
data <- cbind(Y, X)
train.data <- data[1:700, ]
test.data <- data[701:1000, ]

##############################################################

H <- c(1, 2, 5, 10)
M <- 1000
x.grid <- seq(-100, 100, length.out = M)
MSPE <- rep(0, length(H))

Y.matrix <- matrix(0, M, length(H) + 1)
Y.matrix[, 1] <- f.func(x.grid)

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
    approx.func <- approxfun(x.grid, smooth.y, rule = 2) # rule = 2 is needed to extrapolate

    ## Obtain MSPE from test data
    MSPE.tmp <- 0
    for(i in 1:dim(test.data)[1]){
        MSPE.tmp <- MSPE.tmp + (test.data[i, 1] - approx.func(test.data[i, 2]))^2
    }
    MSPE[k] <- MSPE.tmp/dim(test.data)[1]

    Y.matrix[, k + 1] <- smooth.y
}

colnames(Y.matrix) <- c("true", "h = 1", "h = 2", "h = 5", "h = 10")
plot.df <- data.frame(cbind(x.grid, Y.matrix[, -1]))
plot.df <- plot.df %>% gather(key = "h", value = "value", -1)
point.df <- data.frame(test.y = Y[701:1000], test.x = X.2)

plot <- ggplot(plot.df, aes(x = x.grid, y = value)) + geom_line(aes(color = h)) + theme_classic() + ggtitle("Cross Validation with Different Bandwidths h") + xlab("x") + ylab("Estimated Value") + geom_point(point.df, mapping = aes(x = test.x, y = test.y), alpha = 0.6)
ggsave("cross_validation.png", plot)

###
### A Study of Cross Validation Under Different Scenarios
###

library(mrs)
mrs.load()
mrs.seed()

setwd(paste0(getwd(), "/SDS383D-Schwob/exercises/Exercise 6"))

##
## (B)
##

## Problem Set-up
H <- c(0.1, 0.2, 0.5, 1) # list of h values
ASEP <- matrix(0, length(H), 4) # initialize average squared error in prediction (ASEP) matrix

X.1 <- runif(375, 0, 1) # testing data
X.2 <- runif(125, 0, 1) # training data
low.noise <- rnorm(500, 0, 0.25)
high.noise <- rnorm(500, 0, 0.8)

M <- 1000 # number of "new" points to fit kernel-estimator
x.grid <- seq(0, 1, length.out = M) # "new points"
Y.matrix <- matrix(0, M, length(H) + 1) # save estimated y for each scenario's plot

## "Wiggly" Function
wibbly.func <- function(x){
    result <- cos(10*x)
    return(result)
}

## "Smooth" Function
slide.func <- function(x){
    result <- x^(5/3)
    return(result)
}

## Weight Function
weight.func <- function(x.old, x.new, h){
    int <- (x.old - x.new)/h
    weight <- gauss.kernel(int)/h
    return(weight)
}

## Kernel Function
gauss.kernel <- function(x){
    sol <- exp(-(x)^2/2)/sqrt(2*pi)
    return(sol)
}

## Run Cross Validation in Four Scenarios
for(m in 1:4){

    ## Function Selection
    if(m == 1){
        Y <- wibbly.func(X.1) + high.noise[1:375] # train Y
        Y.test <- wibbly.func(X.2) + high.noise[376:500] # test y
        Y.matrix[, 1] <- wibbly.func(x.grid) # add "true" values for y
    }
    if(m == 2){
        Y <- wibbly.func(X.1) + low.noise[1:375] # train Y
        Y.test <- wibbly.func(X.2) + low.noise[376:500] # test y
        Y.matrix[, 1] <- wibbly.func(x.grid) # add "true" values for y
    }
    if(m == 3){
        Y <- slide.func(X.1) + high.noise[1:375] # train Y
        Y.test <- slide.func(X.2) + high.noise[376:500] # test y
        Y.matrix[, 1] <- slide.func(x.grid) # add "true" values for y
    }
    if(m == 4){
        Y <- slide.func(X.1) + low.noise[1:375] # train Y
        Y.test <- slide.func(X.2) + low.noise[376:500] # test y
        Y.matrix[, 1] <- slide.func(x.grid) # add "true" values for y
    }

    ## Initialize Progress Bar
    pb <- mrs.pb(paste0("I need more validation (", m, ")"), length(H))

    ## For each value of h
    for(k in 1:length(H)){

        pb$tick() # update progress bar

        h <- H[k] # pick h
        smooth.y <- rep(0, M) # initialize vector for smoothed y

        ## Fit kernel-regression estimator to training data
        for(j in 1:length(x.grid)){
            
            weights <- rep(0, length(X.1)) # vector of weights

            ## Remember to normalize the weights
            for(i in 1:length(X.1)){
                weights[i] <- weight.func(X.1[i], x.grid[j], h)
            }
            weights <- weights/sum(weights)

            ## Obtain sum for estimated value
            tmp.sum <- 0
            for(i in 1:length(X.1)){
                tmp.sum <- tmp.sum + weights[i]*Y[i]
            }

            smooth.y[j] <- tmp.sum
        }

        ## Approximate Function
        approx.func <- approxfun(x.grid, smooth.y, rule = 2) # rule = 2 is needed to extrapolate (get rid of NAs)

        ## Obtain ASEP from test data
        ASEP.tmp <- 0
        for(i in 1:length(X.2)){
            ASEP.tmp <- ASEP.tmp + (Y.test[i] - approx.func(X.2[i]))^2
        }
        ASEP[k, m] <- ASEP.tmp/length(X.2)

        Y.matrix[, k + 1] <- smooth.y
    }

    ## Prepare Data for Plotting
    colnames(Y.matrix) <- c("true", "h = 0.1", "h = 0.2", "h = 0.5", "h = 1")
    plot.df <- data.frame(cbind(x.grid, Y.matrix[, -1]))
    plot.df <- plot.df %>% gather(key = "h", value = "value", -1)
    point.df <- data.frame(test.y = Y.test, test.x = X.2)

    ## Plot
    plot <- ggplot(plot.df, aes(x = x.grid, y = value)) + geom_line(aes(color = h)) + theme_classic() + ggtitle("Cross Validation") + xlab("x") + ylab("Estimated Value") + geom_point(point.df, mapping = aes(x = test.x, y = test.y), alpha = 0.6)
    assign(paste0("p", m), plot)

}

## Save Plot
plots <- ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
ggsave("scenarios.png", plots)

## Get ASEP
ASEP
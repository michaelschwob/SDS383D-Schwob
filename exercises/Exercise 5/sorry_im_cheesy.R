###
### Learning Hierarchical Linear Models
### with a Case Study in Price Elasticity of Demand (PED)
###

###
### Set up
###

library(mrs)
mrs.load()
mrs.seed()

data <- read.csv(paste0(getwd(), "/SDS383D-Schwob/data/cheese.csv"))
setwd(paste0(getwd(), "/SDS383D-Schwob/exercises/Exercise 5"))

###
### Structure Data
###

names <- unique(data$store) # store names
n <- length(names) # number of stores
biggest.store <- names(which.max(table(data$store)))
max <- sum(data$store == biggest.store) # largest number of observations per store
Y <- matrix(NA, n, max) # ij matrix
X <- array(0, dim = c(max, 4, n)) # n stores at end so we get subsetable matrix containing ; matches form of X_i in write-up
N.i <- rep(0, n) # initialize vector to determine Ni

for(i in 1:n){ # for each store
    tmp.data <- data[which(data$store==names[i]),] # data set for current store
    Y[i, ] <- log(tmp.data$vol) # stores log(quantity)
    N.i <- dim(tmp.data)[1] # number of observations for store i
    for(j in 1:length(N.i)){
        X[j, , i] <- c(1, log(tmp.data$vol[j]), tmp.data$disp[j], log(tmp.data$vol[j])*tmp.data$disp[j])
    }
}

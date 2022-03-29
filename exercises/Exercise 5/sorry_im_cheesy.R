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
X <- array(NA, dim = c(max, 4, n)) # n stores at end so we get subsetable matrix containing ; matches form of X_i in write-up
N.i <- rep(0, n) # initialize vector to determine Ni

for(i in 1:n){ # for each store
    tmp.data <- data[which(data$store==names[i]),] # data set for current store
    N.i[i] <- dim(tmp.data)[1] # number of observations for store i
    Y[i, 1:(N.i[i])] <- log(tmp.data$vol) # stores log(quantity)
    for(j in 1:N.i[i]){
        X[j, , i] <- c(1, log(tmp.data$vol[j]), tmp.data$disp[j], log(tmp.data$vol[j])*tmp.data$disp[j])
    }
}

###
### Set Starting Values and Priors
###

a.0 <- rgamma(1, 3, rate = 1)
b.0 <- rgamma(1, 3, rate = 1)
sp2.0 <- rinvgamma(4, a.0, b.0)
Sigma.0 <- diag(sp2.0)
mu.beta.0 <- rnorm(4, 0, 2)
beta.i.0 <- rmvnorm(1, mu.beta.0, Sigma.0)
s2i.0 <- rinvgamma(n, 0.5, 0.5)

###
### Save Matrices
###

n.mcmc <- 10000

a.save <- b.save <- rep(0, n.mcmc)
sig2.save <- matrix(0, n, n.mcmc)
s2.save <- mu.beta.save <- matrix(0, 4, n.mcmc)
beta.save <- array(0, dim = c(n, 4, n.mcmc))

a.save[1] <- a.0
b.save[1] <- b.0
sig2.save[, 1] <- s2i.0
s2.save[, 1] <- sp2.0
mu.beta.save[, 1] <- mu.beta.0
for(i in 1:n){
    beta.save[i, , 1] <- beta.i.0
}  

###
### Gibbs Sampler
###
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
tmp.Sig <- diag(s2.save[, 1])
mu.beta.save[, 1] <- mu.beta.0
for(i in 1:n){
    beta.save[i, , 1] <- beta.i.0
}  

pb <- mrs.pb("Making that Chesse: ", n.mcmc)

###
### Gibbs Sampler
###

for(k in 2:n.mcmc){

    pb$tick()

    ###
    ### Update beta
    ###
    
    for(i in 1:n){
        Sig.star.inv <- t(X[1:(N.i[i]), , i])%*%X[1:(N.i[i]), , i]/sig2.save[i, k-1] + solve(tmp.Sig)
        Sig.star <- solve(Sig.star.inv)

        mu.star <- Sig.star%*%( t(X[1:(N.i[i]), , i])%*%Y[i, 1:(N.i[i])]/sig2.save[i, k-1] + solve(tmp.Sig)%*%mu.beta.save[, k-1]  ) # removed Y tranpose !!

        beta.save[i, , k] <- rmvnorm(1, mu.star, Sig.star)
    }
    

    ###
    ### Update mu.beta
    ###

    mu.beta.save[, k] <- rmvnorm(1, rep(mean(beta.save[, , k]), 4), tmp.Sig/n) # does this seem right??

    ###
    ### Update s2.p
    ###

    for(p in 1:4){
        tmp.a <- n/2 + 1/2
        tmp.sum <- 0
        for(i in 1:n){
            tmp.sum <- tmp.sum + (beta.save[i, p, k] - mu.beta.save[p, k])^2 # tranposed beta.save
        }
        tmp.b <- 1/2*(1 + tmp.sum)
        s2.save[p, k] <- rinvgamma(1, tmp.a, tmp.b)
    }
    tmp.Sig <- diag(s2.save[, k])

    ###
    ### Update sigma2.i
    ###

    for(i in 1:n){
        tmp.a <- a.save[k-1]/2 + 1/2
        tmp.mat <- matrix(Y[i, 1:N.i[i]]) - X[1:N.i[i], , i]%*%beta.save[i, , k]
        tmp.b <- 1/2*(b.save[k-1] + t(tmp.mat)%*%tmp.mat) # tranposed beta and Y

        sig2.save[i, k] <- rinvgamma(1, tmp.a, tmp.b)
    }

    ###
    ### Update a
    ###

    a.star <- runif(1, 1, 10)

    mh1 <- sum(dinvgamma(sig2.save[, k], a.star/2, b.save[k-1]/2, log = TRUE)) + dgamma(a.star, 3, 1, log = TRUE) + dunif(a.save[k-1], 0, 10, log = TRUE)
    mh2 <- sum(dinvgamma(sig2.save[, k], a.save[k-1]/2, b.save[k-1]/2, log = TRUE)) + dgamma(a.save[k-1], 3, 1, log = TRUE) + dunif(a.star, 0, 10, log = TRUE)
    mh <- exp(mh1 - mh2)
    if(mh > runif(1)){
        a.save[k] <- a.star # use proposed value
    }else{
        a.save[k] <- a.save[k-1] # retain value
    }

    ### 
    ### Update b
    ###

    b.star <- runif(1, 1, 10)

    mh1 <- sum(dinvgamma(sig2.save[, k], a.save[k]/2, b.star/2, log = TRUE)) + dgamma(b.star, 3, 1, log = TRUE) + dunif(b.save[k-1], 0, 10, log = TRUE)
    mh2 <- sum(dinvgamma(sig2.save[, k], a.save[k]/2, b.save[k-1]/2, log = TRUE)) + dgamma(b.save[k-1], 3, 1, log = TRUE) + dunif(b.star, 0, 10, log = TRUE)
    mh <- exp(mh1 - mh2)
    if(mh > runif(1)){
        b.save[k] <- b.star # use proposed value
    }else{
        b.save[k] <- b.save[k-1] # retain value
    }

}

###
### Inference
###

sig2.df <- data.frame()
###
### Learning Hierarchical Linear Models
### with a Case Study in Price Elasticity of Demand (PED)
###

## Note: To run this, you need the mrs R package,  which can be installed with the following steps:
##       1. install and load the "devtools" R package
##       2. run "install_github("michaelschwob/mrs")"
##       3. load mrs library as you usually would

###
### Set up
###

library(mrs)
library(ggpubr)
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
        X[j, , i] <- c(1, log(tmp.data$price[j]), tmp.data$disp[j], log(tmp.data$price[j])*tmp.data$disp[j])
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
s2i.0 <- rep(1, n) # added

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

pb <- mrs.pb("Making that Cheese: ", n.mcmc)

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

    #mu.beta.save[, k] <- rmvnorm(1, rep(mean(beta.save[, , k]), 4), tmp.Sig/n) # does this seem right??
    
    mu.vec <- rep(0, 4)
    for(i in 1:4){
        mu.vec[i] <- mean(beta.save[, i, k])
    }

    mu.beta.save[, k] <- rmvnorm(1, mu.vec, tmp.Sig/n)

    ###
    ### Update s2.p
    ###

    for(p in 1:4){
        tmp.a <- n/2 + 1/2
        tmp.sum <- 0
        for(i in 1:n){
            tmp.sum <- tmp.sum + (beta.save[i, p, k] - mu.beta.save[p, k])^2
        }
        tmp.b <- 1/2*(1 + tmp.sum)
        s2.save[p, k] <- rinvgamma(1, tmp.a, tmp.b)
    }
    tmp.Sig <- diag(s2.save[, k])

    ###
    ### Update sigma2.i
    ###

    for(i in 1:n){
        tmp.a <- a.save[k-1]/2 + N.i[i]/2 # swich to N.i[i] ; was + 1/2
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

    a.save[k] = 0.5

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

    b.save[k] = 0.5

    #  # update a
    # a_grid <- seq(1,10,length.out=500)
    # pa <- sapply(a_grid, FUN = function(x) {sum((-x/2-1)*log(sig2.save[, k]))} ) 
    # pa <- pa - mean(pa) # to avoid things in exp() going to Inf
    # pa <- exp(pa)
    # a.save[k] <- sample(a_grid, 1, prob = pa)
  
    # # update b
    # b_grid <- seq(0.1,10,length.out=500)
    # pb <- sapply(b_grid,FUN= function(x) {exp(-sum(x/sig2.save[, k]/2))} ) 
    # b.save[k] <- sample(b_grid, 1, prob = pb)

}

###
### Trace Plots for a, b, mu.beta,  and s2_p
###

a.trace <- mrs.trace(a.save, "a")
b.trace <- mrs.trace(b.save, "b")

plots <- ggarrange(a.trace, b.trace, ncol = 1)
ggsave("traces.png", plots)

mu.trace.1 <- mrs.trace(mu.beta.save[1, ], "mu[beta][1]") 
mu.trace.2 <- mrs.trace(mu.beta.save[2, ], "mu[beta][2]") 
mu.trace.3 <- mrs.trace(mu.beta.save[3, ], "mu[beta][3]") 
mu.trace.4 <- mrs.trace(mu.beta.save[4, ], "mu[beta][4]")

sp.trace.1 <- mrs.trace(s2.save[1, ], "s[1]^2")
sp.trace.2 <- mrs.trace(s2.save[2, ], "s[2]^2")
sp.trace.3 <- mrs.trace(s2.save[3, ], "s[3]^2")
sp.trace.4 <- mrs.trace(s2.save[4, ], "s[4]^2")

plots2 <- ggarrange(mu.trace.1, mu.trace.2, mu.trace.3, mu.trace.4, sp.trace.1, sp.trace.2, sp.trace.3, sp.trace.4, nrow = 4, ncol = 2)
ggsave("traces2.png", plots2)

###
### Inference on sig2
###

sig2.post <- rep(0, n)

for(i in 1:n){
    tmp.vec <- sig2.save[i, ]
    outliers <- boxplot(tmp.vec, plot=FALSE)$out
    tmp.vec <- tmp.vec[-which(tmp.vec %in% outliers)]
    sig2.post[i] <- mean(tmp.vec)
}

sig2.df <- data.frame(sig2 = sig2.post, n = 1:n)
p2 <- ggplot(sig2.df, aes(x = n, y = sig2)) + geom_point(color = "red") + theme_classic() + ggtitle("Log-normal Variance Per Store") + xlab("Store") + ylab(expression(sigma[i]^2))
ggsave("sig2_perstore.png", p2)

###
### Histogram of Beta values
###

for(p in 1:4){
    tmp.df <- data.frame(beta = c(beta.save[, p, (.3*n.mcmc):n.mcmc]))
    tmp.plot <- ggplot(tmp.df, aes(x = beta)) + geom_histogram(color = "white", fill = "#BF5700") + theme_classic() + ggtitle("Histogram of Posterior Values") +xlab(paste0("beta", p-1))
    assign(paste0("plot", p), tmp.plot)
}
betas <- ggarrange(plot1, plot2, plot3, plot4, nrow = 2, ncol = 2)
ggsave("beta_histograms.png", betas)

###
### Effect on Demand Curve Due to Advertising
###

## beta 2

beta2.post <- rep(0, n)

for(i in 1:n){
    tmp.vec <- beta.save[i, 3, ]
    outliers <- boxplot(tmp.vec, plot=FALSE)$out
    tmp.vec <- tmp.vec[-which(tmp.vec %in% outliers)]
    beta2.post[i] <- mean(tmp.vec)
}

beta2.df <- data.frame(beta2 = beta2.post, n = 1:n)
p3 <- ggplot(beta2.df, aes(x = n, y = beta2)) + geom_point(color = "red") + theme_classic() + ggtitle("Beta[2] Per Store") + xlab("Store") + ylab(expression(beta[2]))
ggsave("beta2_perstore.png", p3)


## beta 3

###
### Convergence Diagnostics
###

library(coda)
raftery.diag(a.save)

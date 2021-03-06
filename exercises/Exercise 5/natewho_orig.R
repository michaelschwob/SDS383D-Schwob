###
### A Hierarchical Probit Model
### demonstrated on data from the 1988 US presidential election
###

###
### Set up
###

library(mrs)
mrs.load()
library(LaplacesDemon)
library(truncnorm)
mrs.seed()

data <- read.csv("/home/mikel/Desktop/Code/SDS383D-Schwob/data/polls.csv")
setwd(paste0(getwd(), "/SDS383D-Schwob/exercises/Exercise 5"))

###
### Clean Data
###

# only keep desired rows and omit observations that have NAs
data <- data %>% select(-c(org, year, survey, weight)) %>% na.omit()

# add indicator columns
data <- data %>% mutate(Bacc = ifelse(data$edu == "Bacc", 1, 0), HS = ifelse(data$edu == "HS", 1, 0), SomeColl = ifelse(data$edu == "SomeColl", 1, 0), NoHS = ifelse(data$edu == "NoHS", 1, 0)) %>% select(-edu)
data <- data %>% mutate(YA = ifelse(data$age == "18to29", 1, 0), A = ifelse(data$age == "30to44", 1, 0), MA = ifelse(data$age == "45to64", 1, 0), E = ifelse(data$age == "65plus", 1, 0)) %>% select(-age)

###
### Structure Data
###

names <- unique(data$state)
n <- length(names)
biggest.sample <- names(which.max(table(data$state)))
max <- sum(data$state == biggest.sample) # largest number of observations per store
P <- 8 # was 10

Y <- matrix(NA, n, max)
X <- array(NA, dim = c(max, P+1, n))
N.i <- rep(0, n)

for(i in 1:n){
    tmp.data <- data[which(data$state==names[i]),] # data set for current state
    N.i[i] <- dim(tmp.data)[1] # number of observations for store i
    Y[i, 1:(N.i[i])] <- tmp.data$bush # states voters
    for(j in 1:N.i[i]){
        X[j, , i] <- c(1, as.numeric(tmp.data[j, c(3, 4, 5, 6, 7, 9, 10, 11)]))
    }
}

###
### Set Starting Values
###

beta.star <- matrix(rnorm((P+1)*n, 0, 10^4), P+1, n)
B.star <- diag(P+1)
alpha <- matrix(0, P+1, n)
for(i in 1:n){
    alpha[, i] <- rmvnorm(1, beta.star[, i], B.star)
}
Z <- matrix(NA, n, max)
for(i in 1:n){
    for(j in 1:N.i[i]){
        Z[i, j] <- ifelse(Y[i, j] == 1, rtruncnorm(1, a = 0, mean = 0, sd = 1), rtruncnorm(1, b = 0, mean = 0, sd = 1))
    }
}

nu.0 <- P + 2 # check if this is okay

###
### Save Matrices
###

n.mcmc <- 10000
beta.save <- array(NA, dim = c(P+1, n, n.mcmc))
alpha.save <- array(NA, dim = c(P+1, n, n.mcmc))
Z.save <- array(NA, dim = c(n, max, n.mcmc))

beta.save[, , 1] <- beta.star
alpha.save[, , 1] <- alpha
Z.save[, , 1] <- Z

pb <- mrs.pb("Nate Silver has nothing on this: ", n.mcmc)

###
### MCMC
###

for(k in 2:n.mcmc){

    pb$tick()

    ###
    ### Update alpha
    ###

    for(i in 1:n){
        B.tilde <- solve(solve(B.star) + t(X[1:N.i[i], , i])%*%X[1:N.i[i], , i])
        beta.tilde <- B.tilde%*%(  solve(B.star)%*%beta.save[, i, k-1] +  t(X[1:N.i[i], , i])%*%Z.save[i, 1:N.i[i], k-1] )

        alpha.save[, i, k] <- rmvnorm(1, beta.tilde, B.tilde, checkSymmetry = FALSE) # added last argument
    }

    ###
    ### Update Z
    ###

    for(i in 1:n){
        for(j in 1:N.i[i]){
            if(Y[i, j] == 1){
                Z.save[i, j, k] <- rtruncnorm(1, a = 0, mean = X[j, , i]%*%alpha.save[, i, k], sd = 1)
            }
            if(Y[i, j] == 0){
                Z.save[i, j, k] <- rtruncnorm(1, b = 0, mean = X[j, , i]%*%alpha.save[, i, k], sd = 1)
            }
        }
    }

    ###
    ### Update beta.star
    ###

    for(i in 1:n){
        A.inv <- solve( solve(B.star) + solve(10^4*diag(P+1))  ) 
        tmp.mean <- A.inv%*%( solve(B.star)%*%alpha.save[, i, k] )
        beta.save[, i, k] <- rmvnorm(1, tmp.mean, A.inv)
    }

    ###
    ### Update B.star
    ###

    nu <- n + nu.0
    tmp.sum <- 0
    for(i in 1:n){
        tmp.sum <- tmp.sum + (alpha.save[, i, k] - beta.save[, i, k])%*%t(alpha.save[, i, k] - beta.save[, i, k])
    }
    S <- diag(P + 1) + tmp.sum

    B.star <- rinvwishart(nu, S)

}

###
### Trace Plots for alphas for each state
###

n.burn <- .3*n.mcmc
for(i in 1:n){
    tmp.df <- data.frame(iter = n.burn:n.mcmc, mu = alpha.save[1, i, n.burn:n.mcmc], beta1 = alpha.save[2, i, n.burn:n.mcmc], beta2 = alpha.save[3, i, n.burn:n.mcmc], beta3 = alpha.save[4, i, n.burn:n.mcmc], beta4 = alpha.save[5, i, n.burn:n.mcmc], beta5 = alpha.save[6, i, n.burn:n.mcmc], beta6 = alpha.save[7, i, n.burn:n.mcmc], beta7 = alpha.save[8, i, n.burn:n.mcmc])
    plot.df <- tmp.df %>% gather(key = "variable", value = "value", -1)
    traces <- ggplot(plot.df, aes(x = iter, y = value)) + geom_line(aes(color = variable)) + theme_classic() + ggtitle("Trace Plots") + xlab("Iteration") + ylab("Value")
    assign(paste0("t", i), traces)
}

#2, 4, 8, 9, 14
ggsave("texas.png", t8)
ggsave("newyork.png", t4)

###
### Plot Alphas Per State
###

n.burn <- .3*n.mcmc

tmp.mean <- rep(0, n)

p <- 1
ind <- rep(0, n)
for(i in 1:n){
    tmp.mean[i] <- mean(alpha.save[p, i, n.burn:n.mcmc]) # obtain posterior mean for each state
if(abs(tmp.mean[i]) < 10){
    ind[i] <- 1 # we will plot this one
    }
}
tmp.names <- names[which(ind == 1)]
tmp.mean <- tmp.mean[which(ind == 1)]
tmp.df <- data.frame(State = tmp.names, beta = tmp.mean)
pos.points <- data.frame(State = tmp.names[which(tmp.mean>0)], beta = tmp.mean[which(tmp.mean>0)])
neg.points <- data.frame(State = tmp.names[which(tmp.mean <= 0)], beta = tmp.mean[which(tmp.mean <= 0)])
tmp.plot <- ggplot(tmp.df, aes(x = State)) + geom_point(aes(y = beta)) + theme_classic() + ylab(bquote(mu)) + geom_hline(yintercept = 0, linetype = "dashed") + geom_segment(pos.points, mapping = aes(x = State, xend = State, y = 0, yend = beta), color = "red") + geom_segment(neg.points, mapping = aes(x = State, xend = State, y = 0, yend = beta), color = "blue")
assign(paste0("p", p), tmp.plot) + theme(axis.text.x = element_text(angle = 90))

for(p in 2:(P+1)){ # for each parameter
    ind <- rep(0, n)
    for(i in 1:n){
        tmp.mean[i] <- mean(alpha.save[p, i, n.burn:n.mcmc]) # obtain posterior mean for each state
        if(abs(tmp.mean[i]) < 10){
            ind[i] <- 1 # we will plot this one
        }
    }
    tmp.names <- names[which(ind == 1)]
    tmp.mean <- tmp.mean[which(ind == 1)]
    tmp.df <- data.frame(State = tmp.names, beta = tmp.mean)
    pos.points <- data.frame(State = tmp.names[which(tmp.mean>0)], beta = tmp.mean[which(tmp.mean>0)])
    neg.points <- data.frame(State = tmp.names[which(tmp.mean <= 0)], beta = tmp.mean[which(tmp.mean <= 0)])
    tmp.plot <- ggplot(tmp.df, aes(x = State)) + geom_point(aes(y = beta)) + theme_classic() + ylab(bquote(beta[.(p)])) + geom_hline(yintercept = 0, linetype = "dashed") + geom_segment(pos.points, mapping = aes(x = State, xend = State, y = 0, yend = beta), color = "red") + geom_segment(neg.points, mapping = aes(x = State, xend = State, y = 0, yend = beta), color = "blue") + theme(axis.text.x = element_text(angle = 90))
    assign(paste0("p", p), tmp.plot)
}
plots <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3, ncol = 3) + ggtitle("Posterior Means of Coefficients Per State")
ggsave("coefficient_estimates.png", plots)

###
### Posterior Means Across States
###

post.betas <- rep(0, 9)
for(i in 1:9){
    post.betas[i] <- mean(beta.save[i, , n.burn:n.mcmc])
}
post.betas
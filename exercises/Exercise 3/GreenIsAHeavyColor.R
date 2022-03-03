## title: A Heavy-tailed Error Model for Green Buildings
## author: Michael R. Schwob

###
### Set-up
###

set.seed(702)
library(ggplot2)
library(mvtnorm)
library(progress)
library(fields)
library(bayestestR)
M <- 5000 # number of iterations

###
### Get Data
###

setwd(paste0(getwd(), "/SDS383D-Schwob/data"))
data <- read.csv("greenbuildings.csv")
n <- dim(data)[1]

## Response  =  revenue*ft^2/100
y <- data$Rent*data$leasing_rate/100

## design matrix
X <- cbind(rep(1, n), data$green_rating, data$City_Market_Rent, data$age, data$class_a,  data$class_b)
p <- dim(X)[2] # number of parameters

###
### Initializations
###

## Priors
K <- rep(0.001, p)
K <- diag(K)
m <- rep(1, p)
d <- 1
eta <- 1
h <- 1 # was 2

## Initial Values
lambda.vec <- rgamma(n, h/2, h/2)
omega <- rgamma(1, d/2, eta/2)
beta <- rmvnorm(1, m, solve(omega*K))

## Save Matrices
lambda.save <- matrix(0, n, M)
lambda.save[, 1] <- lambda.vec
omega.save <- rep(0, M)
omega.save[1] <- omega
beta.save <- matrix(0, p, M)
beta.save[, 1] <- beta

## Computations
nu.star <- n + d
eta.star <- eta + t(y)%*%(lambda.vec%d*%y) + t(m)%*%K%*%m - (t(y)%*%(lambda.vec%d*%X) + t(m)%*%K)%*%solve(t(X)%*%(lambda.vec%d*%X) + K)%*%t(t(y)%*%(lambda.vec%d*%X) + t(m)%*%K)

## Progress Bar
pb <- progress_bar$new(format = " impressing James(?) [:bar] :percent eta: :eta", total = M, clear = FALSE)

omega.lambda <- matrix(0, n, M)

###
### Gibbs Sampler
###

for(i in 2:M){
    
    pb$tick() # update progress bar

    ###
    ### Sample beta
    ###

    tmp.Sig <- solve(omega.save[i-1]*t(X)%*%(lambda.vec%d*%X) + omega.save[i-1]*K)
    tmp.mn <- tmp.Sig%*%(omega.save[i-1]*t(X)%*%(lambda.vec%d*%y) + omega.save[i-1]*K%*%m)
    beta.save[, i] <- rmvnorm(1, tmp.mn, tmp.Sig)

    ###
    ### Sample omega
    ###

    omega.save[i] <- rgamma(1, nu.star/2, rate = eta.star/2)

    ###
    ### Sample lambda_i
    ###

    for(j in 1:n){
        lambda.save[j, i] <- rgamma(1, (h+1)/2, rate = (omega.save[i]*(y[i] - X[i, ]%*%beta.save[, i])^2 + h)/2)
    }

    ## New computations
    lambda.vec <- lambda.save[, i]
    eta.star <- eta + t(y)%*%(lambda.vec%d*%y) + t(m)%*%K%*%m - (t(y)%*%(lambda.vec%d*%X) + t(m)%*%K)%*%solve(t(X)%*%(lambda.vec%d*%X) + K)%*%t(t(y)%*%(lambda.vec%d*%X) + t(m)%*%K)

    omega.lambda[, i] <- omega.save[i]*lambda.save[, i]
}

###
### "Interesting Plot"
###

post.means <- apply(lambda.save, 1, mean)
df <- data.frame(cbind(y, post = 1/post.means))
ggplot(df, aes(x = y, y = post)) + theme_classic() + geom_point(alpha = 0.25) + ggtitle("Relative Variance for Each Data Point") + ylab(expression(1/lambda[i]))
ggsave("RelativeVariance.png")

###
### Compare 95% Intervals
###

## Use other file
setwd("/home/mikel/Desktop/Code/SDS383D-Schwob/exercises/Exercise 3")
source("GoGreen.R")
setwd("/home/mikel/Desktop/Code/SDS383D-Schwob/exercises/Exercise 3")

## Compare CIs
for(i in 1:p){
    cat(paste("Old Model with ", expression(beta), i, ": "))
    print(ci(betas[, i], method = "HDI"))
    cat("\n")
    cat(paste("NewModel with ", expression(beta), i, ": "))
    print(ci(beta.save[i, ], method = "HDI"))
    cat("\n", "\n")
}

###
### Residual Analysis Comparison
###

res.old <- y - X%*%apply(betas, 2, mean)
res.new <- y - X%*%apply(beta.save, 1, mean)

png("compare_histogram.png")
hist(res.old, main = "Histogram of Residuals (Comparison)", breaks = 50, col = rgb(0, 0, 1, 0.25))
hist(res.new, breaks = 50, col = rgb(1, 0, 0, 0.25), add = TRUE)
#abline(v = 0, col = "red")
legend("topright", legend = c("Homoskedastic", "Heteroskedastic"), col = c("red", "blue"), lty = 1)
dev.off()

sigmai.mean <- apply(omega.lambda, 1, mean)
sigmai.mean <- 1/sqrt(sigmai.mean)
hist(drop(res.new)/sigmai.mean, breaks = 30)

post.lambdas <- apply(lambda.save, 1, mean)

hist(post.lambdas)

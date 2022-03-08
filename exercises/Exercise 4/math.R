## title: Introduction to Hierarchical Models (pt. 1)
## author: Michael R. Schwob

###
### Set-up
###

library(mrs)
mrs.load()
mrs.seed()

setwd(paste0(getwd(), "/SDS383D-Schwob/data"))
data  <- read.csv("mathtest.csv")
setwd(paste0("/home/mikel/Desktop/Code/SDS383D-Schwob/exercises/Exercise 4"))

###
### (b) Plot number of students v. school average
###

P <- length(unique(data[, 1]))
n.stud <- y.avg <- rep(0, P)
for(i in 1:P){
    n.stud[i] <- nrow(filter(data, school == i))
    y.avg[i] <- mean(filter(data, school == i)[, 2]) 
}

df <- data.frame(cbind(n.stud, y.avg))
p1 <- ggplot(df, aes(x = n.stud, y = y.avg)) + theme_classic() + geom_point() + xlab("School Size") + ylab("School-level Average") + ggtitle("School-level Averages v. School Size")
ggsave("math_b.png", p1)

###
### (c) Gibbs Sampling
###

## Initializations
M <- 1000
s2.save <- t2.save <- mu.save <- rep(0, M)
theta.save <- kappa.save <- matrix(0, P, M)
s2.save[1] <- t2.save[1] <- 1
mu.save[1] <- 0
theta.save[, 1] <- rep(1, P)
kappa.save[, 1] <- 1/(t2.save[1]*n.stud + 1)

pb <- mrs.pb("stat stuff", M)

for(m in 2:M){

    pb$tick()

    ###
    ### Sample mu
    ###

    mu.save[m] <- rnorm(1, sum(theta.save[, m-1])/P, sqrt(t2.save[m-1]*s2.save[m-1]/P))

    ###
    ### Sample sigma^2
    ###

    tmp.shape <- (sum(n.stud) + P)/2
    tmp.sum1 <- tmp.sum2 <- 0
    for(i in 1:P){
        for(j in 1:n.stud[i]){
            yij <- data$mathscore[which(data$school == i)[j]]
            tmp.sum1 <- tmp.sum1 + (yij - theta.save[i, m-1])^2/2
        }
        tmp.sum2 <- tmp.sum2 + (theta.save[i, m-1] - mu.save[m])^2/(2*t2.save[m-1])
    }
    tmp.scale <- tmp.sum1 + tmp.sum2

    s2.save[m] <- rinvgamma(1, tmp.shape, tmp.scale)

    ###
    ### Sample tau^2
    ###

    tmp.shape <- (P + 1)/2
    tmp.sum1 <- 0
    for(i in 1:P){
        tmp.sum1 <- tmp.sum1 + (theta.save[i, m-1] - mu.save[m])^2/(2*s2.save[m])
    }
    tmp.scale <- tmp.sum1 + 1/2

    t2.save[m] <- rinvgamma(1, tmp.shape, tmp.scale)

    ###
    ### Sample theta
    ###

    tmp.var.vec <- 1/(n.stud/s2.save[m] + 1/(t2.save[m]*s2.save[m]))
    tmp.b <- rep(0, P)
    for(i in 1:P){
        tmp.b[i] <- sum(filter(data, school == i)[, 2])/s2.save[m] + mu.save[m]/(t2.save[m]*s2.save[m])
    }
    tmp.mn <- tmp.var.vec %d*% tmp.b

    theta.save[, m] <- rmvnorm(1, tmp.mn, diag(tmp.var.vec))

    ###
    ### (d) Compute kappa
    ###

    kappa.save[, m] <- 1/(t2.save[m]*n.stud + 1)
}

###
### Diagnostics
###

n.burn <- M/10

tp1 <- mrs.trace(mu.save, "mu", 1, n.burn)
tp2 <- mrs.trace(t2.save, "tau^2", 1, n.burn)
tp3 <- mrs.trace(s2.save, "sigma^2", 1, n.burn)
arr <- ggarrange(tp1, tp2, tp3, ncol = 1)
ggsave("traces.png", arr)


###
### (d) Plot shrinkage coefficient kappa
###

kappa.post <- apply(kappa.save[, n.burn:M], 1, mean)

df <- data.frame(cbind(n.stud, kappa.post))
p2 <- ggplot(df, aes(x = n.stud, y = kappa.post)) + theme_classic() + geom_point() + xlab("School Size") + ylab("Shrinkage Coefficient") + ggtitle("Shrinkage Coefficient v. School Size")
ggsave("kappa_graph.png", p2)
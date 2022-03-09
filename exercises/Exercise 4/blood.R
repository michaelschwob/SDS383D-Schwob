## title: Introduction to Hierarchical Models (pt. 2)
## author: Michael R. Schwob

###
### Set-up
###

library(mrs)
mrs.load()
mrs.seed()

setwd(paste0(getwd(), "/SDS383D-Schwob/data"))
data  <- read.csv("bloodpressure.csv")
setwd(paste0("/home/mikel/Desktop/Code/SDS383D-Schwob/exercises/Exercise 4"))

###
### (a) t-test for difference in means ; pooled
###

group1 <- data[, 2][which(data[, 4] == 1)]
group2 <- data[, 2][which(data[, 4] == 2)]
testa <- t.test(group1, group2)
testa
testa$stderr

###
### (b) t-test for difference in means ; average
###

ind.means <- data %>% group_by(subject) %>% summarise(sys = mean(systolic), treatment = mean(treatment))

group1 <- ind.means[which(ind.means[,3] == 1), 2]
group2 <- ind.means[which(ind.means[,3] == 2), 2]

testb <- t.test(group1, group2)
testb
testb$stderr

###
### (c) Fit Hierarchical Model via Gibbs Sampler
###

## Initializations
P <- 20 # number of people
M <- 1000
s2.save <- t2.save <- mu.save <- beta.save <- rep(0, M)
theta.save <- matrix(0, P, M)
s2.save[1] <- t2.save[1] <- beta.save[1] <- 1

s2.beta <- 1 # are these...
mu.beta <- 0 # ... good?

x <- c(rep(0, 10), rep(1, 10))

n.stud <- rep(0, P)
for(i in 1:P){
    n.stud[i] <- nrow(filter(data, subject == i))
    theta.save[i, 1] <- mean(filter(data, subject == i)[, 2]) 
}
mu.save[1] <- mean(theta.save[, 1])

pb <- mrs.pb("drawing blood", M)

for(m in 2:M){

    pb$tick()

    ###
    ### Sample mu ; updated
    ###

    mu.save[m] <- rnorm(1, (sum(theta.save[, m-1]) - beta.save[m-1]*sum(x))/P, sqrt(t2.save[m-1]*s2.save[m-1]/P))

    ###
    ### Sample sigma^2 ; updated
    ###

    tmp.shape <- (sum(n.stud) + P)/2 # same
    tmp.sum1 <- tmp.sum2 <- 0
    for(i in 1:P){
        for(j in 1:n.stud[i]){
            yij <- data$systolic[which(data$subject == i)[j]] # check
            #cat(yij,  " ")
            tmp.sum1 <- tmp.sum1 + (yij - theta.save[i, m-1])^2/2
        }
        tmp.sum2 <- tmp.sum2 + (theta.save[i, m-1] - mu.save[m] - beta.save[m-1]*x[i])^2/(2*t2.save[m-1])
    }
    tmp.scale <- tmp.sum1 + tmp.sum2

    s2.save[m] <- rinvgamma(1, tmp.shape, tmp.scale)

    ###
    ### Sample tau^2 ; updated
    ###

    tmp.shape <- (P + 1)/2 # same
    tmp.sum1 <- 0
    for(i in 1:P){
        tmp.sum1 <- tmp.sum1 + (theta.save[i, m-1] - mu.save[m] - beta.save[m-1]*x[i])^2/(2*s2.save[m])
    }
    tmp.scale <- tmp.sum1 + 1/2

    t2.save[m] <- rinvgamma(1, tmp.shape, tmp.scale)

    ###
    ### Sample theta ; updated
    ###

    tmp.var.vec <- 1/(n.stud/s2.save[m] + 1/(t2.save[m]*s2.save[m]))
    tmp.b <- rep(0, P)
    for(i in 1:P){
        tmp.b[i] <- sum(filter(data, subject == i)[, 2])/s2.save[m] + (mu.save[m] + beta.save[m-1]*x[i])/(t2.save[m]*s2.save[m])
    }
    tmp.mn <- tmp.var.vec %d*% tmp.b

    theta.save[, m] <- rmvnorm(1, tmp.mn, diag(tmp.var.vec))

    ###
    ### Sample beta ; added
    ###

    tmp.var <- 1/(sum(x^2)/(t2.save[m]*s2.save[m]) + 1/s2.beta)
    tmp.sum <- 0
    for(i in 1:P){
        tmp.sum <- tmp.sum + (theta.save[i, m] - mu.save[m])*x[i]/(t2.save[m]*s2.save[m])
    }
    tmp.mn <- tmp.var*(tmp.sum + mu.beta/s2.beta)

    beta.save[m] <- rnorm(1, tmp.mn, sqrt(tmp.var))

}

###
### Diagnostics
###

n.burn <- M/10

tp1 <- mrs.trace(mu.save, "mu", 1, n.burn)
tp2 <- mrs.trace(t2.save, "tau^2", 1, n.burn)
tp3 <- mrs.trace(s2.save, "sigma^2", 1, n.burn)
tp4 <- mrs.trace(beta.save, "beta", 1, n.burn)
arr <- ggarrange(tp1, tp2, tp3, tp4, ncol = 1)
ggsave("traces_blood.png", arr)


###
### Beta Inference
###

png("beta_histogram.png")
hist(beta.save, main = expression(paste("Histogram of Posterior Distribution of ", beta)), xlab = expression(beta))
dev.off()

mean(beta.save[n.burn:M])
sd(beta.save[n.burn:M])

###
### (d) Are measurements within individuals independent?
###

dfd <- data.frame(data)
dfd$date <- format(as.Date(dfd$date,format="%Y-%m-%d"), format = "%d")

for(i in 1:P){
    ind <- ifelse(i <= 10, 1, 2)
    tmp <- ggplot(dfd, aes(x = date, y = systolic, group = 1)) + geom_line(dfd %>% filter(subject==i), mapping = aes(x = date, y = systolic), col = ind) + theme_classic() + xlab("Date") + ylab("Systolic Measurement") + ggtitle(paste("Subject ", i)) + theme(plot.title = element_text(hjust = 0.5))
    assign(paste0("p", i), tmp)
}

template <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, ncol = 4, nrow = 5)
ggsave("plotsd.png", template)
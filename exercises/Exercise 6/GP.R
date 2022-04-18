###
### An Exploration of Gaussian Processes
###

setwd(paste0(getwd(), "/SDS383D-Schwob/exercises/Exercise 6"))

library(mrs)
library(ggpubr)
mrs.seed()
mrs.load()

## Create covariance matrix with Matern function (squared exponential)
matern.func <- function(x, b, tau1sq, tau2sq) {
	eucDist = as.matrix(dist(x,diag=T,upper=T))
	kron.delta = diag(nrow=length(x))
	tau1sq*exp(-.5*(eucDist/b)^2) + tau2sq*kron.delta
}

## Matern function with parameter 5/2
matern2.func <- function(x, b, tau1sq, tau2sq) {
	eucDist = as.matrix(dist(x,diag=T,upper=T))
	kron.delta = diag(nrow=length(x))

    tau1sq*(1 + sqrt(5)*eucDist/b + 5*eucDist^2/(3*b^2))*exp(-sqrt(5)*eucDist/b) + tau2sq*kron.delta
}

# X values
n <- 500
x.grid <- sort(runif(n, 0, 1))

t2 <- 0
b <- seq(0.0001, 1, length.out = 10)
t1 <- seq(0, 2, length.out = 10)

pb <- mrs.pb("Gaussian Process fitting:", length(b)*length(t1))

for(i in 1:length(b)){ # across different values of b

    Y.matrix <- matrix(0, nrow = n, ncol = 1) # initialize

    for(j in 1:length(t1)){ # across different values of tau1^2

        pb$tick()

        tmp.Sig <- matern.func(x.grid, b[i], t1[j], t2) # switch this function to use different matern covariance function
        tmp.vec <- rmvnorm(1, rep(0, n), tmp.Sig)
        Y.matrix <- cbind(Y.matrix, t(as.matrix(tmp.vec)))

    }

    Y.matrix <- Y.matrix[, -1] # remove the initial

    name.vec <- 0
    for(ki in 1:length(t1)){
        tmp <- paste0("t1 = ", t1[ki])
        name.vec <- c(name.vec, tmp)
    }
    name.vec <- name.vec[-1]
    colnames(Y.matrix) <- name.vec

    plot.df <- data.frame(x.grid, Y.matrix)
    plot.df <- plot.df %>% gather(key = "pars", value = "value", -1)

    tmp.plot <- ggplot(plot.df, aes(x = x.grid, y = value)) + geom_line(aes(color = pars)) + theme_classic() + ggtitle(paste0("Gaussian Processes (b = ", b[i], ")")) + xlab("x") + ylab("Estimated Value")

    assign(paste0("p", i), tmp.plot)
}

plot.arr <- ggarrange(p1, p2, p3, p4, p5, p6, nrow = 3, ncol = 2)
plot.arr
#ggsave("sqexp.png", plot.arr)
#ggsave("matern52.png", plot.arr)
###
### Fitting a Kernel Smoother
###

library(mrs)
mrs.load()
mrs.seed()

## Kernel Function
gauss.kernel <- function(x){
    sol <- exp(-x^2/2)/sqrt(2*pi)
    return(sol)
}

## f(x) Function
f.func <- function(x){
    return(x^2+4*x^3-12)
}

## Weighting Functions
weight.ind <- function(x.old, x.new, h){
    int <- (x.old - x.new)/h
    weight <- gauss.kernel(int)/h
    return(weight)
}

## Simulate data
N <- 1000
x <- rnorm(N, 0, 10)
noise <- rnorm(N, 0, 2)
y <- f.func(x) + noise

## Center data
x.c <- x - mean(x)
y.c <- y - mean(y)

H <- c(1, 2, 5, 10)
M <- 1000
x.grid <- seq(-10, 10, length.out = M)

Y.matrix <- matrix(0, M, length(H) + 1)
Y.matrix[, 1] <- f.func(x.grid)

pb <- mrs.pb("Smooth like Jagger ", length(H))

for(k in 1:length(H)){

    pb$tick()

    h <- H[k]
    smooth.y <- rep(0, M)

    for(j in 1:length(x.grid)){
        tmp.sum <- 0
        for(i in 1:N){
            tmp.sum <- tmp.sum + weight.ind(x.c[i], x.grid[j], h)*y.c[i]
        }
        smooth.y[j] <- tmp.sum
    }
    Y.matrix[, k + 1] <- smooth.y
}

colnames(Y.matrix) <- c("true", "h = 1", "h = 2", "h = 5", "h = 10")
plot.df <- data.frame(cbind(x.grid, Y.matrix[, -1]))
plot.df <- plot.df %>% gather(key = "h", value = "value", -1)

plot <- ggplot(plot.df, aes(x = x.grid, y = value)) + geom_line(aes(color = h)) + theme_classic() + ggtitle("Estimated Functions with Various h") + xlab("x") + ylab("Estimated Value")
ggsave("kernel.png", plot)

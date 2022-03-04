## title: Introduction to Hierarchical Models
## author: Michael R. Schwob

###
### Set-up
###

library(dplyr)
library(ggplot2)

setwd(paste0(getwd(), "/SDS383D-Schwob/data"))
data  <- read.csv("mathtest.csv")
setwd(paste0("/home/mikel/Desktop/Code/SDS383D-Schwob/exercises/Exercise 4"))

###
### (b) Plot number of students v. school average
###

n.stud <- y.avg <- rep(0, length(unique(data[, 1])))
for(i in 1:length(n.stud)){
    n.stud[i] <- nrow(filter(data, school == i))
    y.avg[i] <- mean(filter(data, school == i)[, 2]) 
}

df <- data.frame(cbind(n.stud, y.avg))
ggplot(df, aes(x = n.stud, y = y.avg)) + theme_classic() + geom_point() + xlab("School Size") + ylab("School-level Average") + ggtitle("School-level Averages v. School Size")
ggsave("math_b.png")

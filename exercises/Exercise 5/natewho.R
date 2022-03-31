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
P <- 10

Y <- matrix(NA, n, max)
X <- array(NA, dim = c(max, P+1, n))
N.i <- rep(0, n)

for(i in 1:n){
    tmp.data <- data[which(data$state==names[i]),] # data set for current state
    N.i[i] <- dim(tmp.data)[1] # number of observations for store i
    Y[i, 1:(N.i[i])] <- tmp.data$bush # states voters
    for(j in 1:N.i[i]){
        X[j, , i] <- c(1, as.numeric(tmp.data[j, 3:12]))
    }
}

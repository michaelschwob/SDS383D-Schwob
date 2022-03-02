## title: A Heavy-tailed Error Model for Green Buildings
## author: Michael R. Schwob

###
### Set-up
###

set.seed(702)
library(ggplot2)

###
### Get Data
###

setwd(paste0(getwd(), "/SDS383D-Schwob/data"))
data <- read.csv("greenbuildings.csv")
n <- dim(data)[1]

## Response  =  revenue*ft^2/100
y <- data$Rent*data$leasing_rate/100

## design matrix
X <- cbind(rep(1:n), data$green_rating, data$City_Market_Rent, data$age, data$class_a,  data$class_b)
p <- dim(X)[2] # number of parameters
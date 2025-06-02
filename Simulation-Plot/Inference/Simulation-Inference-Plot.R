library(ggplot2)
library(softImpute)
library(latex2exp)

setwd("./Result data/")

## Here is the data to produce Figure 2. Produce Figures 4-5, 9-11 by changing to other result RData.  ##

load("density31.RData") 
load("density32.RData")
load("density33.RData")
load("density34.RData")

load("density41.RData") 
load("density42.RData")
load("density43.RData")
load("density44.RData")

Iteration <- 1000


## gamma=1/3 ##
value1df <- data.frame(value1 = density31[1:Iteration])
fig1 <- ggplot(data = value1df, aes(x = value1)) + geom_histogram(aes(y = ..density..), fill="#377eb8", color="black", alpha=0.6, stat = "bin", bins = 20) +
  stat_function(fun = dnorm, colour = "red", size = 1)+xlab("") + theme_classic() 
fig1

ggsave("density31.pdf", fig1, units = "px", width = 1027, height = 856, device="pdf")

value1df <- data.frame(value1 = density32[1:Iteration])
fig2 <- ggplot(data = value1df, aes(x = value1)) + geom_histogram(aes(y = ..density..), fill="#377eb8", color="black", alpha=0.6, stat = "bin", bins = 20) +
  stat_function(fun = dnorm, colour = "red", size = 1)+xlab("") + theme_classic() 
fig2

ggsave("density32.pdf", fig2, units = "px", width = 1027, height = 856, device="pdf")

value1df <- data.frame(value1 = density33[1:Iteration])
fig3 <- ggplot(data = value1df, aes(x = value1)) + geom_histogram(aes(y = ..density..), fill="#377eb8", color="black", alpha=0.6, stat = "bin", bins = 20) +
  stat_function(fun = dnorm, colour = "red", size = 1)+xlab("") + theme_classic() 
fig3

ggsave("density33.pdf", fig3, units = "px", width = 1027, height = 856, device="pdf")

value1df <- data.frame(value1 = density34[1:Iteration])
fig4 <- ggplot(data = value1df, aes(x = value1)) + geom_histogram(aes(y = ..density..), fill="#377eb8", color="black", alpha=0.6, stat = "bin", bins = 20) +
  stat_function(fun = dnorm, colour = "red", size = 1)+xlab("") + theme_classic() 
fig4

ggsave("density34.pdf", fig4, units = "px", width = 1027, height = 856, device="pdf")


## gamma=1/4 ##
value1df <- data.frame(value1 = density41[1:Iteration])
fig1 <- ggplot(data = value1df, aes(x = value1)) + geom_histogram(aes(y = ..density..), fill="#377eb8", color="black", alpha=0.6, stat = "bin", bins = 20) +
  stat_function(fun = dnorm, colour = "red", size = 1)+xlab("") + theme_classic() 
fig1

ggsave("density41.pdf", fig1, units = "px", width = 1027, height = 856, device="pdf")

value1df <- data.frame(value1 = density42[1:Iteration])
fig2 <- ggplot(data = value1df, aes(x = value1)) + geom_histogram(aes(y = ..density..), fill="#377eb8", color="black", alpha=0.6, stat = "bin", bins = 20) +
  stat_function(fun = dnorm, colour = "red", size = 1)+xlab("") + theme_classic() 
fig2

ggsave("density42.pdf", fig2, units = "px", width = 1027, height = 856, device="pdf")

value1df <- data.frame(value1 = density43[1:Iteration])
fig3 <- ggplot(data = value1df, aes(x = value1)) + geom_histogram(aes(y = ..density..), fill="#377eb8", color="black", alpha=0.6, stat = "bin", bins = 20) +
  stat_function(fun = dnorm, colour = "red", size = 1)+xlab("") + theme_classic() 
fig3

ggsave("density43.pdf", fig3, units = "px", width = 1027, height = 856, device="pdf")

value1df <- data.frame(value1 = density44[1:Iteration])
fig4 <- ggplot(data = value1df, aes(x = value1)) + geom_histogram(aes(y = ..density..), fill="#377eb8", color="black", alpha=0.6, stat = "bin", bins = 20) +
  stat_function(fun = dnorm, colour = "red", size = 1)+xlab("") + theme_classic() 
fig4

ggsave("density44.pdf", fig4, units = "px", width = 1027, height = 856, device="pdf")

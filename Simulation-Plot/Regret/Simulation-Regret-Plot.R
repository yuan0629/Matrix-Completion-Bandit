library(ggplot2)
library(softImpute)
library(latex2exp)

setwd("./Result data/")

load("regret3.RData") #gamma=1/3
load("regret4.RData") #gamma=1/4

## regret versus T ##

#gamma=1/3
Timeall <- seq(from = 40000, to = 80000, by = 5000)
regret3new <- apply(regret3, 2, mean)
fig1 <- ggplot() + geom_line(aes(x = Timeall^(2/3), y = regret3new), color = "#377eb8", size = 1) + 
  theme_classic() + ylab("Regret") + xlab(TeX("$T^{2/3}$")) + geom_point(aes(x = Timeall^(2/3), y = regret3new), color = "#377eb8", shape = 17, size = 2.5) + labs(tag = "A")
fig1
ggsave("regret3.pdf", fig1, units = "px", width = 1084, height = 640, device="pdf")

#gamma=1/4
Timeall <- seq(from = 20000, to = 60000, by = 5000)
regret4new <- apply(regret4, 2, mean)
fig2 <- ggplot() + geom_line(aes(x = Timeall^(3/4), y = regret4new), color = "#377eb8", size = 1) + 
  theme_classic() + ylab("Regret") + xlab(TeX("$T^{3/4}$")) + geom_point(aes(x = Timeall^(3/4), y = regret4new), color = "#377eb8", shape = 17, size = 2.5) + labs(tag = "B")
fig2
ggsave("regret4.pdf", fig1, units = "px", width = 1084, height = 640, device="pdf")



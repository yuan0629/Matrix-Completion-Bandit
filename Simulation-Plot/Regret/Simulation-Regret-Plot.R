library(ggplot2)
library(softImpute)
library(latex2exp)




regretnew <- apply(regretall[1:4,], 2, mean)
fig1 <- ggplot() + geom_line(aes(x = Timeall^(2/3), y = regretnew), color = "#377eb8", size = 1) + 
  theme_classic() + ylab("Regret") + xlab(TeX("$T^{2/3}$")) + geom_point(aes(x = Timeall^(2/3), y = regret3new), color = "#377eb8", shape = 17, size = 2.5) + labs(tag = "A")
fig1
ggsave("regret3.pdf", fig1, units = "px", width = 1084, height = 640, device="pdf")

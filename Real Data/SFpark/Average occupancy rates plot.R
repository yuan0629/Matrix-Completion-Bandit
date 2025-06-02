library(ggplot2)
load("total.RData")

## weekday of 02ND ST 200

perioddf <- data.frame(matrix(0, nrow = 55, ncol = 3))
colnames(perioddf) <- c("occupancy", "period", "time")
perioddf$period <- rep(c(1:5), 11)
perioddf$time <- rep(c(7:17), 5)
perioddf$timenum <- rep(c(1, 13, 3, 15, 5, 17, 7, 19, 9, 21, 11), 5)
for (i in 1:55) {
  perioddf$occupancy[i] <- mean(Downtown2[(Downtown2$Tinum == perioddf$timenum[i]) & (Downtown2$Nnumnew == 2) & (Downtown2$period == perioddf$period[i]), ]$occupancy)
}
perioddf$period <- as.factor(perioddf$period)
ggplot(data = perioddf) + geom_line(aes(x = time, y = occupancy, color = period), size = 0.8) + scale_colour_brewer(palette = "Paired") + theme_classic() +
  scale_x_discrete(name ="Time", limits=c(7:17)) + geom_point(aes(x = time, y = occupancy, color = period), shape = 17, size = 2.5,)
fig1 <- ggplot(data = perioddf) + geom_line(aes(x = time, y = occupancy, color = period), size = 0.8) + scale_colour_brewer(palette = "Paired") + theme_classic() +
  scale_x_discrete(name ="Time", limits=c(7:17)) + geom_point(aes(x = time, y = occupancy, color = period), shape = 17, size = 2.5,) + ylab("Occupancy Rate")
fig1


## weekend of BATTERY ST 400

perioddf <- data.frame(matrix(0, nrow = 55, ncol = 3))
colnames(perioddf) <- c("occupancy", "period", "time")
perioddf$period <- rep(c(1:5), 11)
perioddf$time <- rep(c(7:17), 5)
perioddf$timenum <- rep(c(12, 2, 14, 4, 16, 6, 18, 8, 20, 10, 22), 5)
for (i in 1:55) {
  perioddf$occupancy[i] <- mean(Downtown2[(Downtown2$Tinum == perioddf$timenum[i]) & (Downtown2$Nnumnew == 3) & (Downtown2$period == perioddf$period[i]), ]$occupancy)
}
perioddf$period <- as.factor(perioddf$period)
ggplot(data = perioddf) + geom_line(aes(x = time, y = occupancy, color = period), size = 0.8) + scale_colour_brewer(palette = "Paired") + theme_classic() +
  scale_x_discrete(name ="Time", limits=c(7:17)) + geom_point(aes(x = time, y = occupancy, color = period), shape = 17, size = 2.5,)
fig2 <- ggplot(data = perioddf) + geom_line(aes(x = time, y = occupancy, color = period), size = 0.8) + scale_colour_brewer(palette = "Paired") + theme_classic() +
  scale_x_discrete(name ="Time", limits=c(7:17)) + geom_point(aes(x = time, y = occupancy, color = period), shape = 17, size = 2.5,) + ylab("Occupancy Rate")
fig2

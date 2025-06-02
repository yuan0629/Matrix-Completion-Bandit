load("total.RData")
library(ggplot2)
library(nloptr)

reward <- function(x) {
  if (is.na(x)) {
    y <- NA
  } else if (x < 0.6) {
    y <- x - 0.6
  } else if (x > 0.8){
    y <- 0.8 - x
  } else {
    y <- 0
  }
  return(y)
}

Nnew <- length(table(Downtown2$Nnumnew))
Ti <- 22
MDowntown21 <- MDowntown20 <- matrix(0, nrow = Nnew, ncol = Ti)
for (i in 1:Nnew) {
  for (j in 1:Ti) {
    MDowntown20[i, j] <- reward(mean(Downtown2[(Downtown2$action == 0) & (Downtown2$Tinum == j) & (Downtown2$Nnumnew == i), ]$occupancy))
    MDowntown21[i, j] <- reward(mean(Downtown2[(Downtown2$action == 1) & (Downtown2$Tinum == j) & (Downtown2$Nnumnew == i), ]$occupancy))
  }
}


Time <- nrow(Downtown2)
a <- 1/3
T0 <- 8 * Time^(1-a)
r <- 5
lmax <- max(svd(M10)$d[1], svd(M00)$d[1])

U0 <- svd(M00)$u[,1:r] %*% diag(c(sqrt(svd(M00)$d)[1:r]))   
V0 <- svd(M00)$v[,1:r] %*% diag(c(sqrt(svd(M00)$d)[1:r]))
U1 <- svd(M10)$u[,1:r] %*% diag(c(sqrt(svd(M10)$d)[1:r]))   
V1 <- svd(M10)$v[,1:r] %*% diag(c(sqrt(svd(M10)$d)[1:r]))


U <- list(U1, U0)
V <- list(V1, V0)

M1hat <- M10
M0hat <- M00

sum <- rep(0, 5)
sumday <- rep(0, 5)
timeyes <- rep(0, Time)
sumoff <- rep(0, 5)

set.seed(1124)
### Matrix Completion Bandit
for (t in 1:Time){
  i <- Downtown2$Nnumnew[t]
  j <- Downtown2$Tinum[t]
  
  if (t <= T0) {
    epsilon <- 0.3
    eta <- 0.02 * Nnew *Ti*log(Nnew)/lmax/Time^(1-a)
  } else {
    epsilon <-  t^(-a)
    eta <- 0.02 * epsilon*Nnew*Ti*log(Nnew)/lmax/Time^(1-a)
  }
  
  # decision
  pit <- (1-epsilon/2) * (M1hat[i, j]> M0hat[i, j]) + epsilon/2 * (M1hat[i, j]<= M0hat[i, j])
  ac <- sample(c(0,1), prob = c(1-pit, pit), 1)
  de <- ac * pit + (1-ac) * (1-pit)
  X <- matrix(0, nrow = Nnew, ncol = Ti)
  X[i, j] <- 1
  
  ##bandit
  if (ac == 0) {
    yhat <- M0hat[j]
    if (!is.na(MDowntown20[i, j])) {
      timeyes[t] <- 1
      if (Downtown2$action[t] == 0) {
        sampley <- Downtown2$occupancy[t]
      } else if (length(Downtown2[(Downtown2$action == 0) & (Downtown2$Tinum == j) & (Downtown2$Nnumnew == i) & (Downtown2$date == Downtown2$date[t]),]$occupancy > 0)){
        sampley <- Downtown2[(Downtown2$action == 0) & (Downtown2$Tinum == j) & (Downtown2$Nnumnew == i) & (Downtown2$date == Downtown2$date[t]), "occupancy"]
      } else {
        temp <- Downtown2[(Downtown2$action == 0) & (Downtown2$Tinum == j) & (Downtown2$Nnumnew == i), ]
        tempdate <- sort(c(temp$date, Downtown2$date[t]))
        repl <- which(tempdate == Downtown2$date[t])
        if (repl >= nrow(temp)) {
          sampley <- temp$occupancy[repl - 1] 
        } else {
          sampley <- temp$occupancy[repl + 1] 
        }
      }
      y <- reward(sampley)
      if (as.numeric(y) == 0) {
        sumday[Downtown2$period[t]] <- sumday[Downtown2$period[t]] + 1
      }
    }else {
      sum[Downtown2$period[t]] <- sum[Downtown2$period[t]] + 1
      next
    }
  }else {
    yhat <- M1hat[j]
    if (!is.na(MDowntown21[i, j])) {
      timeyes[t] <- 1
      if (Downtown2$action[t] == 1) {
        sampley <- Downtown2$occupancy[t]
      } else if (length(Downtown2[(Downtown2$action == 1) & (Downtown2$Tinum == j) & (Downtown2$Nnumnew == i) & (Downtown2$date == Downtown2$date[t]),]$occupancy > 0)){
        sampley <- Downtown2[(Downtown2$action == 1) & (Downtown2$Tinum == j) & (Downtown2$Nnumnew == i) & (Downtown2$date == Downtown2$date[t]), "occupancy"]
      } else {
        temp <- Downtown2[(Downtown2$action == 1) & (Downtown2$Tinum == j) & (Downtown2$Nnumnew == i), ]
        tempdate <- sort(c(temp$date, Downtown2$date[t]))
        repl <- which(tempdate == Downtown2$date[t])
        if (repl >= nrow(temp)) {
          sampley <- temp$occupancy[repl - 1] 
        } else {
          sampley <- temp$occupancy[repl + 1] 
        }
      }
      y <- reward(sampley)
      if (as.numeric(y) == 0) {
        sumday[Downtown2$period[t]] <- sumday[Downtown2$period[t]] + 1
      }
    } else {
      sum[Downtown2$period[t]] <- sum[Downtown2$period[t]] + 1
      next
    }
  }
  
  
  
  # update U and V
  RU <- svd(t(U[[2-ac]]) %*% U[[2-ac]])$u
  DU <- diag(svd(t(U[[2-ac]]) %*% U[[2-ac]])$d)
  RV <- svd(t(V[[2-ac]]) %*% V[[2-ac]])$u
  DV <- diag(svd(t(V[[2-ac]]) %*% V[[2-ac]])$d)
  QU <- svd(sqrt(DU) %*% t(RU) %*% RV %*% sqrt(DV))$u
  #D <- svd(sqrt(DU) %*% t(RU) %*% RV %*% sqrt(DV))$d
  QV <- svd(sqrt(DU) %*% t(RU) %*% RV %*% sqrt(DV))$v
  
  U[[2-ac]] <- U[[2-ac]] - eta * 1/de * (yhat - y) * X %*% V[[2-ac]] %*% RV %*% solve(sqrt(DV)) %*% QV %*% t(QU) %*% sqrt(DU) %*% t(RU)
  V[[2-ac]] <- V[[2-ac]] - eta * 1/de * (yhat - y) * t(X) %*% U[[2-ac]] %*% RU %*% solve(sqrt(DU)) %*% QU %*% t(QV) %*% sqrt(DV) %*% t(RV)
  
  M1hat <- U[[1]] %*% t(V[[1]])
  M0hat <- U[[2]] %*% t(V[[2]])
}



Downtown2yes <- Downtown2[timeyes == 1,]
result2 <- data.frame(matrix(0, nrow = 5, ncol = 2))
colnames(result2) <- c("SFPark", "MBandit")
result2$MBandit <- sumday/(table(Downtown2$period) - sum)

### SFpark policy
for (i in 1:5) {
  result2$SFPark[i] <- sum((Downtown2yes$occupancy >= 0.6) & (Downtown2yes$occupancy <= 0.8) & (Downtown2yes$period == i))/sum((Downtown2yes$period == i))
}






### offline neighborhood method
coefficient <- matrix(0, nrow = Nnew*Ti, ncol = Nnew+1)
for (j in 1:Ti) {
  price <- matrix(0, nrow = 5, ncol = Nnew)
  for (i in 1:Nnew) {
    for (k in 1:5) {
      price[k, i] <- Downtown2[(Downtown2$period == k) & (Downtown2$Nnumnew == i) & (Downtown2$Tinum == j), ]$action[1]
    }
  }
  price[is.na(price)] <- 0
  y <- rep(5, 0)
  for (i in 1:Nnew) {
    for (k in 1:5) {
      y[k] <- mean(Downtown2[(Downtown2$period == k) & (Downtown2$Nnumnew == i) & (Downtown2$Tinum == j), ]$occupancy) 
    }
    if (all(is.na(y))) {
      next
    } else if (any(is.na(y))) {
      modellm <- lm(y[-is.na(y)]~price[-which(is.na(y)), ])
      coefficient[(j-1)*Nnew+i, ] <- modellm$coefficients 
    } else {
      modellm <- lm(y~price)
      coefficient[(j-1)*Nnew+i, ] <- modellm$coefficients 
    }
  }
}

coefficient[is.na(coefficient)] <- 0
# fn <- function(x) {
#   ot <- rep(0.8, Nnew)
#   ob <- rep(0, Nnew)
#   xplus <- c(1,x)
#   for (i in 1:Nnew) {
#     ob[i] <- sum(xplus*coefficient[(j-1)*Nnew+i, ])
#   }
#   return(sum((ot-ob)^2))
# }
# hin <- function(x) c(x, -x+1)
# auglag(price[5,], fn, gr = NULL, hin = hin)$par


actionoffline <- rep(0, Nnew*Ti)
for (j in 1:Ti) {
  price <- rep(0, 5)
  for (i in 1:Nnew) {
    price[i] <- Downtown2[(Downtown2$period == 2) & (Downtown2$Nnumnew == i) & (Downtown2$Tinum == j), ]$action[1]
  }
  price[is.na(price)] <- 0
  fn <- function(x) {
    ot <- rep(0.7, Nnew)
    ob <- rep(0, Nnew)
    xplus <- c(1,x)
    for (i in 1:Nnew) {
      ob[i] <- sum(xplus*coefficient[(j-1)*Nnew+i, ])
    }
    return(sum((ot-ob)^2))
  }
  hin <- function(x) c(x, -x+1)
  actionoffline[((j-1)*Nnew+1) : (j*Nnew)] <- (auglag(price, fn, gr = NULL, hin = hin)$par)>= 0.5
}

sumoff <- rep(0, 5)
sumdayoff <- rep(0, 5)
timeyes <- rep(0, Time)

T0 <- 8 * Time^(1-a)
Nnew <- length(table(Downtown2$Nnumnew))

r <- 5
lmax <- max(svd(M10)$d[1], svd(M00)$d[1])

U0 <- svd(M00)$u[,1:r] %*% diag(c(sqrt(svd(M00)$d)[1:r]))   
V0 <- svd(M00)$v[,1:r] %*% diag(c(sqrt(svd(M00)$d)[1:r]))
U1 <- svd(M10)$u[,1:r] %*% diag(c(sqrt(svd(M10)$d)[1:r]))   
V1 <- svd(M10)$v[,1:r] %*% diag(c(sqrt(svd(M10)$d)[1:r]))


U <- list(U1, U0)
V <- list(V1, V0)

M1hat <- M10
M0hat <- M00

set.seed(1124)
for (t in 1:Time){
  
  
  i <- Downtown2$Nnumnew[t]
  j <- Downtown2$Tinum[t]
  
  
  # decision
  ac <- actionoffline[(j-1)*Nnew+i]
  
  ##bandit
  if (ac == 0) {
    yhat <- M0hat[j]
    if (!is.na(MDowntown20[i, j])) {
      timeyes[t] <- 1
      if (Downtown2$action[t] == 0) {
        sampley <- Downtown2$occupancy[t]
      } else if (length(Downtown2[(Downtown2$action == 0) & (Downtown2$Tinum == j) & (Downtown2$Nnumnew == i) & (Downtown2$date == Downtown2$date[t]),]$occupancy > 0)){
        sampley <- Downtown2[(Downtown2$action == 0) & (Downtown2$Tinum == j) & (Downtown2$Nnumnew == i) & (Downtown2$date == Downtown2$date[t]), "occupancy"]
      } else {
        temp <- Downtown2[(Downtown2$action == 0) & (Downtown2$Tinum == j) & (Downtown2$Nnumnew == i), ]
        tempdate <- sort(c(temp$date, Downtown2$date[t]))
        repl <- which(tempdate == Downtown2$date[t])
        if (repl >= nrow(temp)) {
          sampley <- temp$occupancy[repl - 1] 
        } else {
          sampley <- temp$occupancy[repl + 1] 
        }
      }
      y <- reward(sampley)
      if (as.numeric(y) == 0) {
        sumdayoff[Downtown2$period[t]] <- sumdayoff[Downtown2$period[t]] + 1
      }
    }else {
      sumoff[Downtown2$period[t]] <- sumoff[Downtown2$period[t]] + 1
      next
    }
  }else {
    yhat <- M1hat[j]
    if (!is.na(MDowntown21[i, j])) {
      timeyes[t] <- 1
      if (Downtown2$action[t] == 1) {
        sampley <- Downtown2$occupancy[t]
      } else if (length(Downtown2[(Downtown2$action == 1) & (Downtown2$Tinum == j) & (Downtown2$Nnumnew == i) & (Downtown2$date == Downtown2$date[t]),]$occupancy > 0)){
        sampley <- Downtown2[(Downtown2$action == 1) & (Downtown2$Tinum == j) & (Downtown2$Nnumnew == i) & (Downtown2$date == Downtown2$date[t]), "occupancy"]
      } else {
        temp <- Downtown2[(Downtown2$action == 1) & (Downtown2$Tinum == j) & (Downtown2$Nnumnew == i), ]
        tempdate <- sort(c(temp$date, Downtown2$date[t]))
        repl <- which(tempdate == Downtown2$date[t])
        if (repl >= nrow(temp)) {
          sampley <- temp$occupancy[repl - 1] 
        } else {
          sampley <- temp$occupancy[repl + 1] 
        }
      }
      y <- reward(sampley)
      if (as.numeric(y) == 0) {
        sumdayoff[Downtown2$period[t]] <- sumdayoff[Downtown2$period[t]] + 1
      }
    } else {
      sumoff[Downtown2$period[t]] <- sumoff[Downtown2$period[t]] + 1
      next
    }
  }
  
}

result2$offline <- sumdayoff/(table(Downtown2$period) - sumoff)


## Figure 8
df2 <- data.frame(matrix(0, ncol = 3, nrow = 15))
colnames(df2) <- c("period", "policy", "occupancy")
df2$period <- rep(c(1:5), 3)
df2$occupancy <- c(result2$MBandit, result2$SFPark, result2$offline)
df2$policy <- c(rep("MCBandit", 5), rep("SFPark", 5),  rep("Neighbor", 5))
fig1 <- ggplot(df2) + geom_line(aes(x = period, y = occupancy, color = policy), linewidth = 1) + 
  geom_point(aes(x = period, y = occupancy, color = policy), size = 2.5, shape = 17) + 
  theme_classic() + ylab("Percentage") + xlab("Period")+ scale_color_manual(values = c("#386cb0", "#fdb462", "#7fc97f"), labels = c("MCBandit", "Neighborhood", "SFPark"))
fig1
ggsave("experiment2.pdf", fig1, units = "px", width = 1440, height = 570, device="pdf")


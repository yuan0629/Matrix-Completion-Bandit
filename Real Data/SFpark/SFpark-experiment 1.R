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


#scree plot Figure 6
Nnew <- length(table(Downtown2$Nnumnew))
Ti <- 22

pca1 <- prcomp(M10, scale = TRUE)
eigenvalue1 <- pca1$sdev^2
dataeigenvalue <- data.frame(pc = c(1:Ti), value = eigenvalue1)
fig <- ggplot(data = dataeigenvalue, aes(x = pc, y = value)) + geom_line(color = "#386cb0") + geom_point(color = "#386cb0", shape = 1) + xlab("Factor Numbers") + ylab("Eigenvalues of Principal Factors") + theme_classic()
fig
ggsave(plot = fig, filename = "sfparkscree1.pdf", units = "px", width = 1163, height = 763, device="pdf")

pca0 <- prcomp(M00, scale = TRUE)
eigenvalue1 <- pca0$sdev^2
dataeigenvalue <- data.frame(pc = c(1:Ti), value = eigenvalue1)
fig <- ggplot(data = dataeigenvalue, aes(x = pc, y = value)) + geom_line(color = "#386cb0") + geom_point(color = "#386cb0", shape = 1) + xlab("Factor Numbers") + ylab("Eigenvalues of Principal Factors") + theme_classic()
fig
ggsave(plot = fig, filename = "sfparkscree0.pdf", units = "px", width = 1163, height = 763, device="pdf")



## Experiment1
Time <- nrow(Downtown2)
a <- 1/3
T0 <- 10*Time^(1-a)
r <- 5
lmax <- max(svd(M10)$d[1], svd(M00)$d[1])

U0 <- svd(M00)$u[,1:r] %*% diag(c(sqrt(svd(M00)$d)[1:r]))   
V0 <- svd(M00)$v[,1:r] %*% diag(c(sqrt(svd(M00)$d)[1:r]))
U1 <- svd(M10)$u[,1:r] %*% diag(c(sqrt(svd(M10)$d)[1:r]))   
V1 <- svd(M10)$v[,1:r] %*% diag(c(sqrt(svd(M10)$d)[1:r]))

result1 <- matrix(0, ncol = 8, nrow = 22)
result1[,7] <- c(rep(2,11), rep(3, 11))
result1[,8] <- c(1, 13, 3, 15, 5, 17, 7, 19, 9, 21, 11, 12, 2, 14, 4, 16, 6, 18, 8, 20, 10, 22)




M1hat <- M10
M0hat <- M00

U <- list(U1, U0)
V <- list(V1, V0)


T0new <- 0
Timenew <- 0

sigma1 <- sigma0 <- 0
S1 <- S0 <- 0
M1unbs <- matrix(0, nrow = Nnew, ncol = Ti)
M0unbs <- matrix(0, nrow = Nnew, ncol = Ti)

# algorithm
set.seed(1124)
for (t in 1:Time){
  i <- Downtown2$Nnumnew[t]
  j <- Downtown2$Tinum[t]
  
  if (t <= T0) {
    epsilon <- 0.2
    eta <- 0.025 * Nnew *Ti*log(Nnew)/lmax/Time^(1-a)
  } else {
    epsilon <-  t^(-a)
    eta <- 0.025 * epsilon*Nnew*Ti*log(Nnew)/lmax/Time^(1-a)
  }
  
  # decision
  pit <- (1-epsilon/2) * (M1hat[i, j]> M0hat[i, j]) + epsilon/2 * (M1hat[i, j]<= M0hat[i, j])
  ac <- sample(c(0,1), prob = c(1-pit, pit), 1)
  de <- ac * pit + (1-ac) * (1-pit)
  X <- matrix(0, nrow = Nnew, ncol = Ti)
  X[i, j] <- 1
  if (ac != Downtown2$action[t]) {
    next
  } else {
    Timenew <- Timenew + 1
    if (ac == 0) {
      yhat <- M0hat[j]
      y <- reward(Downtown2$occupancy[t])
    }else {
      yhat <- M1hat[j]
      y <- reward(Downtown2$occupancy[t])
    }
  }
  
  
  # inference way1
  if (t > T0) {
    T0new <- T0new + 1
    
    M1unbs <- M1unbs + M1hat
    M0unbs <- M0unbs + M0hat
    
    #  debias
    if (ac == 1) {
      sigma1 <- sigma1 + 1 / pit * (y - yhat)^2
      
      # debias
      M1unbs <- M1unbs + Nnew*Ti/pit*(y - yhat) * X
    } else {
      sigma0 <- sigma0 + 1 / (1-pit) * (y - yhat)^2
      
      # debias
      M0unbs <- M0unbs + Nnew*Ti/(1-pit)*(y - yhat) * X
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

b <- Timenew / T0new

# projection
M1unbs <- M1unbs / T0new
L1 <- svd(M1unbs)$u[,1:r]
R1 <- svd(M1unbs)$v[,1:r]
M1proj <- L1 %*% t(L1) %*% M1unbs %*% R1 %*% t(R1)
Lproj1 <- svd(M1proj)$u[,1:r]
Rproj1 <- svd(M1proj)$v[,1:r]
Lproj1c <- svd(M1proj)$u[,(r+1):Ti]
Rproj1c <- svd(M1proj)$v[,(r+1):Ti]

M0unbs <- M0unbs / T0new
L0 <- svd(M0unbs)$u[,1:r]
R0 <- svd(M0unbs)$v[,1:r]
M0proj <- L0 %*% t(L0) %*% M0unbs %*% R0 %*% t(R0)
Lproj0 <- svd(M0proj)$u[,1:r]
Rproj0 <- svd(M0proj)$v[,1:r]
Lproj0c <- svd(M0proj)$u[,(r+1):Ti]
Rproj0c <- svd(M0proj)$v[,(r+1):Ti]

sigma1 <- sigma1 / T0new
sigma0 <- sigma0 / T0new

for (iter in 1:nrow(result1)) {  
  
  print(iter)
  
  w1 <- result1[iter,7]
  w2 <- result1[iter,8]
  W <- matrix(0, nrow = Nnew, ncol = Ti)
  W[w1,w2] <- 1
  
  
  PMw <- Lproj1 %*% t(Lproj1) %*% W %*% Rproj1 %*% t(Rproj1) + Lproj1 %*% t(Lproj1) %*% W %*% Rproj1c %*% t(Rproj1c) + 
    Lproj1c %*% t(Lproj1c) %*% W %*% Rproj1 %*% t(Rproj1)
  S1 <- b*(1/(1+a)*2*sum((PMw[M1hat < M0hat])^2) + 1/Timenew^a*sum((PMw[M1hat > M0hat])^2))
  
  
  PMw <- Lproj0 %*% t(Lproj0) %*% W %*% Rproj0 %*% t(Rproj0) + Lproj0 %*% t(Lproj0) %*% W %*% Rproj0c %*% t(Rproj0c) + 
    Lproj0c %*% t(Lproj0c) %*% W %*% Rproj0 %*% t(Rproj0)
  S0 <- b*(1/(1+a)*2*sum((PMw[M1hat > M0hat])^2) + 1/Timenew^a*sum((PMw[M1hat < M0hat])^2))
  
  
  # result
  result1[iter, 1] <- 1 - pnorm((M1proj-M0proj)[w1,w2], 0, sqrt((sigma1*S1 + sigma0*S0) * Nnew * Ti /Timenew^(1-a)))
  result1[iter, 2] <- (M1proj-M0proj)[w1,w2]
  result1[iter, 3] <- S1
  result1[iter, 4] <- S0
  result1[iter, 5] <- (M1proj-M0proj)[w1,w2] - qnorm(0.975) * sqrt((sigma1*S1 + sigma0*S0) * Nnew * Ti /Timenew^(1-a))
  result1[iter, 6] <- (M1proj-M0proj)[w1,w2] + qnorm(0.975) * sqrt((sigma1*S1 + sigma0*S0) * Nnew * Ti /Timenew^(1-a))
}


### compare with offline method
M1unbsoff <- M10
num1 <- sum(Downtown2$action == 1)
sigma1off <- 0
for (i in 1:nrow(Downtown2)) {
  if (Downtown2$action[i] == 1) {
    w1 <- Downtown2$Nnumnew[i]
    w2 <- Downtown2$Tinum[i]
    X <- matrix(0, nrow = Nnew, ncol = Ti)
    X[w1,w2] <- 1
    y <- reward(Downtown2$occupancy[i])
    M1unbsoff <- M1unbsoff + Nnew * Ti /num1 * (y - M10[w1,w2]) * X
    sigma1off <- sigma1off + 1 / num1* (y - M10[w1,w2])^2
  }
}


M0unbsoff <- M00
num0 <- sum(Downtown2$action == 0)
sigma0off <- 0
for (i in 1:nrow(Downtown2)) {
  if (Downtown2$action[i] == 0) {
    w1 <- Downtown2$Nnumnew[i]
    w2 <- Downtown2$Tinum[i]
    X <- matrix(0, nrow = Nnew, ncol = Ti)
    X[w1,w2] <- 1
    y <- reward(Downtown2$occupancy[i])
    M0unbsoff <- M0unbsoff + Nnew * Ti /num0 * (y - M00[w1,w2]) * X 
    sigma0off <- sigma0off + 1 / num0 * (y - M00[w1,w2])^2
  }
}

resultoff <- matrix(0, ncol = 8, nrow = 22)
resultoff[,7] <- c(rep(2,11), rep(3, 11))
resultoff[,8] <- c(1, 13, 3, 15, 5, 17, 7, 19, 9, 21, 11, 12, 2, 14, 4, 16, 6, 18, 8, 20, 10, 22)


for (iter in 1:nrow(resultoff)) {
  
  print(iter)
  
  w1 <- resultoff[iter,7]
  w2 <- resultoff[iter,8]
  W <- matrix(0, nrow = Nnew, ncol = Ti)
  W[w1,w2] <- 1
  
  
  L1off <- svd(M1unbsoff)$u[,1:r]
  R1off <- svd(M1unbsoff)$v[,1:r]
  M1projoff <- L1off %*% t(L1off) %*% M1unbsoff %*% R1off %*% t(R1off)
  Lproj1off <- svd(M1projoff)$u[,1:r]
  Rproj1off <- svd(M1projoff)$v[,1:r]
  Lproj1coff <- svd(M1projoff)$u[,(r+1):Ti]
  Rproj1coff <- svd(M1projoff)$v[,(r+1):Ti]
  
  L0off <- svd(M0unbsoff)$u[,1:r]
  R0off <- svd(M0unbsoff)$v[,1:r]
  M0projoff <- L0off %*% t(L0off) %*% M0unbsoff %*% R0off %*% t(R0off)
  Lproj0off <- svd(M0projoff)$u[,1:r]
  Rproj0off <- svd(M0projoff)$v[,1:r]
  Lproj0coff <- svd(M0projoff)$u[,(r+1):Ti]
  Rproj0coff <- svd(M0projoff)$v[,(r+1):Ti]
  
  PMwoff <- Lproj1off %*% t(Lproj1off) %*% W %*% Rproj1off %*% t(Rproj1off) + Lproj1off %*% t(Lproj1off) %*% W %*% Rproj1coff %*% t(Rproj1coff) + 
    Lproj1coff %*% t(Lproj1coff) %*% W %*% Rproj1off %*% t(Rproj1off)
  S1off <- sum(PMwoff^2)
  
  
  PMwoff <- Lproj0off %*% t(Lproj0off) %*% W %*% Rproj0off %*% t(Rproj0off) + Lproj0off %*% t(Lproj0off) %*% W %*% Rproj0coff %*% t(Rproj0coff) + 
    Lproj0coff %*% t(Lproj0coff) %*% W %*% Rproj0off %*% t(Rproj0off)
  S0off <- sum(PMwoff^2)
  
  resultoff[iter,1] <- 1 - pnorm((M1projoff-M0projoff)[w1,w2], 0, sqrt((sigma0off*S0off/num0 + sigma1off*S1off/num1) * Nnew * Ti))
  resultoff[iter,2] <- M1projoff[w1,w2] - M0projoff[w1,w2]
  resultoff[iter,3] <- S1off
  resultoff[iter,4] <- S0off
  resultoff[iter,5] <- M1projoff[w1,w2] - M0projoff[w1,w2] - qnorm(0.975) * sqrt((sigma0off*S0off/num0 + sigma1off*S1off/num1) * Nnew * Ti)
  resultoff[iter,6] <- M1projoff[w1,w2] - M0projoff[w1,w2] + qnorm(0.975) * sqrt((sigma0off*S0off/num0 + sigma1off*S1off/num1) * Nnew * Ti)
}

## table 1
table1 <- rbind(c(1-result1[1:5,1], result1[6:11,1]), c(1-resultoff[1:5,1], resultoff[6:11,1]))

## table 2
table2 <- rbind(1-result1[12:22,1], 1-resultoff[12:22,1])


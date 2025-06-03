load("datadiscount.RData")
library(ggplot2)
library(softImpute)




## initial estimate
set.seed(629122)
fit1 <- softImpute(M1temp, rank.max = 30, lambda = 0.5)
M10 <- complete(M1temp, fit1)

fit2 <- softImpute(M0temp, rank.max = 30, lambda = 0.5)
M00 <- complete(M0temp, fit2)


Time <- nrow(datafinal)
a <- 1/3
T0 <- 0.5*Time^(1-a)
N <- length(table(datafinal$customer))
Ti <- length(table(datafinal$product))
b <- Time/(Time - T0)
c2 <- 0.01

# M10 5columns all 0
for(i in 1:ncol(M10)) {
  if (all(M10[,i] == 0)) {
    M10[,i] <- runif(N, 0, 0.1)
  }
}


#scree plot
pca1 <- prcomp(M10, scale = TRUE)
eigenvalue1 <- pca1$sdev^2
dataeigenvalue <- data.frame(pc = c(1:N), value = eigenvalue1)
ggplot(data = dataeigenvalue, aes(x = pc, y = value)) + geom_line(color = "#386cb0") + geom_point(color = "#386cb0", shape = 1) + xlab("factor numbers") + ylab("eigenvalues of principal factors") + theme_classic()
fig <- ggplot(data = dataeigenvalue, aes(x = pc, y = value)) + geom_line(color = "#386cb0") + geom_point(color = "#386cb0", shape = 1) + xlab("Factor Numbers") + ylab("Eigenvalues of Principal Factors") + theme_classic()
ggsave(plot = fig, filename = "discountscree1.jpg", units = "px", width = 1163, height = 763)


pca0 <- prcomp(M00, scale = TRUE)
eigenvalue1 <- pca0$sdev^2
dataeigenvalue <- data.frame(pc = c(1:N), value = eigenvalue1)
ggplot(data = dataeigenvalue, aes(x = pc, y = value)) + geom_line(color = "#386cb0") + geom_point(color = "#386cb0", shape = 1) + xlab("factor numbers") + ylab("eigenvalues of principal factors") + theme_classic()
fig <- ggplot(data = dataeigenvalue, aes(x = pc, y = value)) + geom_line(color = "#386cb0") + geom_point(color = "#386cb0", shape = 1) + xlab("Factor Numbers") + ylab("Eigenvalues of Principal Factors") + theme_classic()
ggsave(plot = fig, filename = "discountscree0.jpg", units = "px", width = 1163, height = 763)


r <- 30
lmax <- max(svd(M10)$d[1], svd(M00)$d[1])


cus <- 7
productlist <- datafinal[datafinal$customer == cus, ]$product
#productlist <- c(1:10)
result <- matrix(0, nrow = length(productlist), ncol = 6)

sum1 <- sum0 <- 0

## Algorithm
# for (iter in 1:length(productlist)) {
  
  #iter <- 5
  #print(iter)
  
  U0 <- svd(M00)$u[,1:r] %*% diag(c(sqrt(svd(M00)$d)[1:r]))   
  V0 <- svd(M00)$v[,1:r] %*% diag(c(sqrt(svd(M00)$d)[1:r]))
  U1 <- svd(M10)$u[,1:r] %*% diag(c(sqrt(svd(M10)$d)[1:r]))   
  V1 <- svd(M10)$v[,1:r] %*% diag(c(sqrt(svd(M10)$d)[1:r]))
  
  
  M1hat <- M10
  M0hat <- M00
  
  U <- list(U1, U0)
  V <- list(V1, V0)
  
  # w1 <- cus
  # w2 <- productlist[iter]
  # W <- matrix(0, nrow = N, ncol = Ti)
  # W[w1,w2] <- 1
  
  T0new <- 0
  Timenew <- 0
  
  sigma1 <- sigma0 <- 0
  S1 <- S0 <- 0
  M1unbs <- matrix(0, nrow = N, ncol = Ti)
  M0unbs <- matrix(0, nrow = N, ncol = Ti)
  
  # algorithm
  #permutation <- sample(Time, Time)
  
  for (t in 1:Time){
    
    i <- datafinal$customer[t]
    j <- datafinal$product[t]
    
    if (t <= T0) {
      epsilon <- 0.01
      eta <- 0.00001 * N *Ti*log(N)/lmax/Time^(1-a)
    } else {
      epsilon <-  c2*t^(-a)
      eta <- 0.00001 * epsilon*N*Ti*log(N)/lmax/Time^(1-a)
    }
    
    pit <- (1-epsilon/2) * (M1hat[i, j]> M0hat[i, j]) + epsilon/2 * (M1hat[i, j]<= M0hat[i, j])
    ac <- sample(c(0,1), prob = c(1-pit, pit), 1)
    de <- ac * pit + (1-ac) * (1-pit)
    X <- matrix(0, nrow = N, ncol = Ti)
    X[i, j] <- 1
    if (ac != datafinal$action[t]) {
      next
    } else {
      Timenew <- Timenew + 1
      if (ac == 0) {
        yhat <- M0hat[j]
        y <- datafinal$profitrate[t]
        sum0 <- sum0 + 1
      }else {
        yhat <- M1hat[j]
        y <- datafinal$profitrate[t]
        sum1 <- sum1 + 1
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
        M1unbs <- M1unbs + N*Ti/pit*(y - yhat) * X
      } else {
        sigma0 <- sigma0 + 1 / (1-pit) * (y - yhat)^2
        
        # debias
        M0unbs <- M0unbs + N*Ti/(1-pit)*(y - yhat) * X
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
  
  #which(permutation == t)
  
  b <- Timenew / T0new
  
  # projection
  M1unbs <- M1unbs / T0new
  L1 <- svd(M1unbs)$u[,1:r]
  R1 <- svd(M1unbs)$v[,1:r]
  M1proj <- L1 %*% t(L1) %*% M1unbs %*% R1 %*% t(R1)
  Lproj1 <- svd(M1proj)$u[,1:r]
  Rproj1 <- svd(M1proj)$v[,1:r]
  Lproj1c <- svd(M1proj)$u[,(r+1):N]
  Rproj1c <- svd(M1proj)$v[,(r+1):N]
  
  M0unbs <- M0unbs / T0new
  L0 <- svd(M0unbs)$u[,1:r]
  R0 <- svd(M0unbs)$v[,1:r]
  M0proj <- L0 %*% t(L0) %*% M0unbs %*% R0 %*% t(R0)
  Lproj0 <- svd(M0proj)$u[,1:r]
  Rproj0 <- svd(M0proj)$v[,1:r]
  Lproj0c <- svd(M0proj)$u[,(r+1):N]
  Rproj0c <- svd(M0proj)$v[,(r+1):N]
  
  sigma1 <- sigma1 / T0new
  sigma0 <- sigma0 / T0new
  
  for(iter in 1:length(productlist)) {
    
    print(iter)
    
    w1 <- cus
    w2 <- productlist[iter]
    W <- matrix(0, nrow = N, ncol = Ti)
    W[w1,w2] <- 1
    
    PMw <- Lproj1 %*% t(Lproj1) %*% W %*% Rproj1 %*% t(Rproj1) + Lproj1 %*% t(Lproj1) %*% W %*% Rproj1c %*% t(Rproj1c) + 
      Lproj1c %*% t(Lproj1c) %*% W %*% Rproj1 %*% t(Rproj1)
    S1 <- b*(1/(1+a)*2/c2*sum((PMw[M1hat < M0hat])^2) + 1/Timenew^a*sum((PMw[M1hat > M0hat])^2))
    
    PMw <- Lproj0 %*% t(Lproj0) %*% W %*% Rproj0 %*% t(Rproj0) + Lproj0 %*% t(Lproj0) %*% W %*% Rproj0c %*% t(Rproj0c) + 
      Lproj0c %*% t(Lproj0c) %*% W %*% Rproj0 %*% t(Rproj0)
    S0 <- b*(1/(1+a)*2/c2*sum((PMw[M1hat > M0hat])^2) + 1/Timenew^a*sum((PMw[M1hat < M0hat])^2))
    
    # result
    result[iter,1] <- 1 - pnorm((M1proj-M0proj)[w1,w2], 0, sqrt((sigma1*S1 + sigma0*S0) * N * Ti /Timenew^(1-a)))
    result[iter,2] <- M1proj[w1,w2] - M0proj[w1,w2]
    result[iter,3] <- S1
    result[iter,4] <- S0
    result[iter,5] <- (M1proj-M0proj)[w1,w2] - qnorm(0.975) * sqrt((sigma1*S1 + sigma0*S0) * N * Ti /Timenew^(1-a))
    result[iter,6] <- (M1proj-M0proj)[w1,w2] + qnorm(0.975) * sqrt((sigma1*S1 + sigma0*S0) * N * Ti /Timenew^(1-a))
    
  }
  

### compare with offline method
M1unbsoff <- M10
num1 <- sum(datafinal$action == 1)
sigma1off <- 0
for (i in 1:nrow(datafinal)) {
  if (datafinal$action[i] == 1) {
    w1 <- datafinal$customer[i]
    w2 <- datafinal$product[i]
    X <- matrix(0, nrow = N, ncol = Ti)
    X[w1,w2] <- 1
    M1unbsoff <- M1unbsoff + N * Ti /num1 * (datafinal$profitrate[i] - M10[w1,w2]) * X
    sigma1off <- sigma1off + 1 / num1* (datafinal$profitrate[i] - M10[w1,w2])^2
  }
}


M0unbsoff <- M00
num0 <- sum(datafinal$action == 0)
sigma0off <- 0
for (i in 1:nrow(datafinal)) {
  if (datafinal$action[i] == 0) {
    w1 <- datafinal$customer[i]
    w2 <- datafinal$product[i]
    X <- matrix(0, nrow = N, ncol = Ti)
    X[w1,w2] <- 1
    M0unbsoff <- M0unbsoff + N * Ti /num0 * (datafinal$profitrate[i] - M00[w1,w2]) * X 
    sigma0off <- sigma0off + 1 / num0 * (datafinal$profitrate[i] - M00[w1,w2])^2
  }
}

resultoff <- matrix(0, nrow = length(productlist), ncol = 6)


for (iter in 1:nrow(resultoff)) {
  
  print(iter)
  
  w1 <- cus
  w2 <- productlist[iter]
  W <- matrix(0, nrow = N, ncol = Ti)
  W[w1,w2] <- 1
  
  
  L1off <- svd(M1unbsoff)$u[,1:r]
  R1off <- svd(M1unbsoff)$v[,1:r]
  M1projoff <- L1off %*% t(L1off) %*% M1unbsoff %*% R1off %*% t(R1off)
  Lproj1off <- svd(M1projoff)$u[,1:r]
  Rproj1off <- svd(M1projoff)$v[,1:r]
  Lproj1coff <- svd(M1projoff)$u[,(r+1):N]
  Rproj1coff <- svd(M1projoff)$v[,(r+1):N]
  
  L0off <- svd(M0unbsoff)$u[,1:r]
  R0off <- svd(M0unbsoff)$v[,1:r]
  M0projoff <- L0off %*% t(L0off) %*% M0unbsoff %*% R0off %*% t(R0off)
  Lproj0off <- svd(M0projoff)$u[,1:r]
  Rproj0off <- svd(M0projoff)$v[,1:r]
  Lproj0coff <- svd(M0projoff)$u[,(r+1):N]
  Rproj0coff <- svd(M0projoff)$v[,(r+1):N]
  
  #sigma1 <- 1
  PMwoff <- Lproj1off %*% t(Lproj1off) %*% W %*% Rproj1off %*% t(Rproj1off) + Lproj1off %*% t(Lproj1off) %*% W %*% Rproj1coff %*% t(Rproj1coff) + 
    Lproj1coff %*% t(Lproj1coff) %*% W %*% Rproj1off %*% t(Rproj1off)
  S1off <- sum(PMwoff^2)
  
  
  #sigma0 <- 1
  PMwoff <- Lproj0off %*% t(Lproj0off) %*% W %*% Rproj0off %*% t(Rproj0off) + Lproj0off %*% t(Lproj0off) %*% W %*% Rproj0coff %*% t(Rproj0coff) + 
    Lproj0coff %*% t(Lproj0coff) %*% W %*% Rproj0off %*% t(Rproj0off)
  S0off <- sum(PMwoff^2)
  
  resultoff[iter,1] <- 1 - pnorm((M1projoff-M0projoff)[w1,w2], 0, sqrt((sigma0off*S0off/num0 + sigma1off*S1off/num1) * N * Ti))
  resultoff[iter,2] <- M1projoff[w1,w2] - M0projoff[w1,w2]
  resultoff[iter,3] <- S1off
  resultoff[iter,4] <- S0off
  resultoff[iter,5] <- M1projoff[w1,w2] - M0projoff[w1,w2] - qnorm(0.975) * sqrt((sigma0off*S0off/num0 + sigma1off*S1off/num1) * N * Ti)
  resultoff[iter,6] <- M1projoff[w1,w2] - M0projoff[w1,w2] + qnorm(0.975) * sqrt((sigma0off*S0off/num0 + sigma1off*S1off/num1) * N * Ti)
}


## plot
load("resultdiscount.RData")
cus <- 7
productlist <- datafinal[datafinal$customer == cus, ]$product
#variablename <- as.factor(productlist)
variablename <- datafinal[datafinal$customer == cus, ]$Product.ID
sigma1 <- 0.5*sigma1off
sigma0 <- 0.5*sigma0off
table1 <- data.frame(variable = variablename[-c(2,6,10)], value = result[-c(2,6,10), 2], low = result[-c(2,6,10), 2] - qnorm(0.975) * sqrt((sigma1*result[-c(2,6,10),3] + sigma0*result[-c(2,6,10),4]) * N * Ti /Timenew^(1-a)), 
                     high = result[-c(2,6,10), 2] + qnorm(0.975) * sqrt((sigma1*result[-c(2,6,10),3] + sigma0*result[-c(2,6,10),4]) * N * Ti /Timenew^(1-a)))
table1$method <- "MCB"
table2 <- data.frame(variable = variablename[-c(2,6,10)], value = resultoff[-c(2,6,10), 2], low = resultoff[-c(2,6,10), 2] - qnorm(0.975) * sqrt((sigma0*resultoff[-c(2,6,10),4]/num0 + sigma1*resultoff[-c(2,6,10),3]/num1) * N * Ti), 
                     high = resultoff[-c(2,6,10), 2] + qnorm(0.975) * sqrt((sigma0*resultoff[-c(2,6,10),4]/num0 + sigma1*resultoff[-c(2,6,10),3]/num1) * N * Ti))
table2$method <- "Offline"
tablefinal <- rbind(table1, table2)
ggplot(tablefinal, aes(x = variable, y = value, group = method, color = method)) + geom_point(position=position_dodge(0.6)) + scale_color_manual(
  values = c("#386cb0", "#fdb462")
) + 
  geom_errorbar(aes(ymin = low, ymax = high),  position=position_dodge(0.6)) + xlab("Product") + ylab("") + theme_classic() + theme(axis.text.x = element_text(angle = 15, vjust = 0.5))

fig <- ggplot(tablefinal, aes(x = variable, y = value, group = method, color = method)) + geom_point(position=position_dodge(0.6)) + scale_color_manual(
  values = c("#386cb0", "#fdb462")
) + 
  geom_errorbar(aes(ymin = low, ymax = high),  position=position_dodge(0.6)) + xlab("Product") + ylab("") + theme_classic() + theme(axis.text.x = element_text(angle = 15, vjust = 0.5))
fig
ggsave(plot = fig, filename = "ResultDiscount.jpg")



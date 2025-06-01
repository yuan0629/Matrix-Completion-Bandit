library(ggplot2)
library(softImpute)
library(latex2exp)



### Regret ###
Iteration <- 100
Timeall <- seq(from = 40000, to = 80000, by = 5000)
a <- 1/3 # exploration probability for phase two \varepsilon_t=c_2t^{-a}
N <- 300 # number of rows
Ti <- 300 # number of columns
r <- 2 # rank
c2 <- 10 # exploration probability for phase two\varepsilon_t=c_2t^{-a}


set.seed(1124)
M1 <- matrix(runif(N*Ti, -100, 100), nrow = N, ncol = Ti)
M0 <- M1 + matrix(runif(N*Ti, -2, 2), nrow = N, ncol = Ti)
regretall <- matrix(0, nrow = Iteration, ncol = length(Timeall))
regretall1 <- matrix(0, nrow = Iteration, ncol = length(Timeall))
U1 <- svd(M1)$u[,1:r] %*% diag(c(sqrt(svd(M1)$d)[1:r]))   
V1 <- svd(M1)$v[,1:r] %*% diag(c(sqrt(svd(M1)$d)[1:r]))
U0 <- svd(M0)$u[,1:r] %*% diag(c(sqrt(svd(M0)$d)[1:r]))   
V0 <- svd(M0)$v[,1:r] %*% diag(c(sqrt(svd(M0)$d)[1:r]))
M1 <- U1 %*% t(V1)
M0 <- U0 %*% t(V0)
lmax <- max(svd(M1)$d[1], svd(M0)$d[1])

# initial estimate
set.seed(1124)
M10 <- M00 <- matrix(NA, nrow = N, ncol = Ti)

ran <- sample(c(1:(N*Ti)), 50*N, replace = FALSE)
M10[ran] <- M1[ran]
fit1 <- softImpute(M10, rank.max = r, lambda = 0.5)
M10 <- complete(M10, fit1)


ran <- sample(c(1:(N*Ti)), 50*N, replace = FALSE)
M00[ran] <- M0[ran]
fit2 <- softImpute(M00, rank.max = r, lambda = 0.5)
M00 <- complete(M00, fit2)



regretall <- matrix(0, nrow = Iteration, ncol = length(Timeall))
for (iter in 1:Iteration) {
  print(iter)
  for (i in 1:length(Timeall)) {
    #print(i)
    Time <- Timeall[i]
    M1hat <- M10
    M0hat <- M00
    regret <- 0
    T0 <- 13.5*Time^(1-a) 
    
    U0 <- svd(M00)$u[,1:r] %*% diag(c(sqrt(svd(M00)$d)[1:r]))   
    V0 <- svd(M00)$v[,1:r] %*% diag(c(sqrt(svd(M00)$d)[1:r]))
    U1 <- svd(M10)$u[,1:r] %*% diag(c(sqrt(svd(M10)$d)[1:r]))   
    V1 <- svd(M10)$v[,1:r] %*% diag(c(sqrt(svd(M10)$d)[1:r]))
    
    
    U <- list(U1, U0)
    V <- list(V1, V0)
    set.seed(iter)
    for (t in 1:Time){
      j <- sample(c(1:(N*Ti)), 1)
      
      if (t <= T0) {
        epsilon <- 0.6 # exploration probability for phase one
        eta <- 0.025*N*Ti*log(N)/lmax/Time^(1-a) # stepsize
      } else {
        epsilon <- c2*t^(-a)
        eta <- 0.025*epsilon*N*Ti*log(N)/lmax/Time^(1-a) # stepsize
      }
      
      # decision
      pit <- (1-epsilon/2) * (M1hat[j]> M0hat[j]) + epsilon/2 * (M1hat[j]<= M0hat[j])
      ac <- sample(c(0,1), prob = c(1-pit, pit), 1)
      de <- ac * pit + (1-ac) * (1-pit)
      X <- matrix(0, nrow = N, ncol = Ti)
      X[j] <- 1
      if (ac == 0) {
        yhat <- M0hat[j]
        y <- M0[j] + rnorm(1, 0, 1)
      }else {
        yhat <- M1hat[j]
        y <- M1[j] + rnorm(1, 0, 1)
      }
      
      # regret
      if (ac > (M1[j]>M0[j])) {
        regret <- regret + M0[j] - M1[j]
      } else if (ac < (M1[j]>M0[j])) {
        regret <- regret + M1[j] - M0[j]
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
    regretall[iter,i] <- regret
  }
}

regretnew <- apply(regretall[1:4,], 2, mean)
fig1 <- ggplot() + geom_line(aes(x = Timeall^(2/3), y = regretnew), color = "#377eb8", size = 1) + 
  theme_classic() + ylab("Regret") + xlab(TeX("$T^{2/3}$")) + geom_point(aes(x = Timeall^(2/3), y = regret3new), color = "#377eb8", shape = 17, size = 2.5) + labs(tag = "A")
fig1
ggsave("regret3.pdf", fig1, units = "px", width = 1084, height = 640, device="pdf")

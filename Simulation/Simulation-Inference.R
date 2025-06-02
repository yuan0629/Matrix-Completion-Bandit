library(ggplot2)
library(softImpute)
library(latex2exp)


### Inference ###

Iteration <- 1000
a <- 1/3 # exploration probability for phase two \varepsilon_t=c_2t^{-a}
N <- 300 # number of rows
Ti <- 300 # number of columns
rreal <- 2 # real rank 
Time <- 60000 # Time horizon
T0 <- 20000 # Length of phase one T0 
c2 <- 10 # exploration probability for phase two\varepsilon_t=c_2t^{-a}
r <- 2 # rank used in practice



set.seed(1124)
M1 <- matrix(runif(N*Ti, -100, 100), nrow = N, ncol = Ti)
M0 <- M1 + matrix(runif(N*Ti, -2, 2), nrow = N, ncol = Ti)

U1 <- svd(M1)$u[,1:rreal] %*% diag(c(sqrt(svd(M1)$d)[1:rreal]))   
V1 <- svd(M1)$v[,1:rreal] %*% diag(c(sqrt(svd(M1)$d)[1:rreal]))
U0 <- svd(M0)$u[,1:rreal] %*% diag(c(sqrt(svd(M0)$d)[1:rreal]))   
V0 <- svd(M0)$v[,1:rreal] %*% diag(c(sqrt(svd(M0)$d)[1:rreal]))
M1 <- U1 %*% t(V1)
M0 <- U0 %*% t(V0)
lmax <- max(svd(M1)$d[1], svd(M0)$d[1])



# initial estimate
M10 <- M00 <- matrix(NA, nrow = N, ncol = Ti)

ran <- sample(c(1:(N*Ti)), 50*N, replace = FALSE)
M10[ran] <- M1[ran]
fit1 <- softImpute(M10, rank.max = rreal, lambda = 0.5)
M10 <- complete(M10, fit1)
sum((M1-M10)^2)

ran <- sample(c(1:(N*Ti)), 50*N, replace = FALSE)
M00[ran] <- M0[ran]
fit2 <- softImpute(M00, rank.max = rreal, lambda = 0.5)
M00 <- complete(M00, fit2)
sum((M0-M00)^2)




value31 <- rep(0, Iteration)
value33 <- rep(0, Iteration)
value32 <- value34 <- rep(0, Iteration)

b <- Time/(Time - T0)


for (iter in 1:Iteration) {
  print(iter)

  
  M1hat <- M10
  M0hat <- M00
  
  
  U0 <- svd(M00)$u[,1:r] %*% diag(c(sqrt(svd(M00)$d)[1:r]))   
  V0 <- svd(M00)$v[,1:r] %*% diag(c(sqrt(svd(M00)$d)[1:r]))
  U1 <- svd(M10)$u[,1:r] %*% diag(c(sqrt(svd(M10)$d)[1:r]))   
  V1 <- svd(M10)$v[,1:r] %*% diag(c(sqrt(svd(M10)$d)[1:r]))
  
  
  U <- list(U1, U0)
  V <- list(V1, V0)
  
  sigma1 <- sigma0 <- 0
  S1 <- S0 <- 0
  M1unbs <- matrix(0, nrow = N, ncol = Ti)
  M0unbs <- matrix(0, nrow = N, ncol = Ti)
  
  # algorithm
  set.seed(iter)
  
  for (t in 1:Time){
    j <- sample(c(1:(N*Ti)), 1)
    
    if (t <= T0) {
      epsilon <- 0.6  # exploration probability for phase one
      eta <- 0.025 * N*Ti*log(N)/lmax/Time^(1-a)  # stepsize
    } else {
      epsilon <- c2 * t^(-a)
      eta <- 0.025 * epsilon*N*Ti*log(N)/lmax/Time^(1-a) # stepsize
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
    
    # debiasing
    if (t > T0) {
      
      M1unbs <- M1unbs + M1hat 
      M0unbs <- M0unbs + M0hat
      
  
      if (ac == 1) {
        sigma1 <- sigma1 + 1 / pit * (y - yhat)^2

        # debias
        M1unbs <- M1unbs + N*Ti/pit*(y - yhat) * X
      }else {
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

  # projection
  M1unbs <- M1unbs / (Time - T0)
  L1 <- svd(M1unbs)$u[,1:r]
  R1 <- svd(M1unbs)$v[,1:r]
  M1proj <- L1 %*% t(L1) %*% M1unbs %*% R1 %*% t(R1)
  Lproj1 <- svd(M1proj)$u[,1:r]
  Rproj1 <- svd(M1proj)$v[,1:r]
  Lproj1c <- svd(M1proj)$u[,(r+1):Ti]
  Rproj1c <- svd(M1proj)$v[,(r+1):Ti]
  
  M0unbs <- M0unbs / (Time - T0)
  L0 <- svd(M0unbs)$u[,1:r]
  R0 <- svd(M0unbs)$v[,1:r]
  M0proj <- L0 %*% t(L0) %*% M0unbs %*% R0 %*% t(R0)
  Lproj0 <- svd(M0proj)$u[,1:r]
  Rproj0 <- svd(M0proj)$v[,1:r]
  Lproj0c <- svd(M0proj)$u[,(r+1):Ti]
  Rproj0c <- svd(M0proj)$v[,(r+1):Ti]
  
  
   
  sigma1 <- sigma1 / (Time - T0)
  sigma0 <- sigma0 / (Time - T0)
  
  
  # first W
  W <- matrix(0, nrow = N, ncol = Ti)
  W[w1,w2] <- 1
  #W[w3,w4] <- -1
  
  PMw <- Lproj1 %*% t(Lproj1) %*% W %*% Rproj1 %*% t(Rproj1) + Lproj1 %*% t(Lproj1) %*% W %*% Rproj1c %*% t(Rproj1c) +
    Lproj1c %*% t(Lproj1c) %*% W %*% Rproj1 %*% t(Rproj1)
  S1 <- (1/(1+a)*2/c2*sum((PMw[M1hat < M0hat])^2) + 1/Time^a*sum((PMw[M1hat > M0hat])^2))*b
  
  
  PMw <- Lproj0 %*% t(Lproj0) %*% W %*% Rproj0 %*% t(Rproj0) + Lproj0 %*% t(Lproj0) %*% W %*% Rproj0c %*% t(Rproj0c) +
    Lproj0c %*% t(Lproj0c) %*% W %*% Rproj0 %*% t(Rproj0)
  S0 <- (1/(1+a)*2/c2*sum((PMw[M1hat > M0hat])^2) + 1/Time^a*sum((PMw[M1hat < M0hat])^2))*b
  

  # result
  value31[iter] <- (M1proj[w1,w2] - M1[w1,w2])/sqrt(sigma1*S1*N*Ti/Time^(1-a))
  value32[iter] <- (M0proj[w1,w2] - M0[w1,w2])/sqrt(sigma0*S0*N*Ti/Time^(1-a))
  value33[iter] <- (M0proj[w1,w2] - M0[w1,w2] - M1proj[w1,w2] + M1[w1,w2])/sqrt((sigma1*S1+sigma0*S0)*N*Ti/Time^(1-a))
  
  # second W
  W <- matrix(0, nrow = N, ncol = Ti)
  W[w1,w2] <- 1
  W[w3,w4] <- -1
  
  
  PMw <- Lproj0 %*% t(Lproj0) %*% W %*% Rproj0 %*% t(Rproj0) + Lproj0 %*% t(Lproj0) %*% W %*% Rproj0c %*% t(Rproj0c) +
    Lproj0c %*% t(Lproj0c) %*% W %*% Rproj0 %*% t(Rproj0)
  S0 <- (1/(1+a)*2/c2*sum((PMw[M1hat > M0hat])^2) + 1/Time^a*sum((PMw[M1hat < M0hat])^2))*b
  
  value34[iter] <- (M0proj[w1,w2] - M0proj[w3,w4] - M0[w1,w2] + M0[w3,w4])/sqrt(sigma0*S0*N*Ti/Time^(1-a))
}    

save(value31, file = "density31.RData") # \inp{M_1}{e_1e_5^{\top}}
save(value32, file = "density32.RData") # \inp{M_0}{e_1e_5^{\top}}
save(value33, file = "density33.RData") # \inp{M_0-M_1}{e_1e_5^{\top}}
save(value34, file = "density34.RData") # \inp{M_0}{e_1e_5^{\top} - e_2e_2^{\top}}





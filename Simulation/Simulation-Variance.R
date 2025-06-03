library(ggplot2)
library(softImpute)
library(latex2exp)


### Inference ###

Iteration <- 1
a <- 1/3 # exploration probability for phase two \varepsilon_t=c_2t^{-a}
N <- 300 # number of rows
Ti <- 300 # number of columns
rreal <- 2 # real rank 
Time <- 60000 # Time horizon
T0 <- 20000 # Length of phase one T0


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
fit1 <- softImpute(M10, rank.max = 10, lambda = 0.5)
M10 <- complete(M10, fit1)
sum((M1-M10)^2)

ran <- sample(c(1:(N*Ti)), 50*N, replace = FALSE)
M00[ran] <- M0[ran]
fit2 <- softImpute(M00, rank.max = 10, lambda = 0.5)
M00 <- complete(M00, fit2)
sum((M0-M00)^2)


W <- matrix(0, nrow = N, ncol = Ti)
w1 <- 1
w2 <- 5
W[w1,w2] <- 1
#W[w3,w4] <- -1

sampleall <- matrix(0, nrow = Time - T0, ncol = 6)

b <- Time/(Time - T0)
c2 <- 1 # exploration probability for phase two\varepsilon_t=c_2t^{-a}
r <- 2 # rank used in practice



  
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
  set.seed(629122)
  
  for (t in 1:Time){
    j <- sample(c(1:(N*Ti)), 1)
    
    if (t <= T0) {
      epsilon <- 0.6
      eta <- 0.025 * N*Ti*log(N)/lmax/Time^(1-a)
    } else {
      epsilon <- c2 * t^(-a)
      eta <- 0.025 * epsilon*N*Ti*log(N)/lmax/Time^(1-a)
    }
       
    # decision
    pit <- (1-epsilon/2) * (M1hat[j]> M0hat[j]) + epsilon/2 * (M1hat[j]<= M0hat[j])
    ac <- sample(c(0,1), prob = c(1-pit, pit), 1)
    de <- ac * pit + (1-ac) * (1-pit)
    X <- matrix(0, nrow = N, ncol = Ti)
    X[j] <- 1
    reala <- M1[j] > M0[j]
    if (ac == 0) {
      yhat <- M0hat[j]
      xi <- rnorm(1, 0, 1)
      y <- M0[j] + xi
    }else {
      yhat <- M1hat[j]
      xi <- rnorm(1, 0, 1)
      y <- M1[j] + xi
    }
    
    
    # inference way1
    if (t > T0) {
      
      M1unbs <- M1unbs + M1hat 
      M0unbs <- M0unbs + M0hat
      
      sampleall[t-T0, 2] <- j
      
      
      # way2 & debias
      if (ac == 1) {
        sigma1 <- sigma1 + 1 / pit * (y - yhat)^2

        # debias
        M1unbs <- M1unbs + N*Ti/pit*(y - yhat) * X
        
        sampleall[t-T0, 1] <- 1/pit*(y - yhat)
        sampleall[t-T0, 5] <- 1/pit*xi
        
        if(reala == 1) {
          sampleall[t-T0, 3] <- 1
        } else {
          sampleall[t-T0, 3] <- 2
        }
        
      }else {
        sigma0 <- sigma0 + 1 / (1-pit) * (y - yhat)^2

        # debias
        M0unbs <- M0unbs + N*Ti/(1-pit)*(y - yhat) * X
        
        sampleall[t-T0, 1] <- 1/(1-pit)*(y - yhat)
        sampleall[t-T0, 5] <- 1/(1-pit)*xi
        
        if(reala == 1) {
          sampleall[t-T0, 3] <- 3
        } else {
          sampleall[t-T0, 3] <- 4
        }
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
  #M1proj <- L1 %*% t(L1) %*% M1unbs %*% R1 %*% t(R1)
 
  M0unbs <- M0unbs / (Time - T0)
  L0 <- svd(M0unbs)$u[,1:r]
  R0 <- svd(M0unbs)$v[,1:r]
  #M0proj <- L0 %*% t(L0) %*% M0unbs %*% R0 %*% t(R0)
 
for (i in 1:nrow(sampleall)) {
  
  X <- matrix(0, nrow = N, ncol = Ti)
  X[sampleall[i, 2]] <- 1
  
  if ((sampleall[i, 3] == 1) | (sampleall[i, 3] == 2)) {
    sampleall[i, 4] <-  sampleall[i, 1] * (L1 %*% t(L1) %*% X %*% R1 %*% t(R1))[w1,w2]
    sampleall[i, 6] <-  sampleall[i, 5] * (L1 %*% t(L1) %*% X %*% R1 %*% t(R1))[w1,w2]
  } else {
    sampleall[i, 4] <-  sampleall[i, 1] * (L0 %*% t(L0) %*% X %*% R0 %*% t(R0))[w1,w2]
    sampleall[i, 6] <-  sampleall[i, 5] * (L0 %*% t(L0) %*% X %*% R0 %*% t(R0))[w1,w2]
  }
}  



## Figure 1 left
df <- data.frame(value = sampleall[(sampleall[, 3] == 1) | (sampleall[, 3] == 2),6], class = as.factor(sampleall[(sampleall[, 3] == 1) | (sampleall[, 3] == 2),3]))  
#ggplot(data = df) + geom_boxplot(aes(x = class, y = value, color = class)) + scale_color_manual(values = c("#386cb0", "#fdb462"), labels = c(TeX(r"($\Omega_1$)"), "sample0")) + theme_classic()
fig1 <- ggplot(data = df) + geom_boxplot(aes(x = class, y = value, color = class)) + 
  scale_color_manual(values = c("#386cb0", "#fdb462"), labels = c(TeX(r"( $\Omega_1$)"), TeX(r"( $\Omega_0$)"))) + theme_classic() + xlab("") + ylab("") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + labs(color = "optimal set") 
fig1
ggsave(plot = fig1, filename = "variance1.pdf", units = "px", width = 1163, height = 831, device="pdf")


## Figure 1 right
df <- data.frame(value = sampleall[(sampleall[, 3] == 3) | (sampleall[, 3] == 4),6], class = as.factor(sampleall[(sampleall[, 3] == 3) | (sampleall[, 3] == 4),3]))  
#ggplot(data = df) + geom_boxplot(aes(x = class, y = value, color = class)) + scale_color_manual(values = c("#386cb0", "#fdb462"), labels = c(TeX(r"($\Omega_1$)"), TeX(r"($\Omega_0$)"))) + theme_classic()
fig1 <- ggplot(data = df) + geom_boxplot(aes(x = class, y = value, color = class)) + 
  scale_color_manual(values = c("#386cb0", "#fdb462"), labels = c(TeX(r"( $\Omega_1$)"), TeX(r"( $\Omega_0$)"))) + theme_classic()  + ylab("") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + labs(color = "optimal set")
fig1
ggsave(plot = fig1, filename = "variance0.pdf", units = "px", width = 1163, height = 831, device="pdf")














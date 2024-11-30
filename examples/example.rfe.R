library(lavaan)
data("HolzingerSwineford1939")
X <- scale(as.matrix(HolzingerSwineford1939[,7:15]))
Y <- matrix(0, nrow(X), 4)
Y[,1] <- ifelse(HolzingerSwineford1939$sex == 1, 1, 0) #sex
Y[,2] <- ifelse(HolzingerSwineford1939$sex == 2, 1, 0)
Y[,3] <- ifelse(HolzingerSwineford1939$school == "Pasteur", 1, 0)
Y[,4] <- ifelse(HolzingerSwineford1939$school == "Grant-White", 1, 0)
aics <- c()
lambdas <- 10^(seq(-3, 1, by = 0.5))
res.list <- list()
idx <- 1
n.starts <- 10
for(lambda in lambdas){
  res.list.tmp <- list()
  f.list.tmp <- c()
  for(i in 1:n.starts){
    res.list.tmp[[i]] <- mdfa.rfe(X, Y, 3, lambda)
    f.list.tmp[i] <- min(res.list.tmp[[i]]$history)
  }
  res.list[[idx]] <- res.list.tmp[[which.min(f.list.tmp)]]
  aics[idx] <- res.list[[idx]]$AIC
  idx <- idx + 1
}
res.list[[which.min(aics)]]

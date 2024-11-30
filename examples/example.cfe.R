X <- read.csv("jobImage.csv")
X <- scale(as.matrix(X[,-1])) #Useful, Good, Firm, Quick, Noisy, Busy
X <- X[,c(2,3,10,7,8,12)]
aics <- c()
lambdas <- 10^(seq(-3, 1, by = 0.5))
res.list <- list()
idx <- 1
n.starts <- 3
for(lambda in lambdas){
  res.list.tmp <- list()
  f.list.tmp <- c()
  for(i in 1:n.starts){
    res.list.tmp[[i]] <- mdfa.cfe(X, 2, 3, lambda)
    f.list.tmp[i] <- min(res.list.tmp[[i]]$history)
  }
  res.list[[idx]] <- res.list.tmp[[which.min(f.list.tmp)]]
  aics[idx] <- res.list[[idx]]$AIC
  idx <- idx + 1
}
plot(res.list[[which.min(aics)]]$Fr)
res.list[[which.min(aics)]]$M

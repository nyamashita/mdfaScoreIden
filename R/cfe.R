#' Common factor score identification in matrix decomposition factor analysis (MDFA)
#'   ; Clustering version
#'
#' Specifies common/unique factor scores in MDFA by solving the following minimization problem
#' min. ||X - (FA' + UD)||^2 + lambda||MC - F||^2 over F, A, U, D, M, and C.\cr
#' \cr
#'  n: number of observations\cr
#'  p: number of observed variables\cr
#'  r: number of common factors\cr
#'  k: number of clusters\cr
#'  X: n * p data matrix\cr
#'  F: n * r common factor score matrix\cr
#'  U: n * p unique factor score matrix\cr
#'  A: p * r factor loading matrix\cr
#'  D: diagonal matrix of square root of uniqueness\cr
#'  M: n * k membership matrix\cr
#'  C: k * p cluster centroid matrix\cr
#'  lambda : tuning parameter\cr
#' @param X A matrix. n * p data matrix
#' @param r An integer. Number of common factors.
#' @param k An integer. Number of clusters.
#' @param lambda A numeric value. Tuning parameter > 0
#' @param eps A numeric value. Convergence criterion.
#' @param itemax An integer. Maximum number of iteration

#' @return List of the following objects.\cr
#' F: A matrix. Common factor score.\cr
#' Fr: A matrix. Common factor score; rotated.\cr
#' U: A matrix. Unique factor score.\cr
#' A: A matrix. Factor loadings.\cr
#' Ar: A matrix. Factor loadings; rotated.\cr
#' C: A matrix. Cluster centroid.\cr
#' Cr: A matrix. Cluster centroid; rotated.\cr
#' M: A matrix. Membership.\cr
#' D: A matrix. Square root of uniqueness.\cr
#' int.factor.corr: A matrix. Inter-factor correlation.\cr
#' fitted.cov: A matrix. Covariance matrix fitted to X.\cr
#' history: A numeric vector. Iteration history.\cr
#' AIC: A numeric value. AIC value.\cr
#' @examples
#' X <- read.csv("jobImage.csv")
#' X <- scale(as.matrix(X[,-1])) #Useful, Good, Firm, Quick, Noisy, Busy
#' X <- X[,c(2,3,10,7,8,12)]
#' aics <- c()
#' lambdas <- 10^(seq(-3, 1, by = 0.5))
#' res.list <- list()
#' idx <- 1
#' n.starts <- 3
#' for(lambda in lambdas){
#'   res.list.tmp <- list()
#'   f.list.tmp <- c()
#'   for(i in 1:n.starts){
#'     res.list.tmp[[i]] <- mdfa.cfe(X, 2, 3, lambda)
#'     f.list.tmp[i] <- min(res.list.tmp[[i]]$history)
#'   }
#'   res.list[[idx]] <- res.list.tmp[[which.min(f.list.tmp)]]
#'   aics[idx] <- res.list[[idx]]$AIC
#'   idx <- idx + 1
#' }
#' plot(res.list[[which.min(aics)]]$Fr)
#' res.list[[which.min(aics)]]$M
#' @author Naoto Yamashita, \email{nyam@kansai-u.ac.jp}
#' @export
mdfa.cfe <- function(X, #data matrix
                     r, #number of factors
                     k, #number of clusters
                     lambda, #tuning parameter
                     eps = 1e-6,
                     itemax = 2000
)
{
  #preparation
  n <- nrow(X)
  p <- ncol(X)
  #loss function
  lf <- function(Z, B, M, C){
    sum((X - Z %*% t(B))^2) + lambda*sum((M %*% C - Z[,1:r])^2)
  }
  #iteration
  iter.flg <- 0
  while(iter.flg == 0){
    #initialization
    history <- c()
    B <- cbind(matrix(runif(p*r, -1, 1), p, r), diag(runif(p)))
    M <- matrix(0, n, k)
    for(row in 1:n){
      M[row,sample(sample(1:3, 1))] <- 1
    }
    C <- matrix(rnorm(k*r), k, r)
    for(ite in 1:itemax){
      #Z-step
      Wtild <- as.matrix(Matrix::bdiag(sqrt(lambda)*diag(r), matrix(0, p, p)))
      Ct <- cbind(t(B), Wtild)
      Xtild <- cbind(X, sqrt(lambda)*M %*% C, matrix(0, n, p))
      res.svd <- svd(Xtild %*% t(Ct) / sqrt(n))
      if(min(res.svd$d) < 1e-8){
        warning("Factor scores are not uniquely determined due to the rank of Xtild.")
        break
      }
      Z <- sqrt(n)*res.svd$u %*% t(res.svd$v)
      #Bstep
      B <- t(X) %*% Z / n
      B[,(r+1):(p+r)] <- diag(diag(B[,(r+1):(p+r)]))
      #M-step
      Mnew <- M
      for(row in 1:n){
        fvec <- c()
        for(col in 1:k){
          Mnew[row,] <- 0
          Mnew[row, col] <- 1
          fvec[col] <- lf(Z, B, Mnew, C)
        }
        M[row,] <- 0
        M[row, which.min(fvec)] <- 1
      }
      #C-step
      C <- MASS::ginv(t(M) %*% M) %*% t(M) %*% Z[,1:r]
      #check convergence
      history[ite] <- lf(Z, B, M, C)
      if(ite > 1){
        if((history[ite-1] - history[ite]) < eps){
          iter.flg <- 1
          break
        }
      }
    }
  }
  fitted.cov = B %*% t(B)
  #rotation
  F <- Z[,1:r]
  A <- B[,1:r]
  res.rot <- GPArotation::geominQ(A)
  T <- res.rot$Th
  Fr <- F %*% T
  Ar <- A %*% t(solve(T))
  Cr <- C %*% T
  int.fac.cor <- t(T) %*% T
  #compute AIC
  Sx <- t(X)%*%X / n
  neglogL <- -log(det(solve(fitted.cov))) - log(det(Sx)) + sum(diag(Sx %*% solve(fitted.cov)))
  AIC <- 2*neglogL + 2*(p*r + p*p)
  #return values
  list(
    F = Z[,1:r],                #common factor score
    Fr = Fr,                    #common score; rotated
    U = Z[,(r+1):(p+r)],        #unique factor score
    A = B[,1:r],                #loadings
    Ar = Ar,                    #loadings; rotated
    C = C,                      #cluster-centroid
    Cr = Cr,                    #cluster-centroid; rotated
    M = M,                      #membership
    D = B[(r+1):(p+r)],         #uniqueness
    int.fac.cor = int.fac.cor,  #inter-factor correlation
    fitted.cov = fitted.cov,    #covariance matrix fitted to X
    history = history,          #iteration history
    AIC = AIC                   #AIC value
  )
}

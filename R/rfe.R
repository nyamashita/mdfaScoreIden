#' Common factor score identification in matrix decomposition factor analysis (MDFA)
#'   ; Regression version
#'
#' min. ||X - (FA' + UD)||^2 + lambda||Y - FW||^2 over F, A, U, D, and W.\cr
#' \cr
#'  n: number of observations\cr
#'  p: number of observed variables\cr
#'  q: number of external criterion\cr
#'  r: number of common factors\cr
#'  X: n * p data matrix\cr
#'  F: n * r common factor score matrix\cr
#'  U: n * p unique factor score matrix\cr
#'  A: p * r factor loading matrix\cr
#'  D: diagonal matrix of square root of uniqueness\cr
#'  W: r * q regression coefficient matrix\cr
#'  lambda : tuning parameter\cr
#' @param X A matrix. n * p data matrix
#' @param Y A matrix. n * q data matrix of external criterion
#' @param r An integer. Number of common factors.
#' @param lambda A numeric value. Tuning parameter > 0
#' @param eps A numeric value. Convergence criterion.
#' @param itemax An integer. Maximum number of iteration

#' @return List of the following objects.\cr
#' F: A matrix. Common factor score.\cr
#' Fr: A matrix. Common factor score; rotated.\cr
#' U: A matrix. Unique factor score.\cr
#' A: A matrix. Factor loadings.\cr
#' Ar: A matrix. Factor loadings; rotated.\cr
#' W: A matrix. Regression coefficients\cr
#' D: A matrix. Square root of uniqueness.\cr
#' int.factor.corr: A matrix. Inter-factor correlation.\cr
#' fitted.cov: A matrix. Covariance matrix fitted to X.\cr
#' history: A numeric vector. Iteration history.\cr
#' AIC: A numeric value. AIC value.\cr
#' @examples
#' library(lavaan)
# 'data("HolzingerSwineford1939")
# 'X <- scale(as.matrix(HolzingerSwineford1939[,7:15]))
# 'Y <- matrix(0, nrow(X), 4)
# 'Y[,1] <- ifelse(HolzingerSwineford1939$sex == 1, 1, 0) #sex
# 'Y[,2] <- ifelse(HolzingerSwineford1939$sex == 2, 1, 0)
# 'Y[,3] <- ifelse(HolzingerSwineford1939$school == "Pasteur", 1, 0)
# 'Y[,4] <- ifelse(HolzingerSwineford1939$school == "Grant-White", 1, 0)
#' aics <- c()
#' lambdas <- 10^(seq(-3, 1, by = 0.5))
#' res.list <- list()
#' idx <- 1
#' n.starts <- 10
#' for(lambda in lambdas){
#'   res.list.tmp <- list()
#'   f.list.tmp <- c()
#'   for(i in 1:n.starts){
#'     res.list.tmp[[i]] <- mdfa.rfe(X, Y, 3, lambda)
#'     f.list.tmp[i] <- min(res.list.tmp[[i]]$history)
#'   }
#'   res.list[[idx]] <- res.list.tmp[[which.min(f.list.tmp)]]
#'   aics[idx] <- res.list[[idx]]$AIC
#'   idx <- idx + 1
#' }
#' res.list[[which.min(aics)]]

#' @export
mdfa.rfe <- function(X, #data matrix
                     Y, #matrix of external criterion
                     r, #number of factors
                     lambda, #tuning parameter
                     eps = 1e-6,
                     itemax = 2000
                     )
{
  #preparation
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  Xtild <- cbind(X, sqrt(lambda)*Y, matrix(0, n, p))
  if(q < r){
    warning("Number of dependent var.s is fewer than the factors. Factor scores are not uniquely determined.")
  }
  #loss function
  lf <- function(Z, B, W){
    sum((X - Z %*% t(B))^2) + lambda*sum((Y - Z[,1:r] %*% W)^2)
  }
  #iteration
  iter.flg <- 0
  while(iter.flg == 0){
    #initialization
    history <- c()
    B <- cbind(matrix(runif(p*r, -1, 1), p, r), diag(runif(p)))
    W <- matrix(runif(r*q, -1, 1), r, q)
    for(ite in 1:itemax){
      #Z-step
      Wtild <- as.matrix(Matrix::bdiag(sqrt(lambda)*W, matrix(0, p, p)))
      C <- cbind(t(B), Wtild)
      res.svd <- svd(Xtild %*% t(C) / sqrt(n))
      if(min(res.svd$d) < 1e-8){
        warning("Factor scores are not uniquely determined due to the rank of Xtild.")
        break
      }
      Z <- sqrt(n)*res.svd$u %*% t(res.svd$v)
      #Bstep
      B <- t(X) %*% Z / n
      B[,(r+1):(p+r)] <- diag(diag(B[,(r+1):(p+r)]))
      #W-step
      W <- t(Z[,1:r]) %*% Y / n
      #check convergence
      history[ite] <- lf(Z, B, W)
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
  U <- Z[,(r+1):(p+r)]
  A <- B[,1:r]
  res.rot <- GPArotation::geominQ(A)
  T <- res.rot$Th
  Fr <- F %*% T
  Ar <- A %*% t(solve(T))
  int.factor.corr <- t(T) %*% T
  Wr <- solve(T) %*% W
  #compute AIC
  Sx <- t(X)%*%X / n
  neglogL <- -log(det(solve(fitted.cov))) - log(det(Sx)) + sum(diag(Sx %*% solve(fitted.cov)))
  AIC <- 2*neglogL + 2*(p*r + p*p)
  #return values
  list(
    F = Z[,1:r],                            #common factor score
    Fr = Fr,                                #common factor score; rotated
    U = Z[,(r+1):(p+r)],                    #unique factor score
    A = B[,1:r],                            #loadings
    Ar = Ar,                                #loadings; rotated
    D = B[,(r+1):(p+r)],                    #uniqueness
    W = W,                                  #regression coefficient
    Wr = Wr,                                #regression coefficient; rotated
    int.factor.corr = int.factor.corr,      #inter factor correlation
    fitted.cov = fitted.cov,                #covariance matrix fitted to X
    history = history,                      #iteration history
    AIC = AIC                               #AIC value
  )
}

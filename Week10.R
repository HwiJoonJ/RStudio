x1 <- c(1:10)
x2 <- c(3,5,4,6,5,7,6,8,7,9)
x3 <- c(4,5,6,1,2,3,7,8,9,4)

x <- cbind(x1, x2, x3)
y <- c(2,6,4,7,3,5,4,4,8,9)

cor(x)

VIFs <- function(ind1, ind2, dep) {
  reg <- lm(dep ~ ind1 + ind2)
  Rsq <- summary(reg)$r.squared
  1/(1-Rsq)
}

VIF1 = VIFs(x2, x3, x1)
VIF2 = VIFs(x1, x3, x2)
VIF3 = VIFs(x1, x2, x3)

#Data Normalization
x_sc <- scale(x)
#Eigen Value
eig <- eigen(cov(x_sc))
#Sort Eigen Value
idx <- order(eig$values, decreasing=T)
eigvec <- eig$vectors[,idx]
eigvec
#convolution eigen vector to PCA
x_pca <- x_sc %*% eigvec
x_pca

y_sc <- scale(y)
y_init <- y_sc
x_init <- x_sc
tb <- tp <- 0

pls <- matrix(ncol=ncol(x_sc), nrow=nrow(x_sc))
for(i in 1:ncol(x_sc)) {
  #초기값
  y_init <- y_init - tb
  x_init <- x_init - tp
  
  #공분산을 최대화하는 선형조합 도출
  a <- t(x_init) %*% y_init
  a <- a/sqrt(sum(a^2))
  
  #선형조합 기반 잠재변수 도출
  pls[,i] <- x_init %*% a #pls -> 자료 14p t에 대한 식.
  
  #선형조합에 대한 회귀계수 산출
  p <- pls[,i] %*% x_init/c(pls[,i] %*% pls[,i])
  b <- pls[,i] %*% y_init/c(pls[,i] %*% pls[,i])
  
  tb <- pls[,i] %*% b
  tp <- pls[,i] %*% p
}




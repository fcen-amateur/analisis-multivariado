# Calcula grilla de valores para forma cuadratica
# provee z para que `contour` dibuje elipses
ellipse_grid <- function(rango_x1, rango_x2, A, b) {
  cuadratica <- function(x, b, A) { as.numeric(t(x - b) %*% A %*% (x - b)) }
  z <- matrix(nrow = length(rango_x1), ncol = length(rango_x2))
  for (i in seq_along(rango_x1)) {
    for (j in seq_along(rango_x2)) {
      x <- c(rango_x1[i], rango_x2[j])
      z[i, j] <- cuadratica(x, b, A)
    }
  }
  return(z)
}

# Test de hotelling univariado
hotelling.test <- function(X, C = diag(dim(X)[2]), mu0 = rep(0, dim(C)[1]),
                           int.conf = NULL, int.A = C, int.alfa = 0.05) {
  n <- dim(X)[1]
  p <- dim(C)[1]
  X <- as.matrix(X)
  x_ <- colMeans(X)
  S <- cov(X)
  y <- C %*% x_ - mu0
  W <- C %*% S %*% t(C)
  T0 <- n * as.numeric(t(y) %*% solve(W) %*% y)
  T0f <- T0 * (n-p) / ((n-1) * p)
  pval <- 1 - pf(T0f, p, n-1)
  
  retval <- list(
    statistic = T0,
    p.value = pval,
    fstat = T0f,
    d = p,
    m = n - 1
  )
  
  if (!is.null(int.conf)) {
    A <- int.A
    q <- qr(A)$rank # cant. ICs LI
    r <- dim(A)[1] # cant. ICs total
    
    if (int.conf=="bonferroni") {
      k <- qt(1 - int.alfa/(2*r), n - 1)
    } else if (int.conf=="scheffe") {
      k <- sqrt(q * (n - 1)/(n - q) *qf(1 - int.alfa, q, n - q))
    } else { stop("Tipo IC no reconocido") }
    ICs <- matrix(nrow = r, ncol = 2)
    for (i in seq.int(r)) {
      sigma_ax <- sqrt( as.numeric(t(A[i,]) %*% S %*% A[i,]) / n )
      delta <- k * sigma_ax
      centro <- as.numeric(t(A[i,]) %*% x_ )
      ICs[i,] <- c(centro, delta)
    }
    retval$intervals <- ICs
    retval$int.conf <- int.conf
    retval$int.params <- c("alfa" = int.alfa, "k" = k, "q" = q, "r" = r)
  }
  return(retval)
}

# Plot de distancias robustas para chequear outliers
check_robdist <- function(X, tol = 2) {
  X <- as.matrix(X)
  n <- dim(X)[1]
  p <- dim(X)[2]
  covR <- covRob(X)
  qq <- qqplot(qchisq(ppoints(n), df = p), covR$dist)
  qqline(covR$dist, distribution = function(x) qchisq(x, df = p))
  outs <- abs(qq$y - qq$x) > tol
  text(qq$x[outs] -0.25, qq$y[outs], labels = seq_along(qq$x)[outs])
}

# Devuelve vector de norma 1 con igual direccion y sentido que x
normalizar <- function(x) {x / sqrt(sum(x^2))}

# Calcula la correlacion muestral a partir de la covarianza muestral
corr_matrix <- function(S) {
  rows <- nrow(S)
  cols <- ncol(S)
  R <- matrix(nrow = rows, ncol = cols)
  for (i in seq.int(rows)) {
    for (j in seq.int(cols)) {
      R[i, j] <- S[i,j] / sqrt(S[i,i]*S[j,j])
    }
  }
  return(R)
}

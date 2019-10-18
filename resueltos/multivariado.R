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

hotelling.test <- function(
  # X: muestra x_i de n VA IID Normal(mu, Sigma)
  # C, mu0: definen la hipotesis nula a testear, H_0 : C %*% mu = mu0
  # int.conf: Cuando no es NULL, genera IC de forma a'x_ +- k * sqrt(a'Sa/n)
  #   para cada  fila `a` de la matriz de CLs A, con prob conjunta 1 - int.alfa
  #   `k` es funcion de alfa y el tipo de IC, según Seber[84] pp. 81-83
  X, C = diag(dim(X)[2]), mu0 = rep(0, dim(C)[1]),
  int.conf = c(NA, "bonferroni", "scheffe"), int.A = C, int.alfa = 0.05) {

  # Calculo las VA necesarias para el test de Hotelling: y, W

  n <- dim(X)[1] # Tamaño muestral
  p <- dim(C)[1] # C. lineales de mu a testear
  X <- as.matrix(X)
  # "x raya" ~ N(mu, Sigma/n), promedio observaciones (est. insesgado de mu)
  x_ <- colMeans(X)
  # S es est. insesgado de Sigma
  S <- cov(X)

  # y ~ N(0, (C' Sigma C)/n) 
  # W = cov(Y) = C %*% cov(X) %*% t(C) ~ Wishart(Sigma, n-1, p)
  # (y, W) independientes
  y <- C %*% x_ - mu0
  W <- C %*% S %*% t(C)

  # "T0 cuadrado" ~ Hotelling(p, n-1) 
  T0sq <- n * as.numeric(t(y) %*% solve(W) %*% y) 
  # F0 ~ F(p, n - p)
  F0 <- T0sq * (n-p) / ((n-1) * p)
  pval <- 1 - pf(F0, p, n - p)

  retval <- list(
    statistic = T0sq,
    p.value = pval,
    fstat = F0,
    d = p,
    m = n - 1
  )

  int.conf <- match.arg(int.conf)
  if (!is.na(int.conf)) {
    A <- int.A
    q <- qr(A)$rank # cant. ICs LI
    r <- dim(A)[1] # cant. ICs total

    if (int.conf=="bonferroni") {
      k <- qt(1 - int.alfa/(2*r), n - 1)
    } else if (int.conf=="scheffe") {
      k <- sqrt(q * (n - 1)/(n - q) *qf(1 - int.alfa, q, n - q))
    }
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

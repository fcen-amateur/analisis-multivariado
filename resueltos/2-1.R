library(readr)
library(dplyr)
library(RobStatTM)
library(mvShapiroTest)

source("resueltos/multivariado.R")
# 2.1.1
# Usa multivariado.ellipse_grip
A1 <- matrix(c(2, 0, 0, 3), nrow = 2)
b1 <- c(0, 0)
c1 <- c(1, 1.2, 1.5)
rango_x1 <- seq(b1[1] - 2, b1[1] + 2, 0.1)
rango_x2 <- seq(b1[2] - 1, b1[2] + 1, 0.1)

z1 <- ellipse_grid(rango_x1, rango_x2, A1, b1)
contour(rango_x1, rango_x2, z1, levels = c1^2, asp = 1)

A2 <- matrix(c(2.5, -.5, -.5, 2.5), nrow = 2)
b2 <- c(1, 3)
c2 <- c(2, 2.1, 2.2)
rango_x1 <- seq(b2[1] - 1.5, b2[1] + 1.5, 0.1)
rango_x2 <- seq(b2[2] - 1.5, b2[2] + 1.5, 0.1)

z2 <- ellipse_grid(rango_x1, rango_x2, A2, b2)
contour(rango_x1, rango_x2, z2, levels = c2^2, asp = 1)

# 2.1.2
tipo_cols1 <- cols(
  Peso = col_double(),
  Altura = col_double()
)
df212 <- read_tsv("datos/ej2pra2-2017.txt", col_types = tipo_cols1) %>%
  mutate(Altura = Altura / 1000)

check_robdist(df212, tol = 3)
outs <- c(39)
mu0 <- c(63.64, 1.61538)
hotelling.test(df212, mu0 = mu0) # con outlier
hotelling.test(df212[-outs,], mu0 = mu0) # sin outlier

rango_x1 <- seq(40, 80, 1)
rango_x2 <- seq(1.4, 2, 0.025)
covR <- covRob(df212[-outs,])
z0 <- ellipse_grid(
  rango_x1,
  rango_x2,
  A = solve(covR$cov),
  b = covR$center)

plot(df212)

contour(
  rango_x1,
  rango_x2,
  z0,
  levels = qchisq(0.95, 2),
  add = TRUE)
text(x = df212$Peso[outs]-1,
     y = df212$Altura[outs],
     labels = seq_along(df212$Peso)[outs])

mvShapiro.Test(as.matrix(df212))$p.value
# .05, rechazo normalidad
mvShapiro.Test(as.matrix(df212)[-outs,])$p.value
# .30, no rechazo normalidad

# 2.1.8
tipo_cols2 <- cols(
  N = col_integer(),
  E = col_integer(),
  S = col_integer(),
  O = col_integer()
)
df218 <- read_csv("datos/2-1-8.csv", col_types = tipo_cols2)

# Testear normalidad multivariada
mvShapiro.Test(as.matrix(df218)) # p-valor 0.03, rechazo

check_robdist(df218, tol = 2)
mvoutlier::pcout(df218, makeplot = T)
outs <- mvoutlier::pcout(df218)$wfinal01 == 0

mvShapiro.Test(as.matrix(df218[!outs,])) # p-valor 0.34, NO rechazo normalidad
# este criterio es medio violento, sigue encontrando outliers...
mvoutlier::pcout(df218[!outs,], makeplot = T) 
check_robdist(df218[!outs,], tol = 1)

# Me quedo sin los outliers de pcout en la priemra vuelta y continuo
# Testear que mu = mu0
mu0 <- c(45, 42, 45, 42)
df218_ <- df218[!outs,]
colMeans(df218_)
hotelling.test(df218_, mu0 = colMeans(df218_)) # p-valor de 1, es coherente
hotelling.test(df218_, mu0 = mu0) # p-valor de 0.34, no rechazo

# Testear que todos los mu_i son iguales
C1 <- rbind(
  c(1, -1, 0, 0),
  c(0, 1, -1, 0),
  c(0, 0, 1, -1))
hotelling.test(df218_, C1) # p-valor de 0.0007, rechazo todas iguales
hotelling.test(df218, C1) # p-valro de 0.002, rechazo todas iguales

# contraste N-S, E-O
C2 <- rbind(
  c(1, 0, -1, 0),
  c(0, 1, 0, -1))
hotelling.test(df218_, C2)$p.value # p-v 0.95, no rechazo, razonable!

# Ambos sets de ICs contienen al 0, pero los de bonferroni son más ajustados
hotelling.test(df218_, C2, int.conf = "bonferroni", int.alfa = 0.05)$intervals
hotelling.test(df218_, C2, int.conf = "scheffe", int.alfa = 0.05)$intervals

# Contraste de a pares
# Se pueden construir ICs para todo contraste tal que los c_i sumen 0,
# los construimos para las comparacioenes de medias 1:1
C3 <- rbind(
  c(1, -1, 0, 0),
  c(1, 0, -1, 0),
  c(1, 0, 0, -1),
  c(0, 1, -1, 0),
  c(0, 1, 0, -1),
  c(0, 0, 1, -1))

hotelling.test(df218_, int.A = C3, int.conf = "bonferroni", int.alfa = 0.05)$intervals
hotelling.test(df218_, int.A = C3, int.conf = "scheffe", int.alfa = 0.05)$intervals
# N-O y S-O no incluyen al cero, y son distintas e/si. N-E y S-E, casi
hotelling.test(df218_, int.A = C3, int.conf = "bonferroni", int.alfa = 0.1)$intervals

# Para consistencia, repito comparaciones de Seber p. 85
## Bonferroni para mu_1, ..., mu_4, alfa = 0.05
hotelling.test(df218, int.A = diag(4), int.conf = "bonferroni", int.alfa = 0.05)$intervals
## Scheffe para difs e/ medias a 0.01
hotelling.test(df218, int.A = C3, int.conf = "scheffe", int.alfa = 0.01)$intervals


# 2.1.9
# Necesito version alternativa de la funcón
hotelling.test.alt <- function(x_, S, n, C = diag(length(x_)), mu0 = rep(0, dim(C)[1]),
                               int.conf = NULL, int.A = C, int.alfa = 0.05) {
  p = dim(C)[1]
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
n <- 20
mu0 <- c(12, 4, 2)
x_ <- c(11.5, 4.3, 1.2)
S <- rbind(
  c(10.53, 4.21, -5.26),
  c(4.21, 12.63, -3.16),
  c(-5.26, -3.16, 4.21))

solve(S)
K1 <- matrix(c(1, 0, 0), 1, 3)
vars <- c("p.value", "intervals") # El resto no me interesa
# Testeo de a un canonico por vez. Con q = r = 1, Bonferroni y Scheffe deben dar el mismo intervalo
hotelling.test.alt(x_, S, n, C = K1, mu0 = mu0[1], int.conf = "scheffe", int.alfa = 0.05)[vars]
hotelling.test.alt(x_, S, n, C = K1, mu0 = mu0[1], int.conf = "bonferroni", int.alfa = 0.05)[vars]
# p-valor 0.499, no rechazo. Se ve que los ICs son identicos y contiene a mu_{0,1}=12

# Repito para los otros dos elementos. En todos los casos, el p-valor no es suf
# para rechazar H_0i, y el IC 0.95 incluye al mu_{0,i}
hotelling.test.alt(x_, S, n, C = matrix(c(0, 1, 0), 1, 3), mu0 = mu0[2],
                   int.conf = "scheffe", int.alfa = 0.05)[vars]
hotelling.test.alt(x_, S, n, C = matrix(c(0, 0, 1), 1, 3), mu0 = mu0[3],
                   int.conf = "scheffe", int.alfa = 0.05)[vars]

# Cuando testeo conjuntamente, p-valor es 0.0176: rechazo globalmente,
# pero todos los ICs incluyen a mu_0. La direccion de rechazo no es un canonico
hotelling.test.alt(x_, S, n, mu0 = mu0, int.conf = "scheffe", int.alfa = 0.05)[vars]
hotelling.test.alt(x_, S, n, mu0 = mu0, int.conf = "bonferroni", int.alfa = 0.05)[vars]
hotelling.test.alt(x_, S, n, mu0 = mu0, int.conf = "bonferroni", int.alfa = 0.05*3)[vars]
rech_dir <- t(solve(S) %*% x_)
rd <- normalizar(rech_dir)
rech_mu0 <- rd %*% mu0
rd_plus2 <- rbind(rech_dir, c(0, 1, 0), c(0, 0, 1))
hotelling.test.alt(x_, S, n, C = rd_plus2, mu0 = rd_plus2 %*% mu0,
                   int.conf = "scheffe", int.alfa = 0.05)

R <- corr_matrix(S)
x_ - mu0
# De R, veo que Cuando x_1 está por debajo  la media, 
# - x_2 tambien deberia estar por debajo de la media ("menos"), y
# - x_3 debería estar por encima de la media
# PERO x_ tiene x_1 por debajo, x_2 por encima, y x_3 por debajo.
# PARA ESTA MATRIZ DE CORR, x_ está bastante lejos de mu0

# Ejercicio 1.b
Sigma <- rbind(
  c(1, 2,   4,     1),
  c(2, 5,   8.5,   4),
  c(4, 8.5, 17.25, 5),
  c(1, 4,   5,     6)
)
L <- t(chol(Sigma))
A <- solve(L, tol = 1e-20)

# Ejercicio 2
C1 <- cbind(diag(10), -diag(10))
I5 <- diag(5)
O5 <- diag(0, 5, 5)
C2a <- cbind(I5, -I5, O5, O5)
C2b <- cbind(O5, O5, I5, -I5)
C2 <- rbind(C2a, C2b)
C3 <- rbind(C1, C2a)
qr(C1)$rank == 10
qr(C2)$rank == 10
qr(C3)$rank == 15

# Ejercicio 3
library(readr)
library(tidyverse)
df3 <- read_delim("../Multivariado2/tortugas-3var.txt", delim = " ", col_names = F)
df3 <- list(
  "M" = df3[1:3],
  "F" = df3[4:6]
) %>%
  map(setNames, c("X1", "X2", "X3")) %>%
  enframe("sexo", "muestra") %>%
  mutate(muestra = map(muestra, as.matrix))


Q <- function(X) {
  X <- as.matrix(X)
  # Xrulo: X centrada
  Xrulo <- sweep(X, 2, colMeans(X), "-")
  return(t(Xrulo) %*% Xrulo) 
}
term_prod <- function(X, n) {det(X/n)^(n/2)}
gl <- 6
df3 <- df3 %>%
  mutate(
    ni = map(muestra, dim) %>% map_int(1),
    Qi = map(muestra, Q),
    muhat = map(muestra, colMeans))
    
h1_test <- df3 %>%
  mutate(term = map2_dbl(Qi, ni, term_prod)) %>%
  summarise(
    num = prod(term),
    denom = det(Reduce("+", Qi)/sum(ni))^(sum(ni)/2),
    gama1 = num/denom,
    logstat = -2*log(gama1),
    pval = pchisq(logstat, df = gl, lower.tail = F)
  )

h1_test

# P-VALOR 0.609, NO rechazo igualdad de sigmas
muestras <- deframe(df3[c("sexo", "muestra")])
Qs <- deframe(df3[c("sexo", "Qi")])
ns <-deframe(df3[c("sexo", "ni")])
muhats <- deframe(df3[c("sexo", "muhat")])

difhat <- muhats[["M"]] - muhats[["F"]]
S <- (Qs[["M"]] + Qs[["F"]]) / (n[["M"]] + n[["F"]] - 2)

lambda <- (ns[["M"]] * ns["F"])/(ns[["M"]] + ns[["F"]])
T0sq <- lambda * as.numeric(t(difhat) %*% solve(S) %*% difhat)
F0 <- T0sq/46 * 44/3
pval <- pf(F0, 3, 44, lower.tail = F)

# d: supuestos
library(mvShapiroTest)
map(muestras, mvShapiro.Test)


# Ejercicio 4
# a
df4 <- read_delim("../Multivariado2/craneos-egipcios-3grupos.txt", delim = " ", col_names = F)
df4 <- list(
  "G1" = df4[1:4],
  "G2" = df4[5:8],
  "G3" = df4[9:12]
) %>%
  map(setNames, c("X1", "X2", "X3", "X4")) %>%
  enframe("grupo", "muestra") %>%
  mutate(muestra = map(muestra, as.matrix))

gl <- 20
df4 <- df4 %>%
  mutate(
    ni = map(muestra, dim) %>% map_int(1),
    Qi = map(muestra, Q),
    muhat = map(muestra, colMeans))

h1_test <- df4 %>%
  mutate(term = map2_dbl(Qi, ni, term_prod)) %>%
  summarise(
    num = prod(term),
    denom = det(Reduce("+", Qi)/sum(ni))^(sum(ni)/2),
    gama1 = num/denom,
    logstat = -2*log(gama1),
    pval = pchisq(logstat, df = gl, lower.tail = F)
  )

h1_test



muestras <- list()
for (i in 1:3) {
  muestras[[i]] <- df4$muestra[[i]]
}

Qmats <- map(muestras, Q)
U <- matrix(0, 4, 4)
for (q in Qmats) {
  U <- U + q
}
H <- matrix(0, 4, 4)
n <- 0
x_ <- colMeans(rbind(muestras[[1]], muestras[[2]], muestras[[3]]))
for (m in muestras) {
  ni <- dim(m)[1]
  xi_ <- colMeans(m)
  Hi <- ni * (xi_ - x_) %*% t(xi_ - x_)
  
  H <- H + Hi
  n <- n + ni
}

ratio_dets <- det(U)/det(U+H)
# Siguiendo a Boente y Wilks'
# Hago los calculos de factores a mano
F0 <- (1-sqrt(ratio_dets))/sqrt(ratio_dets) * 84/4
pf(F0, 8, 168, lower.tail = F) # No rechazo independencia

# U, H, n, k del ejemplo previo
k <- 3
D <- chol(U)
Bhat <- t(solve(D)) %*% H %*% solve(D)
eigen_hat <- eigen(Bhat)
# Solo tomo las coordenadas discriminantes que importan
# Los tomo con el signo invertido para que la primer coordenada sea positiva
T1hat <- -eigen_hat$vectors[,eigen_hat$values > 1e-10]
# Escribo 
A1hat <- sqrt(n-k) * solve(D) %*% T1hat

muestras_canonicas <- list()
for (i in seq.int(length(muestras))) {
  X <- as.matrix(muestras[[i]])
  muestras_canonicas[[i]] <- X %*% A1hat
}
muestras_canonicas

df <- map(muestras_canonicas, as_tibble) %>% bind_rows(.id = "id")
medias <- df %>%
  group_by(id) %>%
  summarise_all(mean)
g1 <- df %>%
  ggplot(aes(V1, V2, color = id)) +
  geom_point() +
  geom_point(data = medias, shape = "square", size = 3) +
  coord_fixed() +
  labs(title = "Primeras 2 coord. discr.",
       color = "Grupo",
       x = "Z1", y = "Z2")

x0 <- c(136, 143, 100, 54)
t(A1hat) %*% x0

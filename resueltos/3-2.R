library(tidyverse)

# 3-2-5
df325 <- read_delim("datos/P3-2-ej5-2019.txt", " ",col_names = F, trim_ws = T)

Y <- df325 %>%
  transmute(
    Y1 = X1 + X3,
    Y2 = X2 + X4,
    Y3 = X1 - X3,
    Y4 = X2 - X4
  ) %>%
  as.matrix()

n <- dim(Y)[1]
(n-1)*cov(Y)

Q <- function(X) {
  X <- as.matrix(X)
  # Xrulo: X centrada
  Xrulo <- sweep(X, 2, colMeans(X), "-")
  return(t(Xrulo) %*% Xrulo) 
}

# Test CMV
g1 <- c(T, T, F, F)
Q00 <- Q(Y)
Q11 <- Q00[g1,g1]
Q22 <- Q00[!g1,!g1]
gama <- det(Q00)/(det(Q11)*det(Q22)) # es un Wilks (N=n-1=24, r=2, p=2)
# Usando distr exacta porque p = r = 2
F0 <- (1-sqrt(gama))/sqrt(gama) * (n-4)/2
pf(F0, 4, 2*n-8, lower.tail = F) # No rechazo independencia
# Usano dist asintotica de CMVs,
L0 <- -2*log(gama)
pchisq(L0, 4, lower.tail = F)
# Tal vez usando Bartlett, D13 Seber

# Test Union-Interseccion
Q12 <- Q(Y)[g1,!g1]
H <- t(Q12) %*% solve(Q11) %*% Q12
U <- Q22 - H
det(U)/det(U+H)


# titaMax <- max(eigen(H%*%solve(U+H))$values)
titaMax <- max(eigen(t(Q12)%*%solve(Q11)%*%Q12%*%solve(Q22))$values)
# Hay que ir a buscarlo a D14 (Seber) con v1 y v2

# Ejemplo diapos normal multivariada
library(Flury)
data("apples")
apples
# Quiero testear Sigma1 = Sigma2 = Sigma3
# Hago la cuenta, tiene 20GL el test para igualdad de Sigmas entre 3 pobls
gl <- 20

term_prod <- function(X, n) {det(X/n)^(n/2)}
h1_test <- as_tibble(apples) %>%
  group_by(Rootstock) %>%
  nest() %>%
  mutate(
    ni = map(data, dim) %>% map_int(1),
    Qi = map(data, Q),
    term = map2_dbl(Qi, ni, term_prod) 
  ) %>%
  filter(Rootstock %in% 1:3) %>%
  summarise(
    num = prod(term),
    denom = det(Reduce("+", Qi)/sum(ni))^(sum(ni)/2),
    gama1 = num/denom,
    logstat = -2*log(gama1),
    pval = pchisq(logstat, df = gl, lower.tail = F)
)
h1_test

# No rechazo, entonces planteo test de igual de medias dada igual varianza


muestras <- list()
for (i in 1:3) {
  muestras[[i]] <- apples[apples$Rootstock==i,2:5]
}

Qmats <- map(muestras, Q)
U <- matrix(0, 4, 4)
for (q in Qmats) {
  U <- U + q
}
H <- matrix(0, 4, 4)
n <- 0
x_ <- colMeans(bind_rows(muestras))
for (m in muestras) {
  Hi <- ni * (xi_ - x_) %*% t(xi_ - x_)
  ni <- dim(m)[1]
  xi_ <- colMeans(m)

    H <- H + Hi
  n <- n + ni
}

ratio_dets <- det(U)/det(U+H)
# Siguiendo a Boente y Wilks'
# Hago los calculos de factores a mano
F0 <- (1-sqrt(ratio_dets))/sqrt(ratio_dets) * 18/4
pf(F0, 8, 36, lower.tail = F) # No rechazo independencia

#### Estas otras dos aproximaciones no me quedan tan claras. Mejor obviarlas.
# Usano dist asintotica de CMVs,
L0 <- -2*log(ratio_dets^(n/2))
pchisq(L0, 8, lower.tail = F)
# Usando U de Rao, p.43 Seber
U0 <- (1-ratio_dets)/ratio_dets * 18/3
pf(U0, 3, 18, lower.tail = F)

# U, H, n, k del ejemplo previo
k <- 3
D <- chol(U)
Bhat <- t(solve(D)) %*% H %*% solve(D)
eigen_hat <- eigen(Bhat)
# Solo tomo las coordenadas discriminantes que importan
T1hat <- eigen_hat$vectors[,eigen_hat$values > 1e-10]
#That <- eigen_hat$vectors
A1hat <- sqrt(n-k) * solve(D) %*% T1hat

muestras_canonicas <- list()
for (i in seq.int(length(muestras))) {
  X <- as.matrix(muestras[[i]])
  muestras_canonicas[[i]] <- X %*% A1hat
}
muestras_canonicas
colores <- c("red", "green", "blue")
plot(NA, NA, xli)
map(muestras_canonicas, as_tibble) %>% bind_rows(.id = "id")
bind_rows(muestras_canonicas)
for (i in seq.int(length(muestras_canonicas))) {
  m <- muestras_canonicas[[i]]
  points(m[,1], m[,2], col = colores[i], add = T)
}
df <- map(muestras_canonicas, as_tibble) %>% bind_rows(.id = "id")
medias <- df %>%
  group_by(id) %>%
  summarise_all(mean)
df %>%
  ggplot(aes(-V1, -V2, color = id)) +
  geom_point() +
  geom_point(data = medias, shape = "square", size = 3) +
  coord_fixed()

setwd("~/maestria/analisis-multivariado/practicas")
library(tidyverse)

# 1-1-6
library(mvtnorm)
nsims <- 100
mu <- c(0, 5)
Sigma <- matrix (c(
  4, 3,
  3, 9), byrow = T, ncol = 2)

# Genero observaciones de la distribucion marginal de X_1, N(0, 4)
set.seed(42)
x1 <- rnorm(nsims, mu[1], sqrt(Sigma[1,1]))
# Genero observaciones de X2|X1
mucond <- mu[2] + Sigma[1,2]/Sigma[1,1]*(x1-mu[1])
sigma2cond <- Sigma[2,2] + Sigma[1,2]^2/Sigma[1,1]
x2_x1 <- map_dbl(mucond, ~rnorm(1, ., sd = sqrt(sigma2cond)))
metodo_a <- matrix(c(x1, x2_x1), ncol = 2, byrow = F)
apply(metodo_a, 2, mean)
cov(metodo_a)

set.seed(42)
directo <- rmvnorm(nsims, mean = mu, sigma = Sigma)
apply(directo, 2, mean)
cov(directo)

set.seed(42)
z <- rmvnorm(nsims, mean = rep(0, 2))
L <- chol(Sigma)
metodo_b <- matrix(mu, byrow = T, ncol = 2, nrow = nsims) + z %*% L
apply(metodo_b, 2, mean)
cov(metodo_b)

# 1-1-7
pnorm(-3.84)
pnorm((2-5.263)/sqrt(.72))
pnorm(2, 5.263, sqrt(.72))

# 1-1-8
tibble(
  df = 1:16,
  p = map_dbl(df, ~pchisq(1, .))) %>%
  ggplot(aes(df, p)) +
  geom_line()

# 1-1-9
library(mvShapiroTest)
alfa <- 0.1
nsims <- 10
rango_n <- c(2, 5, 12, 20, 50, 100) # mvShapiro.Test pide 12 <= n <=5000
rango_distrs <- list(
  uni_std = list(mu=0, Sigma=diag(1)),
  bi_std = list(mu = c(0, 0), Sigma = diag(2)),
  bi_cust = list(mu = c(0, 5), Sigma = rbind(c(4, 3), c(3, 9))),
  tri_std = list(mu = c(0, 0, 0), Sigma = diag(3))
)
possibly_shapiro_stats <- possibly(shapiro.test, list("p.value" = NA_real_))
#ayudante_shapiro_ext <- function(muestra) {mvShapiro.Test(t(muestra))}
possibly_shapiro_ext <- possibly(mvShapiro.Test, list("p.value" = NA_real_))
ayudante_rmvnorm <- function(distr, n) {
  rmvnorm(n, mean = distr$mu, sigma = distr$Sigma)
}
simulacion_1_1_9 <- function(rango_distrs, rango_n, nsims, alfa) {
  enframe(rango_distrs, "distr", "params") %>%
    crossing(n = rango_n,
             nsim = seq.int(nsims)) %>%
    mutate(
      muestra = map2(params, n, ayudante_rmvnorm),
      shapiro_externo = map(muestra, possibly_shapiro_ext),
      shapiro_stats = map(muestra, possibly_shapiro_stats),
      pv_externo = map_dbl(shapiro_externo, "p.value"),
      pv_stats = map_dbl(shapiro_stats, "p.value")) %>%
    select(distr, n, pv_externo, pv_stats) %>%
    gather("test", "p_valor", -distr, -n) %>%
    group_by(distr, n, test) %>%
    summarise(
      alfa_emp = mean(p_valor < alfa, na.rm = T),
      sim_ok = sum(!is.na(p_valor)))
}

res119 <- simulacion_1_1_9(rango_distrs, rango_n, nsims, alfa)

grafico119 <- res119 %>%
  ggplot(aes(factor(n), alfa_emp, color = test, shape = test)) +
  geom_point() +
  geom_hline(yintercept = alfa) +
  scale_y_log10() +
  scale_shape_manual(values = c(3,4)) +
  facet_wrap(~ distr, scales = "free_y")
# ggsave("grafico119.png", grafico119)

str(shapiro.test(ayudante_rmvnorm(rango_distrs$c, 20)))
str(mvShapiroTest::mvShapiro.Test(ayudante_rmvnorm(rango_distrs$c, 20)))
# 1-1-11
mu <- c(-1, 1, 0)
a <- c(1, 2, -3)
Sigma <- matrix(
  c(1, 0, 1,
    0, 3, 1,
    1, 1, 2),
  ncol = 3, byrow = TRUE)

mu_ <- t(a)%*%mu
sigma_ <- t(a)%*%Sigma%*%a

k <- -0.4
C <- matrix(
  c(1, 0 ,0,
    1, k, -1), byrow = T, ncol = 3)
C%*%Sigma%*%t(C)

# 1-2-1.6
# rWishart ya existe!
n <- 100
m <- 20
Sigma <- diag(3)
wisharts <- rWishart(n, m, Sigma)
w_media <- apply(wisharts, c(1,2), mean)
# La media de las wisharts(m, Sigma) tiende en prob  a su esperanza, m*Sigma
abs(w_media - m*Sigma)

# 1-2-2.1
dhotelling <- function(x, m, Sigma, mu=rep(0, ncol(Sigma))) {
  p <- length(mu)
  df1 <- p
  df2 <- m-p+1
  c <- df2/(df1*m) # == (m-p+1)/(p*m)
  ncp <- t(mu) %*% solve(Sigma) %*% mu
  #print(c(p, df1, df2, c, ncp))
  return(c*df(c*x, df1, df2, ncp))
}

# 1-2-2.2
nsims(1000)
rhotelling <- function(n, m, Sigma, mu=rep(0, ncol(Sigma))) {
  retval <- vector("numeric", n)
  X <- rmvnorm(n, mu, Sigma)
  W <- rWishart(n, m, Sigma)
  for (i in seq.int(n)) {
    retval[i] <- m*(t(X[i,]) %*% solve(W[,,i]) %*% X[i,])
  }
  return(retval)
}

grafico_1222 <- function(nsims, m, Sigma, mu=rep(0, ncol(Sigma))) {
  tibble(
    valor = rhotelling(nsims, m, Sigma, mu),
    dens = dhotelling(valor, m, Sigma, mu)) %>%
    ggplot(aes(valor, ..density..)) +
    geom_histogram(bins = 100, alpha=0.3) +
    geom_density(aes(valor), color = "blue") +
    geom_line(aes(valor, dens), color = "red") +
    labs(
      title = "Densidad empirica y estimada de Hotelling",
      subtitle = paste(
        "Histograma de probabilidades en gris, ",
        "densidad empirica en azul (n =", nsims, ") y ",
        "densidad real en rojo", sep=""))
}

nsims <- 10000
m <- 20
Sigma <- rbind(
  c(5, 0, 2),
  c(0, 5, 1),
  c(2, 1, 2))
mu <- c(-1, 2, -3)

grafico_1222_central <- grafico_1222(nsims, m, Sigma)
# ggsave("grafico_1222_central.png", grafico_1222_central, width = 9, height = 5)
grafico_1222_no_central <- grafico_1222(nsims, m, Sigma, mu) 
# ggsave("grafico_1222_no_central.png", grafico_1222_central, width = 9, height = 5)

# 2-2-7
nsims <- 100000
Sigma <- rbind(
  c(5, 0, 2),
  c(0, 5, 1),
  c(2, 1, 2))
mu <- c(-1, 2, -3)
k <- 3
rmvt_ <- function(n, k, mu=rep(0, nrow(Sigma)), Sigma=diag(length(mu))) {
  p <- length(mu)
  C <- chol(Sigma) # devuelve la triangular SUPERIOR : C'C=Sigma
  v <- k/rchisq(n, df = k)
  Z <- rmvnorm(n, mean = rep(0, p), sigma = diag(p))
  return(rep(1, n)%*%t(mu) # matrix de n*p con todas las filas = mu
         + sqrt(v) * Z %*% C) # multiplica c/fila de Z por c/ elemento de sqrt(v)
}

rmvt_2 <- function(n, k, mu=rep(0, nrow(Sigma)), Sigma=diag(length(mu))) {
  v <- k/rchisq(n, df = k)
  Z <- rmvnorm(n, sigma = Sigma)
  return(rep(1, n)%*%t(mu)
         + Z * sqrt(v)) # multiplica c/fila de Z por c/ elemento de sqrt(v)
}
# rep(1, n)%*%t(vector) + matriz == sweep(matriz, 2, vector, "+")

ts <- list(
  t1 = rmvt(n=nsims, sigma=Sigma, df=k, delta=mu, type="Kshirsagar"),
  t2 = rmvt(n=nsims, sigma=Sigma, df=k, delta=mu, type="shifted"),
  t3 = rmvt_(n=nsims, k=k, mu=mu, Sigma=Sigma),
  t4 = rmvt_2(n=nsims, k=k, mu=mu, Sigma=Sigma)
)

map(ts, cov)
map(ts, ~apply(.,2, mean))
# OK!

map(list(t1,t2,t3,t4), cov)
_rmvtstd <- function(n, p, k) {
  v <- k/rchisq(n, df = k)
  Z <- rmvnorm(n, mean = rep(0, p))
  return(Z * sqrt(v)) # multiplica c/fila de Z por c/ elemento de sqrt(v)
}

rmv
rtaupk(10, 4, 100) %>% max


#####################################################
# GENERACION DE NORMALES MULTIVARIADAS: EJEMPLOS
#####################################################

####################################
# Usando la libreria mvtnorm
####################################

library(mvtnorm)

p=3
n=100

mu = rep(0,p)
sigma = diag(p)
  
set.seed(12345)
X = rmvnorm(n = n, mean = mu, sigma = sigma)
X[1:4,]

#####################################################
# OTRO EJEMPLO EN DIM 2
#####################################################

p=2
n=500

sigma <- matrix(c(4,2,2,3), ncol=2)
mu=c(1,2)

set.seed(12345)
x <- rmvnorm(n=n, mean=mu, sigma=sigma,method="chol") #USA LA DESCOMPOSICION DE CHOLESKY
x[1:4,]

plot(x[,1],x[,2],xlab=expression(x[1]),ylab=expression(x[2]),cex.lab=1.2)

colMeans(x)

var(x)
cov(x)

#############################################################################
# USANDO QUE SI Z~ N(0,I) Y SIGMA=C C' ENTONCES C Z + mu ~ N(mu, SIGMA)
# TOMO COMO C LA TRIANGULAR INFERIOR DADA EN LA DESCOMPOSICION DE CHOLESKY
############################################################################

D=chol(sigma)

t(D)%*%D #VERIFICO QUE DA SIGMA, LUEGO D=C'

C=t(D)
  
R <- chol(sigma, pivot = TRUE)

t(R)%*%R # EN ESTE CASO D y R dan lo mismo

xnew= matrix(NA, nrow=n, ncol=2)
set.seed(12345)
for (i in 1:n){
y1=rnorm(1)
y2=rnorm(1)
xnew[i,]=C %*% c(y1,y2)+mu
}
xnew[1:4,]

plot(x[,1],x[,2],xlab=expression(x[1]),ylab=expression(x[2]),cex.lab=1.2)
points(xnew[,1],xnew[,2],col="red",pch=20)

#######################################
# AMBOS METODOS DAN LO MISMO CON ESA MATRIZ SIGMA en particular,

# pero no necesariamente con todas las matrices de covarianza.
# más abajo veremos otro ejemplo donde 
# rmvnorm(...,method="chol") NO COINCIDE con 
# C Z + mu cuando C es la triangular inferior de CHOLESKY
#######################################

#############################################################################
# POR OTRO LADO, si usamos el default en mvtnorm que corresponde a method="eigen"
# es decir escribimos a SIGMA= T LAMBDA^(1/2) LAMBDA^(1/2) T', o sea, C=T LAMBDA^(1/2)
# no obtenemos los mismo datos:
#############################################################################

set.seed(12345)
x <- rmvnorm(n=n, mean=mu, sigma=sigma) #USA EIGEN 
x[1:4,]

points(x[,1],x[,2],xlab=expression(x[1]),col="blue",pch=24)

#####################################################
# COMO LA TRANSFORMACION LINEAL ES DISTINTA NO OBTENGO LOS MISMOS DATOS
#####################################################

###############################################################################
# OTRO EJEMPLO EN DIM 2, 
# donde los dos primeros métodos simulados (arriba) no coinciden.
# Veremos cuál es la descomposición de Choslesky para la cual ambos métodos coinciden.
###############################################################################

p=2
n=500

sigma <- matrix(c(4,3,3,9), ncol=2)
mu=c(1,2)

set.seed(12345)
x <- rmvnorm(n=n, mean=mu, sigma=sigma,method="chol") 
x[1:4,]

plot(x[,1],x[,2],xlab=expression(x[1]),ylab=expression(x[2]),cex.lab=1.2)

############################################################################
# USANDO QUE Z~ N(0,I) Y SIGMA=C C' ENTONCES C Z + mu ~ N(mu, SIGMA)
# TOMO COMO C LA TRIANGULAR INFERIOR DADA EN LA DESCOMPOSICION DE CHOLESKY
############################################################################

D=chol(sigma)

t(D)%*%D #VERIFICO QUE DA SIGMA, LUEGO D=C'

C=t(D)

C

#########################################
# SEA R OTRA DESCOMPOSICION DE CHOLESKY:
#########################################

R <- chol(sigma, pivot = TRUE)

t(R) # AHORA C y t(R) NO SON IGUALES

t(R)%*%R # OBSERVAR QUE INVIERTE EL LUGAR DE SIGMA_{11} y SIGMA_{22}

pivot <- attr(R, "pivot") 
oo <- order(pivot)
C2=t(R[,oo]) # TENGO QUE REORDENARLA
C2

C2%*%t(C2)   # DA SIGMA

xnew=xnew2=matrix(NA, nrow=n, ncol=2)
set.seed(12345)
for (i in 1:n){
y1=rnorm(1)
y2=rnorm(1)
xnew[i,]=C %*% c(y1,y2)+mu
xnew2[i,]=C2%*%c(y1,y2)+mu  #ESTA ES LA QUE USA rmvnorm CON "chol"
}
x[1:4,] # primer método

xnew2[1:4,] # con la descomposición de Cholesky dada por pivot = TRUE

xnew[1:4,] # con la descomposición de Cholesky habitual

plot(x[,1],x[,2],xlab=expression(x[1]),ylab=expression(x[2]),cex.lab=1.2)

points(xnew2[,1],xnew2[,2],col="blue",pch=20) # DA LO MISMO que x 

points(xnew[,1],xnew[,2],col="red",pch=20) # NO DA LO MISMO

##################################
# OBSERVACIÓN:
# con la primera matriz sigma <- matrix(c(4,2,2,3), ncol=2)
# ambas descomposiciones de Cholesky coinciden:

#chol(sigma, pivot = FALSE)
#chol(sigma, pivot = TRUE)
##################################

#####################################################################
# OTRA FORMA DE OBTENER LO MISMO QUE x y xnew2
# TRANSFORMO UNA NORMAL N(0,Identidad) generada con rmvnorm
#####################################################################

set.seed(12345)
y <- rmvnorm(n=n, mean=c(0,0), sigma=diag(p),method="chol")   

z=scale(y%*%t(C2), center=-mu, scale=F)
z[1:4,]

points(z[,1],z[,2],col="green",pch=20) 

###############################################
# DA LO MISMO QUE x =rmvnorm(n=n, mean=mu, sigma=sigma,method="chol")
###############################################

# CUIDADO!!, NO DA LO MISMO QUE:

z2<-y%*%t(C2)+mu
z2[1:4,]

# PERO SÍ DA LO MISMO QUE:

z3<-y%*%t(C2)+ matrix(rep(mu,n),n,2,byrow = T)
z3[1:4,]

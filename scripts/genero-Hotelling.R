library(mvtnorm)
######################################
# GENERO UNA WISHART W(sigma, p, n)
##############################################
p=3
n=10 

sigma = diag(p)
  
set.seed(12345)
y = rmvnorm(n = n, mean =rep(0,p), sigma = sigma)

###########################################################
# a mano
############################################################
w=matrix(0,3,3)
 
for (i in 1:n){
w=w+y[i,]%*%t(y[i,])
 

}

w

###########################################################
#OTRA FORMA ES USA QUE x TIENE A LOS INDIVIDUOS COMO FILAS
# ES DECIR x ES NUESTRA MATRIZ X POR LO TANTO W=X' X
##############################################################

W=t(y)%*%y

W
 
###################################################
# GENERO HOTELLING CENTRAL
# NECESITO UNA NORMAL INDEPENDIENTE DE W CON LA MISMA SIGMA
############################################################


mu = c(1,2,3)
 


x = rmvnorm(n = 1, mean = mu, sigma = sigma)

x.menos.mu=x-mu


T2p.n= n* x.menos.mu %*%solve(W)%*%t(x.menos.mu)

T2p.n

F=((n-p+1)/(n*p)) *T2p.n
F



eme=1000 
F=H=rep(NA,eme)
 
set.seed(12345)


for (i in 1:eme){

y = rmvnorm(n = n, mean =rep(0,p), sigma = sigma)

 #################WISHART##################

W=t(y)%*%y

##############NORMAL INDEPENDIENTE

x = rmvnorm(n = 1, mean =rep(0,p), sigma = sigma) 
 
T2p.n= n* x %*%solve(W)%*%t(x)

H[i]=T2p.n

F[i]=((n-p+1)/(n*p)) *T2p.n
}
 

tes=seq(0.05,max(F),length=1000)
densidad=df(tes,  p, n-p+1,  ncp=0,log = FALSE) 
 

plot(density(F, from=0.05, to=max(F)),lwd=4,ylim=c(0,max(densidad)),main=" ")
lines(tes,densidad,col="red",lwd=3)
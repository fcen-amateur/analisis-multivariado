library(mvtnorm)

p=3
n=10 

mu = rep(0,p)
 
sigma = diag(p)
  
set.seed(12345)
x = rmvnorm(n = n, mean = mu, sigma = sigma)

###########################################################
# a mano
############################################################
w=matrix(0,3,3)
 
for (i in 1:n){
w=w+x[i,]%*%t(x[i,])
 

}

w

###########################################################
#OTRA FORMA ES USA QUE x TIENE A LOS INDIVIDUOS COMO FILAS
# ES DECIR x ES NUESTRA MATRIZ X POR LO TANTO W=X' X
##############################################################

W=t(x)%*%x

 

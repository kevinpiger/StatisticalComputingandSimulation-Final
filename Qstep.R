
# input data
mu <- matrix(c(4,60,3,70),ncol = 2, byrow = T)
sigma <- matrix(c(1,0,0,1,1,0,0,1),ncol = 4, byrow = T)
alpha <- c(0.5,0.5)
iter_num  = 10


# Qstep
Qf <- function(data,mu,sigma,alpha,iter_num ){
  require(dplyr)
  require(mvtnorm)
newalpha <- alpha
newmu <- mu
newsigma1 <- sigma
for(j in 1:iter_num){
  p = matrix(0,nrow(faithful),length(alpha))
  temp_mu = matrix(0,nrow(faithful),2*length(alpha))
  temp_sigma = matrix(0,nrow(faithful),4*length(alpha))
  
  for( i in 1:nrow(faithful)){
    x <- faithful[i,] %>% as.numeric()
    temp = 0
    for(k in 1:length(alpha)){
      temp = temp + alpha[k]*dmvnorm(x,mu[k,],matrix(sigma[k,],2,2))
    }
    for(k in 1:length(alpha)){
      p[i,k] <- alpha[k]*dmvnorm(x,mu[k,],matrix(sigma[k,],2,2))/temp
      temp_mu[i,(2*k-1):(2*k)] <- x*p[i,k]
      temp_sigma[i,(4*k-3):(4*k)] <- (x - mu[k,])%*%t(x - mu[k,])*p[i,k]
    }}
  
   for(k in 1:length(alpha)){
    alpha[k] <- mean(p[,k])
    temp_mu1 <- temp_mu[,(2*k-1):(2*k)]
    mu[k,] <- apply(temp_mu1,2,sum)/sum(p[,k]) 
    sigma1 <- temp_sigma[,(4*k-3):(4*k)]
    sigma[k,] <- apply(sigma1,2,sum)/sum(p[,k])
   }
}
colnames(mu) <- colnames(data)
temp = list()
n1 <- NULL
n2 <- NULL
for (w in 1:nrow(mu)){
  temp[1] <- list(c(alpha))
  temp[w+1] = list(mu[w,])
  temp[(w+nrow(mu))+1] = list(matrix(sigma[1,],2,2))
  n1 <- c(n1,paste0("mu",w))
  n2 <- c(n2,paste0("sigma",w))
}
names(temp) <- c("alpha",n1,n2)
return(temp)
}

Qf <- Qf(data,mu,sigma,alpha,iter_num )



# Plot --------------------------------------------------------------------
# refer from other code
# refer from https://commons.wikimedia.org/wiki/File:Em_old_faithful.gif
xpts <- seq(from=1,to=6,length.out=100)
ypts <- seq(from=40,to=100,length.out=100)

#function to plot current data
plot.em <- function(theta){
  mixture.contour <- outer(xpts,ypts,function(x,y) {
    theta$alpha[1]*dmvnorm(cbind(x,y),mean=theta$mu1,sigma=theta$sigma1) + theta$alpha[2]*dmvnorm(cbind(x,y),mean=theta$mu2,sigma=theta$sigma2)
  })
  contour(xpts,ypts,mixture.contour,nlevels=5,drawlabel=FALSE,col="red",xlab="Eruption time",ylab="Waiting time",main="Waiting time vs Eruption time of the Old Faithful geyser")
  points(faithful)
}

plot.em(Qf)

## Describition of old faithful dataset
data = datasets::faithful
plot(data)
waiting = data$waiting
eruption = data$eruptions
hist(waiting)
hist(eruption)
s1 = sort(waiting)
s2 = sort(eruption)
table(s1) 
## from here we observe that the freq. between " 77~79"are larger than 10 and 
## both sides of them are smaller than 10,which means that the maximum pts maybe located around them.
## from here we observe that the freq. between " 49~55" are larger than 5 and 
## both sides of them are smaller than 5,which means that the maximum pts maybe located around them.

table(s2)
sd1 = sd(waiting);sd1
sd2 = sd(eruption);sd2



# Hierarchical Tree -------------------------------------------------------

# Ward
data = datasets::faithful
library(cluster)
x = daisy(data, stand=T) 
agn = agnes(x,metric="euclidean",method="ward")
plot(agn,which.plots = 1) ## two clusters
plot(agn,which.plots = 2) ## two clusters
ac = agn$ac ## 0.9975374

# Diana

di = diana(x,metric="euclidean")
plot(di, which.plots = 1) ## aslo shows that there are two clusters
plot(di, which.plots = 2)
dc = di$dc ## 0.981531


# EM ----------------------------------------------------------------------
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

# K-means -----------------------------------------------------------------


k = 2
num <- sample(1:nrow(faithful),k,replace = F)
old.center <- faithful[num,]
inter.num = 30

kms <- function(data, old, inter.num, k){
  for(z in 1:inter.num){
    new.cc <- matrix(0,nrow = k,ncol = ncol(data))
    dddd <- matrix(0,nrow=nrow(data),k)
    d.min <- matrix(0,nrow=nrow(data),1)
    for(i in 1:nrow(data)){
      for(j in 1:k){
        dddd[i,j] <- dist(rbind(data[i,],old[j,])) %>% as.numeric()}
      d.min[i] <- which.min(dddd[i,])}
    for(i in 1:k){
      A <- which(d.min == i)
      new.cc[i,] <- apply(data[A,],2,mean)
      if(new.cc[i,] == "NaN"){new.cc[i,] =0 }}
    old <- new.cc}
  
  d.min <- t(d.min)
  colnames(d.min) <- rownames(data)
  rownames(d.min) <- "Cluster"
  colnames(new.cc) <- colnames(data)
  rownames(new.cc) <- c(1:k)
  temp <- list(center = new.cc, cluster = d.min)
  return(temp)
}


kms(faithful, old.center, inter.num, k)



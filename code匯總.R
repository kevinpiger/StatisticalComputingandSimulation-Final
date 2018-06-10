
# Old faithful fountain ---------------------------------------------------

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



# Hierarchical Tree

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


# EM 
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



# Plot
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

# K-means 


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



# NIP ---------------------------------------------------------------------

library(dplyr)
library(textstem)

# Worldcloud 

library(SnowballC)
library(tm)
library(RColorBrewer)
library(wordcloud)
library(rvest)
options(scipen=999)
temp <- read.csv("~/Desktop/新統模/NIPS_1987-2015.csv")

data <- temp[,-1]
rownames(data) <- temp[,1]
nrow(data)
data[1:6,1:6]
xxxx <- apply(data,2,sum)

# 清掉全0的文章
c <- which(xxxx== 0) %>% as.numeric
data <- data[,-c]

data <- data %>% as.matrix
# 清除時態
words <- temp[,1] %>% as.character()
a <- lemmatize_words(words)
xx <- unique(a)


temp2 <- matrix(0,nrow =length(xx), ncol = ncol(data) )
for(i in 1:length(xx)){
  t <- data[which(a == xx[i]),]
  if(length(which(a == xx[i]))>1){
    temp2[i,] <- apply(t,2,sum)}else{
      temp2[i,] <- data[which(a == xx[i]),]}
}
rownames(temp2) <- xx
colnames(temp2) <- colnames(data)
nrow(temp2)
# temp2 為去除時態的資料


##########################################################################

# 去除ly，但不完全
#ly <- grep("ly",rownames(temp2))
#for(i in 2:length(rownames(temp2))){
#if(substr(rownames(temp2)[i],nchar(rownames(temp2)[i])-1,nchar(rownames(temp2)[i])) == "ly")
#  {if(substr(rownames(temp2)[i],1,nchar(rownames(temp2)[i])-2) == rownames(temp2)[i-1] ){
#  rownames(temp2)[i] = substr(rownames(temp2)[i],1,nchar(rownames(temp2)[i])-2)}}}


###########################################################################

###########################################################################
# 最常出現的字數
dataaa <- apply(temp2,1,sum)
data_o <- temp2[order(dataaa,decreasing = T),]
# 前100筆資料
data.100 <- data_o[1:100,]
data.100.sum <- apply(data.100,1,sum)
barplot(data.100.sum[1:10])
data.100.sum[1:10]

data.100.t <- t(data.100)

# km
# 用這100字對文本分類
x <- data.100.t
km<-kmeans(x,3,20)
km

pca.newseed<-princomp(scale(x,scale=TRUE,center=TRUE),cor=FALSE)
pcs.newseed<-predict(pca.newseed)

plot(pcs.newseed[,1:2], type="n")
text(pcs.newseed,as.character(km$cluster),col=km$cluster,cex=0.9)

# wordclude
oo <- as.matrix(data.100)
rownames(oo) <- rownames(data.100)
m <- oo
v <- sort(rowSums(m),decreasing=TRUE)
d<- data.frame(word=names(v), freq=v)

wordcloud(words=d$word,freq=d$freq,min.freq=2,
          max.words=50,random.order=FALSE,rot.per=0.5,
          colors=brewer.pal(5,"Dark2"))

###########################################################################

# tf
data_tf <- matrix(0,nrow = nrow(data),ncol=ncol(data))
temp1 <- apply(data,2,sum) %>% as.numeric()
for(i in 1:ncol(data)){
  data_tf[,i] <- data[,i]/temp1[i]
}

# iDF
temp1 <- ((data) > 0)*1
temp.idf <- apply(temp1,1,sum) %>% as.numeric()
data_idf <- log(ncol(data)/temp.idf)

# tf X iDF
data_w <- matrix(0,nrow = nrow(data),ncol=ncol(data))
for(i in 1:nrow(data)){
  data_w[i,] <- data_tf[i,]*data_idf[i]
}

rownames(data_w) <- rownames(data)
colnames(data_w) <- colnames(data)
data_w <- data_w %>% round(.,digits = 4)

data_i <- apply(data_w, 1,sum)
data_i.o <- data_i[order(data_i,decreasing = T)]

data_i.o[1:10]
barplot(data_i.o[1:10])

# 透過tf-idf找出前100個字
data.ti100 <- data[order(data_i,decreasing = T),]
rownames(data.ti100) <- rownames(data_w)[order(data_i,decreasing = T)]
data.ti100 <- data.ti100[1:100,]


# tf X idf 前100最大的字
temp3 <- data_w[order(data_i,decreasing = T),]
data_ii <- temp3[1:100,]
data_ii[1:10,1:10]

# wordclude
oo <- as.matrix(data_ii)
rownames(oo) <- rownames(data_ii)
m <- oo
v <- sort(rowSums(m),decreasing=TRUE)
d<- data.frame(word=names(v), freq=v)

wordcloud(words=d$word,freq=d$freq,min.freq=0.35,
          max.words=20,random.order=FALSE,rot.per=0.5,
          colors=brewer.pal(5,"Dark2"))



# 用這100字對文本分類
x <- t(data.ti100)
km<-kmeans(x,5,20)

pca.newseed<-princomp(scale(x,scale=TRUE,center=TRUE),cor=FALSE)
pcs.newseed<-predict(pca.newseed)

plot(pcs.newseed[,1:2], type="n")
text(pcs.newseed,as.character(km$cluster),col=km$cluster,cex=0.9)


# 合併年份
# tf/idf
library(stringr)

# 字切開了
a <- str_split(string =colnames(data.ti100),pattern = "X", n = 2)
temp4 <- matrix(unlist(a),ncol = 2,byrow = T)[,2]
xxx <- str_split(string =temp4,pattern = "_", n = 2)
temp4 <- matrix(unlist(xxx),ncol = 2,byrow = T)[,1] %>% as.numeric()

# 每年的文章數
table(temp4)
data.ti100.1<- data.ti100
colnames(data.ti100.1) <- temp4


temp5 <- matrix(0,nrow(data.ti100.1),ncol = length(unique(temp4)))
for(i in 1:length(unique(temp4))){
  ys <- which(colnames(data.ti100.1) == unique(temp4)[i])
  temp5[,i] <- apply(data.ti100.1[,ys],1,sum) %>% as.numeric()
}
colnames(temp5) <- unique(temp4)
rownames(temp5) <- rownames(data.ti100.1)

# 每個字每年的使用次數
temp5 %>% View()


# tf
data_tf100 <- matrix(0,nrow = nrow(temp5),ncol=ncol(temp5))
temp1 <- apply(temp5,2,sum) %>% as.numeric()
for(i in 1:ncol(temp5)){
  data_tf100[,i] <- temp5[,i]/temp1[i]
}
rownames(data_tf100) <- rownames(temp5)
colnames(data_tf100) <- colnames(temp5)

data.tf100 <- data_tf100 %>% round(., digits = 4)
data.tf100[1:10,1:10]


# 畫出每個字對於年份的趨勢 
par(mfrow = c(2,5))
for(j in 1:10){
  for(i in 1:10){
    plot(data_tf100[(i+10*(j-1)),], main = rownames(data_tf100)[(i+10*(j-1))], ylab = "frequence")
  }}

# 用年份將字群分群

km2 <- kmeans(temp5,3,20)
pca.newseed<-princomp(scale(temp5,scale=TRUE,center=TRUE),cor=FALSE)
pcs.newseed<-predict(pca.newseed)

plot(pcs.newseed[,1:2], type="n")
text(pcs.newseed,as.character(km2$cluster),col=km2$cluster,cex=0.9)

km2$cluster



# km handmade -------------------------------------------------------------

k = 3
a <- sample(1:nrow(temp5),k,replace = F)
aaaa <- temp5[a,]
inter.num = 50

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
  
  temp2 <- c(d.min)
  rownames(d.min) <- rownames(data)
  colnames(d.min) <- "Cluster"
  colnames(new.cc) <- colnames(data)
  rownames(new.cc) <- c(1:k)
  temp <- list(center = new.cc, cluster = d.min)
  return(temp)
}
km3 <- km.ll(temp5, aaaa, inter.num, k)


pca.newseed<-princomp(scale(temp5,scale=TRUE,center=TRUE),cor=FALSE)
pcs.newseed<-predict(pca.newseed)

plot(pcs.newseed[,1:2], type="n")
text(pcs.newseed,as.character(km3$cluster),col=km3$cluster,cex=0.9)

km3$cluster



# K-means by ourhand ------------------------------------------------------

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


# EM by ourhand -----------------------------------------------------------


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




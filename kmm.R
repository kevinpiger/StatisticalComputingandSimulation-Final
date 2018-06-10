# KM

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

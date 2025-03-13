source("./function.R")
library(bayesm)

## seed 1-500
set.seed(seed)



# load data
count <- read.csv("../data/count_order.csv")
count <- count[,-1]

dist <- read.csv("../data/dist_order.csv")
dist <- as.matrix(dist[,-1])



init_dataset_acc_sim <- function(m,n,count,dist,k,neighbor = 5, sigma = 1){
  d <- nrow(count)
  O_list <- list()
  w_list <- list()
  X_list <- list()
  
  A <- matrix(0,d,d)
  for(i in 1:d){
    idx <- order(dist[i,],decreasing = F)[2:(neighbor+1)]
    A[i,idx] <- exp(-dist[i,idx]/sigma)
    A[idx,i] <- exp(-dist[i,idx]/sigma)
  }
  D1 <- diag(rowSums(A))
  L <- D1-A
  svd_res <- svd(L)
  U <- svd_res$u
  U <- U[,(ncol(U)-k+1):(ncol(U))]

  
  w_list <- matrix(0,m,d)
  
  for(i in 1:m){
    idx1 <- sample(ncol(count),n)
    X0 <- as.matrix(count[,idx1])
    X <- sapply(1:n,function(j)rdirichlet(X0[,j]+0.1))
    X_lib <- X%*%diag(sample(10000:15000,n,replace = T))
    
    
    w_space <- U
    weight <- 1-2*runif(k)
    
    w <- (w_space%*%as.matrix(weight))[,1]
    w <- (w-min(w)+0.05)
    w <- w/max(w)
    
    
    O_list[i] <- list(diag(w)%*%X_lib)
    w_list[i,] <- w
    X_list[i] <- list(X)
  }
  return(list(O_list, w_list, X_list))
}

ce_accuracy <- function(O_list,w_list,w_list_t){
  m <- nrow(w_list)
  d <- ncol(w_list)
  w_list_r <- w_list/mean(w_list)*mean(w_list_t)
  w_diff <- sapply(1:m,function(i)abs(w_list_r[i,]-w_list_t[i,]))
  acc <- rowMeans(w_diff)
  O <- do.call("cbind",O_list)
  abundance_rank <- rank(rowSums(O))
  res <- data.frame("acc"=acc, "rank" = abundance_rank)
  res <- res[order(res$rank,decreasing=F),]
  return(res)
}


m <- 2
n <- 100
alpha <- 0.01
gamma <- 10
r <- 113
k <- 1

data1 <- init_dataset_acc_sim(m,n,count,dist,k=10)
w_list_t <- data1[[2]]
O_list <- data1[[1]]

# result after stage 1
metadict_res0 <- metadict(O_list,0,0,0,dist,meta.list=NULL,neighbor=3,sigma=1000)
w_list_raw <- metadict_res0$w
res0 <-  ce_accuracy(O_list,w_list_raw,w_list_t)
#filename = paste0("./table/ACC_beta_raw/","metadict", "_simuiter", seed, ".csv")
#write.csv(res0, file = filename)

# result after stage 1 and stage 2
metadict_res1 <- metadict(O_list,alpha,1,gamma,dist,meta.list= NULL, neighbor=3, sigma=1000)
w_list1 <- metadict_res1$w
res2 <-  ce_accuracy(O_list,w_list1,w_list_t)
#filename = paste0("./table/ACC_beta1/","metadict", "_simuiter", seed, ".csv")
#write.csv(res2, file = filename)

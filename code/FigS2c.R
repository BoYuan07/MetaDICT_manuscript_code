source("./function.R")
library(bayesm)

## seed 1-500

#set.seed(seed)


# load order level data
count <- read.csv("../data/count_order.csv")
count <- count[,-1]

dist <- read.csv("../data/dist_order.csv")
dist <- as.matrix(dist[,-1])



init_dataset_acc_sim <- function(m,nmax,count,dist,k,neighbor = 5, sigma = 1, balance = T){
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
    if(balance){
      n <- nmax
    }else{
      n <- sample(50:nmax,1)
    }
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

balance.list <- c(T,F)
m <- 4
n <- 100
res <- matrix(0,2,1)
alpha <- 0.01
beta <- 1
gamma <- 10

for(j in 1:2){
  balance <- balance.list[j]
  data1 <- init_dataset_acc_sim(m,n,count,dist,k=10,neighbor=5,balance=balance)
  w_list_t <- data1[[2]]
  
  
  O_list <- data1[[1]]
  metadict_res <- metadict(O_list,alpha,beta,gamma,dist,meta.list=NULL,neighbor=3,sigma=1000)
  w_list <- metadict_res$w

  w_all <- c(w_list)
  w_all_t <- c(w_list_t)
  res[j,1] <- mean(cor(w_all,w_all_t))
}

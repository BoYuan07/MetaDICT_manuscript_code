source("./function.R")
library(bayesm)

## seed 1-500
set.seed(seed)


# load order level data
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



beta.list <- c(0,1,5)
m <- 4
n <- 100
res <- matrix(0,length(beta.list),2)

alpha <- 0.01
gamma <- 10
r <- 113
k <- 1

data1 <- init_dataset_acc_sim(m,n,count,dist,k=10)
w_list_t <- data1[[2]]
O_list <- data1[[1]]

metadict_res0 <- metadict(O_list,0,0,0,dist,meta.list=NULL,neighbor=3,sigma=1000)
w_list0 <- metadict_res0$w

metadict_res <- metadict(O_list,alpha,1,gamma,dist,meta.list=NULL,neighbor=3,sigma=1000)
w_list <- metadict_res$w

w_all0 <- c(w_list0)
w_all_t <- c(w_list_t)
w_all1 <- c(w_list)

# pearson correlation after stage 1
res[1,2] <- mean(cor(w_all0,w_all_t))
# pearson correlation after stage 1 and stage 2
res[2,2] <- mean(cor(w_all1,w_all_t))



#filename = paste0("./table/ACC_refine/","metadict", "_simuiter", seed, ".rds")
#saveRDS(res, file = filename)

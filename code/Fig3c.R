library(ConQuR)
library(doParallel)
library(MMUPHin)
source("./function.R")

## seed 1-500
#jobid <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#simu.iter = jobid
#set.seed(jobid)

count = read.csv("../data/count_order.csv")
count = count[,-1]

dist = read.csv("../data/dist_order.csv")
dist = as.matrix(dist[,-1])


library(bayesm)
init_dataset = function(m,n,count,dist,k,neighbor = 5, sigma = 1, ez = 10, lib_var = T, eta = c(1/2,1/2)){
  d = nrow(count)
  O_list = list()
  X_list = list()
  meta.list = list()
  meta.list1 = list()
  
  prevalence = rowSums(count!=0)
  d1 = sample(1:d,40)
  
  A = matrix(0,d,d)
  for(i in 1:d){
    idx = order(dist[i,],decreasing = F)[2:(neighbor+1)]
    A[i,idx] = exp(-dist[i,idx]/sigma)
    A[idx,i] = exp(-dist[i,idx]/sigma)
  }
  D1 = diag(rowSums(A))
  L = D1-A
  svd_res = svd(L)
  U = svd_res$u
  U = U[,(ncol(U)-k+1):(ncol(U))]
  
  w_list = matrix(0,m,d)
  
  for(i in 1:m){
    idx1 = sample(ncol(count),n)
    X0 = as.matrix(count[,idx1])
    X0 = t(t(X0)/colSums(X0))
    Y = sample(1:2,size=n,replace=TRUE,prob=c(eta[i],1-eta[i]))
    X0[d1,which(Y==1)] = (X0[d1,which(Y==1)]+0.1)*ez
    X = sapply(1:n,function(j)rdirichlet(X0[,j]))

    if(lib_var){
      if(i%%2==1){
        X_lib = X_lib%*%diag(sample(1000:2000,n,replace = T))
      }else{
        X_lib = X_lib%*%diag(sample(4000:5000,n,replace = T))
        }
    }else{
      X = X%*%diag(sample(10000:11000,n,replace = T))
      X_lib = X
    }
    
    w_space = U
    weight = 1-2*runif(k)
    w = (w_space%*%as.matrix(weight))[,1]
    w = (w-min(w)+0.05)
    w = w/max(w)
    
    Y2 = as.factor(rbinom(n,1,1/2))
    O_list[i] = list(floor(diag(w)%*%X_lib))
    w_list[i,] = w
    X_list[i] = list(floor(X_lib))
    meta.list[[i]] = data.frame("Y" = sapply(Y,function(x)paste("Group",x)),"Y2" = Y2)
    meta.list1[[i]] = data.frame("Y2" = Y2)
  }
  return(list(O_list, w_list, X_list, meta.list, meta.list1))
}

eta.list = list(c(1/2,1/2),c(1/4,3/4),c(1/6,5/6))
alpha = 1
beta = 0.1
gamma = 10
m = 2
n = 100

d = nrow(count)
res = data.frame("R2" = numeric(),"Type" = character(),"Method" = character(),"Level" = character(),"Cov" = character())
levels = c("None","Medium","High")

for(i in 1:3){
  eta = eta.list[[i]]
  data = init_dataset(m,n,count,dist,k=10, ez=10, lib_var = F, eta = eta)
  
  O_list = data[[1]]
  w_list_t = data[[2]]
  X_list = data[[3]]
  meta.list = data[[4]]
  meta.list.sub = data[[5]]
  
  meta = do.call("rbind",meta.list)
  batchid = as.factor(do.call(c,lapply(1:m, function(x)rep(x,n))))
  
  O = do.call(cbind,O_list)
  X_t = do.call(cbind,X_list)
  
  meta$batch =as.factor(do.call(c,lapply(1:m, function(x)rep(paste("Batch",x),n))))
  meta_sub = meta[,-1]
  dataset_info = meta$batch
  sample_info = meta$Y


  metadict_res_1 = metadict(X_list,alpha,beta,gamma,dist,meta.list = meta.list.sub)
  metadict_res_2 = metadict(X_list,alpha,beta,gamma,dist,meta.list = meta.list)

  res_ComBatSeq_1 = sva::ComBat_seq(as.matrix(X_t),dataset_info,covar_mod = as.data.frame(meta_sub[,-2]))
  res_ComBatSeq_2 = sva::ComBat_seq(as.matrix(X_t),dataset_info, covar_mod = meta[,-3])

  colnames(X_t) = sapply(1:ncol(X_t),function(x)paste("Sample",x))
  rownames(meta) = sapply(1:ncol(X_t),function(x)paste("Sample",x))
  meta$batch = as.factor(meta$batch)
  res_mmuphin_1 = adjust_batch(feature_abd = X_t,
                                batch = "batch",
                              covariates = "Y2",
                                data = meta)$feature_abd_adj
  res_mmuphin_2 = adjust_batch(feature_abd = X_t,
                              batch = "batch",
                              covariates = c("Y","Y2"),
                              data = meta)$feature_abd_adj

  tax_tab = t(X_t)
  res_ConQuR_1 = ConQuR(tax_tab,batchid,batch_ref = 1,covariates = as.data.frame(meta_sub[,-2]))
  meta$Y = as.factor(meta$Y)
  res_ConQuR_2 = ConQuR(tax_tab,batchid,batch_ref = 1,covariates = meta[,-3])

  res = rbind(res,c(permanova(O,batchid),"Batch","Unprocessed",levels[i],"Not observed"))
  res = rbind(res,c(permanova(O,batchid),"Batch","Unprocessed",levels[i],"Observed"))
  res = rbind(res,c(permanova(X_t,batchid),"Batch","Truth",levels[i],"Not observed"))
  res = rbind(res,c(permanova(X_t,batchid),"Batch","Truth",levels[i],"Observed"))
  res = rbind(res,c(permanova(O,sample_info),"Sample","Unprocessed",levels[i],"Not observed"))
  res = rbind(res,c(permanova(O,sample_info),"Sample","Unprocessed",levels[i],"Observed"))
   res = rbind(res,c(permanova(X_t,sample_info),"Sample","Truth",levels[i],"Not observed"))
  res = rbind(res,c(permanova(X_t,sample_info),"Sample","Truth",levels[i],"Observed"))

  res = rbind(res,c(permanova(metadict_res_1$X,batchid),"Batch","MetaDICT",levels[i],"Not observed"))
  res = rbind(res,c(permanova(metadict_res_2$X,batchid),"Batch","MetaDICT",levels[i],"Observed"))
  res = rbind(res,c(permanova(metadict_res_1$X,sample_info),"Sample","MetaDICT",levels[i],"Not observed"))
  res = rbind(res,c(permanova(metadict_res_2$X,sample_info),"Sample","MetaDICT",levels[i],"Observed"))

  res = rbind(res,c(permanova(res_ComBatSeq_1,batchid),"Batch","ComBatSeq",levels[i],"Not observed"))
  res = rbind(res,c(permanova(res_ComBatSeq_2,batchid),"Batch","ComBatSeq",levels[i],"Observed"))
  res = rbind(res,c(permanova(res_ComBatSeq_1,sample_info),"Sample","ComBatSeq",levels[i],"Not observed"))
  res = rbind(res,c(permanova(res_ComBatSeq_2,sample_info),"Sample","ComBatSeq",levels[i],"Observed"))

  res = rbind(res,c(permanova(t(res_ConQuR_1),batchid),"Batch","ConQuR",levels[i],"Not observed"))
  res = rbind(res,c(permanova(t(res_ConQuR_2),batchid),"Batch","ConQuR",levels[i],"Observed"))
  res = rbind(res,c(permanova(t(res_ConQuR_1),sample_info),"Sample","ConQuR",levels[i],"Not observed"))
  res = rbind(res,c(permanova(t(res_ConQuR_2),sample_info),"Sample","ConQuR",levels[i],"Observed"))

  res = rbind(res,c(permanova(res_mmuphin_1,batchid),"Batch","MMUPHin",levels[i],"Not observed"))
  res = rbind(res,c(permanova(res_mmuphin_2,batchid),"Batch","MMUPHin",levels[i],"Observed"))
  res = rbind(res,c(permanova(res_mmuphin_1,sample_info),"Sample","MMUPHin",levels[i],"Not observed"))
  res = rbind(res,c(permanova(res_mmuphin_2,sample_info),"Sample","MMUPHin",levels[i],"Observed"))
}



#filename = paste0("./table/AA_null/","sim_", simu.iter, ".csv")
#write.csv(res, file = filename)



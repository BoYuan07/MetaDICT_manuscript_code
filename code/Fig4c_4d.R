library(ConQuR)
library(doParallel)
library(MMUPHin)
library(bayesm)
source("./function.R")

## seed 1-500
#jobid <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#simu.iter = jobid
#set.seed(jobid)

count = read.csv("../data/count_order.csv")
count = count[,-1]

dist = read.csv("../data/dist_order.csv")
dist = as.matrix(dist[,-1])


init_dataset = function(m,n,count,dist,k,neighbor = 5, sigma = 1, ez = 10, lib_var = T, p = list(c(1/4,1/4,1/4,1/4),c(1/4,1/4,1/4,1/4),c(1/4,1/4,1/4,1/4),c(1/4,1/4,1/4,1/4),c(1/4,1/4,1/4,1/4))){
  d = nrow(count)
  O_list = list()
  X_list = list()
  meta.list = list()
  meta.list1 = list()
  
  taxagroup = sample(1:5,size=d,replace=TRUE,prob=c(1/5,1/5,1/5,1/5,1/5))
  
  for(i in 1:d){
    idx = order(dist[i,],decreasing = F)[2:(neighbor+1)]
    A[i,idx] = exp(-dist[i,idx]/sigma)
    A[idx,i] = exp(-dist[i,idx]/sigma)
  }
  
  A = matrix(0,d,d)
  D1 = diag(rowSums(A))
  L = D1-A
  svd_res = svd(L)
  U = svd_res$u
  U = U[,(ncol(U)-k+1):(ncol(U))]
  
  w_list = matrix(0,m,d)
  
  for(i in 1:m){
    idx1 = sample(ncol(count),n)
    X = as.matrix(count[,idx1])
    Y = sample(1:4,size=n,replace=TRUE,prob=p[[i]])
    
    X[which(taxagroup == 2),which(Y==1)] = (X[which(taxagroup == 2),which(Y==1)]+1)*ez
    X[which(taxagroup == 3),which(Y==2)] = (X[which(taxagroup == 3),which(Y==2)]+1)*ez
    X[which(taxagroup == 4),which(Y==3)] = (X[which(taxagroup == 4),which(Y==3)]+1)*ez
    X[which(taxagroup == 5),which(Y==4)] = (X[which(taxagroup == 5),which(Y==4)]+1)*ez

    w_space = U
    weight = 1-2*runif(k)
    w = (w_space%*%as.matrix(weight))[,1]
    w = (w-min(w)+0.05)
    w = w/max(w)

    if(lib_var){
      if(i%%2==1){
        X_lib = X%*%diag(sample(100:500,n,replace = T))
      }else{
        X_lib = X%*%diag(sample(10:50,n,replace = T))
        }
    }
    
    Y_f = as.factor(rbinom(n,1,1/2))
    O_list[i] = list(floor(diag(w)%*%X_lib))
    w_list[i,] = w
    X_list[i] = list(X)
    meta.list[[i]] = data.frame("Y" = sapply(Y,function(x)paste("Group",x)),"Yf" = Y_f)
    meta.list1[[i]] = data.frame("Yf" = Y_f)
  }
  return(list(O_list, w_list, X_list, meta.list, taxagroup, meta.list1))
}

m = 5
n = 50
p1 = list(c(1/4,1/4,1/4,1/4),c(1/4,1/4,1/4,1/4),c(1/4,1/4,1/4,1/4),c(1/4,1/4,1/4,1/4),c(1/4,1/4,1/4,1/4))
p2 = list(c(1/8,1/2,1/8,1/4),c(1/2,1/8,1/4,1/8),c(1/4,1/4,1/4,1/4),c(1/4,1/2,1/8,1/8),c(1/8,1/8,1/2,1/4))
p.list = list(p1,p2)

level = c("None","Medium")
res = data.frame("ARI" = numeric(),"Method" = character(),"level" = character(), "Cov" = character(),"Type" = character())

for(i in 1:2){
  p = p.list[[i]]
  data = data = init_dataset(m,n,count,dist,k=10, ez=50, p=p)
  
  O_list = data[[1]]
  w_list_t = data[[2]]
  X_t = do.call(cbind,data[[3]])
  meta.list = data[[4]]
  diff_set = data[[5]]
  meta.list1 = data[[6]]
  
  meta = do.call("rbind",meta.list)
  meta$Y = as.factor(meta$Y)
  
  meta1 = do.call("rbind",meta.list1)
  batchid = as.factor(do.call(c,lapply(1:m, function(x)rep(x,n))))
  O = do.call(cbind,O_list)
  
  #######unprocessed data
  res = rbind(res,c(adjustedRandIndex(community_detection(O,K=10)$cluster,diff_set),"Unprocessed",level[i],"With","Taxa"))
  res = rbind(res,c(adjustedRandIndex(community_detection(O,K=10)$cluster,diff_set),"Unprocessed",level[i],"Without","Taxa"))
  res = rbind(res,c(adjustedRandIndex(community_detection(t(O),K=50)$cluster,meta$Y),"Unprocessed",level[i],"With","Sample"))
  res = rbind(res,c(adjustedRandIndex(community_detection(t(O),K=50)$cluster,meta$Y),"Unprocessed",level[i],"Without","Sample"))



  #######single dataset
  ari = c()
  avg_silwidth = c()
  for(j in 1:m){
      membership = community_detection(O_list[[j]],K=10)$cluster
      if(length(unique(membership))!=1){
            avg_silwidth= c(avg_silwidth,mean(cluster::silhouette(membership, dist(O_list[[j]], method = "euclidean"))[,3], na.rm = TRUE))
            ari = c(ari,adjustedRandIndex(membership,diff_set))
        }
  }
  res = rbind(res,c(ari[which.max(avg_silwidth)],"Single Dataset",level[i],"With","Taxa"))
  res = rbind(res,c(ari[which.max(avg_silwidth)],"Single Dataset",level[i],"Without","Taxa"))
        


  #######Metadict
  alpha = 1
  beta = 0.1
  gamma = 1
  metadict_res = metadict(O_list,alpha,beta,gamma,dist,meta.list)
  D = metadict_res$D
  R = do.call(cbind,metadict_res$R)
  X = metadict_res$X
  c_metadict = community_detection(D[,1:50],K=10)$cluster
  c_metadict_s = community_detection(t(R),K=50)$cluster
  res = rbind(res,c(adjustedRandIndex(c_metadict,diff_set),"MetaDICT",level[i],"With","Taxa"))
  res = rbind(res,c(adjustedRandIndex(c_metadict_s,meta$Y),"MetaDICT",level[i],"With","Sample"))

  metadict_res1 = metadict(O_list,alpha,beta,gamma,dist,meta.list1)
  D1 = metadict_res1$D
  R1 = do.call(cbind,metadict_res1$R)
  X1 = metadict_res1$X
  c_metadict1 = community_detection(D1[,1:50],K=10)$cluster
  c_metadict_s1 = community_detection(t(R1),K=50)$cluster
  res = rbind(res,c(adjustedRandIndex(c_metadict1,diff_set),"MetaDICT",level[i],"Without","Taxa"))
  res = rbind(res,c(adjustedRandIndex(c_metadict_s1,meta$Y),"MetaDICT",level[i],"Without","Sample"))

  
  #######ComBatSeq
  res_ComBatSeq = sva::ComBat_seq(as.matrix(O),batchid, covar_mod = meta)
  c_combatseq = community_detection(res_ComBatSeq,K=10)$cluster
  c_combatseq_s = community_detection(t(res_ComBatSeq),K=50)$cluster
  res = rbind(res,c(adjustedRandIndex(c_combatseq,diff_set),"ComBatSeq",level[i],"With","Taxa"))
  res = rbind(res,c(adjustedRandIndex(c_combatseq_s,meta$Y),"ComBatSeq",level[i],"With","Sample"))

  res_ComBatSeq1 = sva::ComBat_seq(as.matrix(O),batchid, covar_mod = meta1)
  c_combatseq1 = community_detection(res_ComBatSeq1,K=10)$cluster
  c_combatseq_s1 = community_detection(t(res_ComBatSeq1),K=50)$cluster
  res = rbind(res,c(adjustedRandIndex(c_combatseq1,diff_set),"ComBatSeq",level[i],"Without","Taxa"))
  res = rbind(res,c(adjustedRandIndex(c_combatseq_s1,meta$Y),"ComBatSeq",level[i],"Without","Sample"))


  ########MMUPHin
  colnames(O) = sapply(1:ncol(O),function(x)paste("Sample",x))
  rownames(meta) = sapply(1:ncol(O),function(x)paste("Sample",x))
  meta$batch = as.factor(batchid)
  res_mmuphin = adjust_batch(feature_abd = O,
                                batch = "batch",
                              covariates = c("Y","Yf"),
                                data = meta)$feature_abd_adj
  c_mmuphin = community_detection(res_mmuphin,K=10)$cluster
  c_mmuphin_s = community_detection(t(res_mmuphin),K=50)$cluster
  res = rbind(res,c(adjustedRandIndex(c_mmuphin,diff_set),"MMUPHin",level[i],"With","Taxa"))
  res = rbind(res,c(adjustedRandIndex(c_mmuphin_s,meta$Y),"MMUPHin",level[i],"With","Sample"))

  
  res_mmuphin1 = adjust_batch(feature_abd = O,
                                batch = "batch",
                              covariates = "Yf",
                                data = meta)$feature_abd_adj
  c_mmuphin1 = community_detection(res_mmuphin1,K=10)$cluster
  c_mmuphin_s1 = community_detection(t(res_mmuphin1),K=50)$cluster
  res = rbind(res,c(adjustedRandIndex(c_mmuphin1,diff_set),"MMUPHin",level[i],"Without","Taxa"))
  res = rbind(res,c(adjustedRandIndex(c_mmuphin_s1,meta$Y),"MMUPHin",level[i],"Without","Sample"))
  

  ###########ConQUR
  tax_tab = t(O)
  res_conqur = ConQuR(tax_tab,batchid,batch_ref = 1,covariates = meta[,-3])
  c_conqur = community_detection(t(res_conqur),K=10)$cluster
  c_conqur_s = community_detection(res_conqur,K=50)$cluster
  res = rbind(res,c(adjustedRandIndex(c_conqur,diff_set),"ConQuR",level[i],"With","Taxa"))
  res = rbind(res,c(adjustedRandIndex(c_conqur_s,meta$Y),"ConQuR",level[i],"With","Sample"))

  res_conqur1 = ConQuR(tax_tab,batchid,batch_ref = 1,covariates = meta1)
  c_conqur1 = community_detection(t(res_conqur1),K=10)$cluster
  c_conqur_s1 = community_detection(res_conqur1,K=50)$cluster
  res = rbind(res,c(adjustedRandIndex(c_conqur1,diff_set),"ConQuR",level[i],"Without","Taxa"))
  res = rbind(res,c(adjustedRandIndex(c_conqur_s1,meta$Y),"ConQuR",level[i],"Without","Sample"))             
}



#filename = paste0("./table/joint_cluster_cov/","sim_", simu.iter, ".csv")
#write.csv(res, file = filename)



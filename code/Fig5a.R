library(ConQuR)
library(doParallel)
library(MMUPHin)
library(effsize)
source("./function.R")

## seed 1-500
#jobid <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#simu.iter = jobid
#set.seed(jobid)

count = read.csv("../data/count_order.csv")
count = count[,-1]

dist = read.csv("../data/dist_order.csv")
dist = as.matrix(dist[,-1])

misclassify_times = function(bias,res){
  final_res = rep(0,length(res))
  final_res[which(res==T)] = 1
  return(list("False" = final_res,"Bias" = bias))
}

init_dataset_da = function(m,n,count,dist,k,k1,p,ez=5,lib_var=T){
  d = nrow(count)
  O_list = list()
  w_list = list()
  X_list = list()
  meta_list = list()
  
  prevalence = rowSums(count!=0)
  diff_idx = order(prevalence,decreasing = T)[k:1]
  idx1 = sample(k,k1)
  seq_idx1 = diff_idx[idx1]
  seq_idx2 = diff_idx[-idx1]

  A = exp(-dist/0.1)
  L = diag(colSums(A))-A
  svd_res = svd(L)
  U = svd_res$u
  U = U[,(ncol(U)-100):(ncol(U))]
  
  control_sample_list = list()
  
  for(i in 1:m){
    idx = sample(ncol(count),n)
    X = as.matrix(count[,idx])
    Y = sample(0:1,n,T,prob = c(p[i],1-p[i]))
   
    X[seq_idx1,which(Y==0)] = (X[seq_idx1,which(Y==0)])/ez
    X[seq_idx2,which(Y==0)] = (X[seq_idx2,which(Y==0)])*ez
    
    index = sample(100,30)
    w_space = U[,index]
    weight = 1-2*runif(30)
    
    w = (w_space%*%as.matrix(weight))[,1]
    w = (w-min(w)+1e-3)
    w = w/max(w)
    w = w/mean(w)*0.5

  if(lib_var){
      if(i%%2 == 1){
        print(i)
        X_lib = X%*%diag(sample(100:500,n,replace = T))
      }else{
        X_lib = X%*%diag(sample(10:100,n,replace = T))
        }
    }else{
      X_lib = X
    }
   
    
    meta_list[[i]] = data.frame("Y" = Y) 
    O_list[i] = list(floor(diag(w)%*%X_lib))
    w_list[i] = list(w)
    X_list[i] = list(floor(X_lib))
    diff_idx = c(seq_idx1,seq_idx2)
  }
  return(list(O_list, w_list, X_list, diff_idx, meta_list))
}

t_test_taxa = function(X,Y,Clade,alpha,method = "BH"){
  d = nrow(X)
  data = cbind(t(X),Y)
  t_p <- sapply(1:d, function(i)t.test(data[, i] ~ Y, data = data)$p.value)
  p <- p.adjust(t_p, method=method)
  p.res = which(p<alpha)
  if(length(res)==0){
    fdr = 0
  }else{
    fdr = mean(!p.res %in% Clade)
  }
  res = rep(F,d)
  res[which(p<alpha)] = T
  return(as.factor(res))
}

m = 2
n = 50
alpha = 0.01
beta = 0.01
gamma = 1
alpha = 0.1
d = nrow(count)

res = init_dataset_da(m,n,count,dist,floor(0.05*d),floor(0.02*d),ez=5, lib_var = T, p = c(1/4,3/4))
O_list = res[[1]]
X_t = do.call(cbind,res[[3]])
diff_seq = res[[4]]
meta_list = res[[5]]
meta = do.call(rbind,meta_list)
O = do.call(cbind,O_list)
batchid = as.factor(do.call(c,lapply(1:m, function(x)rep(x,n))))
meta$batch = batchid

Y = as.factor(meta$Y)

bias = sapply(1:nrow(O),function(i)cohen.d(O[i,]~batchid)$estimate)

X = metadict(O_list,alpha,beta,gamma,dist,meta_list)$X
res.ConQuR = ConQuR(t(O),batchid,batch_ref = 1,covariates = Y)
res.ComBatSeq = sva::ComBat_seq(as.matrix(O),batchid,Y)
colnames(O) = sapply(1:ncol(O),function(x)paste("Sample",x))
rownames(meta) = sapply(1:ncol(O),function(x)paste("Sample",x))
res.mmuphin = adjust_batch(feature_abd = O,
                              batch = "batch",
                              covariates = "Y",
                              data = meta)$feature_abd_adj

res1 = t_test_taxa(O,Y,diff_seq,alpha = 0.1,method = "BH")
res2 = t_test_taxa(X,Y,diff_seq,alpha = 0.1,method = "BH")
res3 = t_test_taxa(res.ComBatSeq,Y,diff_seq,alpha = 0.1,method = "BH")
res4 = t_test_taxa(res.mmuphin,Y,diff_seq,alpha = 0.1,method = "BH")
res5 = t_test_taxa(t(res.ConQuR),Y,diff_seq,alpha = 0.1,method = "BH")

m1 = misclassify_times(bias[-diff_seq],res1[-diff_seq])
m2 = misclassify_times(bias[-diff_seq],res2[-diff_seq])
m3 = misclassify_times(bias[-diff_seq],res3[-diff_seq])
m4 = misclassify_times(bias[-diff_seq],res4[-diff_seq])
m5 = misclassify_times(bias[-diff_seq],res5[-diff_seq])

table = data.frame("False" = c(m1[[1]],m2[[1]],m3[[1]],m4[[1]],m5[[1]]),"Bias" = c(m1[[2]],m2[[2]],m3[[2]],m4[[2]],m5[[2]]),"Abundance" = c(m1[[3]],m2[[3]],m3[[3]],m4[[3]],m5[[3]]),"Rank" = rep(c(1:length(m1[[1]])),5),
"Method" = c(rep("Unprocessed",length(m1[[1]])),rep("MetaDICT",length(m1[[1]])),rep("ComBatSeq",length(m1[[1]])),rep("MMUPHin",length(m1[[1]])),rep("ConQuR",length(m1[[1]]))))

#filename = paste0("./table/DA_taxa/","sim_", simu.iter, ".csv")
#write.csv(table, file = filename)



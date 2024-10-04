library(ConQuR)
library(doParallel)
library(MMUPHin)
source("./function.R")

count = read.csv("../data/count_order.csv")
count = count[,-1]

dist = read.csv("../data/dist_order.csv")
dist = as.matrix(dist[,-1])

## seed 1-500
#jobid <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#simu.iter = jobid
#set.seed(jobid)

init_dataset_da = function(m,n,count,dist,lib_var=T,p = c(1/4,1/6)){
  d = nrow(count)
  O_list = list()
  w_list = list()
  X_list = list()
  meta_list = list()

  A = exp(-dist/0.1)
  L = diag(colSums(A))-A
  svd_res = svd(L)
  U = svd_res$u
  U = U[,(ncol(U)-100):(ncol(U))]
  
  for(i in 1:m){
    idx = sample(ncol(count),n)
    X = as.matrix(count[,idx])

    Y = sample(0:1,n,T,prob = c(p[i],1-p[i]))
    
    index = sample(100,30)
    w_space = U[,index]
    weight = 1-2*runif(30)
    
    w = (w_space%*%as.matrix(weight))[,1]
    w = (w-min(w)+1e-3)
    w = w/max(w)
    w = w/mean(w)*0.5

    if(lib_var){
      if(i%%2 == 1){
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
    X_list[i] = list(floor(X))
  }
  return(list(O_list, w_list, X_list, meta_list))
}


t_test_res = function(X,Y,alpha=0.1){
  X = X[rowSums(X)>0,]
  t_p = sapply(1:(nrow(X)), function(i)t.test(as.numeric(X[i,]) ~ Y)$p.value)
  p <- p.adjust(t_p, method='BH')
  return(mean(p<alpha))
}

p1 = c(0.5,0.5,0.5)
p2 = c(0.8,0.2,0.8)
p3 = c(0.99,0.01,0.99)
p_list = list(p1,p2,p3)

m = 3
n = 100
alpha = 0.01
beta = 0.01
gamma = 1
alpha = 0.1
d = nrow(count)

fdr_table = data.frame("FDR" = numeric(),"Method" = character(), "Level" = character())
Level.list = c("None","Moderate","Strong")

for(i in 1:length(p_list)){
  p = p_list[[i]]
  res = init_dataset_da(m,n,count,dist,lib_var = T,p = p)
  O_list = res[[1]]
  meta_list = res[[4]]
  meta = do.call(rbind,meta_list)
  O = do.call(cbind,O_list)
  batchid = as.factor(do.call(c,lapply(1:m, function(x)rep(x,n))))
  meta$batch = batchid

  Y = as.factor(meta$Y)
  
  X = metadict(O_list,alpha,beta,gamma,dist,meta_list)$X
  res.ConQuR = ConQuR(t(O),batchid,batch_ref = 1,covariates = Y)
  res.ComBatSeq = sva::ComBat_seq(as.matrix(O),batchid,Y)
  colnames(O) = sapply(1:ncol(O),function(x)paste("Sample",x))
  rownames(meta) = sapply(1:ncol(O),function(x)paste("Sample",x))
  res.mmuphin = adjust_batch(feature_abd = O,
                                batch = "batch",
                                covariates = "Y",
                                data = meta)$feature_abd_adj
  
  t_raw = t_test_res(O,meta$Y)
  fdr_table = rbind(fdr_table,c(t_raw,"Unprocessed",Level.list[i]))
  
  t_combatseq = t_test_res(res.ComBatSeq,meta$Y)
  fdr_table = rbind(fdr_table,c(t_combatseq,"ComBatSeq",Level.list[i]))
  
  t_mmuphin = t_test_res(res.mmuphin,meta$Y)
  fdr_table = rbind(fdr_table,c(t_mmuphin,"MMUPHin",Level.list[i]))
  
  t_metadict = t_test_res(X,meta$Y)
  fdr_table = rbind(fdr_table,c(t_metadict,"MetaDICT",Level.list[i]))
  
  t_conqur = t_test_res(t(res.ConQuR),meta$Y)
  fdr_table = rbind(fdr_table,c(t_conqur,"ConQuR",Level.list[i]))
}

#filename = paste0("./table/rd_fdr/","sim_", simu.iter, ".csv")
#write.csv(fdr_table, file = filename)

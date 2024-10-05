library(ConQuR)
library(doParallel)
library(MMUPHin)
source("./function.R")

load("../data/CRC_Duvallet.RData")

## seed 1-500
#jobid <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#simu.iter = jobid
#set.seed(jobid)

t_test_res = function(X,Y,alpha=0.1){
  X = X[rowSums(X)>0,]
  t_p = sapply(1:(nrow(X)), function(i)t.test(as.numeric(X[i,]) ~ Y)$p.value)
  p <- p.adjust(t_p, method='BH')
  return(mean(p<alpha))
}

p1 = c(0.5,0.5,0.5)
p2 = c(0.8,0.5,0.15)
p3 = c(0.99,0.5,0.01)
p_list = list(p1,p2,p3)

n_zeller = sum(meta_h$dataset == "Zeller")
n_zack = sum(meta_h$dataset == "Zackular")
n_baxter = sum(meta_h$dataset == "Baxter")

fdr_table = data.frame("FDR" = numeric(),"Method" = character(), "Level" = character())
Level.list = c("None","Moderate","Strong")

alpha = 0.001
beta = 0.01
gamma = 1

for(i in 1:3){
    
  p = p_list[[i]]
  meta_h$Y = 0
  meta_h$Y[which(meta_h$dataset=="Zeller")] = rbinom(n_zeller,1,p[1])
  meta_h$Y[which(meta_h$dataset=="Zackular")] = rbinom(n_zack,1,p[2])
  meta_h$Y[which(meta_h$dataset=="Baxter")] = rbinom(n_baxter,1,p[3])
  meta_h$Y = as.factor(meta_h$Y)
  
  batchid = as.factor(meta_h$dataset)
  
  t_raw = t_test_res(count.h,meta_h$Y)
  fdr_table = rbind(fdr_table,c(t_raw,"Unprocessed",Level.list[i]))
  
  res.ComBatSeq = sva::ComBat_seq(as.matrix(count.h),batchid,group = meta_h$Y)
  t_combatseq = t_test_res(res.ComBatSeq,meta_h$Y)
  fdr_table = rbind(fdr_table,c(t_combatseq,"ComBatSeq",Level.list[i]))
  
  meta_mmuphin = meta_h[,c(4,5)]
  res.mmuphin = adjust_batch(feature_abd = count.h,
                             batch = "dataset",
                             covariates = "Y",
                             data = meta_mmuphin)$feature_abd_adj
  t_mmuphin = t_test_res(res.mmuphin,meta_h$Y)
  fdr_table = rbind(fdr_table,c(t_mmuphin,"MMUPHin",Level.list[i]))
  
  otu_list = list()
  otu_list[[1]] = as.matrix(count.h[,meta_h$dataset=="Zeller"])
  otu_list[[2]] = as.matrix(count.h[,meta_h$dataset=="Baxter"])
  otu_list[[3]] = as.matrix(count.h[,meta_h$dataset=="Zackular"])
  
  meta_h1 = meta_h$Y
  meta_list = list()
  meta_list[[1]] = data.frame("Y" = meta_h1[meta_h$dataset=="Zeller"])
  meta_list[[2]] = data.frame("Y" = meta_h1[meta_h$dataset=="Baxter"])
  meta_list[[3]] = data.frame("Y" = meta_h1[meta_h$dataset=="Zackular"])
  
  metadict.res = metadict(otu_list,alpha,beta,gamma,dist_genus.sub,meta_list)
  X = metadict.res$X
  t_metadict = t_test_res(X,meta_h$Y)
  fdr_table = rbind(fdr_table,c(t_metadict,"MetaDICT",Level.list[i]))
  
  res.ConQuR = t(ConQuR(t(count.h),batchid =batchid, batch_ref = 'Baxter', covariates = meta_h$Y))
  t_conqur = t_test_res(res.ConQuR,meta_h$Y)
  fdr_table = rbind(fdr_table,c(t_conqur,"ConQuR",Level.list[i]))
}

#filename = paste0("./table/rd_fdr/","sim_", simu.iter, ".csv")
#write.csv(fdr_table, file = filename)

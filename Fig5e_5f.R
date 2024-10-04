library(ConQuR)
library(doParallel)
library(MMUPHin)
library(pROC)
source("./function.R")

## seed 1-500
#jobid <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#simu.iter = jobid
#set.seed(jobid)

count = read.csv("../data/count_order.csv") #change
count = count[,-1]

dist = read.csv("../data/dist_order.csv") #change
dist = as.matrix(dist[,-1])

knn_pred = function(train,test,Y_train,Y_test){
  test_pred <- class::knn(train = t(train), test = t(test), cl = Y_train, k=25)
  actual <- Y_test
  roc_object <- roc(actual,as.numeric(test_pred))
  return(auc(roc_object))
}

RF_pred = function(train,test,Y_train,Y_test){
  rf <- randomForest::randomForest(x = t(train), y = Y_train)
  test_pred <- predict(rf, newdata = t(test))
  actual <- Y_test
  roc_object <- roc(actual,as.numeric(test_pred))
  return(auc(roc_object))
}


init_dataset = function(m,n,count,dist,k,k1,k3=100,k4=30,ez=50,lib_var=T,p = c(1/2,1/2,1/2,1/2,1/2)){
  d = nrow(count)
  O_list = list()
  w_list = list()
  X_list = list()
  meta_list = list()
  meta_list1 = list()
  
  prevalence = rowSums(count!=0)
  diff_idx = order(prevalence,decreasing = T)[k:1]
  seq_idx1 = diff_idx
  
  seq_idx2 = sample(1:d,k1,replace = F)

  sample_info = c()

  A = exp(-dist/0.1)
  L = diag(colSums(A))-A
  svd_res = svd(L)
  U = svd_res$u
  U = U[,(ncol(U)-k3):(ncol(U))]
  
  for(i in 1:m){
    idx = sample(ncol(count),n)
    X = as.matrix(count[,idx])
    
    Y = sample(0:1,n,T,prob = c(p[i],1-p[i]))
    X[seq_idx1,which(Y==0)] = (X[seq_idx1,which(Y==0)]+1)*ez**2
    index = sample(k3,k4)
    w_space = U[,index]
    weight = 1-2*runif(k4)
    
    w = (w_space%*%as.matrix(weight))[,1]
    w = (w-min(w)+1e-3)
    w = w/max(w)
    w = w/mean(w)*0.5

    
    if(lib_var){
      if(i%%2==1){
        X_lib = X%*%diag(sample(1:5,n,replace = T))
      }else{
        X_lib = X
        }
    }else{
      X_lib = X
    }

    Y_f = as.factor(rbinom(n,1,0.5))
    O_list[i] = list(floor(diag(w)%*%X_lib))
    w_list[i] = list(w)
    X_list[i] = list(floor(X))
    meta_list[[i]] = data.frame("Y" = as.factor(Y),  "Y_f" = Y_f)
    meta_list1[[i]] = data.frame("Y_f" = Y_f)
    diff_idx = c(seq_idx1, seq_idx2)
  }
  return(list(O_list, w_list, X_list, diff_idx, meta_list, meta_list1))
}


m = 2
n = 100
p = c(1/2,1/2)
ez = 2
d = nrow(count)

alpha = 0.01
beta = 0.01
gamma = 1
table = c()


data = init_dataset(m,n,count,dist,floor(0.05*d),floor(0.05*d),ez=ez,p=p)
O_list = data[[1]]
meta_list = data[[5]]
meta_list1 = data[[6]]
meta = do.call("rbind",meta_list)
meta1 = do.call("rbind",meta_list1)

Y = as.factor(meta$Y)
Y_f = as.factor(meta$Y_f)
O = do.call(cbind,O_list)

split.ind = sample(1:ncol(O),floor(0.75*ncol(O)))
Y_train = Y[split.ind]
Y_test = Y[-split.ind]
Y_f_train = Y_f[split.ind]
Y_f_test = Y_f[-split.ind]

X = do.call("cbind",data[[3]])

batchid = as.factor(do.call(c,lapply(1:m, function(x)rep(x,n))))
res.metadict = metadict(O_list,alpha,beta,gamma,dist,meta_list)
X = res.metadict$X
res.metadict1 = metadict(O_list,alpha,beta,gamma,dist,meta_list1)
X1 = res.metadict1$X

res.ConQuR = ConQuR(t(O),batchid,batch_ref = 1,covariates = meta)
res.ConQuR1 = ConQuR(t(O),batchid,batch_ref = 1,covariates = meta1)

res.ComBatSeq = sva::ComBat_seq(as.matrix(O),batchid,covar_mod = meta)
res.ComBatSeq1 = sva::ComBat_seq(as.matrix(O),batchid,covar_mod = meta1)

colnames(O) = sapply(1:ncol(O),function(x)paste("Sample",x))
rownames(meta) = sapply(1:ncol(O),function(x)paste("Sample",x))
meta$batchid = batchid
res.mmuphin = adjust_batch(feature_abd = O,
                              batch = "batchid",
                              covariates = c("Y","Y_f"),
                              data = meta)$feature_abd_adj
res.mmuphin1 = adjust_batch(feature_abd = O,
                              batch = "batchid",
                              covariates = "Y_f",
                              data = meta)$feature_abd_adj

acc1 = RF_pred(X[,split.ind],X[,-split.ind],Y_train,Y_test)
acc12 = RF_pred(X1[,split.ind],X1[,-split.ind],Y_train,Y_test)
acc11 = RF_pred(X[,split.ind],X[,-split.ind],Y_f_train,Y_f_test)

table = rbind(table,c(acc1,"Random Forest","Sample","MetaDICT","With",ez))
table = rbind(table,c(acc11,"Random Forest","Negative Control","MetaDICT","With",ez))
table = rbind(table,c(acc12,"Random Forest","Sample","MetaDICT","Without",ez))

acc2 = RF_pred(O[,split.ind],O[,-split.ind],Y_train,Y_test)
acc21 = RF_pred(O[,split.ind],O[,-split.ind],Y_f_train,Y_f_test)

table = rbind(table,c(acc2,"Random Forest","Sample","Unprocessed","With",ez))
table = rbind(table,c(acc21,"Random Forest","Negative Control","Unprocessed","With",ez))
table = rbind(table,c(acc2,"Random Forest","Sample","Unprocessed","Without",ez))

acc3 = RF_pred(res.ComBatSeq[,split.ind],res.ComBatSeq[,-split.ind],Y_train,Y_test)
acc31 = RF_pred(res.ComBatSeq[,split.ind],res.ComBatSeq[,-split.ind],Y_f_train,Y_f_test)
acc32 = RF_pred(res.ComBatSeq1[,split.ind],res.ComBatSeq1[,-split.ind],Y_train,Y_test)

table = rbind(table,c(acc3,"Random Forest","Sample","ComBatSeq","With",ez))
table = rbind(table,c(acc31,"Random Forest","Negative Control","ComBatSeq","With",ez))
table = rbind(table,c(acc32,"Random Forest","Sample","ComBatSeq","Without",ez))


acc4 = RF_pred(res.mmuphin[,split.ind],res.mmuphin[,-split.ind],Y_train,Y_test)
acc41 = RF_pred(res.mmuphin[,split.ind],res.mmuphin[,-split.ind],Y_f_train,Y_f_test)
acc42 = RF_pred(res.mmuphin1[,split.ind],res.mmuphin1[,-split.ind],Y_train,Y_test)

table = rbind(table,c(acc4,"Random Forest","Sample","MMUPHin","With",ez))
table = rbind(table,c(acc41,"Random Forest","Negative Control","MMUPHin","With",ez))
table = rbind(table,c(acc42,"Random Forest","Sample","MMUPHin","Without",ez))

acc5 = RF_pred(t(res.ConQuR)[,split.ind],t(res.ConQuR)[,-split.ind],Y_train,Y_test)
acc51 = RF_pred(t(res.ConQuR)[,split.ind],t(res.ConQuR)[,-split.ind],Y_f_train,Y_f_test)
acc52 = RF_pred(t(res.ConQuR1)[,split.ind],t(res.ConQuR1)[,-split.ind],Y_train,Y_test)

table = rbind(table,c(acc5,"Random Forest","Sample","ConQuR","With",ez))
table = rbind(table,c(acc51,"Random Forest","Negative Control","ConQuR","With",ez))
table = rbind(table,c(acc52,"Random Forest","Sample","ConQuR","Without",ez))

acc1 = knn_pred(X[,split.ind],X[,-split.ind],Y_train,Y_test)
acc12 = knn_pred(X1[,split.ind],X1[,-split.ind],Y_train,Y_test)
acc11 = knn_pred(X[,split.ind],X[,-split.ind],Y_f_train,Y_f_test)

table = rbind(table,c(acc1,"K-NN","Sample","MetaDICT","With",ez))
table = rbind(table,c(acc11,"K-NN","Negative Control","MetaDICT","With",ez))
table = rbind(table,c(acc12,"K-NN","Sample","MetaDICT","Without",ez))

acc2 = knn_pred(O[,split.ind],O[,-split.ind],Y_train,Y_test)
acc21 = knn_pred(O[,split.ind],O[,-split.ind],Y_f_train,Y_f_test)

table = rbind(table,c(acc2,"K-NN","Sample","Unprocessed","With",ez))
table = rbind(table,c(acc21,"K-NN","Negative Control","Unprocessed","With",ez))
table = rbind(table,c(acc2,"K-NN","Sample","Unprocessed","Without",ez))

acc3 = knn_pred(res.ComBatSeq[,split.ind],res.ComBatSeq[,-split.ind],Y_train,Y_test)
acc31 = knn_pred(res.ComBatSeq[,split.ind],res.ComBatSeq[,-split.ind],Y_f_train,Y_f_test)
acc32 = knn_pred(res.ComBatSeq1[,split.ind],res.ComBatSeq1[,-split.ind],Y_train,Y_test)

table = rbind(table,c(acc3,"K-NN","Sample","ComBatSeq","With",ez))
table = rbind(table,c(acc31,"K-NN","Negative Control","ComBatSeq","With",ez))
table = rbind(table,c(acc32,"K-NN","Sample","ComBatSeq","Without",ez))


acc4 = knn_pred(res.mmuphin[,split.ind],res.mmuphin[,-split.ind],Y_train,Y_test)
acc41 = knn_pred(res.mmuphin[,split.ind],res.mmuphin[,-split.ind],Y_f_train,Y_f_test)
acc42 = knn_pred(res.mmuphin1[,split.ind],res.mmuphin1[,-split.ind],Y_train,Y_test)

table = rbind(table,c(acc4,"K-NN","Sample","MMUPHin","With",ez))
table = rbind(table,c(acc41,"K-NN","Negative Control","MMUPHin","With",ez))
table = rbind(table,c(acc42,"K-NN","Sample","MMUPHin","Without",ez))

acc5 = knn_pred(t(res.ConQuR)[,split.ind],t(res.ConQuR)[,-split.ind],Y_train,Y_test)
acc51 = knn_pred(t(res.ConQuR)[,split.ind],t(res.ConQuR)[,-split.ind],Y_f_train,Y_f_test)
acc52 = knn_pred(t(res.ConQuR1)[,split.ind],t(res.ConQuR1)[,-split.ind],Y_train,Y_test)

table = rbind(table,c(acc5,"K-NN","Sample","ConQuR","With",ez))
table = rbind(table,c(acc51,"K-NN","Negative Control","ConQuR","With",ez))
table = rbind(table,c(acc52,"K-NN","Sample","ConQuR","Without",ez))




#filename = paste0("./table/pred_sim1/","pred", "_simuiter", simu.iter, ".csv")
#write.csv(table, file = filename)

library(stats)
library(ecodist)
library(vegan)
library(viridis)
library(mclust)
library(ggplot2)

target_func <- function(x,O_list,alpha,beta,gamma,m,d,r,sample_num,L){
  res <- convert_from_vec(x,m,r,d,sample_num)
  w_list <- res[[1]]
  D <- res[[2]]
  R_list <- res[[3]]
  target <- 0
  for(i in 1:m){
    w <- w_list[i,]
    W <- diag(w)
    O_diff <- O_list[[i]]-W%*%D%*%R_list[[i]]
    W_diff <- t(w)%*%L%*%w
    target <- target+gamma*norm(O_diff,"F")**2+beta*W_diff/(d*d)+norm(R_list[[i]],"F")**2*alpha/(2*r*sample_num[i])
    
  }
  target <- target+norm(D,"F")**2*alpha/(2*d*r)
  return(target)
}

gradient_func <- function(x,O_list,alpha,beta,gamma,m,d,r,sample_num,L){
  res <- convert_from_vec(x,m,r,d,sample_num)
  w_list <- res[[1]]
  D <- res[[2]]
  R_list <- res[[3]]
  gradw <- matrix(0,m,d)
  gradD <- matrix(0,d,r)
  gradR <- list()
  for(i in 1:m){
    w <- w_list[i,]
    W <- diag(w)
    O_diff <- O_list[[i]]-W%*%D%*%R_list[[i]]
    W_diff <- t(w)%*%L%*%w
    gradD <- gradD+gamma*(-2*W%*%O_diff%*%t(R_list[[i]]))
    gradR[[i]] <- gamma*(-2)*t(D)%*%W%*%O_diff+R_list[[i]]*alpha/(r*sample_num[i])
    B <- gamma*(-2)*(O_diff)%*%t(R_list[[i]])%*%t(D)
    gradw[i,] <- diag(B)+2*t(w)%*%L*beta/(d*d)
  }
  gradD <- gradD+alpha*D/(d*r)
  grad <- convert_to_vec(m,gradw,gradD,gradR)
  return(grad)
}


metadict <- function(O_list,alpha,beta,gamma,dist,meta.list,r = NULL, normalization = "uq", max_iter=10000, neighbor=5, sigma=10, trace = 0){
  m <- length(O_list)
  d <- nrow(dist)

  if(is.null(r)){
      O <- do.call("cbind",O_list)
      r <- sum(svd(O)$d>1e-3)
  }
  
  sample_num <- sapply(O_list,ncol)
  if(normalization == "uq"){
      O_norm <- lapply(O_list,function(x)uq(x)$P)
  }else if(normalization == "rsim"){
      O_norm <- lapply(O_list,function(x)rsim(x)$P)
      O_norm <- laaply(O_norm, function(x)x/max(x))
  }else if(normalization == "tss"){
      O_norm <- lapply(O_list,function(x)t(t(x)/colSums(x)))
  }
  else if(!normalization){
      O_norm <- O_list
  }
  
  scale <- max(unlist(O_norm))
  O_list_test <- lapply(O_norm,function(x)x/scale)
  A <- matrix(0,d,d)
  for(i in 1:d){
    idx <- order(dist[i,],decreasing = F)[2:(neighbor+1)]
    A[i,idx] <- exp(-dist[i,idx]/sigma)
    A[idx,i] <- exp(-dist[i,idx]/sigma)
  }
  D1 <- diag(rowSums(A))
  L <- D1-A
  initial <- init_algo(m,d,O_list_test,r,meta.list)

  x0 <- convert_to_vec(m,initial[[1]],initial[[2]],initial[[3]])
  lower <- c(rep(0,m*d),rep(-Inf,d*r+r*sum(sample_num)))
  upper <- c(rep(1,m*d),rep(Inf,d*r+r*sum(sample_num)))
  optim.res <- optim(x0, fn = (function(x) target_func(x,O_list_test,alpha,beta,gamma,m,d,r,sample_num,L)), gr = (function(x) gradient_func(x,O_list_test,alpha,beta,gamma,m,d,r,sample_num,L)), method = "L-BFGS-B",
      lower = lower, upper = upper, control = list(maxit = max_iter, trace = trace))
  print(optim.res$convergence)
  x <- optim.res$par
  para <- convert_from_vec(x,m,r,d,sample_num)
  w_list <- para[[1]]
  D <- para[[2]]
  R_list <- para[[3]]
  X_list <- list()
  for(i in 1:m){
    X_list[[i]] <- D%*%R_list[[i]]*scale
  }
  X <- do.call(cbind,X_list)
  X[X<0] <- 0
  X[O==0] <- 0
  res <- list()
  res[["X"]] <- X
  
  var_res <- varimax(D, normalize = FALSE)
  D_rot <- var_res$loadings
  T_mat <- var_res$rotmat
  res[["D"]] <- D_rot
  res[["R"]] <- lapply(R_list,function(x)t(T_mat)%*%x*scale)
  res[["w"]] <- w_list
  effective_r <- effective_rank(D_rot)
  error_each <- sapply(1:length(O_norm), function(i) norm(O_norm[[i]]-diag(w_list[i,])%*%X_list[[i]], "F")^2/norm(O_norm[[i]],"F")^2)
  cat("Maximum Relative error:", max(error_each), "\n", "Effective rank of D", effective_r)
  return(res)
}

init_algo <- function(m,d,O_list,r,meta.list){
  w_list <- matrix(0,nrow = m,ncol = d)
  w_list[1,] <- 1
  for(i in 2:m){
    O1 <- O_list[[1]]
    O2 <- O_list[[i]]
    batchnum <- as.factor(c(rep(1,ncol(O1)),rep(i,ncol(O2))))
    if(is.null(meta.list)){
      meta <- data.frame("lib" = c(colSums(O1),colSums(O2)),"batch" = batchnum)
    }else{
      meta <- rbind(meta.list[[1]],meta.list[[i]])
      meta$batch <- batchnum
      meta$lib <- c(colSums(O1),colSums(O2))
    }
    mylogit <- glm(batch ~ ., data = meta, family = "binomial")
    meta$psvalue <- predict(mylogit, type="response")
    meta$weight <- ifelse(meta$batch==i,1/meta$psvalue,1/(1-meta$psvalue))
    O_adj_1 <- t(t(O1)*meta$weight[which(meta$batch==1)])
    O_adj_2 <- t(t(O2)*meta$weight[which(meta$batch==i)])
    w_list[i,] <- (rowMeans(O_adj_2)+1e-6)/(rowMeans(O_adj_1)+1e-6)
  }
  w_list <- t(t(w_list)/colSums(w_list))
  O_list1 <- lapply(1:m,function(i)diag(1/w_list[i,])%*%O_list[[i]])
  O <- do.call(cbind,O_list1)
  svd.res <- svd(O)
  D <- svd.res$u[,1:r]
  R_list <- lapply(1:m,function(i)t(D)%*%O_list1[[i]])
  initial_val <- list()
  initial_val[[1]] <- w_list
  initial_val[[2]] <- D
  initial_val[[3]] <- R_list
  return(initial_val)
}

convert_from_vec <- function(x,m,r,d,sample_num){
  x.w_list <- x[1:(d*m)]
  x.D <- x[(d*m+1):(d*m+d*r)]
  R_list <- list()
  n0 <- d*m+d*r
  for(i in 1:length(sample_num)){
    num <- sample_num[[i]]
    R_list[[i]] <- matrix(x[(n0+1):(n0+num*r)],r,num)
    n0 <- n0+num*r
  }
  
  D <- matrix(x.D,d,r)
  w_list <- matrix(x.w_list,m,d)
  return(list(w_list,D,R_list))
}

convert_to_vec <- function(m,w_list,D,R_list){
  x <- c(c(w_list),c(D),unlist(R_list))
  return(x)
}



effective_rank <- function(A) {
  svd_vals <- svd(A)$d  # Compute singular values
  
  # Normalize singular values to get probabilities
  p <- svd_vals / sum(svd_vals)
  
  # Compute entropy
  entropy <- -sum(p * log(p), na.rm = TRUE)
  
  # Compute effective rank
  exp(entropy)
}


##################################################Normalization

uq<-function(X){
    dds <- edgeR::calcNormFactors(as.matrix(X+1), method = "upperquartile")
    upperQ <- dds*colSums(X)
    upq.res <- scale(X,center=FALSE,scale=upperQ)
  return(list('P' = upq.res, 'sf' = upperQ))
}

CStat <- function(X){
  d <- nrow(X)
  R <- X
  S1 <- apply(R,1,order)
  S1 <- S1 - colMeans(S1);
  S1 <- S1 / sqrt(colSums(S1^2));
  corr_s <- crossprod(S1)
  med <- as.data.frame(matrixStats::colMedians(corr_s))
  return(as.numeric(med[,1]))

}

rsim <- function(X, eta = 0.1){
    d <- nrow(X)
    v <- CStat(X)
    I0.1 <- which(v>0.8)
    X0 <- X[I0.1,]
    v0 <- replicate(3,CStat(X0[sample(1:nrow(X0),0.5*nrow(X0)),]))
    w <- v[v>0.8]
    f1 <- sapply(w,function(x)mean(v>x))
    f0 <- sapply(w,function(x)mean(v0>x))
    pi <- sum(f1*f0)/sum(f0^2)
    vord <- order(v,decreasing = T)
    res <- sapply(1:length(vord),function(x)(1-pi*length(vord)*mean(v0>v[x])/(which(vord==x))))
    lowerx <- max(which(res[vord]<eta))
    ref <- vord[1:lowerx]
    tc.cn <- apply(X,2,function(x)sum(x[ref]))
    f.cn <- tc.cn/(mean(tc.cn))
    f.cn <- ifelse(f.cn==0,1,f.cn)
    cn.res <- scale(X,center=FALSE,scale=f.cn)
  return(list('P' = cn.res, 'I0' = ref, 'pi0'= pi, 'sf'=f.cn))
}

#################################################### help functions

########### Community Detection

cluster_core <- function(D,k,resolution,method = "Louvain"){
    knn.info <- RANN::nn2(D, k=k)
    knn <- knn.info$nn.idx
    adj <- matrix(0, nrow(D), nrow(D))
    for(i in seq_len(nrow(D))) {
        adj[i,knn[i,2:k]] <- 1
    }
    g <- igraph::graph.adjacency(adj, mode="undirected")
    g <- igraph::simplify(g)
    
    if(method == "Louvain"){
        km <- igraph::cluster_louvain(g,resolution = resolution)
    }else if(method == "Walktrap"){
        km <- igraph::cluster_walktrap(g)
    }
    return(list("membership"=km$membership,"graph"=g))
}

community_detection <- function(D, K = 50, resolution=1, tol = 1e-3, method = "Louvain", min_k = 5){
  avg_silwidth <- numeric()
  k.list <- c()
  for(k in min_k:K){
   membership <- cluster_core(D,k,resolution,method)$membership
   if(length(unique(membership))!=1){
        avg_silwidth<- c(avg_silwidth,mean(cluster::silhouette(membership, dist(D, method = "euclidean"))[,3], na.rm = TRUE))
        k.list <- c(k.list,k)
   }
  }
  best_k <- k.list[which.max(avg_silwidth)]
  best_res <- cluster_core(D,best_k,resolution,method)
  return(list("cluster" = best_res$membership, "graph" = best_res$graph))
}


############# PCoA plot
# return R2
permanova_pcoa <- function(distP, Y) {
  df.Y <- as.data.frame(Y)
  Re <- adonis2(distP~Y, data = df.Y)
  return(Re$R2[1])
}

tss <- function(X){
    return(t(t(X)/colSums(X)))
}

pcoa.plot.discrete <- function(A, gl, main, distance = "bray-curtis", colorset = "Set1", pointsize = 1){
    if(distance == "bray-curtis"){
        dist_matrix = bcdist(t(A))
    }else if(distance == "euclidean"){
        dist_matrix = dist(t(A))
    }else{
        stop("Selected distance is not supported.")
    }
  mds.stuff <- cmdscale(dist_matrix, eig=T, x.ret=T)
  mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100,1)
  mds.values <- mds.stuff$points
  mds.data <- data.frame( X=mds.values[,1],
                         Y=mds.values[,2],
                         Sample = gl)

  r2 <- permanova_pcoa(dist_matrix,gl)
  ggplot(data = mds.data, aes(x = X, y = Y, color = Sample)) +
  geom_point(size = pointsize) +
  xlab(paste("PCoA1: ", mds.var.per[1], '%', sep = "")) +
  ylab(paste("PCoA2: ", mds.var.per[2], '%', sep = "")) +
  labs(title = main,
       subtitle = paste("R2 =", round(r2, digits = 2))) +
  scale_color_brewer(palette = colorset) +
  theme_bw() +
  theme(
    #legend.key.height=unit(0.9,"cm"),
    legend.title = element_text(size = 14),  # Increase legend title size
    legend.text  = element_text(size = 14)
  )
}


pcoa.plot.continuous = function(A, gl, main, distance = "bray-curtis", pointsize = 1){
    if(distance == "bray-curtis"){
        dist_matrix = bcdist(t(A))
    }else if(distance == "euclidean"){
        dist_matrix = dist(t(A))
    }else{
        stop("Selected distance is not supported.")
    }
  mds.stuff <- cmdscale(dist_matrix, eig=T, x.ret=T)
  mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100,1)
  mds.values <- mds.stuff$points
  mds.data <- data.frame(X=mds.values[,1],
                         Y=mds.values[,2],
                         Signal = gl)
  r2 <- permanova_pcoa(dist_matrix,gl)
  ggplot(data=mds.data, aes(x=X, y=Y,col=Signal))+
    geom_point(size=pointsize)+
    scale_colour_gradientn(colors = viridis(10))+
    theme(legend.title=element_blank()) +
    xlab(paste("PCoA1: ", mds.var.per[1], '%', sep=""))+
    ylab(paste("PCoA2: ", mds.var.per[2], '%', sep=""))+
    labs(title = main,
              subtitle = paste("R2 =",round(r2, digits = 2)))+
    theme_bw()+
    theme(legend.key.height=unit(0.5,"cm"),
          legend.title = element_text(size = 14),  # Increase legend title size
          legend.text  = element_text(size = 12))
}

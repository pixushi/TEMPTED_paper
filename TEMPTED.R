library(ggplot2)
library(np)


#' Temporal Functional SVD
#' @param datlist A list of n matrices data, each matrix represents a subject;
#' the first row represents the sampling timestamp; 
#' following rows represent the feature measurements.
#' @param r CP-rank/Number of principle components. Default: 3.
#' @param interval The range of time points. Default: range of all observed data. 
#' @param resolution Resolution for the output eigencurves. Default: 100.
#' @param smooth Smoothing parameter for RKHS norm. Larger means smoother functions. Default: 1e-8.
#' @param maxiter Maximum number of itereation. Default: 20.
#' @param epsilon Convergence criteria for a and b. Default: 1e-4.
#' @return The estimations of the principle components in suject/feature/time;
#'       Var.prop: Explained variances.
tempted <- function(datlist, interval = NULL, r = 3, resolution = 251, smooth=1e-6,
                  maxiter=20, epsilon=1e-4){
  n = length(datlist)
  p = nrow(datlist[[1]])-1
  
  Lambda = rep(0, r)
  A = matrix(0, n, r)
  B = matrix(0, p, r)
  Phi = matrix(0, resolution, r)
  PCname <- paste('Component', 1:r)
  colnames(A) = PCname
  colnames(B) = PCname
  colnames(Phi) = PCname
  rownames(A) = names(datlist)
  rownames(B) = rownames(datlist[[1]])[-1]
  
  
  # Calculate range.
  timestamps.all = do.call(c,lapply(datlist, FUN=function(u){u[1,]}))
  
  timestamps.all = sort(unique(timestamps.all))
  if (is.null(interval)){
    interval = c(timestamps.all[1], timestamps.all[length(timestamps.all)])
  }
  
  # rescale the time to 0-1.
  input.time.range = c(timestamps.all[1], timestamps.all[length(timestamps.all)])
  for (i in 1:n){
    datlist[[i]][1,] = (datlist[[i]][1,] - input.time.range[1]) / (input.time.range[2] - input.time.range[1])
  }
  interval = (interval - input.time.range[1]) / (input.time.range[2] - input.time.range[1])
  
  res <- NULL
  Lambda <- rep(0, r)
  X <- NULL
  y0 <- NULL
  Rsq <- accumRsq <- rep(0, r)
  
  ti <- vector(mode='list', length=n)
  for (i in 1:n){
    temp <- 1 + round((resolution-1) * (datlist[[i]][1,] - interval[1]) / (interval[2] - interval[1]))
    temp[which(temp<=0 | temp>resolution)] <- 0
    ti[[i]] <- temp
  }
  
  tipos <- vector(mode='list', length=n)
  for (i in 1:n){
    keep <- ti[[i]]>0
    tipos[[i]] <- keep
    y0 <- c(y0, as.vector(t(datlist[[i]][2:(p+1),keep])))
  }
  
  Lt <- list()
  ind_vec <- NULL
  for (i in 1:n){
    Lt <- c(Lt, list(datlist[[i]][1,]))
    ind_vec <- c(ind_vec, rep(i,length(Lt[[i]])))
  }
  
  tm <- unlist(Lt)
  Kmat <- bernoulli_kernel(tm, tm)
  Kmat_output <- bernoulli_kernel(seq(interval[1],interval[2],length.out = resolution), tm)
  
  # calculate rank-1 component sequentially.  
  for (s in 1:r){ 
    # Step 1: initialization.
    print(sprintf("Calculate the %dth Component", s))
    
    # intialization of b
    data.unfold = NULL
    y <- NULL
    for (i in 1:n){
      data.unfold = cbind(data.unfold, datlist[[i]][2:(p+1),])
      y <- c(y, as.vector(t(datlist[[i]][2:(p+1),tipos[[i]]])))
    }
    b.initials <- svd(data.unfold, nu=r, nv=r)$u
    b.hat <- b.initials[,1]
    # initialization of a
    a.hat <- rep(1,n)/sqrt(n)
    
    # iteratively update a,b,phi
    t <- 0
    dif <- 1
    while(t<=maxiter & dif>epsilon){
      # update phi:
      Ly <- list()
      for (i in 1:n){
        Ly <- c(Ly, list(a.hat[i]*as.numeric(b.hat%*%datlist[[i]][2:(p+1),])))
      }
      phi.hat <- freg_rkhs(Ly, a.hat, ind_vec, Kmat, Kmat_output, smooth=smooth)
      phi.hat <- phi.hat / sqrt(sum(phi.hat^2))
      
      # update a:
      a.tilde <- rep(0,n)
      for (i in 1:n){
        t.temp <- tipos[[i]]
        a.tilde[i] <- b.hat %*% datlist[[i]][2:(p+1),t.temp] %*% phi.hat[ti[[i]][t.temp]] 
        a.tilde[i] <- a.tilde[i] / sum((phi.hat[ti[[i]][t.temp]])^2)
      }
      a.new <- a.tilde / sqrt(sum(a.tilde^2))
      dif <- sum((a.hat - a.new)^2)
      a.hat <- a.new
      
      # update b:
      temp.num <- matrix(0,p,n)
      temp.denom <- rep(0,n)
      for (i in 1:n){
        t.temp <- tipos[[i]]
        temp.num[,i] <- datlist[[i]][2:(p+1),t.temp] %*% phi.hat[ti[[i]][t.temp]]
        temp.denom[i] <-sum((phi.hat[ti[[i]][t.temp]])^2)
      }
      b.tilde <- as.numeric(temp.num%*%a.hat) / as.numeric(temp.denom%*%(a.hat^2))
      b.new <- b.tilde / sqrt(sum(b.tilde^2))
      dif <- max(dif, sum((b.hat - b.new)^2))
      b.hat <- b.new
      
      t <- t+1
    }
    
    # calculate lambda
    x <- NULL
    for (i in 1:n){
      t.temp <- ti[[i]]
      t.temp <- t.temp[t.temp>0]
      x <- c(x,as.vector(t(a.hat[i]*b.hat%o%phi.hat[t.temp])))
    }
    X <- cbind(X, x)
    l.fit <- lm(y~x-1)
    lambda <- as.numeric(l.fit$coefficients)
    A[,s] <- a.hat
    B[,s] <- b.hat
    Phi[,s] <- t(phi.hat)
    Lambda[s] <- lambda
    Rsq[s] <- summary(l.fit)$r.squared
    accumRsq[s] <- summary(lm(y0~X-1))$r.squared
    
    # update datlist
    for (i in 1:n){
      temp <- tipos[[i]]
      datlist[[i]][2:(p+1),which(temp)] <- datlist[[i]][2:(p+1),which(temp)] - 
        Lambda[s] * A[i,s] * (B[,s] %*% t(Phi[ti[[i]][temp],s])) 
    }
    print(paste0("Convergence reached at dif=", dif, ', iter=', t))
  }
  l.fit <- lm(y0~X-1)
  Lambda <- as.numeric(l.fit$coefficients)
  
  # revise the sign of Lambda
  for (r in 1:length(Lambda)){
    if (Lambda[r]<0){
      Lambda[r] <- -Lambda[r]
      A[,r] <- -A[,r]
    }
  }
  # revise the signs to make sure summation of phi is nonnegative
  sgn_phi <- sign(colSums(Phi))
  sgn_phi[sgn_phi==0] <- 1
  for (r in 1:ncol(Phi)){
    Phi[,r] <- sgn_phi[r]*Phi[,r]
    A[,r] <- sgn_phi[r]*A[,r]
  }
  # revise the signs to make sure summation of B is nonnegative
  sgn_B <- sign(colSums(B))
  sgn_B[sgn_B==0] <- 1
  for (r in 1:ncol(Phi)){
    B[,r] <- sgn_B[r]*B[,r]
    A[,r] <- sgn_B[r]*A[,r]
  }
  time.return <- seq(interval[1],interval[2],length.out = resolution)
  time.return <- time.return * (input.time.range[2] - input.time.range[1]) + input.time.range[1]
  results <- list("A.hat" = A, "B.hat" = B, 
                 "Phi.hat" = Phi, "time" = time.return,
                 "Lambda" = Lambda, "r.square" = Rsq, "accum.r.square" = accumRsq)
  return(results)
}



svd_centralize <- function(datlist, r = 1){
  n <- length(datlist)
  p <- nrow(datlist[[1]])-1
  mean_mat <- matrix(0,n,p)
  for (i in 1:length(datlist)){
    mean_mat[i,] <- apply(datlist[[i]][-1,], 1, mean)
  }
  mean_mat.svd <- svd(mean_mat, nu=r, nv=r)
  mean_mat.svd1 <- mean_mat.svd$u %*% t(mean_mat.svd$v * mean_mat.svd$d[1:r])
  mf.new <- datlist
  for (i in 1:length(datlist)){
    mf.new[[i]][-1,] <- datlist[[i]][-1,] - mean_mat.svd1[i,]
  }
  results <- list("datlist" = mf.new, "A.tilde" = mean_mat.svd$u, 
                  "B.tilde" = mean_mat.svd$v, "lambda.tilde" = mean_mat.svd$d[1:r])
  return(results)
}


#' Format data table into input of tempted
#' @param taxon_table A table of read counts, with n rows for samples and p columns for taxa.
#' @param time_point The time stamp of each sample, relative to the start of the study. 
#' A length n vector.
#' @param subjectID The subject ID of each sample. A length n vector.
#' @param threshold A threshold for taxon filtering. 
#' Taxa with zero counts percentage >= threshold will be excluded.
#' @param pseudo_count A small number to add to all the counts before 
#' normalizing into proportions and log transformation.
#' @return Input for tempted. A list of matrices, each representing the 
format_tempted <- function(taxon_table, time_point, subjectID, threshold=0.95, 
                         feature_names=NULL, 
                         pseudo_count=0.5, transform="clr"){
  # format data table into a list as input for tempted(), 
  # read counts all have 1/2 added, before being normalized into proportions and log transformation
  # check length of subID and time_point
  ntm <- which(table(subjectID)==1)
  if(length(ntm)>0)
    stop(paste('Please remove these subjects with only one time point:', 
               paste(names(ntm), collapse=', ')))
  if (length(subjectID)!=nrow(taxon_table)) 
    stop('length of subjectID does not match taxon_table!')
  if (length(time_point)!=nrow(taxon_table)) 
    stop('length of time_point does not match taxon_table!')
  # keep taxon that has non-zeros in >=1-threshold samples
  if (is.null(feature_names)){
    taxon_table <- taxon_table[,colMeans(taxon_table==0)<threshold]
  }else{
    taxon_table <- taxon_table[,feature_names]
  }
  if(transform=='log_comp'){
    taxon_table <- taxon_table+pseudo_count
    taxon_table <- t(log(taxon_table/rowSums(taxon_table)))
  }else if(transform=='comp'){
    taxon_table <- taxon_table
    taxon_table <- t(taxon_table/rowSums(taxon_table))
  }else if(transform=='ast'){
    taxon_table <- taxon_table
    taxon_table <- t(asin(sqrt(taxon_table/rowSums(taxon_table))))
  }else if(transform=='clr'){
    taxon_table <- taxon_table+pseudo_count
    taxon_table <- log(taxon_table/rowSums(taxon_table))
    taxon_table <- t(taxon_table-rowMeans(taxon_table))
  }else if(transform=='logit'){
    taxon_table <- taxon_table+pseudo_count
    taxon_table <- t(taxon_table/rowSums(taxon_table))
    taxon_table <- log(taxon_table/(1-taxon_table))
  }else if(transform=='none'){
    taxon_table <- t(taxon_table)
  }else{
    print('Input transformation method is wrong! log_comp is applied instead')
    taxon_table <- taxon_table+pseudo_count
    taxon_table <- t(log(taxon_table/rowSums(taxon_table)))
  }
  taxon_table <- rbind(time_point, taxon_table)
  rownames(taxon_table)[1] <- 'time_point'
  subID <- unique(subjectID)
  nsub <- length(subID)
  
  # construct list of data matrices, each element representing one subject
  datlist <- vector("list", length = nsub)
  names(datlist) <- subID
  
  # Each slice represents an individual (unequal sized matrix).
  for (i in 1:nsub){
    # print(i)
    datlist[[i]] <- taxon_table[, subjectID==subID[i]]
    datlist[[i]] <- datlist[[i]][,order(datlist[[i]][1,])]
    datlist[[i]] <- datlist[[i]][,!duplicated(datlist[[i]][1,])]
  }
  return(datlist)
}




##########################
# RKHS functional regression
#########################
bernoulli_kernel <- function(x, y){
  k1.x <- x-0.5
  k1.y <- y-0.5
  k2.x <- 0.5*(k1.x^2-1/12)
  k2.y <- 0.5*(k1.y^2-1/12)
  xy <- abs(x %*% t(rep(1,length(y))) - rep(1,length(x)) %*% t(y))
  k4.xy <- 1/24 * ((xy-0.5)^4 - 0.5*(xy-0.5)^2 + 7/240)
  kern.xy <- k1.x %*% t(k1.y) + k2.x %*% t(k2.y) - k4.xy + 1
  return(kern.xy)
}

freg_rkhs <- function(Ly, a.hat, ind_vec, Kmat, Kmat_output, smooth=1e-8){
  A <- Kmat
  for (i in 1:length(Ly)){
    A[ind_vec==i,] <- A[ind_vec==i,]*a.hat[i]^2
  }
  cvec <- unlist(Ly)
  
  A.temp <- A + smooth*diag(ncol(A))
  beta <- solve(A.temp)%*%cvec
  
  phi.est <- Kmat_output %*% beta
  return(phi.est)
}

####################

#' Multiply loadings from res_tempted into denoised tensor
#' @param res_tempted Output of tempted
#' @param mean_svd Output of svd_centralize
#' @return The denoised functional tensor
tdenoise <- function(res_tempted, mean_svd=NULL){
  n <- nrow(res_tempted$A.hat)
  p <- nrow(res_tempted$B.hat)
  resol <- nrow(res_tempted$Phi.hat)
  tensor.est <- array(0,dim=c(n,p,resol))
  if (!is.null(mean_svd))
    tensor.est <- (mean_svd$A.tilde %*% t(mean_svd$B.tilde * mean_svd$lambda.tilde)) %o%
    rep(1, resol)
  for (i in 1:ncol(res_tempted$A.hat)){
    tensor.est <- tensor.est+res_tempted$A.hat[,i]%o%res_tempted$B.hat[,i]%o%res_tempted$Phi.hat[,i]*res_tempted$Lambda[i]
  }
  dimnames(tensor.est)[[3]] <- res_tempted$time
  return(tensor.est)
}




plot_time_loading <- function(res, r=NULL, ...){
  Phi.data <- res$Phi.hat
  if(is.null(r)) r <- ncol(Phi.data)
  Phi.data <- Phi.data[,1:r]
  ntime <- nrow(Phi.data)
  Phi.data <- data.frame(time=res$time, value=as.vector(Phi.data), 
                         component=as.factor(as.vector(t(matrix(rep(1:r,ntime),r,)))))
  ptime <- ggplot(data=Phi.data, aes(x=time, y=value, color=component)) + geom_line(aes(...))
  return(ptime)
}

# aggregate features using feature loadings
#' @param res_tempted output of tempted()
#' @param res_svd output of svd_centralize()
#' @param datlist output of format_tempted()
#' @param pct the percent of features to aggregate, features ranked by abs(feature loading)
#' @param contrast a matrix choosing how components are combined, 
#' each column is a contrast of length r
#' @param get_contrast a vector denoting which components to run SVD on
#' the loading will become part of the contrast
aggregate_feature <- function(res_tempted, res_svd=NULL, datlist, 
                               pct=1, 
                              contrast=NULL, get_contrast=NULL){
  B.data <- as.data.frame(res_tempted$B.hat)
  r <- ncol(B.data)
  if (!is.null(get_contrast)){
    datlist.agg <- sapply(datlist, function(x){t(B.data)%*%x[-1,]}, simplify=F)
    metafeature.obs <- NULL
    for (i in 1:length(datlist.agg)){
      tmp <- data.frame(value=as.vector(datlist.agg[[i]]), 
                        subID=names(datlist.agg)[i],
                        timepoint=as.vector(t(matrix(datlist[[i]][1,], ncol(datlist[[i]]), ncol(B.data)))),
                        PC=rep(rownames(datlist.agg[[i]]), ncol(datlist.agg[[i]])))
      
      metafeature.obs <- rbind(metafeature.obs, tmp)
    }
    metafeature.obs.wide <- reshape(metafeature.obs, timevar="PC", 
                                    idvar=c("subID", "timepoint"), 
                                    direction="wide")
    colnames(metafeature.obs.wide)[-c(1,2)] <- substring(colnames(metafeature.obs.wide)[-c(1,2)],7)
    scale_mat <- scale(metafeature.obs.wide[-c(1,2)][get_contrast])
    component_svd <- svd(scale_mat)
    contrast_svd <- matrix(0, r, ncol(component_svd$v))
    contrast_svd[get_contrast,] <- component_svd$v
    contrast <- cbind(contrast, contrast_svd)
  }
  if(!is.null(contrast)){
    contrast.data <- res_tempted$B.hat%*%contrast
    colnames(contrast.data) <- paste('Contrast', 1:ncol(contrast))
    B.data <- cbind(B.data, contrast.data)
  }
  toppct <- apply(abs(B.data), 2, function(x){x>quantile(x, 1-pct)})
  datlist.agg <- sapply(datlist, function(x){t(B.data*toppct)%*%x[-1,]}, simplify=F)
  metafeature.obs <- NULL
  for (i in 1:length(datlist.agg)){
    tmp <- data.frame(value=as.vector(datlist.agg[[i]]), 
                      subID=names(datlist.agg)[i],
                      timepoint=as.vector(t(matrix(datlist[[i]][1,], ncol(datlist[[i]]), ncol(B.data)))),
                      PC=rep(rownames(datlist.agg[[i]]), ncol(datlist.agg[[i]])))
    
    metafeature.obs <- rbind(metafeature.obs, tmp)
  }
  metafeature.obs <- metafeature.obs[,c("value", "subID", "timepoint", "PC")]
  metafeature.obs$type <- 'observed'
  
  # estimated
  tensor.est <- tdenoise(res_tempted, res_svd)
  tensor.est.agg <- apply(tensor.est, c(1,3), function(x){(t(B.data*toppct)%*%x)})
  metafeature.est <- NULL
  for (i in 1:r){ 
    tmp <- data.frame(value=as.vector(tensor.est.agg[i,,]),
                      subID=rep(dimnames(tensor.est.agg)[[2]], dim(tensor.est.agg)[3]),
                      timepoint=as.vector(t(matrix(res_tempted$time,length(res_tempted$time),dim(tensor.est.agg)[2]))),
                      PC=colnames(B.data)[i])
    metafeature.est <- rbind(metafeature.est,tmp)
  }
  metafeature.est <- metafeature.est[,c("value", "subID", "timepoint", "PC")]
  metafeature.est$type <- 'estimated'
  
  return(list(metafeature.obs=metafeature.obs, 
              metafeature.est=metafeature.est,
              contrast=contrast,
              toppct=toppct))
}



# take log ratio of the total abundance of top features over bottom features
#' @param res_tempted output of tempted()
#' @param datlist output of format_tempted() with option transform="none"
#' @param pct the percent of features to sum up, features ranked by abs(feature loading)
#' @param contrast a matrix choosing how components are combined, 
#' each column is a contrast of length r
#' @param get_contrast a vector denoting which components to run SVD on
#' the loading will become part of the contrast
ratio_feature <- function(res_tempted, datlist, 
                              pct=0.05, absolute=FALSE, contrast=NULL){
  B.data <- as.data.frame(res_tempted$B.hat)
  if (!is.null(contrast)){
    contrast.data <- res_tempted$B.hat%*%contrast
    colnames(contrast.data) <- paste('Contrast', 1:ncol(contrast))
    B.data <- cbind(B.data, contrast.data)
  }
  if(!absolute){
    toppct <- apply(B.data, 2, function(x){x>quantile(x, 1-pct) & x>0})
    bottompct <- apply(-B.data, 2, function(x){x>quantile(x, 1-pct) & x>0})
  }else{
    toppct <- apply(B.data, 2, function(x){abs(x)>quantile(abs(x), 1-pct) & x>0})
    bottompct <- apply(B.data, 2, function(x){abs(x)>quantile(abs(x), 1-pct) & x<0})
  }
  
  datlist.ratio <- sapply(datlist, 
                          function(x){
                            tt <- (t(toppct)%*%x[-1,])
                            bb <- (t(bottompct)%*%x[-1,])
                            return(log((tt+0.5)/(bb+0.5)))},
                          simplify=F)
  metafeature.ratio <- NULL
  for (i in 1:length(datlist.ratio)){
    tmp <- data.frame(value=as.vector(datlist.ratio[[i]]), 
                      subID=names(datlist.ratio)[i],
                      timepoint=as.vector(t(matrix(datlist[[i]][1,], ncol(datlist[[i]]), ncol(B.data)))),
                      PC=rep(rownames(datlist.ratio[[i]]), ncol(datlist.ratio[[i]])))
    
    metafeature.ratio <- rbind(metafeature.ratio, tmp)
  }
  metafeature.ratio <- metafeature.ratio[,c("value", "subID", "timepoint", "PC")]
  return(list(metafeature.ratio=metafeature.ratio, 
              contrast=contrast,
              toppct=toppct, bottompct=bottompct))
}



# plot nonparametric smoothing lines and error band
#' @param feature_mat sample by feature matrix
#' @param time_vec time vector correspond to samples
#' @param group_vec factor variable to group samples
#' @param nrow number of rows to wrap the figures
plot_feature_summary <- function(feature_mat, time_vec, group_vec, 
                                   coverage=0.95, bws=NULL, nrow=1){
  nfeature <- ncol(feature_mat)
  if(class(group_vec)!='factor') group_vec <- as.factor(group_vec)
  group_level <- levels(group_vec)
  time_all <- NULL
  mean_all <- NULL
  merr_all <- NULL
  feature_all <- NULL
  group_all <- NULL
  if(is.null(colnames(feature_mat))){stop('feature_mat needs to have column names!')}
  CI_length <- -qnorm((1-coverage)/2)
  for (jj in 1:nfeature){
    for (ii in 1:length(group_level)){
      ind <- group_vec==group_level[ii]
      if(is.null(bws)){
        model.np <- npreg(feature_mat[ind,jj]~time_vec[ind],
                          regtyle="ll", bwmethod="cv.aic", gradients=T)
      }else{
        model.np <- npreg(feature_mat[ind,jj]~time_vec[ind], bws=bws,
                          regtyle="ll", bwmethod="cv.aic", gradients=T)
      }
      time_eval <- as.vector(t(model.np$eval))
      mean_eval <- model.np$mean[order(time_eval)]
      merr_eval <- model.np$merr[order(time_eval)]
      time_eval <- sort(time_eval)
      
      time_all <- c(time_all, time_eval)
      mean_all <- c(mean_all, mean_eval)
      merr_all <- c(merr_all, merr_eval)
      feature_all <- c(feature_all, 
                       rep(colnames(feature_mat)[jj], length(time_eval)))
      group_all <- c(group_all,
                     rep(group_level[ii], length(time_eval)))
    }
  }
  group_all <- factor(group_all, levels=group_level)
  tab_summary <- data.frame(time=time_all, mean=mean_all, merr=merr_all,
                            group=group_all, feature=feature_all)
  
  p_summary <- ggplot(data=tab_summary, 
                      aes(x=time, y=mean, group=group, color=group)) + 
    geom_line() + 
    geom_ribbon(aes(ymin=mean-CI_length*merr, ymax=mean+CI_length*merr, 
                    color=group, fill=group), linetype=2, alpha=0.3) + 
    ylab(paste0('mean +/- ', round(CI_length,2), '*se')) + facet_wrap(.~feature, scales="free", nrow=nrow)
  return(p_summary)
}

tabular_feature_est <- function(tensor.est, feature_sel){
  tensor.est <- tensor.est[,feature_sel,]
  tab_est <- NULL
  tm <- as.numeric(dimnames(tensor.est)[[3]])
  for (i in 1:dim(tensor.est)[1]){
    tmp <- data.frame(subID=rownames(tensor.est)[i],
                      time=as.vector(t(matrix(tm, length(tm), ncol(tensor.est)))),
                      feature=rep(colnames(tensor.est), length(time)),
                      value=as.vector(tensor.est[i,,]))
    tab_est <- rbind(tab_est, tmp)
  }
  tab_est$type <- 'estimated'
  return(tab_est)
}


tabular_feature_obs <- function(datlist, feature_sel){
  tab_obs <- NULL
  for (i in 1:length(feature_sel)){
    value <- unlist(sapply(datlist, function(x){x[feature_sel[i],]}, simplify=F))
    time_point <- unlist(sapply(datlist, function(x){x['time_point',]}, simplify=F))
    nobs <- sapply(datlist, function(x){ncol(x)}, simplify=F)
    subID <- unlist(mapply(function(i){rep(names(datlist)[i], nobs[i])}, 
                           1:length(nobs), SIMPLIFY=F))
    tmp <- data.frame(subID=subID, time=time_point, feature=feature_sel[i], value=value)
    rownames(tmp) <- NULL
    tab_obs <- rbind(tab_obs, tmp)
  }
  tab_obs$type <- 'observed'
  return(tab_obs)
}




#' Calculate ROC using logistic regression
#' @param clust the observed labels
#' @param Xdata the predictors
#' @return output of function roc 
roc_logistic <- function(xtrain, ytrain, xtest, ytest){
  dftrain <- data.frame(y=ytrain, x=xtrain)
  dftest <- data.frame(y=ytest, x=xtest)
  fit <- glm(y ~ ., data = dftrain, family = "binomial")
  prob <- predict(fit, newdata=dftest, type = c("response"))
  g <- roc(dftest$y~prob)
  return(g$auc)
}

#' Estimate A of testing data based on B and Phi from training data
#' @param datlist testing data
#' @param res_tempted tempted result from training data
#' @return estimated A of testing data
est_A <- function(datlist, res_tempted, mean_svd=NULL){
  B <- res_tempted$B.hat
  Phi <- res_tempted$Phi.hat
  Lambda <- res_tempted$Lambda
  time.return <- res_tempted$time
  n <- length(datlist)
  p <- nrow(B)
  r <- ncol(B)
  resolution <- length(time.return)
  A_test <- matrix(0,n,r)
  y <- NULL
  ti <- vector(mode = "list", length = n)
  # get the coordinate of observed time points in the returned time grid
  for (i in 1:n){
    ti[[i]] <- sapply(datlist[[i]][1,], function(x){which.min(abs(x-time.return))})
    y <- c(y, as.numeric(t(datlist[[i]][-1,ti[[i]]>0])))
  }
  mf.new <- datlist
  if(!is.null(mean_svd)){
    mean_mat <- matrix(0,n,p)
    for (i in 1:length(datlist)){
      mean_mat[i,] <- apply(datlist[[i]][-1,], 1, mean)
    }
    mean_mat.svd1 <- mean_mat%*%tcrossprod(mean_svd$B.tilde)
    mf.new <- datlist
    for (i in 1:length(datlist)){
      mf.new[[i]][-1,] <- datlist[[i]][-1,] - mean_mat.svd1[i,]
    }
  }

  for (s in 1:r){
    for (i in 1:n){
      t.temp <- ti[[i]]>0
      A_test[i,s] <- B[,s] %*% mf.new[[i]][2:(p+1),t.temp] %*% Phi[ti[[i]][t.temp],s] 
      A_test[i,s] <- A_test[i,s] / sum((Phi[ti[[i]][t.temp],s])^2) / Lambda[s]
      mf.new[[i]][2:(p+1),t.temp] <- mf.new[[i]][2:(p+1),t.temp] - 
        Lambda[s] * A_test[i,s] * (B[,s] %*% t(Phi[ti[[i]][t.temp],s])) 
    }
  }
  rownames(A_test) <- names(mf.new)
  colnames(A_test) <- paste('Component', 1:r)
  return(A_test)
}

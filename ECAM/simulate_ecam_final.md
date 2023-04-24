---
title: "Simulation Based on ECAM Data"
author: "Pixu Shi"
date: "2023-04-24"
output:
  html_document:
    toc: true
    theme: united
---

# Preparation

## Library


```r
rm(list=ls())
# for data
library(qiime2R) # read in Qiime artifacts
library(dplyr) # data formatting
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
# for computing
library(reticulate) # run py codes
library(vegan) # distance matrix
```

```
## Loading required package: permute
```

```
## Loading required package: lattice
```

```
## This is vegan 2.6-4
```

```r
library(PERMANOVA) # permanova
```

```
## Loading required package: Matrix
```

```
## Loading required package: xtable
```

```
## Loading required package: MASS
```

```
## 
## Attaching package: 'MASS'
```

```
## The following object is masked from 'package:dplyr':
## 
##     select
```

```
## Loading required package: scales
```

```
## Loading required package: deldir
```

```
## deldir 1.0-6      Nickname: "Mendacious Cosmonaut"
```

```
## 
##      The syntax of deldir() has had an important change. 
##      The arguments have been re-ordered (the first three 
##      are now "x, y, z") and some arguments have been 
##      eliminated.  The handling of the z ("tags") 
##      argument has been improved.
##  
##      The "dummy points" facility has been removed. 
##      This facility was a historical artefact, was really 
##      of no use to anyone, and had hung around much too 
##      long.  Since there are no longer any "dummy points", 
##      the structure of the value returned by deldir() has 
##      changed slightly.  The arguments of plot.deldir() 
##      have been adjusted accordingly; e.g. the character 
##      string "wpoints" ("which points") has been 
##      replaced by the logical scalar "showpoints". 
##      The user should consult the help files.
```

```
## 
## Attaching package: 'deldir'
```

```
## The following object is masked from 'package:permute':
## 
##     getCol
```

```r
library(randomForest) # random forest
```

```
## randomForest 4.7-1.1
```

```
## Type rfNews() to see new features/changes/bug fixes.
```

```
## 
## Attaching package: 'randomForest'
```

```
## The following object is masked from 'package:dplyr':
## 
##     combine
```

```r
library(PRROC) # roc and pr

# for plotting
library(ggpubr)
```

```
## Loading required package: ggplot2
```

```
## 
## Attaching package: 'ggplot2'
```

```
## The following object is masked from 'package:randomForest':
## 
##     margin
```

```
## 
## Attaching package: 'ggpubr'
```

```
## The following object is masked from 'package:qiime2R':
## 
##     mean_sd
```

```r
library(ggplot2)
library(gridExtra)
```

```
## 
## Attaching package: 'gridExtra'
```

```
## The following object is masked from 'package:randomForest':
## 
##     combine
```

```
## The following object is masked from 'package:dplyr':
## 
##     combine
```

```r
library(RColorBrewer)


# set working directory to be where the current script is located
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
color_RB <- brewer.pal(3,'Set1')[1:2]
source('../TEMPTED.R')
```

```
## Nonparametric Kernel Methods for Mixed Datatypes (version 0.60-16)
## [vignette("np_faq",package="np") provides answers to frequently asked questions]
## [vignette("np",package="np") an overview]
## [vignette("entropy_np",package="np") an overview of entropy-based methods]
```

## Read in the data


```r
count_tab <- read.csv("data/otu_count_cleaned_q2.csv", row.names=1)
meta_tab <- read.csv("data/otu_metadata_cleaned_q2.csv", row.names=1)
taxon_tab <- read.csv("data/otu_taxonomy_cleaned_q2.csv", row.name=1)
table(rownames(count_tab)==rownames(meta_tab))
```

```
## 
## TRUE 
##  694
```

```r
# remove distal gut samples
ind <- meta_tab$qiita_empo_3=='anthropogenic sample' | is.na(meta_tab$qiita_empo_3)
count_tab <- count_tab[ind,]
meta_tab <- meta_tab[ind,]
meta_tab$delivery_ind <- 'Vaginal'==meta_tab$delivery
table(rownames(count_tab)==rownames(meta_tab))
```

```
## 
## TRUE 
##  694
```

```r
metauni <- unique(meta_tab[,c('studyid', 'delivery_ind', 'diet')])
rownames(metauni) <- metauni$studyid
```

## Plot timeline of study


```r
meta_reordered <- meta_tab[order(meta_tab$delivery, meta_tab$studyid, decreasing=T),]
meta_reordered$studyid <- factor(meta_reordered$studyid, 
                                 unique(meta_reordered$studyid))
colnames(meta_reordered)[10] <- "Delivery"
p_timeline <- ggplot(data=meta_reordered, 
       aes(x=day_of_life, y=studyid, 
           group=studyid, color=Delivery, shape=Delivery)) +
  geom_line() + geom_point() + scale_color_manual(values=color_RB) +
  labs(y="Host ID", x="Day of Life") + 
  theme_bw() +
  theme(legend.position="bottom")
p_timeline
```

![plot of chunk timeline](figure/timeline-1.png)

```r
pdf('../figure_table/ECAM_timeline.pdf', width=7.7, height=5)
print(p_timeline)
dev.off()
```

```
## quartz_off_screen 
##                 2
```

# Simulated data by resampling ECAM samples


```r
nkeep <- c(2,3,4,5,6,7,8,9,10)
nsim <- 100
set.seed(0)
```


```r
sampleID <- matrix(list(NULL), nsim, length(nkeep))
for (ss in 1:nsim){
  print(ss)
  for (jj in 1:length(nkeep)){
    sampleID_temp <- NULL
    for (i in 1:nrow(metauni)){
      meta_sub <- dplyr::filter(meta_tab, studyid==metauni$studyid[i])
      meta_sub <- dplyr::arrange(meta_sub, -day_of_life)
      meta_sub <- meta_sub[!duplicated(meta_sub$day_of_life),]
      meta_sub <- dplyr::arrange(meta_sub, day_of_life)
      nsample1 <- sum(meta_sub$day_of_life<=365)
      nsample2 <- sum(meta_sub$day_of_life>365)
      prob_vec <- c(rep(1,nsample1), rep(2,nsample2))
      prob_vec <- prob_vec/sum(prob_vec)
      sel <- sample(rownames(meta_sub), 
                    size=min(nsample1+nsample2, nkeep[jj]),
                    prob=prob_vec)
      sel <- sel[order(meta_sub[sel,'day_of_life'])]
      sampleID_temp <- c(sampleID_temp, sel)
    }
    sampleID[ss,jj] <- list(sampleID_temp)
  }
}
colnames(sampleID) <- nkeep
save(sampleID, file="simdata/realsim_ecam_sampleID.Rdata")
```


```r
load(file.path("simdata", "realsim_ecam_sampleID.Rdata"))
nkeep <- as.numeric(colnames(sampleID))
nsim <- nrow(sampleID)
npc <- 2 
```

# TEMPTED

## Run TEMPTED, save subject loading and trajectory


```r
for (jj in 1:length(nkeep)){
  for (ss in 1:nsim){
    count_sub <- count_tab[sampleID[ss,jj][[1]],]
    meta_sub <- meta_tab[sampleID[ss,jj][[1]],c("studyid", "day_of_life", "month", "delivery")]
    meta_sub$studyid <- as.character(meta_sub$studyid)
    meta_sub$month <- as.numeric(meta_sub$month)
    meta_sub$delivery_ind <- meta_sub$delivery=="Vaginal"
    metauni <- unique(meta_sub[,c("studyid", "delivery_ind")])
    subdata <- format_tempted(count_sub, meta_sub$day_of_life, meta_sub$studyid, 
                            threshold=0.95, pseudo_count=0.5,
                            transform="clr")
    # run tempted with all subjects and run permanova test
    svd_sub <- svd_centralize(subdata)
    res_tempted <- tempted(svd_sub$datlist, r = npc, resolution = 101, smooth=1e-4)
    agg_feat <- aggregate_feature(res_tempted, NULL, subdata)
    aggfeat_mat <- reshape(agg_feat$metafeature.obs[,1:4], 
                           idvar=c("subID","timepoint") , 
                           v.names=c("value"), timevar="PC",
                           direction="wide")
    colnames(aggfeat_mat)[1] <- 'studyid'
    aggfeat_mat <- merge(aggfeat_mat, metauni, by='studyid')
    aggfeat_mat[,2+1:npc] <- apply(aggfeat_mat[,2+1:npc], 2, scale)
    tempted_sub <- as.data.frame(res_tempted$A.hat)
    tempted_sub$`intercept` <- svd_sub$A.tilde
    tempted_sub$studyid <- rownames(tempted_sub)
    tempted_sub <- merge(metauni, tempted_sub)
    write.csv(aggfeat_mat, 
              file=paste0('simresult/tempted_traj_sim',ss,'_ntime',nkeep[jj], '.csv'))
    write.csv(tempted_sub, 
              file=paste0('simresult/tempted_subj_sim',ss,'_ntime',nkeep[jj], '.csv'))
  }
}
```

## Sample-level

### PERMANOVA F value


```r
Fmodel_tempted <- matrix(NA, nsim, length(nkeep))
colnames(Fmodel_tempted) <- paste0("nsample", nkeep)
for (jj in 1:length(nkeep)){
  print(jj)
  for (ss in 1:nsim){
    fname <- paste0('tempted_traj_sim',ss,'_ntime',nkeep[jj], '.csv')
    aggfeat_mat <- read.csv(file.path("simresult", fname), row.names=1)
    # calculate PERMANOVA F
    dist_aggft <- vegdist(aggfeat_mat[,2+1:npc], method='euclidean')
    res_perm <- adonis2(dist_aggft ~ aggfeat_mat$delivery)
    Fmodel_tempted[ss,jj] <- res_perm$F[1]
  }
}

write.csv(Fmodel_tempted, 
          file='result/realsim_ecam_Fvalue_tempted_clr.csv')
```

### In-sample ROC & PR

-   logistic regression
-   random forest


```r
roc_sample_glm_tempted <- matrix(NA, nsim, length(nkeep))
colnames(roc_sample_glm_tempted) <- paste0("nsample", nkeep)
pr_sample_rf_tempted <- roc_sample_rf_tempted <- pr_sample_glm_tempted <- roc_sample_glm_tempted 
for (jj in 1:length(nkeep)){
  print(jj)
  for (ss in 1:nsim){
    fname <- paste0('tempted_traj_sim',ss,'_ntime',nkeep[jj], '.csv')
    aggfeat_mat <- read.csv(file.path("simresult", fname), row.names=1)
    # logistic regression, even though dimension reduction has all training samples, prediction is leave-one-out
    predprob_glm <- rep(NA, nrow(aggfeat_mat))
    for (ii in 1:nrow(aggfeat_mat)){
      glm_fit <- glm(delivery_ind ~ value.Component.1+value.Component.2,
               data = aggfeat_mat[-ii,], family = "binomial")
      predprob_glm[ii] <- predict(glm_fit, newdata=aggfeat_mat[ii,], type = "response")
    }
    roc_sample_glm_tempted[ss,jj] <- roc.curve(predprob_glm[aggfeat_mat$delivery_ind], 
                                        predprob_glm[!aggfeat_mat$delivery_ind])$auc
    pr_sample_glm_tempted[ss,jj] <- pr.curve(predprob_glm[aggfeat_mat$delivery_ind], 
                                      predprob_glm[!aggfeat_mat$delivery_ind])$auc.integral
    # out-of-bag prediction for random forest
    rf_fit <- randomForest(as.factor(delivery_ind) ~ value.Component.1+value.Component.2,
                   data = aggfeat_mat)
    predprob_rf <- predict(rf_fit, type = "prob")[,"TRUE"]
    roc_sample_rf_tempted[ss,jj] <- roc.curve(predprob_rf[aggfeat_mat$delivery_ind], 
                                        predprob_rf[!aggfeat_mat$delivery_ind])$auc
    pr_sample_rf_tempted[ss,jj] <- pr.curve(predprob_rf[aggfeat_mat$delivery_ind], 
                                      predprob_rf[!aggfeat_mat$delivery_ind])$auc.integral
  }
}
write.csv(roc_sample_glm_tempted, 
        file='result/realsim_ecam_roc_sample_glm_tempted_clr.csv')
write.csv(pr_sample_glm_tempted, 
        file='result/realsim_ecam_pr_sample_glm_tempted_clr.csv')
write.csv(roc_sample_rf_tempted, 
        file='result/realsim_ecam_roc_sample_rf_tempted_clr.csv')
write.csv(pr_sample_rf_tempted, 
        file='result/realsim_ecam_pr_sample_rf_tempted_clr.csv')
```

## Subject-level

### In-sample ROC & PR

-   logistic regression
-   random forest


```r
roc_sub_glm_tempted <- matrix(NA, nsim, length(nkeep))
colnames(roc_sub_glm_tempted) <- paste0("nsample", nkeep)
pr_sub_rf_tempted <- roc_sub_rf_tempted <- pr_sub_glm_tempted <- roc_sub_glm_tempted
for (jj in 1:length(nkeep)){
  print(jj)
  for (ss in 1:nsim){
    fname <- paste0('tempted_subj_sim',ss,'_ntime',nkeep[jj], '.csv')
    tempted_sub <- read.csv(file.path("simresult", fname), row.names=1)
    # glm
    predprob_glm <- rep(NA, nrow(tempted_sub))
    for (ii in 1:nrow(tempted_sub)){
      glm_fit <- glm(delivery_ind ~ Component.1+Component.2,
                   data = tempted_sub[-ii,], family = "binomial")
      predprob_glm[ii] <- predict(glm_fit, newdata=tempted_sub[ii,], type = "response")
    }
    roc_sub_glm_tempted[ss,jj] <- roc.curve(predprob_glm[tempted_sub$delivery_ind], 
                                        predprob_glm[!tempted_sub$delivery_ind])$auc
    pr_sub_glm_tempted[ss,jj] <- pr.curve(predprob_glm[tempted_sub$delivery_ind], 
                                      predprob_glm[!tempted_sub$delivery_ind])$auc.integral
    # random forest
    rf_fit <- randomForest(as.factor(delivery_ind) ~ Component.1+Component.2,
                   data = tempted_sub)
    predprob_rf <- predict(rf_fit, type = "prob")[,"TRUE"]
    roc_sub_rf_tempted[ss,jj] <- roc.curve(predprob_rf[tempted_sub$delivery_ind], 
                                        predprob_rf[!tempted_sub$delivery_ind])$auc
    pr_sub_rf_tempted[ss,jj] <- pr.curve(predprob_rf[tempted_sub$delivery_ind], 
                                      predprob_rf[!tempted_sub$delivery_ind])$auc.integral
  }
}

write.csv(roc_sub_glm_tempted, 
          file='result/realsim_ecam_roc_sub_glm_tempted_clr.csv')
write.csv(pr_sub_glm_tempted, 
          file='result/realsim_ecam_pr_sub_glm_tempted_clr.csv')
write.csv(roc_sub_rf_tempted, 
          file='result/realsim_ecam_roc_sub_rf_tempted_clr.csv')
write.csv(pr_sub_rf_tempted, 
          file='result/realsim_ecam_pr_sub_rf_tempted_clr.csv')
```

### Out-of-Sample ROC & PR

-   Logistic regression
-   Random forest


```r
npc <- 2
roc_glm_oos <- matrix(NA, nsim, length(nkeep))
colnames(roc_glm_oos) <- paste0("nsample", nkeep)
pr_rf_oos <- roc_rf_oos <- pr_glm_oos <- roc_glm_oos
for (jj in 1:length(nkeep)){
  print(jj)
  for (ss in 1:nsim){
    count_sub <- count_tab[sampleID[ss,jj][[1]],]
    meta_sub <- meta_tab[sampleID[ss,jj][[1]],]
    subdata <- format_tempted(count_sub, meta_sub$day_of_life, meta_sub$studyid, 
                            threshold=0.95, pseudo_count=0.5, transform='clr')
    metauni_sub <- metauni[names(subdata),]
    # leave one out prediction
    predprob_glm <- predprob_rf <- rep(NA, length(subdata))
    
    for (ii in 1:length(subdata)){
      print(ii)
      svd_train <- svd_centralize(subdata[-ii])
      res_train <- tempted(svd_train$datlist, r = npc, resolution = 101, smooth=1e-4)
      A_test <- est_A(subdata[ii], res_train, svd_train)
      dftrain <- data.frame(y=metauni[-ii,'delivery_ind'], x=res_train$A)
      dftest <- data.frame(y=metauni[ii,'delivery_ind'], x=A_test)
      # logistic regression
      glm_fit <- glm(y ~ ., data = dftrain, family = "binomial")
      predprob_glm[ii] <- predict(glm_fit, newdata=dftest, type = "response")
      # random forest
      rf_fit <- randomForest(as.factor(y) ~ ., data = dftrain)
      predprob_rf[ii] <- predict(rf_fit, newdata=dftest, type = "prob")[,"TRUE"]
    }
    roc_glm_oos[ss,jj] <- roc.curve(predprob_glm[metauni_sub$delivery_ind], 
                                        predprob_glm[!metauni_sub$delivery_ind])$auc
    pr_glm_oos[ss,jj] <- pr.curve(predprob_glm[metauni_sub$delivery_ind], 
                                      predprob_glm[!metauni_sub$delivery_ind])$auc.integral
    roc_rf_oos[ss,jj] <- roc.curve(predprob_rf[metauni_sub$delivery_ind], 
                                        predprob_rf[!metauni_sub$delivery_ind])$auc
    pr_rf_oos[ss,jj] <- pr.curve(predprob_rf[metauni_sub$delivery_ind], 
                                      predprob_rf[!metauni_sub$delivery_ind])$auc.integral
  }
  write.csv(pr_glm_oos, file="result/realsim_ecam_pr_oos_glm_tempted.csv")
  write.csv(roc_glm_oos, file="result/realsim_ecam_roc_oos_glm_tempted.csv")
  write.csv(pr_rf_oos, file="result/realsim_ecam_pr_oos_rf_tempted.csv")
  write.csv(roc_rf_oos, file="result/realsim_ecam_roc_oos_rf_tempted.csv")
}
```

# CTF

## Run CTF, save subject distance matrix and trajectory


```r
for (jj in 1:length(nkeep)){
  ntime <- nkeep[jj]
  for (ss in 1:nsim){
    count_sub <- count_tab[sampleID[ss,jj][[1]],]
    meta_sub <- meta_tab[sampleID[ss,jj][[1]],c("studyid", "day_of_life", "month", "delivery")]
    meta_sub$studyid <- as.character(meta_sub$studyid)
    meta_sub$month <- as.numeric(meta_sub$month)
    meta_sub$delivery_ind <- meta_sub$delivery=="Vaginal"
    metauni <- unique(meta_sub[,c("studyid", "delivery_ind")])
    py_run_file(file="run_ecam_ctf.py",convert=F)
  }
}
```

## Sample-level

### PERMANOVA F value


```r
Fmodel_ctf <- matrix(NA, nsim, length(nkeep))
colnames(Fmodel_ctf) <- paste0("nsample", nkeep)
for (jj in 1:length(nkeep)){
  ntime <- nkeep[jj]
  for (ss in 1:nsim){
    count_sub <- count_tab[sampleID[ss,jj][[1]],]
    meta_sub <- meta_tab[sampleID[ss,jj][[1]],c("studyid", "day_of_life", "month", "delivery")]
    meta_sub$studyid <- as.character(meta_sub$studyid)
    meta_sub$month <- as.numeric(meta_sub$month)
    meta_sub$delivery_ind <- meta_sub$delivery=="Vaginal"
    metauni <- unique(meta_sub[,c("studyid", "delivery_ind")])
    # calculate PERMANOVA F
    fname <- paste0("distance-matrix_sim", ss, "_ntime", ntime, ".qza")
    ctf_dist <- as.matrix(read_qza(file.path("simresult", fname))$data)
    deliv <- meta_sub[rownames(ctf_dist),]$delivery
    res_perm <- adonis2(ctf_dist ~ deliv)
    Fmodel_ctf[ss,jj] <- res_perm$F[1]
  }
}

write.csv(Fmodel_ctf, 
          file='result/realsim_ecam_Fvalue_ctf.csv')
```

### In-sample ROC & PR

-   logistic regression
-   random forest


```r
roc_sample_glm_ctf <- matrix(NA, nsim, length(nkeep))
colnames(roc_sample_glm_ctf) <- paste0("nsample", nkeep)
pr_sample_rf_ctf <- roc_sample_rf_ctf <- pr_sample_glm_ctf <- roc_sample_glm_ctf

for (jj in 1:length(nkeep)){
  print(jj)
  ntime <- nkeep[jj]
  for (ss in 1:nsim){
    count_sub <- count_tab[sampleID[ss,jj][[1]],]
    meta_sub <- meta_tab[sampleID[ss,jj][[1]],c("studyid", "day_of_life", "month", "delivery")]
    meta_sub$studyid <- as.character(meta_sub$studyid)
    meta_sub$month <- as.numeric(meta_sub$month)
    meta_sub$delivery_ind <- meta_sub$delivery=="Vaginal"
    # calculate ROC & PR
    fname <- paste0("distance-matrix_sim", ss, "_ntime", ntime, ".qza")
    ctf_dist <- as.matrix(read_qza(file.path("simresult", fname))$data)
    ctf_mds <- as.data.frame(cmdscale(ctf_dist, k=2))
    colnames(ctf_mds) <- c("PC1", "PC2")
    ctf_mds$delivery_ind <- meta_sub[rownames(ctf_dist),]$delivery_ind
    predprob_glm <- rep(NA, nrow(ctf_mds))
    for(ii in 1:nrow(ctf_mds)){
      glm_fit <- glm(delivery_ind~PC1+PC2, data=ctf_mds[-ii,])
      predprob_glm[ii] <- predict(glm_fit, newdata=ctf_mds[ii,], type = "response")
    }
    roc_sample_glm_ctf[ss,jj] <- roc.curve(predprob_glm[ctf_mds$delivery_ind], 
                                    predprob_glm[!ctf_mds$delivery_ind])$auc
    pr_sample_glm_ctf[ss,jj] <- pr.curve(predprob_glm[ctf_mds$delivery_ind], 
                                  predprob_glm[!ctf_mds$delivery_ind])$auc.integral    
    predprob_rf <- rep(NA, nrow(ctf_mds))
    rf_fit <- randomForest(as.factor(delivery_ind)~PC1+PC2, data=ctf_mds)
    predprob_rf <- predict(rf_fit, type = "prob")[,"TRUE"]
    roc_sample_rf_ctf[ss,jj] <- roc.curve(predprob_rf[ctf_mds$delivery_ind], 
                                    predprob_rf[!ctf_mds$delivery_ind])$auc
    pr_sample_rf_ctf[ss,jj] <- pr.curve(predprob_rf[ctf_mds$delivery_ind], 
                                  predprob_rf[!ctf_mds$delivery_ind])$auc.integral  
  }
}

write.csv(roc_sample_glm_ctf, 
          file='result/realsim_ecam_roc_sample_glm_ctf.csv')
write.csv(pr_sample_glm_ctf, 
          file='result/realsim_ecam_pr_sample_glm_ctf.csv')
write.csv(roc_sample_rf_ctf, 
          file='result/realsim_ecam_roc_sample_rf_ctf.csv')
write.csv(pr_sample_rf_ctf, 
          file='result/realsim_ecam_pr_sample_rf_ctf.csv')
```

## Subject-level

### In-sample ROC & PR

-   logistic regression
-   random forest


```r
roc_sub_glm_ctf <- matrix(NA, nsim, length(nkeep))
colnames(roc_sub_glm_ctf) <- paste0("nsample", nkeep)
pr_sub_rf_ctf <- roc_sub_rf_ctf <- pr_sub_glm_ctf <- roc_sub_glm_ctf
for (jj in 1:length(nkeep)){
  print(jj)
  ntime <- nkeep[jj]
  for (ss in 1:nsim){
    count_sub <- count_tab[sampleID[ss,jj][[1]],]
    meta_sub <- meta_tab[sampleID[ss,jj][[1]],c("studyid", "day_of_life", "month", "delivery")]
    meta_sub$studyid <- as.character(meta_sub$studyid)
    meta_sub$month <- as.numeric(meta_sub$month)
    meta_sub$delivery_ind <- meta_sub$delivery=="Vaginal"
    metauni_sub <- unique(meta_sub[,c('studyid', 'delivery_ind')])
    rownames(metauni_sub) <- metauni_sub$studyid
    fname <- paste0("subject-biplot_sim", ss, "_ntime", ntime, ".qza")
    ctf_sub <- read_qza(file.path("simresult", fname))$data$Vectors
    colnames(ctf_sub)[1] <- "studyid"
    ctf_sub <- merge(ctf_sub, metauni_sub)
    # logistic regression
    predprob_glm <- rep(NA, nrow(ctf_sub))
    for(ii in 1:nrow(ctf_sub)){
      glm_fit <- glm(delivery_ind~PC1+PC2, data=ctf_sub[-ii,])
      predprob_glm[ii] <- predict(glm_fit, newdata=ctf_sub[ii,], type = "response")
    }
    roc_sub_glm_ctf[ss,jj] <- roc.curve(predprob_glm[ctf_sub$delivery_ind], 
                                    predprob_glm[!ctf_sub$delivery_ind])$auc
    pr_sub_glm_ctf[ss,jj] <- pr.curve(predprob_glm[ctf_sub$delivery_ind], 
                                  predprob_glm[!ctf_sub$delivery_ind])$auc.integral
    # random forest
    rf_fit <- randomForest(as.factor(delivery_ind)~PC1+PC2, data=ctf_sub)
    predprob_rf <- predict(rf_fit, type = "prob")[,"TRUE"]
    roc_sub_rf_ctf[ss,jj] <- roc.curve(predprob_rf[ctf_sub$delivery_ind], 
                                    predprob_rf[!ctf_sub$delivery_ind])$auc
    pr_sub_rf_ctf[ss,jj] <- pr.curve(predprob_rf[ctf_sub$delivery_ind], 
                                  predprob_rf[!ctf_sub$delivery_ind])$auc.integral
  }
}
write.csv(roc_sub_glm_ctf, 
          file='result/realsim_ecam_roc_sub_glm_ctf.csv')
write.csv(pr_sub_glm_ctf, 
          file='result/realsim_ecam_pr_sub_glm_ctf.csv')
write.csv(roc_sub_rf_ctf, 
          file='result/realsim_ecam_roc_sub_rf_ctf.csv')
write.csv(pr_sub_rf_ctf, 
          file='result/realsim_ecam_pr_sub_rf_ctf.csv')
```

# PCoA

## Sample-level

### PERMANOVA F value


```r
metric_file <- list.files('data', pattern = "distance_matrix")
metric_file
metric_name <- c("Aitchison", "Bray-Curtis", 
                 "gUniFrac-alpha0", "gUniFrac-alpha1", "gUniFrac-alpha5", 
                 "Jaccard", "unweighted_UniFrac", "Weighted_UniFrac")
distmat_all <- vector(mode="list", length(metric_name))
names(distmat_all) <- metric_name
for (ii in 1:length(metric_name)){
  tmp <- read_qza(paste0("data/", metric_file[ii]))$data
  distmat_all[[ii]] <- as.matrix(tmp)
}

Fmodel <- array(0, dim=c(length(metric_name), nsim, length(nkeep)),
                     dimnames=list(metric_name, paste0('sim',1:nsim), paste0('nsample=',nkeep)))
for (jj in 1:length(nkeep)){
  print(jj)
  for (ss in 1:nsim){
    print(paste0('nsample=', nkeep[jj]))
    ind_sub <- sampleID[ss,jj][[1]]
    meta_sub <- meta_tab[ind_sub,]
    for (ii in 1:length(metric_name)){
      dist_sub <- distmat_all[[ii]][ind_sub,ind_sub]
      res_perm <- adonis2(dist_sub ~ meta_sub$delivery)
      Fmodel[ii,ss,jj] <- res_perm$F[1]
    }
  }
  for (ii in 1:length(metric_name)){
    write.csv(Fmodel[ii,,], 
       file=paste0('result/realsim_ecam_Fvalue_PCoA_', metric_name[ii], '.csv'))
  }
}
```

### In-sample ROC & PR

-   logistic regression
-   random forest


```r
metric_file <- list.files('data', pattern = "distance_matrix")
metric_file
metric_name <- c("Aitchison", "Bray-Curtis", 
                 "gUniFrac-alpha0", "gUniFrac-alpha1", "gUniFrac-alpha5", 
                 "Jaccard", "unweighted_UniFrac", "Weighted_UniFrac")
distmat_all <- vector(mode="list", length(metric_name))
names(distmat_all) <- metric_name
for (ii in 1:length(metric_name)){
  tmp <- read_qza(paste0("data/", metric_file[ii]))$data
  distmat_all[[ii]] <- as.matrix(tmp)
}

pr_glm_array <- roc_glm_array <- 
  pr_rf_array <- roc_rf_array <-
  array(0, dim=c(length(metric_name), nsim, length(nkeep)),
                     dimnames=list(metric_name, paste0('sim',1:nsim), paste0('nsample=',nkeep)))
for (jj in 1:length(nkeep)){
  print(jj)
  for (ss in 1:nsim){
    ind_sub <- sampleID[ss,jj][[1]]
    meta_sub <- meta_tab[ind_sub,]
    for (mm in 1:length(metric_name)){
      pcoa_dist <- distmat_all[[mm]][ind_sub,ind_sub]
      pcoa_mds <- as.data.frame(cmdscale(pcoa_dist, k=2))
      colnames(pcoa_mds) <- c("PC1", "PC2")
      pcoa_mds$delivery_ind <- meta_sub[rownames(pcoa_dist),]$delivery_ind
      # glm
      predprob_glm <- rep(NA, nrow(pcoa_mds))
      for(ii in 1:nrow(pcoa_mds)){
        glm_fit <- glm(delivery_ind~PC1+PC2, data=pcoa_mds[-ii,])
        predprob_glm[ii] <- predict(glm_fit, newdata=pcoa_mds[ii,], type = "response")
      }
      roc_glm_array[mm,ss,jj] <- roc.curve(predprob_glm[pcoa_mds$delivery_ind], 
                                      predprob_glm[!pcoa_mds$delivery_ind])$auc
      pr_glm_array[mm,ss,jj] <- pr.curve(predprob_glm[pcoa_mds$delivery_ind], 
                                    predprob_glm[!pcoa_mds$delivery_ind])$auc.integral
      # random forest
      rf_fit <- randomForest(as.factor(delivery_ind)~PC1+PC2, data=pcoa_mds)
      predprob_rf <- predict(rf_fit, type = "prob")[,"TRUE"]
      roc_rf_array[mm,ss,jj] <- roc.curve(predprob_rf[pcoa_mds$delivery_ind], 
                                      predprob_rf[!pcoa_mds$delivery_ind])$auc
      pr_rf_array[mm,ss,jj] <- pr.curve(predprob_rf[pcoa_mds$delivery_ind], 
                                    predprob_rf[!pcoa_mds$delivery_ind])$auc.integral
    }
  }
  for (mm in 1:length(metric_name)){
    write.csv(roc_glm_array[mm,,], 
       file=paste0('result/realsim_ecam_roc_sample_glm_PCoA_', metric_name[mm], '.csv'))
    write.csv(pr_glm_array[mm,,], 
       file=paste0('result/realsim_ecam_pr_sample_glm_PCoA_', metric_name[mm], '.csv'))
    write.csv(roc_rf_array[mm,,], 
       file=paste0('result/realsim_ecam_roc_sample_rf_PCoA_', metric_name[mm], '.csv'))
    write.csv(pr_rf_array[mm,,], 
       file=paste0('result/realsim_ecam_pr_sample_rf_PCoA_', metric_name[mm], '.csv'))
  }
}
```

# Summarize results

## Summarize subject-level ROC and PR


```r
method  <- c("tempted_clr", "ctf")
measure <- c("roc", "pr")
classify <- c("glm", "rf")

tab_auc_subj <- NULL

for (ms in measure){
  for (cls in classify){ 
    # in sample
    for (mthd in method){
      fname <- paste0("result/realsim_ecam_", ms, "_sub_", cls, "_", mthd, ".csv")
      tab0 <- read.csv(fname, row.names=1)
      tab0 <- 1-as.matrix(tab0)
      tab_tmp <- data.frame(auc=as.vector(tab0),
                      nsample=as.vector(t(matrix(rep(nkeep,nsim),length(nkeep),nsim))),
                      method=toupper(mthd), measure=toupper(ms), classify=cls,
                      type="In-Sample")
      tab_auc_subj <- rbind(tab_auc_subj, tab_tmp)
    }
    # out of sample
    fname <- paste0("result/realsim_ecam_", ms, "_oos_", cls, "_tempted.csv")
    tab0 <- read.csv(fname, row.names=1)
    tab0 <- 1-as.matrix(tab0)
    tab_tmp <- data.frame(auc=as.vector(tab0),
                    nsample=as.vector(t(matrix(rep(nkeep,nsim),length(nkeep),nsim))),
                    method="TEMPTED", measure=toupper(ms), classify=cls,
                    type="Out-of-Sample")
    tab_auc_subj <- rbind(tab_auc_subj, tab_tmp)
  }  
}
tab_auc_subj$classify <- gsub("glm", "Logistic Regression", tab_auc_subj$classify)
tab_auc_subj$classify <- gsub("rf", "Random Forest", tab_auc_subj$classify)
tab_auc_subj$method <- gsub("TEMPTED_CLR", "TEMPTED", tab_auc_subj$method)
```

### Plot all subject-level ROC and PR


```r
tab_auc_subj$nsample <- as.factor(tab_auc_subj$nsample)
tab1 <- aggregate(auc~nsample+measure+method+classify+type, data=tab_auc_subj, 
                  FUN=mean)
names(tab1)[6] <- 'mean'
tab2 <- aggregate(auc~nsample+measure+method+classify+type, data=tab_auc_subj, 
                  FUN=function(x){sd(x)/sqrt(length(x))})
names(tab2)[6] <- 'se'
rownames(tab1) <-rownames(tab2) <- NULL
tab_auc_subj_summary <- merge(tab1, tab2)
tab_auc_subj_summary$nsample <- factor(tab_auc_subj_summary$nsample, 
                                     level=as.character(nkeep))

color_method <- c('#33a02c', #green for CTF
                  '#377eb8') #blue for TEMPTED
p_subj_roc <- ggplot(data=dplyr::filter(tab_auc_subj_summary, measure=="ROC"), 
                           aes(x=nsample, y=mean, group=paste0(type,method), color=method)) + 
  geom_line(size=1, position=position_dodge(0.5), aes(linetype=type)) + 
  geom_point(size=2, position=position_dodge(0.5)) +
  geom_errorbar(aes(ymin=mean-2*se, ymax=mean+2*se, width=0.5), 
                position=position_dodge(0.5), size=1) + 
  facet_wrap(.~classify) +
  labs(y='AUC-ROC error', x='# time points') + 
  scale_color_manual(values=color_method) +
  theme_bw()
```

```
## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
## ℹ Please use `linewidth` instead.
## This warning is displayed once every 8 hours.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.
```

```r
p_subj_roc
```

![plot of chunk plot_sub_auc](figure/plot_sub_auc-1.png)

```r
p_subj_pr <- ggplot(data=dplyr::filter(tab_auc_subj_summary, measure=="PR"), 
                            aes(x=nsample, y=mean, group=paste0(type,method), color=method)) + 
  geom_line(size=1, position=position_dodge(0.5), aes(linetype=type)) + 
  geom_point(size=2, position=position_dodge(0.5)) +
  geom_errorbar(aes(ymin=mean-2*se, ymax=mean+2*se, width=0.5), 
                position=position_dodge(0.5), size=1) + 
  facet_wrap(.~classify) +
  labs(y='AUC-PR error', x='# time points') +
  scale_color_manual(values=color_method) +
  theme_bw()
p_subj_pr
```

![plot of chunk plot_sub_auc](figure/plot_sub_auc-2.png)

```r
pdf('../figure_table/realsim_ecam_pr_roc.pdf', width=7,height=5)
print(p_subj_pr)
print(p_subj_roc)
dev.off()
```

```
## quartz_off_screen 
##                 2
```

## Summarize sample-level PERMANOVA F values


```r
metric_name <- c("Bray-Curtis",
                 "unweighted_UniFrac", "Weighted_UniFrac")
# read in results
Fmodel <- array(0, dim=c(length(metric_name)+2, nsim, length(nkeep)),
                dimnames=list(c(metric_name, "TEMPTED", "CTF"), paste0('sim',1:nsim), paste0('nsample=',nkeep)))
for (ii in 1:length(metric_name)){
  Fmodel[ii,,] <- as.matrix(read.csv(paste0('result/realsim_ecam_Fvalue_PCoA_', metric_name[ii], '.csv'),
                           row.names=1))
}
Fmodel_tempted <- read.csv("result/realsim_ecam_Fvalue_tempted_clr.csv", row.names=1)
Fmodel[ii+1,,] <- as.matrix(Fmodel_tempted)
Fmodel_ctf <- read.csv("result/realsim_ecam_Fvalue_ctf.csv", row.names=1)
Fmodel[ii+2,,] <- as.matrix(Fmodel_ctf)
# make table
tab_Fmodel <- NULL
for (ii in 1:dim(Fmodel)[1]){
  tab_tmp <- data.frame(Fvalue=as.vector(Fmodel[ii,,]),
                        nsample=as.vector(t(matrix(rep(nkeep,nsim),length(nkeep),nsim))),
                        method=dimnames(Fmodel)[[1]][ii])
  tab_Fmodel <- rbind(tab_Fmodel, tab_tmp)
}
tab_Fmodel$method[tab_Fmodel$method=="unweighted_UniFrac"] <- "UniFrac"
tab_Fmodel$method[tab_Fmodel$method=="Weighted_UniFrac"] <- "Weighted UniFrac"
```

### Plot all PERMANOVA F-values


```r
tab1 <- aggregate(Fvalue~nsample+method, data=tab_Fmodel, 
                  FUN=mean)
names(tab1)[3] <- 'mean'
tab2 <- aggregate(Fvalue~nsample+method, data=tab_Fmodel, 
                  FUN=function(x){sd(x)/sqrt(length(x))})
names(tab2)[3] <- 'se'
rownames(tab1) <-rownames(tab2) <- NULL
tab_Fmodel_summary <- merge(tab1, tab2)

color_method <- c('#ff7f00', #orange for Bray
                  '#33a02c', #green for CTF
                  '#377eb8', #blue for TEMPTED
                  '#6a3d9a', #purple for Unifrac
                  '#e31a1c') #red for unweighted Unifrac
p_Fmodel_summary <- ggplot(data=tab_Fmodel_summary, 
                                aes(x=nsample, y=mean, group=method, color=method)) + 
  geom_line(size=1, position=position_dodge(0.1)) + geom_point(size=2, position=position_dodge(0.1)) +
  geom_errorbar(aes(ymin=mean-2*se, ymax=mean+2*se, width=1), 
                position=position_dodge(0.1), size=1) + 
  scale_x_continuous(breaks=2:10) + 
  scale_color_manual(values=color_method) +
  ylab('PERMANOVA F value') + xlab('# time points') + 
  theme_bw() + 
  theme(legend.position = "bottom")
p_Fmodel_summary
```

![plot of chunk plot_Fvalue](figure/plot_Fvalue-1.png)

```r
pdf('../figure_table/realsim_ecam_Fvalue.pdf', width=7.7,height=5)
p_Fmodel_summary
dev.off()
```

```
## quartz_off_screen 
##                 2
```

## Summarize sample-level ROC & PR


```r
method  <- c("tempted_clr", "ctf", "PCoA_Bray-Curtis", "PCoA_unweighted_UniFrac", "PCoA_Weighted_UniFrac")
measure <- c("roc", "pr")
classify <- c("glm", "rf")

tab_auc_sample <- NULL

for (ms in measure){
  for (cls in classify){  
    for (mthd in method){
      fname <- paste0("result/realsim_ecam_", ms, "_sample_", cls, "_", mthd, ".csv")
      tab0 <- read.csv(fname, row.names=1)
      tab0 <- 1-as.matrix(tab0)
      tab_tmp <- data.frame(auc=as.vector(tab0),
                      nsample=as.vector(t(matrix(rep(nkeep,nsim),length(nkeep),nsim))),
                      method=mthd, measure=toupper(ms), classify=cls,
                      type="In-Sample")
      tab_auc_sample <- rbind(tab_auc_sample, tab_tmp)
    }
  }  
}
tab_auc_sample$classify <- gsub("glm", "Logistic Regression", tab_auc_sample$classify)
tab_auc_sample$classify <- gsub("rf", "Random Forest", tab_auc_sample$classify)
tab_auc_sample$method <- gsub("tempted_clr", "TEMPTED", tab_auc_sample$method)
tab_auc_sample$method <- gsub("ctf", "CTF", tab_auc_sample$method)
tab_auc_sample$method <- gsub("PCoA_Bray-Curtis", "Bray-Curtis", tab_auc_sample$method)
tab_auc_sample$method <- gsub("PCoA_unweighted_UniFrac", "UniFrac", tab_auc_sample$method)
tab_auc_sample$method <- gsub("PCoA_Weighted_UniFrac", "Weighted UniFrac", tab_auc_sample$method)
```

### Plot sample level ROC & PR


```r
tab1 <- aggregate(auc~nsample+method+type+measure+classify, data=tab_auc_sample, 
                  FUN=mean)
names(tab1)[6] <- 'mean'
tab2 <- aggregate(auc~nsample+method+type+measure+classify, data=tab_auc_sample, 
                  FUN=function(x){sd(x)/sqrt(length(x))})
names(tab2)[6] <- 'se'
rownames(tab1) <-rownames(tab2) <- NULL
tab_auc_sample_summary <- merge(tab1, tab2)
tab_auc_sample_summary$nsample <- factor(tab_auc_sample_summary$nsample, 
                                         level=as.character(nkeep))

color_method <- c('#ff7f00', #orange for Bray
                  '#33a02c', #green for CTF
                  '#377eb8', #blue for TEMPTED
                  '#6a3d9a', #purple for Unifrac
                  '#e31a1c') #red for Weighted Unifrac
p_sample_roc <- ggplot(data=dplyr::filter(tab_auc_sample_summary, measure=="ROC"), 
                           aes(x=nsample, y=mean, group=paste0(type,method), color=method)) + 
  geom_line(position=position_dodge(0.1), size=1, aes(linetype=type)) + 
  geom_point(position=position_dodge(0.1), size=2) +
  geom_errorbar(aes(ymin=mean-2*se, ymax=mean+2*se, width=0.5), 
                position=position_dodge(0.1), size=1) + 
  facet_wrap(.~classify) +
  labs(y='AUC-ROC error', x='# time points') + 
  scale_color_manual(values=color_method) +
  theme_bw() +
  theme(legend.position="bottom")
p_sample_roc
```

![plot of chunk plot_sample_auc](figure/plot_sample_auc-1.png)

```r
p_sample_pr <- ggplot(data=dplyr::filter(tab_auc_sample_summary, measure=="PR"), 
                           aes(x=nsample, y=mean, group=paste0(type,method), color=method)) + 
  geom_line(position=position_dodge(0.1), size=1, aes(linetype=type)) + 
  geom_point(position=position_dodge(0.1), size=2) +
  geom_errorbar(aes(ymin=mean-2*se, ymax=mean+2*se, width=0.5), 
                position=position_dodge(0.1), size=1) + 
  facet_wrap(.~classify) +
  labs(y='AUC-PR error', x='# time points') + 
  scale_color_manual(values=color_method) +
  theme_bw() +
  theme(legend.position="bottom")
p_sample_pr
```

![plot of chunk plot_sample_auc](figure/plot_sample_auc-2.png)

```r
pdf('../figure_table/realsim_ecam_sample_auc.pdf', width=7.7,height=5)
p_sample_roc
p_sample_pr
dev.off()
```

```
## quartz_off_screen 
##                 2
```

#### get legend color for all


```r
color_method <- c('#ff7f00', #orange for Bray
                  '#33a02c', #green for CTF
                  '#377eb8', #blue for TEMPTED
                  '#6a3d9a', #purple for Unifrac
                  '#e31a1c') #red for unweighted Unifrac
tab_lgd <- data.frame(Method=rep(c("Bray-Curtis", "CTF", 
                               "TEMPTED", "UniFrac","Weighted UniFrac"),4),
                      value=rnorm(20), time=rnorm(20),
                      Type=c(rep("In-Sample",each=10), rep("Out-of-Sample", each=10)))
p_lgd <- ggplot(data=tab_lgd, aes(x=time, y=value, color=Method)) + 
  geom_point(size=2) + geom_line(aes(linetype=Type), size=1) + 
  scale_color_manual(values=color_method) +
  theme(legend.position="bottom")
p_lgd
```

![plot of chunk get_legend](figure/get_legend-1.png)

```r
lgd <- get_legend(p_lgd)
```

## Print plot for manuscript


```r
lay <- rbind(c(1,2,2),c(1,2,2),c(1,2,2),c(1,2,2),c(1,2,2),
             c(3,3,3))
p_all <- grid.arrange(p_Fmodel_summary + theme(legend.position="none"), 
             p_subj_pr + theme(legend.position="none"),
             lgd,
             layout_matrix=lay)
```

![plot of chunk plot_all](figure/plot_all-1.png)

```r
ggsave('../figure_table/realsim_ecam.pdf', width=9, height=4, plot=p_all)
```
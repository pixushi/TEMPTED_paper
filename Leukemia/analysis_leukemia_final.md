---
title: "Analysis of Mice Leukemia Data"
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
library(readr) # read tsv
```

```
## 
## Attaching package: 'readr'
```

```
## The following object is masked from 'package:scales':
## 
##     col_factor
```

```r
library(qiime2R) # read in Qiime artifacts
library(dplyr) # data formatting
library(yaml) # for read_qza() in qiime2R
# for computing
library(reticulate) # run py codes
library(phyloseq) # phyloseq object
```

```
## 
## Attaching package: 'phyloseq'
```

```
## The following object is masked from 'package:reticulate':
## 
##     import
```

```r
library(vegan) # distance matrix
library(PERMANOVA) # permanova
library(randomForest) # random forest
library(PRROC) # roc and pr
# for plotting
library(ggpubr)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(plotly)
```

```
## 
## Attaching package: 'plotly'
```

```
## The following object is masked from 'package:ggplot2':
## 
##     last_plot
```

```
## The following object is masked from 'package:MASS':
## 
##     select
```

```
## The following object is masked from 'package:stats':
## 
##     filter
```

```
## The following object is masked from 'package:graphics':
## 
##     layout
```

```r
# set working directory to be where the current script is located
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source('../TEMPTED.R')
col_group <- c(brewer.pal(6,'Set2')[6], brewer.pal(6,'Set1')[c(1,3)])
col_microbe <- c(brewer.pal(6,'Set1')[c(2,4)], brewer.pal(6,'Set2')[4])
npc <- 3 # number of components
bws <- 8 # for plotting
set.seed(1)
```


## Read in the count data, meta data, distance matrices


```r
meta_all <- read.csv("data/metadata_cleaned.csv", row.names=1)
meta_all$group <- factor(meta_all$group, levels=
                           c("genotype=WT, disease=no",
                             "genotype=Pax5+/-, disease=no",
                             "genotype=Pax5+/-, disease=pB-ALL"))
count_all <- read.csv("data/count_cleaned.csv", row.names=1)
load("data/distance_mat.Rdata")
metauni <- unique(meta_all[,c('hostID', 'genotype', 'diseased', 'group')])
metauni$group <- factor(metauni$group, sort(levels(metauni$group)))
rownames(metauni) <- metauni$hostID
table(metauni$group)
```

```
## 
##     genotype=Pax5+/-, disease=no genotype=Pax5+/-, disease=pB-ALL 
##                               10                               17 
##          genotype=WT, disease=no 
##                               17
```


## Plot timeline of study


```r
meta_ordered <- meta_all[order(meta_all$group, meta_all$host_subject_id, decreasing=T),]
meta_ordered$host_subject_id <- factor(meta_ordered$host_subject_id, 
                                       levels=unique(meta_ordered$host_subject_id))
colnames(meta_ordered)[22] <- "Group"
p_timeline <- ggplot(data=meta_ordered, 
       aes(x=week, y=host_subject_id, 
           group=as.factor(host_subject_id), color=Group, shape=Group)) +
  geom_line() + geom_point() + scale_color_manual(values=col_group) +
  labs(y="Host ID", x="Week") + theme_bw() +
  theme(legend.position="bottom") + 
  scale_x_continuous(breaks=seq(from=0,to=90,by=5))
p_timeline
```

![plot of chunk plot_timeline](figure/plot_timeline-1.png)

```r
pdf('../figure_table/leukemia_timeline.pdf', width=7.7, height=5)
print(p_timeline)
dev.off()
```

```
## quartz_off_screen 
##                 2
```



# TEMPTED

## Run analysis


```r
datlist_all <- format_tempted(count_all, meta_all$week, meta_all$hostID, 
                            threshold=0.95, pseudo_count=0.5, transform='clr')
print(dim(datlist_all[[1]]))
```

```
## [1] 1065   10
```

```r
svd_all <- svd_centralize(datlist_all, 1)
res_tempted_all <- tempted(svd_all$datlist, r = npc, resolution = 51, smooth=1e-4)
```

```
## [1] "Calculate the 1th Component"
## [1] "Convergence reached at dif=1.29075362466149e-05, iter=3"
## [1] "Calculate the 2th Component"
## [1] "Convergence reached at dif=7.67492932197077e-05, iter=9"
## [1] "Calculate the 3th Component"
## [1] "Convergence reached at dif=3.20041217205195e-05, iter=7"
```

```r
save(svd_all, res_tempted_all, datlist_all,
     file='result/res_leukemia_tempted.Rdata')
```


```r
load('result/res_leukemia_tempted.Rdata')

res_tempted <- res_tempted_all
svd_tempted <- svd_all
meta_tab <- meta_all
count_tab <- count_all
datlist <- datlist_all
```

## Plot TEMPTED time loading


```r
p_time <- plot_time_loading(res_tempted, r=3) + scale_color_manual(values=col_microbe) + 
  labs(x='Week', y="Temporal Loading", color="Component") + 
  geom_line(size=1.5) + theme_bw() + 
  theme(legend.position='bottom')
p_time
```

![plot of chunk time_loading](figure/time_loading-1.png)


## Plot TEMPTED subject loadings


```r
A.hat <- metauni
rownames(A.hat) <- A.hat$hostID
table(rownames(A.hat)==rownames(res_tempted$A.hat))
```

```
## 
## TRUE 
##   44
```

```r
A.hat <- cbind(res_tempted$A.hat[,1:3], A.hat)
A.hat$group <- factor(A.hat$group, 
  levels=c("genotype=WT, disease=no", "genotype=Pax5+/-, disease=no",
         "genotype=Pax5+/-, disease=pB-ALL"))
p_sub <- plot_ly(A.hat, x=~`Component 1`, y=~`Component 2`, z=~`Component 3`, 
             color=~group, colors=col_group)
p_sub <- p_sub %>% add_markers(marker=list(size=5))
p_sub <- p_sub %>% layout(scene=list(xaxis=list(title="Component 1", titlefont = list(size = 20)), 
                             yaxis=list(title="Component 2", titlefont = list(size = 20)), 
                             zaxis=list(title="Component 3", titlefont = list(size = 20))),
                          legend = list(font = list(size = 20), orientation = "h"))
p_sub
```

```
## Error in loadNamespace(name): there is no package called 'webshot'
```

```r
htmlwidgets::saveWidget(widget=p_sub, 
                        file="../figure_table/leukemia_sub.html",
                        selfcontained=F)
```

```r
p_sub23 <- ggplot(data=A.hat, aes(x=`Component 2`, y=`Component 3`)) +
  geom_point(aes(color=group)) + 
  scale_color_manual(values=col_group) + 
  theme_bw()+
  theme(legend.position='bottom')
```



## Plot TEMPTED trajectory of log ratio of top features


```r
datlist_raw <- format_tempted(count_all, meta_all$week, meta_all$hostID, 
                            threshold=0.95, transform='none')
ratio_feat <- ratio_feature(res_tempted, datlist_raw, pct=0.005)

## summed up, by individual subject

tab_feat_ratio <- ratio_feat$metafeature.ratio
colnames(tab_feat_ratio)[2] <- 'hostID'
tab_feat_ratio <- merge(tab_feat_ratio, metauni)

## summed up, by mean and sd
reshape_feat_ratio <- reshape(tab_feat_ratio, 
                            idvar=c("hostID","timepoint") , 
                            v.names=c("value"), timevar="PC",
                            direction="wide")
CC <- grep("value", colnames(reshape_feat_ratio))
colnames(reshape_feat_ratio)[CC] <- paste("Component", 1:length(CC))
feature_mat_ratio <- reshape_feat_ratio[,c("Component 2", "Component 3")]

time_vec_ratio <- reshape_feat_ratio$timepoint
group_vec_ratio <- factor(reshape_feat_ratio$group, 
                        levels=c("genotype=WT, disease=no", "genotype=Pax5+/-, disease=no",
                                 "genotype=Pax5+/-, disease=pB-ALL"))
p_feat_ratio_summary <- plot_feature_summary(feature_mat_ratio, 
                                           time_vec_ratio, 
                                           group_vec_ratio, bws=bws, nrow=1) + 
  xlab('Week') + theme_bw() + 
  theme(legend.position='bottom') +
  scale_color_manual(values=col_group) + scale_fill_manual(values=col_group)
p_feat_ratio_summary
```

![plot of chunk ratio_top_feature](figure/ratio_top_feature-1.png)


## Plot top features associated with disease


```r
# the log ratio is constructed using the following features
ftnames <- rownames(datlist_raw[[1]])[-1]
tab_ftnames <- NULL
for (ii in 1:npc){
  tab_ftnames <- rbind(tab_ftnames, 
                       data.frame(sequence=ftnames[ratio_feat$toppct[,ii]],
                                  PC=ii, rank="Top"))
  tab_ftnames <- rbind(tab_ftnames, 
                       data.frame(sequence=ftnames[ratio_feat$bottompct[,ii]],
                                  PC=ii, rank="Bottom"))
}
taxtab <- read.csv("data/gg2-taxonomy.tsv", sep="\t", row.names=1)
tab_ftnames <- cbind(tab_ftnames, taxtab[tab_ftnames$sequence,])
tab_ftnames$dummy <- paste0("ASV", 1:nrow(tab_ftnames))
write.csv(tab_ftnames, file="../figure_table/leukemia_topfeature.csv")
```


```r
tab_ftnames_pc2 <- dplyr::filter(tab_ftnames, PC==2)
topfeat_pc2 <- tab_ftnames_pc2$sequence

prop_all <- (count_all)/rowSums(count_all)
feat_mat_pc2 <- prop_all[meta_all$genotype!="WT",topfeat_pc2]
nfeat <- ncol(feat_mat_pc2)
meta_pax <- meta_all[meta_all$genotype!="WT",]
tab_topfeat <- data.frame(RA=as.vector(as.matrix(feat_mat_pc2)),
                          feature=rep(tab_ftnames_pc2$dummy, each=nrow(feat_mat_pc2)),
                          time=rep(meta_pax$week, nfeat),
                          group=rep(meta_pax$diseased, nfeat),
                          hostID=rep(meta_pax$hostID, nfeat))
p_topfeat <- ggplot(data=tab_topfeat) + 
  geom_point(aes(x=time, y=RA, color=group), size=0.8) + 
  facet_wrap(vars(feature), nrow=2, scales="free_y") + theme_bw() + 
  labs(x="Week", y="Relative Abundance", color="Group")
pdf("../figure_table/leukemia_topfeat.pdf", height=5, width=10)
print(p_topfeat)
dev.off()
```

```
## quartz_off_screen 
##                 2
```



```r
lay <- rbind(c(1,2,3,3), c(1,2,3,3), c(1,2,3,3), c(1,2,3,3), c(1,2,3,3), 
             c(4,5,5,5))
lgd1 <- get_legend(p_time)
lgd2 <- get_legend(p_sub23)
p_tempted <- grid.arrange(p_time+theme(legend.position="none"),
                          p_sub23+theme(legend.position="none"),
                          p_feat_ratio_summary+theme(legend.position="none"),
                          lgd1, lgd2,
                          layout_matrix=lay)
```

![plot of chunk tempted_allplot](figure/tempted_allplot-1.png)

```r
ggsave('../figure_table/leukemia_plot.pdf', width=10, height=3, plot=p_tempted)
```


# CTF

## Run analysis


```r
py_run_file(file="analysis_leukemia_ctf.py",convert=F)
```


## CTF subject loadings


```r
ctf_sub <- read_qza("result/subject-biplot_leukemia.qza")$data$Vectors
ctf_sub$group <- metauni[ctf_sub$SampleID,"group"]
p_sub <- plot_ly(ctf_sub, x=~`PC1`, y=~`PC2`, z=~`PC3`, 
                 color=~group, colors=col_group)
p_sub <- p_sub %>% add_markers(marker=list(size=5))
p_sub <- p_sub %>% layout(scene=list(xaxis=list(title="Component 1", titlefont = list(size = 20)), 
                                     yaxis=list(title="Component 2", titlefont = list(size = 20)), 
                                     zaxis=list(title="Component 3", titlefont = list(size = 20))),
                          legend = list(font = list(size = 20), orientation = "h"))
p_sub
```

```
## Error in loadNamespace(name): there is no package called 'webshot'
```

```r
htmlwidgets::saveWidget(widget=p_sub, 
                        file="../figure_table/leukemia_sub_ctf.html",
                        selfcontained=F)
p_sub12 <- ctf_sub %>% ggplot() + 
  geom_point(aes(x=PC1, y=PC2, color=group)) +
  scale_color_manual(values=col_group) +
  theme_bw() +
  theme(legend.position="bottom")
p_sub13 <- ctf_sub %>% ggplot() + 
  geom_point(aes(x=PC1, y=PC3, color=group)) +
  scale_color_manual(values=col_group) + 
  theme_bw() + 
  theme(legend.position='bottom')
p_sub23 <- ctf_sub %>% ggplot() + 
  geom_point(aes(x=PC2, y=PC3, color=group)) +
  scale_color_manual(values=col_group) +
  theme_bw() + 
  theme(legend.position='bottom')
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend <- g_legend(p_sub12)
```

## CTF trajectories


```r
# function to read in trajectories
read_traj <- function(file){
  tmp <- tempdir()
  rm <- TRUE
  unzip(file, exdir=tmp) 
  unpacked<-unzip(file, exdir=tmp, list=TRUE)
  artifact<-read_yaml(paste0(tmp,"/", paste0(gsub("/..+","", unpacked$Name[1]),"/metadata.yaml")))
  artifact$contents<-data.frame(files=unpacked)
  artifact$contents$size=sapply(paste0(tmp, "/", artifact$contents$files), file.size)
  artifact$version=read.table(paste0(tmp,"/",artifact$uuid, "/VERSION"))
  artifact$data <- read_tsv(paste0(tmp,"/",artifact$uuid,"/data/trajectory.tsv"))
  return(artifact)
}
ctf_traj <- as.data.frame(read_traj("result/state-subject-ordination_leukemia.qza")$data)
```

```
## Rows: 327 Columns: 8
## ── Column specification ─────────────────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (3): #SampleID, subject_id, group
## dbl (5): PC1, PC2, PC3, week_int, week
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
ctf_traj$group <- meta_all[ctf_traj$`#SampleID`, "group"]
p_traj <- plot_feature_summary(ctf_traj[,c("PC1","PC2")], 
                     ctf_traj$week_int, 
                     ctf_traj$group,
                     bws=bws) + 
  xlab('Week') + theme_bw() + 
  theme(legend.position='bottom') +
  scale_color_manual(values=col_group) + 
  scale_fill_manual(values=col_group)
p_traj
```

![plot of chunk ctf_traj](figure/ctf_traj-1.png)


## CTF time loadings


```r
ctf_time <- read_qza("result/state-biplot_leukemia.qza")$data$Vectors
colnames(ctf_time)[1] <- "week_int"
tab_time <- data.frame(time=ctf_time$week_int, value=as.vector(as.matrix(ctf_time[,2:4])), 
                       component=as.factor(as.vector(t(matrix(rep(1:3,nrow(ctf_time)),3,)))))
p_time <- tab_time %>% ggplot(aes(x=time, y=value, color=component)) + 
  geom_line(size=1.5) + 
  scale_color_manual(values=col_microbe) + 
  labs(x='Week', y="Temporal Loading", color="Component") + 
  theme_bw() + 
  theme(legend.position='bottom')
p_time
```

![plot of chunk ctf_time](figure/ctf_time-1.png)

```r
mylegend2 <- g_legend(p_time)
```



```r
lay <- rbind(c(1,2,3,3), c(1,2,3,3), c(1,2,3,3), c(1,2,3,3), c(1,2,3,3), 
             c(4,5,5,5))
lgd1 <- get_legend(p_time)
lgd2 <- get_legend(p_sub23)
p_ctf_all <- grid.arrange(p_time+theme(legend.position="none"),
                          p_sub12+theme(legend.position="none"),
                          p_traj+theme(legend.position="none"),
                          lgd1, lgd2,
                          layout_matrix=lay)
```

![plot of chunk ctf_allplot](figure/ctf_allplot-1.png)

```r
ggsave("../figure_table/leukemia_ctf.pdf", width=10, height=3, plot=p_ctf_all)
```

# PCoA


```r
# get mds results
mds_mat_all <- NULL
for (i in 1:length(dm_list)){
  dist_tab <- dm_list[[i]]
  mds_all <- cmdscale(dist_tab, k=3)
  mds_mat <- data.frame(hostID=meta_all[rownames(mds_all),'hostID'],
                        PC=mds_all)
  mds_mat$timepoint <- meta_all[rownames(mds_all),'week']
  mds_mat <- merge(mds_mat, metauni, by='hostID')
  mds_mat$method <- names(dm_list)[i]
  mds_mat_all <- rbind(mds_mat_all, mds_mat)
}
mds_mat_all$group <- factor(mds_mat_all$group, 
                            levels=c("genotype=WT, disease=no", "genotype=Pax5+/-, disease=no",
                                     "genotype=Pax5+/-, disease=pB-ALL"))
colnames(mds_mat_all)
```

```
## [1] "hostID"    "PC.1"      "PC.2"      "PC.3"      "timepoint" "genotype"  "diseased" 
## [8] "group"     "method"
```

```r
mds_mat_all$method <- as.factor(mds_mat_all$method)
```


## Week 35-70 PERMANOVA for only Pax mice wrt disease


```r
# other methods
mds_mat_Pax <- dplyr::filter(mds_mat_all, genotype!='WT')
# result from tempted
ratio_mat_Pax <- dplyr::filter(reshape_feat_ratio, genotype!='WT')

mth <- levels(mds_mat_all$method)
pval_permanova_late <- rep(NA, length(mth)+2)
names(pval_permanova_late) <- c(mth, "TEMPTED", "CTF")
Fval_permanova_late <- pval_permanova_late
for (j in 1:length(mth)){
  tmp <- dplyr::filter(mds_mat_Pax, 
                       timepoint>=35 & timepoint<=70 & method==mth[j])
  res_adonis <- adonis2(tmp[,2:4] ~ diseased, method='euclidean',
                       data=tmp)
  pval_permanova_late[j] <- res_adonis$`Pr(>F)`[1]
  Fval_permanova_late[j] <- res_adonis$F[1]
}
# for TEMPTED
tmp <- dplyr::filter(ratio_mat_Pax, timepoint>=35 & timepoint<=70)
res_adonis <- adonis2(tmp[,6:8] ~ diseased, method='euclidean',
                     data=tmp)
pval_permanova_late[length(mth)+1] <- res_adonis$`Pr(>F)`[1]
Fval_permanova_late[length(mth)+1] <- res_adonis$F[1]
# for CTF
tmp <- dplyr::filter(ctf_traj, week_int>=35 & week_int<=70 & group!="genotype=WT, disease=no")
res_adonis <- adonis2(tmp[,2:4] ~ tmp$group, method='euclidean')
pval_permanova_late[length(mth)+2] <- res_adonis$`Pr(>F)`[1]
Fval_permanova_late[length(mth)+2] <- res_adonis$F[1]

cbind(pval_permanova_late,
      Fval_permanova_late)
```

```
##          pval_permanova_late Fval_permanova_late
## bray                   0.143           1.9960173
## jaccard                0.143           1.9300087
## unifrac                0.342           0.9555013
## wunifrac               0.798           0.1450785
## TEMPTED                0.001           6.9220221
## CTF                    0.032           3.7071993
```

```r
write.csv(cbind(pval_permanova_late,
                Fval_permanova_late), file="../figure_table/leukemia_permanova.csv")
```

## Trajectory plots for other methods and Wilcoxon rank test of Week 35-70 samples


```r
group_level <- levels(mds_mat_all$group)
nmethod <- length(dm_list)
CI_length <- -qnorm((1-0.95)/2)
tab_summary_all <- NULL
pval_wilcox_late <- Wval_wilcox_late <- matrix(NA, nmethod+2, npc)
rownames(pval_wilcox_late) <- rownames(Wval_wilcox_late) <- c(names(dm_list), "TEMPTED", "CTF")
for (kk in 1:nmethod){
  mthd <- names(dm_list)[kk]
  time_all <- NULL
  mean_all <- NULL
  merr_all <- NULL
  feature_all <- NULL
  group_all <- NULL
  ind_mthd <- mds_mat_all$method==mthd
  feature_mat <- mds_mat_all[ind_mthd,c(1+1:npc)]
  time_vec <- mds_mat_all$timepoint[ind_mthd]
  group_vec <- mds_mat_all$group[ind_mthd]
  
  # Wilcoxon test
  for (jj in 1:ncol(feature_mat)){
    sel <- time_vec>=35 & time_vec<=70 & 
    group_vec!='genotype=WT, disease=no'
    table(group_vec[sel])
    test_tmp <- wilcox.test(feature_mat[sel,jj]~group_vec[sel])
    pval_wilcox_late[kk,jj] <- test_tmp$p.value
    Wval_wilcox_late[kk,jj] <- test_tmp$statistic
  }
  for (jj in 1:ncol(feature_mat)){
    for (ii in 1:length(group_level)){
      ind <- group_vec==group_level[ii]
      model.np <- npreg(feature_mat[ind,jj]~time_vec[ind], bws=bws,
                          regtyle="ll", bwmethod="cv.aic", gradients=T)
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
  tab_summary$method <- names(dm_list)[kk]
  tab_summary_all <- rbind(tab_summary_all, tab_summary)
}
colnames(pval_wilcox_late) <- colnames(Wval_wilcox_late) <- colnames(feature_mat)
pval_wilcox_late <- as.data.frame(pval_wilcox_late)
Wval_wilcox_late <- as.data.frame(Wval_wilcox_late)

# for TEMPTED
sel <- time_vec_ratio>=35 & time_vec_ratio<=70 & 
  group_vec_ratio!='genotype=WT, disease=no'
wilcox_tempted <- wilcox.test(feature_mat_ratio$`Component 2`[sel]~
                             group_vec_ratio[sel])
pval_wilcox_late["TEMPTED","PC.2"] <- wilcox_tempted$p.value
Wval_wilcox_late["TEMPTED","PC.2"] <- wilcox_tempted$statistic
# for CTF
ctf_traj_late <- dplyr::filter(ctf_traj, week_int>=35 & week_int<=70 & 
  group!='genotype=WT, disease=no')
wilcox_ctf <-  wilcox.test(ctf_traj_late$PC1~ctf_traj_late$group)
pval_wilcox_late["CTF","PC.1"] <- wilcox_ctf$p.value
Wval_wilcox_late["CTF","PC.1"] <- wilcox_ctf$statistic

tab_summary_all$method <- factor(tab_summary_all$method, 
                                 levels=unique(tab_summary_all$method))
tab_summary_all$feature <- factor(tab_summary_all$feature, 
                                  levels=unique(tab_summary_all$feature))
p_summary <- ggplot(data=tab_summary_all, 
                    aes(x=time, y=mean, group=group, color=group)) + 
  geom_line() + theme_bw()+ 
  geom_ribbon(aes(ymin=mean-CI_length*merr, ymax=mean+CI_length*merr, 
                  color=group, fill=group), linetype=2, alpha=0.3) + 
  ylab(paste0('mean +/- ', round(CI_length,2), '*se')) + xlab("Week") +
  facet_grid(feature~method) + 
  theme(legend.position="bottom") + 
  scale_color_manual(values=col_group) + 
  scale_fill_manual(values=col_group) + 
  scale_x_continuous(breaks=c(0,25,50,75,100))
p_summary
```

![plot of chunk traj_all](figure/traj_all-1.png)

```r
pdf("../figure_table/leukemia_mds.pdf", width=8, height=5)
print(p_summary)
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
# Wilcoxon rank test
pval_wilcox_late
```

```
##                 PC.1         PC.2      PC.3
## unifrac  0.135431301 0.1310596149 0.8488535
## wunifrac 0.672460669 0.6130066586 0.7591171
## bray     0.006876426 0.6130066586 0.4910435
## jaccard  0.005848879 0.5229909514 0.4303680
## TEMPTED           NA 0.0000156904        NA
## CTF      0.029485574           NA        NA
```

```r
write.csv(pval_wilcox_late, '../figure_table/leukemia_wilcox_pcoa.csv')
```



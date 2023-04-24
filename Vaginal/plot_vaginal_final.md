---
title: "Plot Result from Vaginal Microbiome Study"
author: "Pixu Shi"
date: "2023-04-19"
output: html_document
---

## Preparation


```r
rm(list=ls())
library(gridExtra)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)

# set working directory to be where the current script is located
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


source('../TEMPTED.R')
color_RB <- brewer.pal(3,'Set1')[1:2]
```


## Plot feature loading bar plot


```r
tb_vaginal <- read.csv("result/Vaginal_microbiome_feature_loadings_w_taxonomy.csv", row.names=1, sep=",")

plot(tb_vaginal$Component.1, tb_vaginal$Component.2)
```

![plot of chunk feature_loading](figure/feature_loading-1.png)

```r
sel <- tb_vaginal$Component.1 %>% (function(x){which(abs(x)>=quantile(abs(x),0.95))})
tb_vaginal_sel <- tb_vaginal[sel,] %>% mutate(sign=ifelse(Component.1<0,"Preterm","Term"))
# get species labels, if NA, use family labels
label_s <- substring(tb_vaginal_sel$Species, 4, 30)
label_f <- substring(tb_vaginal_sel$Family, 4, 30)
label_s[which(nchar(label_s)==0)] <- label_f[which(nchar(label_s)==0)]
tb_vaginal_sel$label <- paste0(label_s,
                               " PG", substring(rownames(tb_vaginal_sel),5,6))
tb_vaginal_sel <- tb_vaginal_sel %>% arrange(-Component.1)
tb_vaginal_sel$label <- factor(tb_vaginal_sel$label, levels=unique(tb_vaginal_sel$label))


p_bar <- tb_vaginal_sel %>% 
  ggplot(aes(x=label, y=Component.1)) + 
  geom_bar(stat="identity", colour="black", width = 0.5, aes(fill = sign)) + 
  geom_hline(yintercept = 0) +
  theme_bw() +
  scale_fill_manual(values=color_RB) +
  coord_flip() + 
  labs(y = "Component 1", x = "") + 
  theme(legend.position="none")

p_bar
```

![plot of chunk feature_loading](figure/feature_loading-2.png)




## Plot log ratio trajectory plot


```r
ratio_vaginal <- read.csv("result/Vaginal_microbiome_log_ratios_per_sample_Pixu.csv", row.names=1, sep=",")

ratio_vaginal <- filter(ratio_vaginal, !is.na(Log_ratios))
ratio_vaginal$Group2 <- ifelse(ratio_vaginal$Group=="TB", "Term", "Preterm")
with(ratio_vaginal, table(Group, Group2))
```

```
##      Group2
## Group Preterm Term
##   PTB     126    0
##   TB        0  357
```

```r
time_vec_ratio <- ratio_vaginal$Timepoint
group_vec_ratio <- factor(ratio_vaginal$Group2)
feature_mat_ratio <- ratio_vaginal[,4, drop=FALSE]
p_feat_ratio_summary <- plot_feature_summary(feature_mat_ratio, 
                                             time_vec_ratio, 
                                             group_vec_ratio, bws=8, nrow=1) + 
  xlab('Gestational Week') + theme_bw() + 
  theme(legend.position=c(0.85,0.15), 
        legend.background = element_rect(fill="white", size=0.2, linetype="solid", 
                                         colour ="black")) + 
  scale_fill_manual("Outcome", values=color_RB) +
  scale_color_manual("Outcome", values=color_RB) +
  theme(strip.text = element_blank(), 
        strip.background = element_blank())
p_feat_ratio_summary
```

![plot of chunk log_ratio](figure/log_ratio-1.png)



```r
# plot timeline
meta_ordered <- ratio_vaginal %>%
  arrange(Group2, -Timepoint) %>% 
  mutate(PregID=factor(subjects, levels=unique(subjects)))
p_timeline <- ggplot(data=meta_ordered, 
       aes(x=Timepoint, y=PregID, 
           PregOut=as.factor(PregID), color=Group2, shape=Group2)) +
  geom_line() + 
  geom_point()  + 
  scale_color_manual(values=color_RB) +
  labs(y="Host ID", x="Gestational Week", color="Outcome", shape="Outcome") + 
  theme_bw() +
  theme(legend.position="bottom") + 
  scale_x_continuous(breaks=seq(from=0,to=36,by=4)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
pdf('../figure_table/vaginal_timeline.pdf', width=11.5, height=7.5)
print(p_timeline)
dev.off()
```

```
## quartz_off_screen 
##                 2
```




## plot boxplot of subject loading


```r
sub_vaginal <- read.csv("result/Vaginal_microbiome_subject_loadings.csv",
                        row.names=1)
sub_vaginal$group2 <- ifelse(sub_vaginal$group=="TB", "Term", "Preterm")
with(sub_vaginal, table(group, group2))
```

```
##      group2
## group Preterm Term
##   PTB      35    0
##   TB        0  111
```

```r
pval <- round(wilcox.test(sub_vaginal$Component.1~sub_vaginal$group2)$p.value,4)
pval
```

```
## [1] 0.006
```

```r
p_sub_box <- ggplot(sub_vaginal) + 
  geom_boxplot(aes(x=group2, y=Component.1, fill=group2)) + 
  theme_bw() + scale_fill_manual(values=color_RB) +
  labs(x="Outcome") + 
  theme(legend.position="none") +
  annotate("text", x=1.5, y=0.22, label= paste0("p-value = ", pval))
p_sub_box
```

![plot of chunk subj_loading](figure/subj_loading-1.png)

```r
lay <- matrix(c(1,1,2,2,2,2,3,3,3),1,)
pdf('../figure_table/vaginal_plot.pdf', width=16, height=3.5)
grid.arrange(p_sub_box, p_bar, p_feat_ratio_summary, layout_matrix=lay)
dev.off()
```

```
## quartz_off_screen 
##                 2
```


## plot simulation results for PEMANOVA F


```r
tab_Fmodel <- read.csv("result/perm-res.csv", row.names=1)
tab_Fmodel$metric[tab_Fmodel$metric=="FT-SVD"] <- "TEMPTED"
tab_Fmodel <- dplyr::filter(tab_Fmodel, seq_depth>=50000 & metric!="Jaccard")
tab_Fmodel$seq_depth <- paste0("Seq. Depth = ", tab_Fmodel$seq_depth/1000, "K")
tab_Fmodel$seq_depth <- factor(tab_Fmodel$seq_depth, levels=
                                 c("Seq. Depth = 50K", "Seq. Depth = 100K", "Seq. Depth = 150K"))

tab1 <- aggregate(test.statistic~n_samples+metric+seq_depth, data=tab_Fmodel, 
                  FUN=mean)
names(tab1)[4] <- 'mean'
tab2 <- aggregate(test.statistic~n_samples+metric+seq_depth, data=tab_Fmodel, 
                  FUN=function(x){sd(x)/sqrt(length(x))})
names(tab2)[4] <- 'sd'
rownames(tab1) <-rownames(tab2) <- NULL
tab_Fmodel_summary <- merge(tab1, tab2)

color_method <- c('#ff7f00', #orange for Bray
                  '#33a02c', #green for CTF
                  '#377eb8') #blue for TEMPTED
p_Fmodel_summary <- ggplot(data=tab_Fmodel_summary, 
                           aes(x=n_samples, y=mean, group=metric, color=metric)) + 
  geom_line(size=1) + geom_point(position=position_dodge(0.2), size=2) +
  geom_errorbar(aes(ymin=mean-2*sd, ymax=mean+2*sd, width=0.5), 
                position=position_dodge(0.2), size=1) + 
  facet_grid(~seq_depth) +
  scale_x_continuous(breaks=2:10) + 
  scale_color_manual(values=color_method) +
  ylab('PERMANOVA F value') + xlab('# time points') + 
  theme_bw() + 
  theme(legend.position = "bottom")
p_Fmodel_summary
```

![plot of chunk simulation](figure/simulation-1.png)

```r
pdf("../figure_table/realsim_vaginal_Fvalue.pdf", width=10, height=4)
p_Fmodel_summary
dev.off()
```

```
## quartz_off_screen 
##                 2
```



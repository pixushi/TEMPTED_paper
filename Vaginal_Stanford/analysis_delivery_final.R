rm(list=ls())


library(readr)
library(qiime2R)
library(phyloseq)
library(gridExtra)
library(ggplot2)
library(vegan)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
#library(magrittr)
library(pROC)

# set working directory to be where the current script is located
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source('../TEMPTED.R')
color_RB <- brewer.pal(3,'Set1')[1:2]

##############################################
#### read in the count data and meta data ####

# Get data in 2nd and 3rd trimester
load("data/data_post.RData")
metadata_all <- post.df %>% 
  filter(DayVsDel < 0, LibrarySize >= 4e4, PregOut!="Miscarriage") %>%
  select(Cohort, SampleID, PregID, DayVsConc, DayVsDel, PregOut, GDdel) %>%
  group_by(PregID) %>% filter(n() > 1) %>% ungroup()
head(metadata_all)
metadata <- metadata_all %>% filter(DayVsConc >= 12*7, DayVsConc < 37*7)
ps_filtered <- prune_samples(metadata$SampleID, post.ps)
ps_filtered <- filter_taxa(ps_filtered, function(x) sum(x!=0) >= 5, TRUE)
count_tab <- data.frame(otu_table(ps_filtered))
meta_tab <- data.frame(sample_data(ps_filtered))

#######################
#### plot timeline ####

meta_ordered <- metadata_all %>%
  arrange(PregOut, GDdel) %>% 
  mutate(PregID=factor(PregID, levels=unique(PregID)))
pdf('../figure_table/vaginal_Stanford_timeline.pdf', width=11.5, height=7.5)
ggplot(data=meta_ordered, 
       aes(x=DayVsConc, y=PregID, 
           PregOut=as.factor(PregID), color=PregOut, shape=PregOut)) +
  geom_line() + 
  geom_point()  + 
  scale_color_manual(values=color_RB) +
  labs(y="Host ID", x="Gestational Day", color="Outcome", shape="Outcome") + 
  theme_bw() +
  theme(legend.position="bottom") + 
  scale_x_continuous(breaks=seq(from=84,to=294,by=28)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  geom_vline(xintercept=c(12,37)*7, linetype="dashed")
dev.off()


datlist_conc <- format_tempted(count_tab, meta_tab$DayVsConc, meta_tab$PregID, 
                              threshold=0.95, pseudo_count=0.5, transform='clr')
datlist_raw_conc <- format_tempted(count_tab, meta_tab$DayVsConc, meta_tab$PregID, 
                                  threshold=0.95, transform='none')
npc <- 2
svd_conc <- svd_centralize(datlist_conc, 1)
res_tempted_conc <- tempted(svd_conc$datlist, r = npc, resolution = 51, smooth=1e-3)


save(ps_filtered, 
     svd_conc, res_tempted_conc, datlist_conc, datlist_raw_conc,
     file='result/res_delivery_tempted.Rdata')


load('result/res_delivery_tempted.Rdata')
count_tab <- data.frame(otu_table(ps_filtered))
meta_tab <- data.frame(sample_data(ps_filtered))
res_tempted <- res_tempted_conc
svd_tempted <- svd_conc
datlist <- datlist_conc
datlist_raw <- datlist_raw_conc
npc <- 2


###########################
#### plot time loading ####


p_time <- plot_time_loading(res_tempted, r=npc) +
  labs(x='Week', y="Temporal Loading", title=" ", color="Component") + 
  geom_line(size=1.5) + theme_bw() + 
  theme(legend.position='bottom')
p_time

###############################
#### plot subject loadings ####

metauni <- unique(meta_tab[,c('PregID', 'PregOut')])
rownames(metauni) <- metauni$PregID
A.hat <- metauni
rownames(A.hat) <- A.hat$PregID
table(rownames(A.hat)==rownames(res_tempted$A.hat))
A.hat <- cbind(res_tempted$A.hat[,1:npc], A.hat)
A.hat$PregOut <- as.factor(A.hat$PregOut)

pval <- round(wilcox.test(A.hat$`Component 1`~A.hat$PregOut)$p.value,4)
pval
p_sub_box <- ggplot(A.hat) + 
  geom_boxplot(aes(x=PregOut, y=`Component 1`, fill=PregOut)) + 
  theme_bw() + scale_fill_manual(values=color_RB) + 
  labs(x="Outcome") + 
  theme(legend.position="none") +
  annotate("text", x=1.5, y=0.15, label= paste0("p-value = ", pval))
p_sub_box

res_glm <- with(A.hat, glm(I(PregOut=='Term') ~`Component 1`+`Component 2`,
                           family=binomial(link = "logit")))
auc(res_glm$y, predict.glm(res_glm, type="response"))



#########################################################################
#### plot trajectory of log ratio of total abundance of top features ####

contrast <- NULL
ratio_feat <- ratio_feature(res_tempted, datlist_raw,
                            contrast=contrast, pct=0.05,
                            absolute=TRUE)

## summed up, by individual subject

tab_feat_ratio <- ratio_feat$metafeature.ratio
colnames(tab_feat_ratio)[2] <- 'PregID'
tab_feat_ratio <- merge(tab_feat_ratio, metauni)

## summed up, by mean and sd
reshape_feat_ratio <- reshape(tab_feat_ratio, 
                              idvar=c("PregID","timepoint") , 
                              v.names=c("value"), timevar="PC",
                              direction="wide")
CC <- grep('Component',colnames(reshape_feat_ratio))
colnames(reshape_feat_ratio)[CC] <- paste('Component', 1:npc)
feature_mat_ratio <- reshape_feat_ratio[,CC]
colnames(feature_mat_ratio) <- paste("Component", 1:npc)


time_vec_ratio <- reshape_feat_ratio$timepoint
group_vec_ratio <- factor(reshape_feat_ratio$PregOut)
p_feat_ratio_summary <- plot_feature_summary(feature_mat_ratio, 
                                             time_vec_ratio, 
                                             group_vec_ratio, bws=12, nrow=1)
feature_mat_ratio1 <- feature_mat_ratio[,1,drop=FALSE]
colnames(feature_mat_ratio1) <- "Component 1"
p_feat_ratio_summary <- plot_feature_summary(feature_mat_ratio1, 
                                             time_vec_ratio, 
                                             group_vec_ratio, bws=12, nrow=1) + 
  xlab('Gestational Day') + theme_bw() + 
  theme(legend.position=c(0.85,0.15), 
        legend.background = element_rect(fill="white", size=0.2, linetype="solid", 
                                                                     colour ="black")) + 
  scale_fill_manual("Outcome", values=color_RB) +
  scale_color_manual("Outcome", values=color_RB) +
  scale_x_continuous(breaks=seq(from=84,to=294,by=28)) + ylim(c(-3.5, 7.5))+
  theme(strip.text = element_blank(), 
        strip.background = element_blank())
p_feat_ratio_summary


##########################################
#### bar plot of top feature loading  ####

# taxanomic information of top and bottom bacteria
sum(ratio_feat$toppct[,1])
sum(ratio_feat$bottompct[,1])
tab_taxa <- rbind(tax_table(ps_filtered)[names(which(ratio_feat$toppct[,1])),],
tax_table(ps_filtered)[names(which(ratio_feat$bottompct[,1])),])
tab_taxa <-as.data.frame(tab_taxa)
tab_taxa$loading_PC1 <- res_tempted$B.hat[rownames(tab_taxa),1]
write.csv(tab_taxa, file="../figure_table/vaginal_Stanford_topfeat.csv")
rownames(tab_taxa) <- NULL
tab_taxa

# barplot
tab_taxa <- tab_taxa %>% mutate(outcome=ifelse(loading_PC1<0,"Preterm","Term"))
tab_taxa$label <- paste0(substring(tab_taxa$asvLab, 6, 50), " ",
                               tab_taxa$asvNum)
tab_taxa <- tab_taxa %>% arrange(-loading_PC1)
tab_taxa$label <- factor(tab_taxa$label, levels=unique(tab_taxa$label))

p_bar <- tab_taxa %>% 
  ggplot(aes(x=label, y=loading_PC1)) + 
  geom_bar(stat="identity", colour="black", width = 0.5, aes(fill = outcome)) + 
  geom_hline(yintercept = 0) +
  theme_bw() +
  scale_fill_manual(values=color_RB) +
  coord_flip() + 
  labs(y = "Component 1", x = "") + 
  theme(legend.position="none")

p_bar

lay <- matrix(c(1,1,2,2,2,2,3,3,3),1,)
pdf('../figure_table/vaginal_Stanford_plot.pdf', width=16, height=3.5)
grid.arrange(p_sub_box, p_bar, p_feat_ratio_summary, layout_matrix=lay)
dev.off()

rm(list=ls())
library(dirmult)
library(ggplot2)
library(reticulate)
library(readr)
library(gridExtra)
library(ggpubr)
library(grid)
library(RColorBrewer)
library(tempted)

# set working directory to be where the current script is located
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


col_group <- brewer.pal(7,'Set1')[c(7,5,3)]
col_microbe <- c(brewer.pal(6,'Set1')[c(2,4)], brewer.pal(6,'Set2')[4])
col_microbe <- c(col_microbe[c(1,3,2)], '#999999')

#### function for simulation dirichlet multinomial distribution ####

rDM <- function(size, mu, overdisp){
  # size is a length n vector, mu is n by p matrix, overdisp is scalar
  n <- nrow(mu)
  size <- round(size)
  if (overdisp>0){
    sum_alpha <- 1/overdisp
    alpha <- mu*sum_alpha  # when overdisp goes to 0, it becomes multinomial
    probmat <- t(mapply(function(i){rdirichlet(1, alpha[i,])}, 1:n))
  }else if (overdisp==0){
    probmat <- mu
  }else{
    stop('overdisp needs to be no smaller than 1!')
  }
  #probmat <- apply(alpha, 1:2, function(a){rgamma(1,a)})
  #probmat <- probmat/rowSums(probmat)
  xx <- t(mapply(function(i){rmultinom(1, size[i], probmat[i,])}, 1:n))
  return(xx)
}

#### the three temporal trends to simulate ####

x <- (0:1000)/1000
# oscillating function
wfunc1 <- function(x){
  return(0.5*sin(x*pi*4-pi/2)+0.5)
}
plot(x, wfunc1(x))

wfunc2 <- function(x){
  return(atan((x-0.5)*10)/pi+atan((0.5)*10)/pi)
}
plot(x, wfunc2(x))

wfunc3 <- function(x){
  return(0.5*sin(x*pi*2-pi/2)+0.5)
}
plot(x, wfunc3(x))


#### Simulate data for each group ####

set.seed(0)
nsub <- 20
ntime <- 20
nfeat <- 100

# relative abundance for group 3
pmat_all <- NULL
sub_all <- NULL
time_all <- NULL
for (i in 1:nsub){
  x <- runif(ntime)
  pmat <- matrix(0, ntime, nfeat)
  pmat[,1] <- wfunc1(x)+runif(1,0,0.5)
  pmat[,2] <- wfunc2(x)+runif(1,0,0.25)
  for (j in 3:nfeat){
    pmat[,j] <- runif(ntime, 0, 0.25)
  }
  pmat <- pmat/rowSums(pmat)
  pmat_all <- rbind(pmat_all, pmat)
  time_all <- c(time_all, x)
  sub_all <- c(sub_all, rep(paste0('sub',i), ntime))
}
df_signal1 <- data.frame(microbe=pmat_all, time=time_all, subject=sub_all)
df_signal1$group <- 'Group3'

# relative abundance for group 2
pmat_all <- NULL
sub_all <- NULL
time_all <- NULL
for (i in 1:nsub){
  x <- sort(runif(ntime))
  pmat <- matrix(0, ntime, nfeat)
  pmat[,1] <- wfunc1(x)+runif(1,0,0.5)
  pmat[,2] <- runif(ntime, 0, 0.25)
  for (j in 3:nfeat){
    pmat[,j] <- runif(ntime, 0, 0.25)
  }
  pmat <- pmat/rowSums(pmat)
  pmat_all <- rbind(pmat_all, pmat)
  time_all <- c(time_all, x)
  sub_all <- c(sub_all, rep(paste0('sub',i+nsub), ntime))
}
df_signal2 <- data.frame(microbe=pmat_all, time=time_all, subject=sub_all)
df_signal2$group <- 'Group2'

# relative abundance for group 1
pmat_all <- NULL
sub_all <- NULL
time_all <- NULL
for (i in 1:nsub){
  x <- sort(runif(ntime))
  pmat <- matrix(0, ntime, nfeat)
  pmat[,1] <- wfunc1(x)+runif(1,0,0.5)
  pmat[,2] <- wfunc2(x)+runif(1,0,0.25)
  pmat[,3] <- wfunc3(x)+runif(1,0,0.5)
  for (j in 4:nfeat){
    pmat[,j] <- runif(ntime, 0, 0.25)
  }
  pmat <- pmat/rowSums(pmat)
  pmat_all <- rbind(pmat_all, pmat)
  time_all <- c(time_all, x)
  sub_all <- c(sub_all, rep(paste0('sub',i+2*nsub), ntime))
}
df_signal3 <- data.frame(microbe=pmat_all, time=time_all, subject=sub_all)
df_signal3$group <- 'Group1'

df_signal <- rbind(df_signal1, df_signal2, df_signal3)

# simulating sequencing depth
depth <- round(exp(rnorm(nrow(df_signal), 10, 0.5)))

# simulating counts
df_count <- rDM(size=depth, mu=as.matrix(df_signal[,1:nfeat]), overdisp=0)

# gather metadata
colnames(df_count) <- paste0('Microbe',1:ncol(df_count))
df_count[,c(1,3)] <- df_count[,c(3,1)]
df_meta <- df_signal[,c('time', 'subject', 'group')]

write.csv(df_count, file='toy_example_count.csv')
write.csv(df_meta, file='toy_example_meta.csv')
write.csv(df_signal, file='toy_example_signal.csv')


df_count <- read.csv('toy_example_count.csv', row.names=1)
df_meta <- read.csv('toy_example_meta.csv', row.names=1)
df_signal <- read.csv('toy_example_signal.csv', row.names=1)
# original data
transform_count_tab <- df_count/rowSums(df_count)
meta_tab <- df_meta

tab_ori <- cbind(transform_count_tab[,1:4], meta_tab)
tab_ori <- reshape(tab_ori, direction="long",
                   varying=list(1:4), v.names="abundance",
                   timevar="Microbe")
tab_ori$Microbe <- paste0('Microbe',tab_ori$Microbe)
tab_ori$Microbe[tab_ori$Microbe=="Microbe4"] <- "1 of 97 Others"
tab_ori$Microbe <- factor(tab_ori$Microbe, levels=unique(tab_ori$Microbe))
p_ori <- ggplot(tab_ori) +
  geom_point(aes(x=time, y=abundance, color=Microbe), size=0.3) +
  facet_grid(Microbe~group) +
  scale_y_continuous(breaks=seq(0, 0.1, 0.05)) +
  scale_color_manual(values=col_microbe) +
  labs(x='', y='Observed Relative Abundance', color='Microbe') +
  theme_bw() +
  theme(panel.spacing.y = unit(0, "lines"), panel.spacing.x = unit(1, "lines"),
        legend.position="none", axis.text.x=element_blank(),
        plot.margin=margin(0.2,0.2,0,0.2,'cm'))

p_ori



p_ori2 <- ggplot_gtable(ggplot_build(p_ori))
strips <- which(grepl('strip-', p_ori2$layout$name))

for (i in 1:3) {
  k <- which(grepl('rect', p_ori2$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', p_ori2$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  p_ori2$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- col_group[i]
  #g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$col <- pal[i + 1]
}
for (i in 4:7) {
  k <- which(grepl('rect', p_ori2$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', p_ori2$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  p_ori2$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- col_microbe[i-3]
  #g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$col <- pal[i + 1]
}

plot(p_ori2)

#p_ori_summary

meta_tab2 <- meta_tab
meta_tab2$subject <- as.character(as.numeric(substring(meta_tab2$subject,4))%%20)
meta_tab2$microbe <- "Timeline"

p_timeline <- ggplot(data=meta_tab2,
                     aes(x=time, y=subject,
                         group=subject, color=group)) +
  geom_line(size=0.15) + geom_point(size=0.3) +
  scale_color_manual(values=col_group) +
  labs(y="Subjects",x="Time",color="Group") + facet_grid(microbe~group) +
  theme_bw() +
  theme(axis.text.y=element_blank(), legend.position='none',
        panel.spacing.y = unit(0, "lines"), panel.spacing.x = unit(1, "lines"),
        strip.text.x=element_blank(), strip.background.y = element_rect("#FFFFFF"),
        plot.margin=margin(-0.4,0.2,0,0.88,'cm'),
        panel.grid.major = element_blank(), panel.grid.minor=element_blank())



pdf('../figure_table/toyexample_signal.pdf', width=6, height=5)
lay <- rbind(c(1,1),
             c(1,1),
             c(1,1),
             c(1,1),
             c(2,2))
grid.arrange(p_ori2, p_timeline, nrow=2, layout_matrix=lay)
dev.off()


count_tab <- df_count
meta_tab <- df_meta
datlist <- format_tempted(count_tab, meta_tab$time, meta_tab$subject,
                        threshold=1, pseudo=0.5, transform='clr')
npc <- 3

res_svd <- svd_centralize(datlist, 1)
res_tempted <- tempted(res_svd$datlist, r=npc, resolution=101, smooth=1e-4)
res_tempted$A_hat[,2] <- -res_tempted$A_hat[,2]
res_tempted$B_hat[,2] <- -res_tempted$B_hat[,2]
save(res_tempted, file="toy_example.Rdata")

load("toy_example.Rdata")

# time loading
p_time <- plot_time_loading(res_tempted)+ theme_bw() +
     labs(x="Time", y="Temporal Loading", color="Component") +
     scale_color_manual(values=col_microbe) +
     theme(legend.position="bottom") +
     geom_line(size=1.5)
p_time

# subject loading
metauni <- unique(meta_tab[,c('subject', 'group')])
A_hat <- metauni
rownames(A_hat) <- A_hat$subject
table(rownames(A_hat)==rownames(res_tempted$A_hat))
A_hat <- cbind(res_tempted$A_hat, A_hat)
tab_A <- reshape(A_hat, direction="long", idvar="subject",
                 v.names="value", timevar="Component",
                 varying=list(1:3))
tab_A$subject <- as.numeric(substring(tab_A$subject,4))
tab_A$order <- NA
for (j in 1:3){
  for(i in 1:3){
    id <- tab_A$group==paste0('Group',i) & tab_A$Component==j
    tab_A[id,]$order <- order(order(tab_A[id,]$value))+20*(i-1)
  }
}
tab_A$Component <- paste("Component", tab_A$Component)



p_sub <- ggplot(data=tab_A, aes(y=value,x=order, fill=group)) +
  geom_bar(stat = "identity") + scale_fill_manual(values=col_group) +
  theme_bw() + theme(legend.position="bottom") +
  labs(x="Subjects Ordered by Loading", y="Subject Loading", fill="Group") +
  facet_grid(Component~.)
p_sub

## feature loading
col_feat <- col_microbe[c(1:4,rep(4,96))]
B_hat <- data.frame(res_tempted$B_hat, rownames(res_tempted$B_hat))
colnames(B_hat)[npc+1] <- 'microbe_name'
topfeat <- c(rep(TRUE, 3), rep(FALSE, nfeat-3))
B_hat$microbe_name[!topfeat] <- "Others"
B_hat$rownum <- as.character(1:nrow(B_hat))

tab_B <- reshape(B_hat, direction="long", idvar="rownum",
                 v.names="value", timevar="Component",
                 varying=list(1:3))
tab_B$rownum <- as.numeric(tab_B$rownum)
tab_B$Component <- paste("Component", tab_B$Component)

p_feature <- ggplot(data=tab_B, aes(x=rownum, y=value, fill=microbe_name)) +
  geom_bar(stat="identity") + scale_fill_manual(values=col_microbe) +
  theme_bw() + theme(legend.position="bottom") +
  labs(x="Microbes", y="Feature Loading", fill="Microbe") +
  facet_grid(Component~.)
p_feature

pdf('../figure_table/toyexample_TEMPTED.pdf', width=13, height=4)
grid.arrange(p_sub, p_feature, p_time, nrow=1)
dev.off()


library(writexl)
res_format <- res_tempted
for (i in 1:length(res_format)){
  res_format[[i]] <- as.data.frame(res_format[[i]])
  if (i<=2){
    res_format[[i]] <- cbind(" "=rownames(res_format[[i]]),
                             res_format[[i]])
  }
  if (i>=4){colnames(res_format[[i]]) <- NULL}
}
write_xlsx(res_format, path="../figure_table/toyexample_TEMPTED.xlsx")

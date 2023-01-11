library(readr)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(rtracklayer)
library(stringr)
library(immunedeconv)
library(RColorBrewer)
library(tibble)
library(DelayedArray)
library(circlize)
library(reshape2)
library(ggsci)
library(ggpubr)
cnt_df <- read_delim("data/IPF/IPF_counts.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
summary(cnt_df)
colnames(cnt_df) <- c("GeneName",paste0(rep("IPF", 67),1:67),
                      paste0(rep("CTRL", 63),1:63))
pheno <- data.frame(type=c(rep("IPF", 67),rep("CTRL", 63)),
                       sampleID=c(paste0(rep("IPF", 67),1:67),
                                  paste0(rep("CTRL", 63),1:63)))
#TPM
cnt <- cnt_df[, 2:ncol(cnt_df)]
rownames(cnt) <- cnt_df$GeneName
cnt_mat <- as.matrix(cnt)

load("immunedeconv/effLen.Rdata")

countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}
expr_cnt <- na.omit(left_join(effLen,cnt_df))
Len <- expr_cnt[,c(1:4)]
cnt_matrix <- as.matrix(expr_cnt[5:ncol(expr_cnt)])
rownames(cnt_matrix) <- expr_cnt$GeneName
ipftpms <- countToTpm(counts = cnt_matrix,effLen = Len$Lenkb)
write.table(ipftpms,file="immunedeconv/ipftpmkb.txt",sep="\t",quote=F)

##CIBERSORT
source('immunedeconv/CIBERSORT.R')
ipf_cibersort <- CIBERSORT('immunedeconv/LM22.txt','immunedeconv/ipftpmkb.txt', perm = 100, QN = F) 

write.table(ipf_cibersort, file="immunedeconv/ipf_cibersort.txt", sep="\t", col.names=T, row.names=F, quote=F)
save(ipf_cibersort, file="immunedeconv/ipf_cibersort.Rdata")
load("immunedeconv/ipf_cibersort.Rdata")



group_list=c(rep('IPF',67),rep('control',63))
group_list=factor(group_list) 
group_list
data_m <- melt(ipf_cibersort[,1:22])
head(data_m)
data_m$group <- group_list                  
head(data_m)
colnames(data_m)=c("ID","Celltype","value","group")
head(data_m)

plot_order = data_m[data_m$group == "IPF",] %>% 
  group_by(Celltype) %>% 
  summarise(m = median(value)) %>% 
  arrange(desc(m)) %>% 
  pull(Celltype)
data_m$Celltype = factor(data_m$Celltype,levels = plot_order)
library(ggplot2)
if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"), 
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 45, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  ) }

box_TME <- ggplot(data_m, aes(x = Celltype, y = value))+ 
  labs(y="Cell composition",x= NULL,title = "TME Cell composition")+  
  geom_boxplot(aes(fill = group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
  scale_fill_manual(values = c("#1CB4B8", "#EB7369"))+
  theme_test() + mytheme + 
  stat_compare_means(aes(group = group),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = T)
box_TME
ggsave(box_TME,file="immunedeconv/ipf_TME.pdf",width =13.5,height = 6)

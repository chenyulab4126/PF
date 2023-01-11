library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(stringr)
library(readr)
library(dplyr)

load("data/covres.RData")
cov_MMPs <- str_subset(res_df$GeneName,"^MMP")
cov_TIMPs <- str_subset(res_df$GeneName,"^TIMP")
cov_Genes<- as.data.frame(c(cov_MMPs,cov_TIMPs))
colnames(cov_Genes) <- c("GeneName")

cov_expr_cnt <- left_join(cov_Genes,res_df)

cov_genes <- cov_expr_cnt$GeneName
write.table(cov_genes,file = "cov/genelist.txt",
            row.names = F,col.names = F,quote = F)
cov_deg <- cov_expr_cnt[,c(1,8,12)]
if (! dir.exists('cov/cytoscape')) dir.create('cov/cytoscape')
write.csv(cov_deg,file = "cov/cytoscape/DEG.csv",
            row.names = F,quote = F)

load("data/sleres.RData")
sle_MMPs <- str_subset(res_df$GeneName,"^MMP")
sle_TIMPs <- str_subset(res_df$GeneName,"^TIMP")
sle_Genes<- as.data.frame(c(sle_MMPs,sle_TIMPs))
colnames(sle_Genes) <- c("GeneName")
sle_expr_cnt <- left_join(sle_Genes,res_df)
sle_genes <- sle_expr_cnt$GeneName
write.table(sle_genes,file = "sle/genelist.txt",
            row.names = F,col.names = F,quote = F)
sle_deg <- sle_expr_cnt[,c(1,8,12)]
if (! dir.exists('sle/cytoscape')) dir.create('sle/cytoscape')
write.csv(sle_deg,file = "sle/cytoscape/DEG.csv",
          row.names = F,quote = F)

load("data/sleres.RData")
sle_MMPs <- str_subset(res_df$GeneName,"^MMP")
sle_TIMPs <- str_subset(res_df$GeneName,"^TIMP")
sle_Genes<- as.data.frame(c(sle_MMPs,sle_TIMPs))
colnames(sle_Genes) <- c("GeneName")
sle_expr_cnt <- left_join(sle_Genes,res_df)
sle_genes <- sle_expr_cnt$GeneName
write.table(sle_genes,file = "sle/genelist.txt",
            row.names = F,col.names = F,quote = F)
sle_deg <- sle_expr_cnt[,c(1,8,12)]
if (! dir.exists('sle/cytoscape')) dir.create('sle/cytoscape')
write.csv(sle_deg,file = "sle/cytoscape/DEG.csv",
          row.names = F,quote = F)

load("data/ipfres.RData")
MMPs <- str_subset(res_df3$GeneName,"^MMP")
TIMPs <- str_subset(res_df3$GeneName,"^TIMP")
ipf_genes<- as.data.frame(c(MMPs,TIMPs))
colnames(ipf_genes) <- c("GeneName")
write.table(ipf_genes,file = "ipf/genelist.txt",
            row.names = F,col.names = F,quote = F)

expr_cnt <- left_join(ipf_genes,res_df3)
expr <- expr_cnt[,c(1:132)]
ipf_deg <- expr_cnt[,c(1,134,138)]
if (! dir.exists('ipf/cytoscape')) dir.create('ipf/cytoscape')
write.csv(ipf_deg,file = "ipf/cytoscape/DEG.csv",
          row.names = F,quote = F)


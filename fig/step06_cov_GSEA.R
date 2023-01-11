###CoV GSEA
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(readr)
library(R.utils)
R.utils::setOption("clusterProfiler.download.method",'auto')
cov_deg <- read.csv("cov/DE.csv",header = T,check.names = F)
genelist <- cov_deg[,c("GeneName","log2FoldChange")]
gene <- bitr(genelist$GeneName,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
colnames(gene) <- c("GeneName","ENTREZID")
genelist <- left_join(gene,genelist) %>% arrange(desc(log2FoldChange))
geneList <- genelist$log2FoldChange
names(geneList) <- genelist$ENTREZID
kk <- gseKEGG(geneList     = geneList,
              organism     = 'hsa',
              minGSSize    = 120,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
head(kk)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
gse <- as.data.frame(kk@result) 
write.csv(gse,"cov/GSEA_KEGG.csv")
p1 <- gseaplot2(kk, geneSetID = 24, title = kk$Description[24],
                color = "green",ES_geom = "line")
p1
ggsave(paste0("cov/Neutrophil extracellular trap formation", "_GSEA", ".pdf"), p1, width = 5.2, height = 4.6)

BP <- gseGO(geneList     = geneList,
      OrgDb        = org.Hs.eg.db,
      ont          = "BP",
      minGSSize    = 100,
      maxGSSize    = 500,
      pvalueCutoff = 0.05,
      verbose      = F)
BP <- setReadable(BP, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
BP_gse <- as.data.frame(BP@result)
write.csv(BP_gse,"cov/GSEA_BP.csv")
p3 <- gseaplot2(BP, geneSetID = c(2), title = BP$Description[2],
                color = "green",ES_geom = "line")
p3
ggsave(paste0("cov/regulation_of _cell_growth", "_GSEA", ".pdf"), p3, width = 5.2, height = 4.3)

CC <- gseGO(geneList = geneList,
            OrgDb = org.Hs.eg.db,
            ont = "CC",
            minGSSize = 100,
            maxGSSize = 500,
            pvalueCutoff = 0.05,
            verbose = F)
CC <- setReadable(CC, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
CC_gse <- as.data.frame(CC@result)
write.csv(CC_gse,"cov/GSEA_CC.csv")

MF <- gseGO(geneList = geneList,
            OrgDb = org.Hs.eg.db,
            ont = "MF",
            minGSSize = 100,
            maxGSSize = 500,
            pvalueCutoff = 0.05,
            verbose = F)
MF <- setReadable(MF, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
MF_gse <- as.data.frame(MF@result)
write.csv(MF_gse,"cov/GSEA_MF.csv")

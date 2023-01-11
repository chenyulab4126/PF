###ipf GSEA分析
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(readr)
ipf_deg <- read.csv("ipf/DE.csv",header = T,check.names = F)
genelist <- ipf_deg[,c("GeneName","log2FoldChange")]
gene <- bitr(genelist$GeneName,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
colnames(gene) <- c("GeneName","ENTREZID")
genelist <- left_join(gene,genelist) %>% arrange(desc(log2FoldChange))
geneList <- genelist$log2FoldChange
names(geneList) <- genelist$ENTREZID
kk <- gseKEGG(geneList     = geneList,
              organism     = 'hsa',
              minGSSize    = 10,
              pvalueCutoff = 1,
              verbose      = FALSE)
head(kk)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
gse <- as.data.frame(kk@result) 
write.csv(gse,"ipf/GSEA_KEGG.csv")
p1 <- gseaplot2(kk, geneSetID = 2, title = kk$Description[2],
                color = "green",ES_geom = "line")
p1
ggsave(paste0("ipf/NET", "_GSEA", ".pdf"), p1, width = 5.2, height = 4.6)

BP <- gseGO(geneList     = geneList,
            OrgDb        = org.Hs.eg.db,
            ont          = "BP",
            minGSSize    = 100,
            maxGSSize    = 500,
            pvalueCutoff = 0.05,
            verbose      = F)
BP <- setReadable(BP, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
BP_gse <- as.data.frame(BP@result)
write.csv(BP_gse,"ipf/GSEA_BP.csv")
p2 <- gseaplot2(BP, geneSetID = 12, title = BP$Description[12],
                color = "green",ES_geom = "line")
p2
ggsave(paste0("ipf/ECM", "_GSEA", ".pdf"), p2, width = 5.2, height = 4.3)
CC <- gseGO(geneList = geneList,
            OrgDb = org.Hs.eg.db,
            ont = "CC",
            minGSSize = 100,
            maxGSSize = 500,
            pvalueCutoff = 0.05,
            verbose = F)
CC <- setReadable(CC, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
CC_gse <- as.data.frame(CC@result)
write.csv(CC_gse,"ipf/GSEA_CC.csv")


MF <- gseGO(geneList = geneList,
            OrgDb = org.Hs.eg.db,
            ont = "MF",
            minGSSize = 100,
            maxGSSize = 500,
            pvalueCutoff = 0.05,
            verbose = F)
MF <- setReadable(MF, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
MF_gse <- as.data.frame(MF@result)
write.csv(MF_gse,"ipf/GSEA_MF.csv")
p3 <- gseaplot2(MF, geneSetID = 1, title = MF$Description[1],
                color = "green",ES_geom = "line")
p3
ggsave(paste0("ipf/ECM_st", "_GSEA", ".pdf"), p3, width = 5.2, height = 4.3)








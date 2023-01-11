###SLE GSEA
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(readr)
sle_deg <- read.csv("sle/DE.csv",header = T,check.names = F)
genelist <- sle_deg[,c("GeneName","log2FoldChange")]
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
write.csv(gse,"sle/GSEA_KEGG.csv")
p1 <- gseaplot2(kk, geneSetID = 7, title = kk$Description[7],
          color = "green",ES_geom = "line")
p1
ggsave(paste0("sle/NET", "_GSEA", ".pdf"), p1, width = 5.2, height = 4.6)

BP <- gseGO(geneList     = geneList,
            OrgDb        = org.Hs.eg.db,
            ont          = "BP",
            minGSSize    = 10,
            maxGSSize    = 500,
            pvalueCutoff = 0.05,
            verbose      = F)
BP <- setReadable(BP, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
BP_gse <- as.data.frame(BP@result)
write.csv(BP_gse,"sle/GSEA_BP.csv")

CC <- gseGO(geneList = geneList,
            OrgDb = org.Hs.eg.db,
            ont = "CC",
            minGSSize = 10,
            maxGSSize = 500,
            pvalueCutoff = 0.5,
            verbose = F)
CC <- setReadable(CC, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
CC_gse <- as.data.frame(CC@result)
write.csv(CC_gse,"sle/GSEA_CC.csv")


MF <- gseGO(geneList = geneList,
            OrgDb = org.Hs.eg.db,
            ont = "MF",
            minGSSize = 10,
            maxGSSize = 500,
            pvalueCutoff = 0.5,
            verbose = F)
MF <- setReadable(MF, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
MF_gse <- as.data.frame(MF@result)
write.csv(MF_gse,"sle/GSEA_MF.csv")
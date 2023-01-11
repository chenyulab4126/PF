#SLE and COV comparison analysis
library(readr)
library(dplyr)
library(ggplot2)
library(DelayedArray)
library(VennDiagram)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(readr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tinyarray)
library(venn)        
library(VennDiagram) 

SLE_DE <- read.csv(file.path("sle", "DE.csv"),check.names = F)
COV_DE <- read.csv(file.path("cov", "DE.csv"),check.names = F)
a <- unique(COV_DE$GeneName)
f_out <- "degcmp"
if (!dir.exists(f_out)) dir.create(f_out)
#common up
Up_gene <- list(
  "CoV Up" = COV_DE$GeneName[COV_DE$Tag == "Up"],
  "SLE Up" = SLE_DE$GeneName[SLE_DE$Tag == "Up"]
)

venn( Up_gene,
     zcolor='style', 
     opacity = 0.3,  
     box = F,        
     ilcs =1.5,     
     sncs = 1.5     
     )
down_gene <- list(
  "CoV Down" = COV_DE$GeneName[COV_DE$Tag == "Down"],
  "SLE Down" = SLE_DE$GeneName[SLE_DE$Tag == "Down"]
)
venn( down_gene,
      zcolor='#008844, #4B0091', 
      opacity = 0.3, 
      box = F,        
      ilcs =1.5,  
      sncs = 1.5       
)

DE_gene <- list(
  "SLE Up" = SLE_DE$GeneName[SLE_DE$Tag == "Up"],
  "SLE Down" = SLE_DE$GeneName[SLE_DE$Tag == "Down"],
  "CoV Up" = COV_DE$GeneName[COV_DE$Tag == "Up"],
  "CoV Down" = COV_DE$GeneName[COV_DE$Tag == "Down"]
)

# Plot the venn diagram
venn.diagram(
  DE_gene,
  file.path(f_out, "DEG.intersect.venn.tiff"),
  height = 450, width = 450, resolution = 150,
  fill = brewer.pal(4, "Set1"),
  cat.cex = 0.55
)
multi_size_df <- data.frame(GeneName=sort(unique(c(SLE_DE$GeneName, COV_DE$GeneName))))
multi_size_df$SLEUp <- multi_size_df$GeneName %in% DE_gene[["SLE Up"]]
multi_size_df$SLEDown <- multi_size_df$GeneName %in% DE_gene[["SLE Down"]]
multi_size_df$COVUp <- multi_size_df$GeneName %in% DE_gene[["CoV Up"]]
multi_size_df$COVDown <- multi_size_df$GeneName %in% DE_gene[["CoV Down"]]
multi_size_df <- multi_size_df[rowSums(multi_size_df[, 2:ncol(multi_size_df)]) > 1, ]

write_tsv(multi_size_df, file.path(f_out, "DEG.intersect.SourceData.tsv"))
coup <- filter(multi_size_df,SLEUp & COVUp)
codown <- filter(multi_size_df,SLEDown & COVDown)
codeg <- rbind(coup,codown)
codeg <- left_join(codeg,COV_DE)
codeg <- codeg[,c("GeneName","GeneID","CPM-Ctrl1","CPM-Ctrl2","CPM-COV1","CPM-COV2","Tag")]
codeg <- left_join(codeg,SLE_DE)
codeg <- codeg[,c("GeneName","GeneID","CPM-Ctrl1","CPM-Ctrl2","CPM-COV1","CPM-COV2","CPM-SLE1","CPM-SLE2","Tag")]
# PPI
genelist <- codeg$GeneName
write.table(genelist,file = "degcmp/genelist.txt",
            row.names = F,col.names = F,quote = F)
deg <- codeg[,c(1,9)]
write.csv(deg,file = "degcmp/DEG.csv",
          row.names = F,quote = F)
# heatmap
co_expr_cnt <- codeg[,c(3,4,5,6,7,8)]
rownames(co_expr_cnt) <- codeg$GeneName
co_expr_matrix <- t(scale(t(as.matrix(co_expr_cnt))))
colnames(co_expr_matrix) <- c("CTRL1","CTRL2","CoV1","CoV2","SLE1","SLE2")
co_expr_matrix <- t(co_expr_matrix)
p <- Heatmap(co_expr_matrix,
             name = "Scaled Expr.",
             cluster_columns = FALSE,
             col = colorRamp2(c(-2,0,2),c("blue","white","red")),
             row_names_gp = gpar(fontsize = 6),
             column_names_gp = gpar(fontsize = 3),
             column_names_centered = T,
             row_dend_width = unit(5, "mm"),
             column_dend_height = unit(5, "mm"),
             column_names_rot = 90,
             width = 2,
             heatmap_legend_param = list(
               labels_gp = gpar(fontsize = 4),
               title_gp = gpar(fontsize = 4),
               grid_width = unit(1, "mm"),
               grid_height = unit(1, "mm")))
p
inche_cm <- 2.54
pdf(file.path(f_out,"degcompare.pdf"),
    width= (3 + 0.1*ncol(co_expr_matrix))/inche_cm, 
    height = (3 + 0.15*nrow(co_expr_matrix))/inche_cm)
print(p)
dev.off()


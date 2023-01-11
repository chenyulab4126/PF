#SLE transplotion sample cor analysis
library(ComplexHeatmap)
library(circlize)
library(readr)
library(ggplot2)
library(dplyr)
library(edgeR)
library(gridExtra)
library(FactoMineR)
library(ggrepel)
library(factoextra)

SLE_tpl_count <- read.table("data/SLE lung transplation/sle_counts.txt",
                            header = T,check.names = F)
f_out <- "sle"
if (!dir.exists("sle")) dir.create("sle")
gene_info <- read.table("data/hg38_gencode.v32.info.tsv",header = T)
colnames(gene_info) <- c("GeneID","GeneName")
sle_counts <- left_join(gene_info,SLE_tpl_count)
colnames(sle_counts) <- c("Geneid","GeneName","CTRL1","CTRL2","SLE1","SLE2")
save(sle_counts,file = "data/step01sle.Rdata")
save(gene_info,file = "data/geneinfo.Rdata")

cnt <- sle_counts[,c(3:ncol(sle_counts))]
row.names(cnt) <- sle_counts$Geneid
cnt_mat <- as.matrix(cnt)
expr_cor <- cor(cnt_mat, method = "spearman")
expr_cor

write.csv(format(expr_cor,digits = 4),
          file.path(f_out,"expr_cor.csv"),
          row.names = T,
          quote = F)

tmp_df1 <- data.frame(expr1 = cnt_mat[, 3], expr2 = cnt_mat[, 4])
p1 <- ggplot(tmp_df1, aes(x = log10(expr1), y = log10(expr2))) +
  geom_point(size = 0.1, color = "black") +
  annotate(
    "text",
    x = 0.85 * log10(max(tmp_df1$expr1)),
    y = 0.15 * log10(max(tmp_df1$expr2)),
    label = sprintf("Cor = %.2f", expr_cor[3, 4]),
    size = 2
  ) +
  theme_bw() +
  labs(
    x = sprintf("Gene read counts in %s (log10)", colnames(cnt_mat)[3]),
    y = sprintf("Gene read counts in %s (log10)", colnames(cnt_mat)[4])
  ) +
  theme(
    axis.text = element_text(family = "ArialMT", color = "black", size = 7),
    axis.title = element_text(family = "ArialMT", color = "black", size = 7),
    panel.grid = element_blank(),
  )    

plot_size <- ncol(cnt_mat) - 1
inche_cm <- 2.54
pdf(
  file.path(f_out,"sle_Cor.pdf"),
  width = 5.8*plot_size/inche_cm,
  height = 5.8*plot_size/inche_cm,
  colormodel = "cmyk",
  family = "ArialMT"
)
p1
dev.off()


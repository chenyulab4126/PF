#COVID-19 transplotion sample cor analysis
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

COV_tpl_count <- read.table("data/COVID-19 lung trasplation/cov_counts.tsv",
                            header = T,check.names = F) 

gene_info <- read.table("data/hg38_gencode.v32.info.tsv",header = T)
colnames(gene_info) <- c("Geneid","GeneName")
cov_counts <- left_join(gene_info,COV_tpl_count)
colnames(cov_counts) <- c("Geneid","GeneName","CTRL1","CTRL2","COV1","COV2")
save(cov_counts,file="data/step01cov.Rdata")

cnt <- cov_counts[,c(3:ncol(cov_counts))]
row.names(cnt) <- cov_counts$Geneid
cnt_mat <- as.matrix(cnt)
expr_cor <- cor(cnt_mat, method = "spearman")
expr_cor
f_out <- "cov"
if (!dir.exists("cov")) dir.create("cov")
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
  file.path(f_out,"cov_Cor.pdf"),
  width = 4.8*plot_size/inche_cm,
  height = 4.8*plot_size/inche_cm,
  colormodel = "cmyk",
  family = "ArialMT"
)
p1
dev.off()
tmp_df2 <- data.frame(expr1 = cnt_mat[, 1], expr2 = cnt_mat[, 2])
p2 <- ggplot(tmp_df2, aes(x = log10(expr1), y = log10(expr2))) +
  geom_point(size = 0.1, color = "black") +
  annotate(
    "text",
    x = 0.85 * log10(max(tmp_df2$expr1)),
    y = 0.15 * log10(max(tmp_df2$expr2)),
    label = sprintf("Cor = %.2f", expr_cor[1, 2]),
    size = 2
  ) +
  theme_bw() +
  labs(
    x = sprintf("Gene read counts in %s (log10)", colnames(cnt_mat)[1]),
    y = sprintf("Gene read counts in %s (log10)", colnames(cnt_mat)[2])
  ) +
  theme(
    axis.text = element_text(family = "ArialMT", color = "black", size = 7),
    axis.title = element_text(family = "ArialMT", color = "black", size = 7),
    panel.grid = element_blank(),
  )    

plot_size <- ncol(cnt_mat) - 1
inche_cm <- 2.54
pdf(
  file.path(f_out,"control_Cor.pdf"),
  width = 4.8*plot_size/inche_cm,
  height = 4.8*plot_size/inche_cm,
  colormodel = "cmyk",
  family = "ArialMT"
)
p2
dev.off()


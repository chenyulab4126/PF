#COVID-19 差异基因分析
library(DESeq2)
library(dplyr)
library(ggplot2)
library(readr)
fc_cutoff <- 8
padj_cutoff <- 0.05
f_out <- "cov"
load("data/step01cov.Rdata")
all_cnt_mat <- as.matrix(cov_counts[, c(
  "CTRL1",
  "CTRL2",
  "COV1",
  "COV2"
)])
rownames(all_cnt_mat) <- cov_counts$Geneid
condition <- factor(c(rep("CTRL", 2), rep("COV", 2)), levels = c("CTRL", "COV"))
dds <- DESeqDataSetFromMatrix(all_cnt_mat, DataFrame(condition), ~ condition)
dds <- DESeq(dds)
res <- results(dds)
res <- as.data.frame(res)
res$GeneID <- rownames(res)
expr_res <- na.omit(res)
expr_res$Tag <- "NC"
expr_res$Tag[expr_res$log2FoldChange > log2(fc_cutoff) & expr_res$padj < padj_cutoff] <- "Up"
expr_res$Tag[expr_res$log2FoldChange < (-1 * log2(fc_cutoff)) & expr_res$padj < padj_cutoff] <- "Down"
expr_res$Tag <- factor(expr_res$Tag, levels = c("Up", "NC", "Down"))
res_df <- expr_res[, c("GeneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "Tag")]
#表达矩阵CPM--按列标化
expr_mat <- as.matrix(cov_counts[,3:ncol(cov_counts)])
row.names(expr_mat) <- cov_counts$Geneid
norm_mat <- 1e6 * t(t(expr_mat) / colSums(expr_mat))
norm_expr_df <- as.data.frame(norm_mat)
norm_expr_df$GeneID <- rownames(norm_mat)
gene_info <- read.table("data/hg38_gencode.v32.info.tsv",header = T)
colnames(gene_info) <- c("GeneID", "GeneName")
norm_expr_df <- left_join(gene_info, norm_expr_df)
names(norm_expr_df) <- c("GeneID", "GeneName", "CPM-Ctrl1", "CPM-Ctrl2", "CPM-COV1", "CPM-COV2")
res_df <- right_join(norm_expr_df, res_df)
res_df <- res_df[!duplicated(res_df$GeneName),]

write.csv(res_df,file.path(f_out,"DE.csv"),quote = F,row.names = F)
save(res_df,file = "data/step02_covres.RData")
###火山图
deseq2_info <- res_df %>% group_by(Tag) %>% summarise(Num = n())
deseq2_info$Text <- sprintf("N=%d", deseq2_info$Num)
deseq2_info <- deseq2_info[which(deseq2_info$Tag != "NC"), ]
deseq2_info
deseq2_info$x <- -9 # for label text position
deseq2_info$x[deseq2_info$Tag == "Up"] <- -deseq2_info$x[deseq2_info$Tag == "Up"]
dat  = res_df
col_li <- c("#53868B","#e9e9e9","orange")
names(col_li) <- c("Down","NC","Up")
p1 <- ggplot(dat, aes(
  x = log2FoldChange,
  y = -log10(padj)))  +
  geom_point(size = 1.2,shape=16,aes(color=Tag)) +
  labs(x = "log2 (Fold Change)", y = "Adjusted P-value",title = "COVID-19 differential expression genes") +
  scale_y_continuous(
    breaks = c(0, 10,20,30,40),
    labels = c(
      expression("10" ^ "-0"),
      expression("10" ^ "-10"),
      expression("10" ^ "-20"),
      expression("10" ^ "-30"),
      expression("10" ^ "-40")
    )
  ) +
  scale_x_continuous(
    breaks = c(-8,-4,0,4,8),
    labels = c(
      expression("-8"),
      expression("-4"),
      expression("0"),
      expression("4"),
      expression("8")
    )
  ) +
  scale_color_manual(values = col_li) +
  coord_cartesian(xlim = c(-8, 8), ylim = c(0, 40)) +
  geom_vline(xintercept=c(-3,3),lty=4,col="#787878",lwd=0.3) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="#787878",lwd=0.3) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
    text = element_text(
      family = "ArialMT",
      color = "black",
      size = 12
    ),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.key.size = unit(12, "mm"),
    panel.grid = element_blank()
  )
p1
ggsave(
  file.path(f_out,"volcano_plot.pdf"),
  p1,
  width = 7,
  height =5.5
)
dev.off()
Immune <- read.delim("data/immmune_genelist_unique.txt",
                      sep = "\t",header = T)

Immunelist <- Immune$Symbol

if(T){
  for_label <- dat%>% 
    filter(!Tag == "NC") %>%
    filter(GeneName %in% Immunelist) 
  for_label$group <- "Immune genes" 
}

volcano_plot <- p1 +
  geom_point(size = 1,shape=16, data = for_label,color="red") +
  ggrepel::geom_text_repel(aes(label = GeneName),
                           data = for_label,
                           color="black",
                           size = 3,
                           segment.color = "#cccccc", 
                           segment.size = 0.4,
                          ) 
volcano_plot 
ggsave(
  file.path(f_out,"volcano_plot_immune.pdf"),
  volcano_plot,
  width = 7,
  height =5.5
)
dev.off()


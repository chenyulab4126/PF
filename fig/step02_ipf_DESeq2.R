#ipf差异分析
library(DESeq2)
library(dplyr)
library(readr)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)

fc_cutoff <- 2
padj_cutoff <- 0.05
f_out <- "ipf"
if (!dir.exists("ipf")) dir.create("ipf")

cnt_df <- read_delim("data/IPF/IPF_counts.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
summary(cnt_df)
colnames(cnt_df) <- c("Name",paste0(rep("IPF", 67),1:67),
                      paste0(rep("CTRL", 63),1:63))
cnt <- cnt_df[, 2:ncol(cnt_df)]
rownames(cnt) <- cnt_df$Name
cnt_mat <- as.matrix(cnt)
expr_mat <- cnt_mat
expr_cor <- cor(expr_mat, method = "spearman")

expr_cor
write.table(
  format(expr_cor, digits = 4),
  file.path(f_out, "Cor.SourceData.csv"),
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  sep = "\t")

### Correlation clustering
ht <- Heatmap(
  expr_cor,
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  name = "Expr. Cor.",
  col = colorRamp2(c(0.8, 1), c("white", "#D7191C")),
  row_names_gp = gpar(fontsize = 1),
  column_names_gp = gpar(fontsize = 1),
  row_dend_width = unit(5, "mm"),
  column_dend_height = unit(5, "mm"),
  show_row_dend = FALSE,
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 6),
    title_gp = gpar(fontsize = 6),
    grid_width = unit(2, "mm"),
    grid_height = unit(2, "mm")
  )
)
inche_cm <- 2.54
pdf(file.path(f_out, "Cor.ht.pdf"), width = 7/inche_cm, height = 6/inche_cm)
print(ht)
dev.off()
dim(cnt_mat)
## DEG analysis
condition <- factor(c(rep("IPF", 67), rep("CTRL", 63)), levels = c("CTRL", "IPF"))
dds <- DESeqDataSetFromMatrix(cnt_mat, DataFrame(condition), ~ condition)

dds <- DESeq(dds)
res <- results(dds)
res <- as.data.frame(res)
res$GeneName <- rownames(res)
expr_res <- na.omit(res)
expr_res$Tag <- "NC"
expr_res$Tag[expr_res$log2FoldChange > log2(fc_cutoff) & expr_res$padj < padj_cutoff] <- "Up"
expr_res$Tag[expr_res$log2FoldChange < (-1 * log2(fc_cutoff)) & expr_res$padj < padj_cutoff] <- "Down"
expr_res$Tag <- factor(expr_res$Tag, levels = c("Up", "NC", "Down"))
res_df <- expr_res[, c("GeneName","baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj","Tag")]
dim(res_df)
dim(expr_mat)
norm_mat <- 1e6 * t(t(expr_mat) / colSums(expr_mat))
norm_expr_df <- as.data.frame(norm_mat)
norm_expr_df$GeneName<- rownames(norm_mat)

norm_expr_df <- dplyr::select(norm_expr_df,"GeneName",everything())
dim(norm_expr_df)
names(norm_expr_df) <- c("GeneName",paste0(rep("CPM-IPF", 67),1:67) , paste0(rep("CPM-CTRL", 63),1:63))
res_df <- right_join(norm_expr_df, res_df)
head(res_df)
gene_info <- read.table("data/hg38_gencode.v32.info.tsv",header = T)
colnames(gene_info) <- c("GeneID", "GeneName")
res_df1 <- left_join(res_df,gene_info)
res_df2 <- na.omit(res_df1)
res_df3 <- dplyr::select(res_df2,"GeneID",everything())
write.csv(res_df3,file.path(f_out,"DE.csv"),quote = F,row.names = F)
save(res_df3,file = "data/ipfres.RData")
###火山图
deseq2_info <- res_df3 %>% group_by(Tag) %>% summarise(Num = n())
deseq2_info$Text <- sprintf("N=%d", deseq2_info$Num)
deseq2_info <- deseq2_info[which(deseq2_info$Tag != "NC"), ]
deseq2_info
deseq2_info$x <- -9 # for label text position
deseq2_info$x[deseq2_info$Tag == "Up"] <- -deseq2_info$x[deseq2_info$Tag == "Up"]
dat  = res_df3
col_li <- c("#53868B","#e9e9e9","orange")
names(col_li) <- c("Down","NC","Up")
p1 <- ggplot(dat, aes(
  x = log2FoldChange,
  y = -log10(padj)))  +
  geom_point(size = 1.2,shape=16,aes(color=Tag)) +
  labs(x = "log2 (Fold Change)", y = "Adjusted P-value",title = "IPF differential expression genes") +
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
    breaks = c(-4,-2,0,2,4),
    labels = c(
      expression("-4"),
      expression("-2"),
      expression("0"),
      expression("2"),
      expression("4")
    )
  ) +
  scale_color_manual(values = col_li) +
  coord_cartesian(xlim = c(-4,4), ylim = c(0, 40)) +
  geom_vline(xintercept=c(-1,1),lty=4,col="#787878",lwd=0.3) +
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

if(F) {
  #p鍊兼渶灏忕殑鍓?3涓嬭皟鍜屽墠3涓婅皟
  x1 = dat %>% 
    filter(Tag == "Up") %>% 
    arrange(padj) %>%
    head(6)
  x2 = dat %>% 
    filter(Tag == "Down") %>% 
    arrange(padj) %>%
    head(6)
  for_label = rbind(x1,x2)
}


volcano_plot <- p1 +
  geom_point(size = 1,shape=16, data = for_label,color="red") +
  ggrepel::geom_text_repel(aes(label = GeneName),
                           data = for_label,
                           color="black",
                           size = 3,
                           segment.color = "#cccccc", 
                           segment.size = 0.4) 
volcano_plot 
ggsave(
  file.path(f_out,"volcano_plot_immune.pdf"),
  volcano_plot,
  width = 7,
  height =5.5
)
dev.off()

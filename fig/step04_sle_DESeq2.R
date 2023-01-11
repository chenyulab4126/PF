library(DESeq2)
library(dplyr)
library(ggplot2)
library(readr)

fc_cutoff <- 8
padj_cutoff <- 0.05
f_out <- "sle"
if (!dir.exists("sle")) dir.create("sle")
load("data/step01sle.Rdata")
all_cnt_mat <- as.matrix(sle_counts[, c(
  "CTRL1",
  "CTRL2",
  "SLE1",
  "SLE2"
)])
rownames(all_cnt_mat) <- sle_counts$Geneid
condition <- factor(c(rep("CTRL", 2), rep("SLE", 2)), levels = c("CTRL", "SLE"))
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

expr_mat <- as.matrix(sle_counts[,3:ncol(sle_counts)])
row.names(expr_mat) <- sle_counts$Geneid
norm_mat <- 1e6 * t(t(expr_mat) / colSums(expr_mat))
norm_expr_df <- as.data.frame(norm_mat)
norm_expr_df$GeneID <- rownames(norm_mat)
load("data/geneinfo.Rdata")
norm_expr_df <- left_join(gene_info, norm_expr_df)
names(norm_expr_df) <- c("GeneID", "GeneName", "CPM-Ctrl1", "CPM-Ctrl2", "CPM-SLE1", "CPM-SLE2")
res_df <- right_join(norm_expr_df, res_df)
res_df <- res_df[!duplicated(res_df$GeneName),]

write.csv(res_df,file.path(f_out,"DE.csv"),quote = F,row.names = F)
save(res_df,file = "data/sleres.RData")

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
  labs(x = "log2 (Fold Change)", y = "Adjusted P-value",title = "SLE differential expression genes") +
  scale_y_continuous(
    breaks = c(0, 1,2,3,4),
    labels = c(
      expression("10" ^ "-0"),
      expression("10" ^ "-1"),
      expression("10" ^ "-2"),
      expression("10" ^ "-3"),
      expression("10" ^ "-4")
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
  coord_cartesian(xlim = c(-8,8), ylim = c(0, 4)) +
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

if(F) {
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
                           size = 2,
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

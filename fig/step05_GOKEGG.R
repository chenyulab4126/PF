#GO-term and KEGG enrichment analysis of DEGs
library(VennDiagram)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)

#' GO enrichment analysis by ontology class
goByOntology <- function(genelist, outprefix, ont="BP", pvalcut=0.05, qvalcut=0.2, backgroundgene=NULL) {
  stopifnot(ont %in% c("BP", "MF", "CC"))
  go <- enrichGO(
    genelist,
    keyType = 'ENTREZID', # IDtype
    OrgDb = org.Hs.eg.db,
    ont = ont,
    pAdjustMethod = 'BH',
    pvalueCutoff = pvalcut,
    qvalueCutoff = qvalcut,
    universe = backgroundgene
  )
  go@result <- go@result[order(go@result$pvalue, decreasing = F), ]
  
  fprefix <- paste0(outprefix, '_GO_', ont)
  write.csv(as.data.frame(go), paste0(fprefix, '.csv'), row.names = F)
  
  ggsave(
    paste0(fprefix, '_bar', '.pdf'),
    barplot(go, showCategory = 20, drop = T),
    width = 10,
    height = 10)
  ggsave(
    paste0(fprefix, '_bub', '.pdf'),
    dotplot(go, showCategory = 10),
    width = 10,
    height = 10
  )
  return(as.data.frame(go))
}

# Do KEGG enrichment analysis
keggEnrich <- function(genelist, outprefix, pvalcut=0.05, qvalcut=0.2) {
  kegg <- enrichKEGG(
    genelist,
    organism = 'hsa',
    keyType = 'kegg',
    pvalueCutoff = pvalcut,
    pAdjustMethod = 'BH',
    minGSSize = 10,
    maxGSSize = length(genelist),
    qvalueCutoff = qvalcut,
    use_internal_data = FALSE
  )
  kegg@result = kegg@result[order(kegg@result$pvalue, decreasing = F), ]
  
  fprefix <- paste0(outprefix, '_KEGG')
  write.csv(as.data.frame(kegg), paste0(fprefix, '.csv'), row.names = F)
  
  ggsave(
    paste0(fprefix, '_bar', '.pdf'),
    barplot(kegg, showCategory = 20, drop = T),
    width = 10,
    height = 10)
  ggsave(
    paste0(fprefix, '_bub', '.pdf'),
    dotplot(kegg, showCategory = 30),
    width = 10,
    height = 10
  )
  return(as.data.frame(kegg))
}

goKeggEnrich <- function(proj="cov", group="Down") {
  f_out <- file.path(proj, "TERM")
  if (!dir.exists(f_out)) dir.create(f_out)
  outprefix <- file.path(f_out, group)
  
  f_DE <- file.path(proj, "DE.csv")
  data <- read.csv(f_DE, header = TRUE,check.names = F,stringsAsFactors = FALSE)
  head(data)
  if (group == "Deg") {
    gdf <- data[data$Tag == "Up" | data$Tag == "Down", c("GeneID", "GeneName")]
  } else {
    gdf <- data[data$Tag == group, c("GeneID", "GeneName")]
  }
  
  genes <- bitr(gdf$GeneName, fromType = "SYMBOL", toType = c("ENTREZID", "ENSEMBL"), OrgDb = "org.Hs.eg.db")
  dfcomb <- inner_join(gdf, genes, by = c("GeneName" = "SYMBOL"))
  dfcomb
  genelist <- unique(genes$ENTREZID)
  
  # Use expressed genes as background
  bkgenes <- bitr(data$GeneName, fromType = "SYMBOL", toType = c("ENTREZID", "ENSEMBL"), OrgDb = "org.Hs.eg.db")
  bkgenelist <- unique(bkgenes$ENTREZID)
  
  # Do GO-term enrichment analysis
  BP_df <- goByOntology(genelist, outprefix, ont="BP", backgroundgene=bkgenelist)
  MF_df <- goByOntology(genelist, outprefix, ont="MF")
  CC_df <- goByOntology(genelist, outprefix, ont="CC")
  KEGG_df <- keggEnrich(genelist, outprefix)
  term_all <- rbind(BP_df, MF_df, CC_df, KEGG_df)
  
  # Choose top 10 terms by category to present
  topterms <- function(df, ont = "BP", n = 10) {
    dft <- cbind(df, ont)
    dft <- dft[1:min(n, nrow(dft)), c(2, 5, ncol(dft))]
    colnames(dft) <- c("Description", "pvalue", "Category")
    return(dft)
  }
  BP_dftop <- topterms(BP_df, "BP")
  MF_dftop <- topterms(MF_df, "MF")
  CC_dftop <- topterms(CC_df, "CC")
  KEGG_dftop <- topterms(KEGG_df, "KEGG")
  term_top <- rbind(BP_dftop, MF_dftop, CC_dftop, KEGG_dftop)
  term_top$Description <- factor(term_top$Description, levels = rev(unique(term_top$Description)))
  term_top
  
  # Generate figure for top terms by category
  ymax <- max(-log10(term_top$pvalue))
  ybreaks <- seq(0, ymax, floor(ymax/4.))
  p <- ggplot(term_top, aes(
    x = Description,
    fill = Category,
    y = -log10(pvalue)
  )) +
    theme_bw() +
    geom_bar(stat = "identity", width = 0.6) +
    geom_text(aes(label = Description), size = 2, hjust = 0) +
    scale_fill_brewer(palette = "Dark2") +
    labs(y = "-log10(P-value)") +
    coord_flip() +
    scale_y_continuous(breaks = ybreaks,
                       limits = c(0, ymax*2),
                       labels = sprintf("%.0E", 10 ** ybreaks)) +
    theme(
      axis.text = element_text(
        family = "ArialMT",
        color = "black",
        size = 6
      ),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(),
    )
  ggsave(paste0(outprefix, "_GOKEGG", ".pdf"), p, width = 10, height = 10, units = "cm")
  write.csv(term_top, file=paste0(outprefix, "_GOKEGG_top", ".csv"), row.names = F)
  write.csv(term_all, file=paste0(outprefix, "_GOKEGG_all", ".csv"), row.names = F)
}


for (proj in c("cov", "sle","ipf")) {
  if (!dir.exists(proj)) dir.create(proj)
  for (group in c("Up", "Down", "Deg")) {
    cat(proj, group, "\n")
    goKeggEnrich(proj, group)
  }
}



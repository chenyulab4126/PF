library(ggplot2)
library(dplyr)
library(VennDiagram)
term_top <- read.delim("data/GOKEGGterm/cov_selectGOKEGGterm.txt",header = T)
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
  geom_bar(stat = "identity", width = 0.7) +
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
p
ggsave(paste0("cov/", "select_GOKEGG", ".pdf"), p, width = 10, height = 5.5, units = "cm")


term_top <- read.delim("data/GOKEGGterm/sle_selectGOKEGGterm.txt",header = T)
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
  geom_bar(stat = "identity", width = 0.7) +
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
p
ggsave(paste0("sle/", "select_GOKEGG", ".pdf"), p, width = 10, height = 5.5, units = "cm")



term_top <- read.delim("data/GOKEGGterm/ipf_selectGOKEGGterm.txt",header = T)
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
  geom_bar(stat = "identity", width = 0.7) +
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
p
ggsave(paste0("ipf/", "select_GOKEGG", ".pdf"), p, width = 10, height =5.5, units = "cm")



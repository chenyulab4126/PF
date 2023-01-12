
## load packages
library(pheatmap)
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(ComplexHeatmap)
library(circlize)
library(Hmisc)

## import data
final <- read.table("Data-Figure5_heatmap.csv",header = T,sep = ",",row.names = 1) 
head(final)

table(final$Drug)
## create matrix for heatmap
hmap_bt <- as.matrix(final[,-c(1:2)])
hmap_bt <- as.matrix(t(scale(t(hmap_bt))))

# heatmap color set
col_fun = colorRamp2(c(max(hmap_bt), mean(as.matrix(hmap_bt)), min(hmap_bt)), c("red","white","blue"))

## rownames color setting
col <- c( rep("deeppink3",1), rep("forestgreen",2),rep("grey30",14),rep("royalblue1",2),rep("chocolate2",3),rep("cyan4",2))

## image
png("Figure5_heatmap.png", width = 2.3, height = 4.5, units = 'in', res = 600)
Heatmap <- Heatmap(hmap_bt, 
                   col = col_fun,
                   column_names_rot = 60,
                   cluster_rows = T, 
                   cluster_columns = FALSE,
                   row_names_gp = gpar(col = col, fontsize = 7,fontface = "bold"),
                   column_split = data.frame(c(rep("A_CTRL",2),rep("B_COV",2),rep("C_SLE",2))),
                   column_title_gp = gpar(col = c("white","white","white")),
                   column_names_gp = gpar(col = "black", fontsize = 8, fontface = "bold"),
                   heatmap_legend_param = list(title= "VST", legend_direction = "horizontal",legend_width = unit(2.8, "cm")))
draw(Heatmap, show_heatmap_legend = T, heatmap_legend_side="bottom")
dev.off()








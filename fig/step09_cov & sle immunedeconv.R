library(rtracklayer)
library(tidyr)
library(readr)
library(stringr)
library(dplyr)
library(immunedeconv)
library(tibble)
library(IOBR)
library(remotes)
library(ggplot2)
library(tidyverse)
#gtf
gtf = rtracklayer::import("data/Homo_sapiens.GRCh38.102.chr.gtf")
class(gtf)
gtf = as.data.frame(gtf);dim(gtf)
table(gtf$type)
save(gtf,file='data/gtf.Rdata')
load("data/gle.Rdata")
colnames(gle) <- c("GeneName","length")
load("data/geneinfo.Rdata")
effLen <- na.omit(left_join(gene_info,gle))
effLen <- effLen[!duplicated(effLen$GeneName),]
effLen$Lenkb <- effLen$length/1000
save(effLen,file = "data/effLen.Rdata")
tra = gtf[gtf$type=="transcript",
          c("start","end","gene_name")]
length(unique(tra$gene_name))
glt = mutate(tra, length = end - start) %>%
  arrange(desc(length)) %>% 
  filter(!duplicated(gene_name)) %>% 
  dplyr::select(c(3,4))
head(glt)

load("data/step01cov.Rdata")
load("data/step01sle.Rdata")
expr_cnt <- left_join(cov_counts,sle_counts)
expr_cnt <- expr_cnt[!duplicated(expr_cnt$GeneName),]
expr_matrix <- as.matrix(expr_cnt[,3:8])
rownames(expr_matrix) <- expr_cnt$GeneName

countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}
tpms <- countToTpm(counts = expr_matrix,effLen = effLen$Lenkb)
write.table(tpms,file="data/tpmkb.txt",sep="\t",quote=F)


#CIBERSORT

source('immunedeconv/CIBERSORT.R')

res_cibersort <- CIBERSORT('immunedeconv/LM22.txt','immunedeconv/DATA.txt', perm = 100, QN = F) 
write.table(res_cibersort, file="immunedeconv/res_cibersort.txt", sep="\t", col.names=T, row.names=F, quote=F)
save(res_cibersort, file="immunedeconv/res_cibersort.Rdata")
load("immunedeconv/res_cibersort.Rdata")

colour <- c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB",
            "#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3",
            "#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
            "#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B",
            "#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
            "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347",
            "#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

my_theme <- function(){
  theme(panel.grid = element_blank(),       
        panel.border = element_blank(),    
        legend.position="right",          
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        axis.line = element_line(size=1),   
        text = element_text(family="Times"),
        axis.text.y = element_text(size = 8,face='bold',color='black'),
        axis.text.x = element_text(size = 8,face='bold',color='black',angle=90,hjust=1),        # x轴标签样式，angle=45 倾斜 45 度
        axis.title = element_text(size=10,face="bold"),  
        plot.title = element_text(hjust=0.5,size=10))    
}  

p1 <- res_cibersort[,1:22] %>% reshape2::melt() %>%
  ggplot(aes(x=Var1,y=value,fill=Var2)) +
  geom_bar(stat='identity') +
  coord_flip()  +
  scale_fill_manual(values =colour ) +
  theme_bw()+ theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+my_theme()
p1
pdf("cibersort.pdf",width = 10,height = 4)
p1
dev.off()

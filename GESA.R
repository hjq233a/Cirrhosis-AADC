rm(list = ls())  
options(stringsAsFactors = F)
# library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(ggstatsplot)
library("GSEABase")
library("GSVA")

#geneset 两列，第一列为通路名，第二列为基因
gene_set <- read.csv("geneset.csv")  

#genelist,两列：基因名，变化倍数
need_DEG <-read.csv("FMT-P-geneset5.csv")
geneList <- need_DEG$log2FC
geneList <-as.numeric(geneList)
names(geneList) <- need_DEG$gene
geneList <- geneList[is.finite(geneList)]
gene_list_vector <- sort(geneList, decreasing = T) 

gsea_results <- GSEA(geneList = gene_list_vector,
                            minGSSize = 1, # minimum gene set size, can be adjusted
                            maxGSSize = 5000, # maximum gene set size, can be adjusted
                            pvalueCutoff = 0.9, # can be adjusted based on your preference
                            verbose = TRUE,
                            by = "fgsea",
                            TERM2GENE = gene_set)
#暂时没学会
#gsea_results1 <- gsva(gene_list_vector, gene_set, method = "gsea")

#GESA作图
gseap1 <- gseaplot2(gsea_results,
                    gsea_results$ID[3],#富集的ID编号
                    title = NULL,#标题
                    color = "black", #GSEA线条颜色
                    base_size = 30,#基础字体大小
                    rel_heights = c(1.5, 0.2),#副图的相对高度
                    subplots = c(1,2),   #要显示哪些副图 如subplots=c(1,3) 
                    ES_geom = "line", #enrichment score用线还是用点"dot",line
                    pvalue_table = F) 
gseap1


ggsave(gseap1, filename = 'GSEA_genest5.HEup.pdf', width =10, height =6,dpi = 400)


#批量做GO或KEGG
GO_kk_entrez <- gseGO(geneList     = geneList,
                      ont          = "BP",  # "BP"、"MF"和"CC"或"ALL"
                      OrgDb        = "org.Mm.eg.db",#人类org.Hs.eg.db 鼠org.Mm.eg.db
                      keyType      = "SYMBOL",
                      pvalueCutoff = 0.25)   #实际为padj阈值可调整

GO_kk <- DOSE::setReadable(GO_kk_entrez, 
                           OrgDb="org.Hs.eg.db",
                           keyType="SYMBOL")#转化id 


kk_gse <-gsea_results   #GESA结果：gsea_results，GO_kk
kk_gse_entrez <- GO_kk_entrez

kk_gse_cut <- kk_gse[kk_gse$pvalue<0.05 & kk_gse$p.adjust<0.25 & abs(kk_gse$NES)>1]
kk_gse_cut_down <- kk_gse_cut[kk_gse_cut$NES < 0,]
kk_gse_cut_up <- kk_gse_cut[kk_gse_cut$NES > 0,]
write.csv(kk_gse_cut_up,"kk_gse_cut_up.csv")

#选择展现NES前几个通路，这一步也可以选择自己感兴趣的通路 
down_gsea <- kk_gse_cut_down[tail(order(kk_gse_cut_down$NES,decreasing = T),200),]
up_gsea <- kk_gse_cut_up[head(order(kk_gse_cut_up$NES,decreasing = T),200),]
diff_gsea <- kk_gse_cut[head(order(abs(kk_gse_cut$NES),decreasing = T),200),]

#### 经典的GSEA图 
up_gsea$Description
i=2
gseap1 <- gseaplot2(kk_gse,
                    kk_gse$ID[1],#富集的ID编号
                    title = kk_gse$ID[1],#标题
                    color = "red", #GSEA线条颜色
                    base_size = 20,#基础字体大小
                    rel_heights = c(1.5, 0.5, 1),#副图的相对高度
                    subplots = 1:3,   #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
                    ES_geom = "line", #enrichment score用线还是用点"dot"
                    pvalue_table = T) #显示pvalue等信息
ggsave(gseap1, filename = 'GSEA_up_monoamine.self.pdf', width =10, height =8)


#### 合并 GSEA通路 
gseap2 <- gseaplot2(kk_gse,
                    up_gsea$ID,#富集的ID编号
                    title = "UP_GSEA_all",#标题
                    color = "red",#GSEA线条颜色
                    base_size = 20,#基础字体大小
                    rel_heights = c(1.5, 0.5, 1),#副图的相对高度
                    subplots = 1:3, #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
                    ES_geom = "line",#enrichment score用线还是用点"dot"
                    pvalue_table = T) #显示pvalue等信息
ggsave(gseap2, filename = "GSEA_up_all.pdf",width =12,height =12)


## gsearank plot 绘制出属于特定基因集的基因排序列表
##绘制up_gsea前3个富集通路
library(cowplot)
library(ggplot2)
pp <- lapply(1:3, function(i) {
  anno <- up_gsea[i, c("NES", "pvalue", "p.adjust")]
  lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")
  
  gsearank(kk_gse,
           up_gsea$ID[i], #可以变成想要的geneset吗？
           title = up_gsea$Description[i]) + 
    xlab(NULL) +
    # ylab(NULL) +
    annotate("text", 2500,
             up_gsea[i, "enrichmentScore"] * .75, 
             label = lab, 
             hjust=0, vjust=0)
})
rankp <- plot_grid(plotlist=pp, ncol=1)
ggsave(rankp, filename = "gsearank_up.pdf",width=8,height=10)

 

#
anno <-  kk_gse_cut_up[116, c("NES", "pvalue", "p.adjust")]
lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")
pp<-gsearank(kk_gse,
             kk_gse_cut_up$ID[116], #可以变成想要的geneset吗？
             title =  kk_gse_cut_up$Description[116]) + 
  xlab(NULL) +
  # ylab(NULL) +
  annotate("text", 2500,
           up_gsea[i, "enrichmentScore"] * .75, 
           label = lab, 
           hjust=0, vjust=0)
rankp <- plot_grid(plotlist=pp, ncol=1)
ggsave(rankp, filename = "gsearank_mono.pdf",width=8,height=10)



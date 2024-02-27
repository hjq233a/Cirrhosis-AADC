getwd()
setwd("C:/Users/86156/Desktop/")
library(clusterProfiler)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(pathview)
library("org.Hs.eg.db")
library("pSI")
#BiocManager::install("pSI")
library(vegan)
library(ggplot2)
library(limma)
library(dplyr)
library(org.Mm.eg.db) 

#数据读入与准备
data <- read.csv("filtered_output.txt", sep = "\t",check.names=F, header = TRUE)
data1<-as.data.frame(t(data)) 
colnames(data1) <- data1[1,]
data1<-data1[-1,]
data1$ID<-rownames(data1)
data1 <- data1[, c("ID", setdiff(names(data1), "ID"))] 
data[,1] <- sub("\\..*", "", data[,1])
write.table(data, "bulk.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

cell <- read.csv("cell.txt", sep = "\t",check.names=F, header = TRUE)
meta <- read.csv("meta.txt", sep = "\t",check.names=F, header = TRUE)
x<-merge(meta,data1,by="ID")


#one-way ANOVA test  followed by a post-hoc Tukey HSD test
#函数1
perform_anova_test <- function(agro1,agro2,agro3) {
CYCname<-c(agro1,agro2,agro3)
CYCname <- factor(CYCname, levels = c(agro1,agro2,agro3))
x<-x[x$Group2==CYCname[1]|x$Group2==CYCname[2]|x$Group2==CYCname[3],]

gene_name<-c()
pvalue<-c()
mean_1<-c()
mean_2<-c()
mean_3<-c()

for (i in 5:ncol(x))
{
  g1=colnames(x)[i]
  gene_name=c(gene_name,g1)
  val<-as.numeric(x[,i])
  gro <- x[,4]
  aa <-na.omit(data.frame(gro,val))
  model <- aov(val ~ gro,data=aa)
  pvalue=c(pvalue,summary(model)[[1]][["Pr(>F)"]][1])
  mean <- tapply(aa[,2],aa[,1],mean)
  a<-as.data.frame(mean)
  mean_1=c(mean_1,a[1,])
  mean_2=c(mean_2,a[2,])
  mean_3=c(mean_3,a[3,])
}
anova_results<-data.frame(gene_name,mean_1,mean_2,mean_3,pvalue)
# Filter significant genes based on a threshold (e.g., p-value < 0.05)
significant_genes <- anova_results$gene_name[anova_results$pvalue < 0.05]
significant_genes <- na.omit(significant_genes)


# Create a list to return multiple values
results_list <- list(CYCname = CYCname, significant_genes = significant_genes, anova_results = anova_results)

return(results_list)  # This returns the list
}


# Run Tukey's HSD test on the significant genes
#函数2
perform_tukey_test <- function(Tgro1, Tgro2) {
  list_of_tukey_tests <- list()
  groups<-x[x$Group2== Tgro1|x$Group2==Tgro2,"Group2"]
  
for (gene in significant_genes) {
  response <- unlist(x[x$Group2== Tgro1|x$Group2==Tgro2,gene])
  model <- aov(response ~ groups)
  tukey_results <- TukeyHSD(model)
  
  # Add the results to the list
  list_of_tukey_tests[[gene]] <- tukey_results
}


# Initialize an empty dataframe to store the results
all_results <- data.frame()

for (gene in names(list_of_tukey_tests)) {
  tukey <- list_of_tukey_tests[[gene]]
  
  # Extract the comparisons and p-values; these are stored in a complex structure within the TukeyHSD object
  tukey_summary <- as.data.frame(tukey$groups)  # assuming 'group' is your factor of interest
  
  # Since the results are for one gene, we create a new column to specify the gene name
  tukey_summary$Gene <- gene
  
  # Bind the results to the initialized dataframe
  all_results <- rbind(all_results, tukey_summary)
}

# Write the consolidated results to a text file
all_results %>%
  mutate(change=case_when(
    `p adj` < 0.05 & diff>0 ~ "up",
    `p adj` < 0.05 & diff<0 ~ "down",
    TRUE~ "NOT-CHANGE"
  )) ->dif

selected_genes_up <-  dif$Gene[dif$change == "up"]
selected_genes_down<-  dif$Gene[dif$change == "down"]
selected_genes_unchange <-  dif$Gene[dif$change == "NOT-CHANGE"]
tukey_list <- list(unchange =selected_genes_unchange, down = selected_genes_down, up = selected_genes_up, result=dif)

return(tukey_list)
}


# Main loop to process each condition
#需要手动改变参数1("FMT-H-CON","FMT-H-HE","FMT-H-HE+PDCI")  ("FMT-H-CON","FMT-H-HE","FMT-H-HE+PDCI") ("FMT-W-CON","FMT-W-HE","FMT-W-HE+PDCI")
results <- perform_anova_test("FMT-W-CON","FMT-W-HE","FMT-W-HE+PDCI")
CYCname <- results$CYCname
significant_genes <- results$significant_genes 
anova_results <- results$anova_results
write.csv(anova_results,"anova_results.cell.H.csv")

for (i in 1:3) {
  if (i < 3) {
    Tgro1 <- CYCname[i]
    Tgro2 <- CYCname[i + 1]
  } else {
    Tgro1 <- CYCname[1]
    Tgro2 <- CYCname[3]
  }
  
  file_name <- paste0(Tgro1, '_', Tgro2)
  print(file_name)
  # Perform Tukey tests
  tukey_tests_results <- perform_tukey_test(Tgro1, Tgro2)
  assign(file_name, tukey_tests_results, envir = .GlobalEnv)
  write.csv(tukey_tests_results$result,paste0(file_name,".csv"))
}

#需要手动改变参数2
geneset1<-intersect(`FMT-W-CON_FMT-W-HE`$up,`FMT-W-HE_FMT-W-HE+PDCI`$down)
geneset2<-intersect(`FMT-W-CON_FMT-W-HE`$up,`FMT-W-CON_FMT-W-HE+PDCI`$unchange)
geneset3<-intersect(`FMT-W-CON_FMT-W-HE`$down,`FMT-W-HE_FMT-W-HE+PDCI`$up)
geneset4<-intersect(`FMT-W-CON_FMT-W-HE`$down,`FMT-W-CON_FMT-W-HE+PDCI`$unchange)
geneset5<-`FMT-W-CON_FMT-W-HE`$up
geneset6<-`FMT-W-CON_FMT-W-HE`$down
#tukey_tests<-`FMT-W-CON_FMT-W-HE+PDCI`$result
geneset<-list(geneset1=geneset1,geneset2=geneset2,geneset3=geneset3,geneset4=geneset4,geneset5=geneset5,geneset6=geneset6)
for (gset in names(geneset)) {
write.csv(geneset[[gset]],paste0("FMT-W-",gset,".csv"))
}

#差异基因通路富集
#需要手动改变参数3,散在的4个w
list.geneset <- list(geneset1=geneset1, geneset2=geneset2,geneset3=geneset3,geneset4=geneset4,geneset5=geneset5,geneset6=geneset6)
for (var.name in names(list.geneset)){
  stripped_genes <- sub("\\..*", "", list.geneset[[var.name]])
  symbols <- mapIds(org.Mm.eg.db, 
                  keys=stripped_genes, 
                  column="SYMBOL",
                  keytype="ENSEMBL",
                  multiVals="first")

  #Go enrichment
  ego <- enrichGO(gene         = symbols,
                OrgDb        = org.Mm.eg.db,
                keyType      = 'SYMBOL',
                ont          = "ALL",  # BP,MF,CC
                pAdjustMethod = "BH",  # Benjamini-Hochberg
                qvalueCutoff = 0.05,
                readable     = TRUE)
  
  filename = paste0("P.",var.name, ".go.tif")
  print(filename)
  # Check if 'egg_results' is empty
  if (nrow(ego) == 0) {
    # 'ego_results' is empty; create a blank file
    # Create a blank plot
    blank_plot <- ggplot2::ggplot() + 
      ggplot2::geom_blank() +
      ggplot2::theme_void() # theme_void makes sure no axes or labels are plotted
    
    # Save the blank plot
    ggplot2::ggsave(filename = filename, plot = blank_plot, device = "tiff", width = 16, height = 18, units = "cm", dpi = 600)
    
  } else {
    # 'ego_results' is not empty; proceed as normal
    my_plot <- clusterProfiler::dotplot(ego, showCategory=10)
    ggplot2::ggsave(filename = filename, plot = my_plot, device = "tiff", width = 16, height = 18, units = "cm", dpi = 600)
  }
  

  # Extract the genes for each GO term
  ego_results <- as.data.frame(ego)
  write.csv(ego_results,paste0("P.",var.name, ".go.csv"))

  ##KEGG enrichment
  symbols <- mapIds(org.Mm.eg.db, 
                  keys=stripped_genes, 
                  column="ENTREZID",
                  keytype="ENSEMBL",
                  multiVals="first")
  kegg_results <- enrichKEGG(gene = symbols,
                           organism = 'mmu',  # 'mmu' is the code for Mus musculus (mouse)
                           keyType = 'kegg', 
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05)  # You can adjust these cutoffs as needed

  filename1 = paste0("P.",var.name, ".kegg.tif")
  print(filename1)
  
  # Check if 'kegg_results' is empty
  if (nrow(kegg_results) == 0) {
    # 'kegg_results' is empty; create a blank file
    # Create a blank plot
    blank_plot <- ggplot2::ggplot() + 
      ggplot2::geom_blank() +
      ggplot2::theme_void() # theme_void makes sure no axes or labels are plotted
    
    # Save the blank plot
    ggplot2::ggsave(filename = filename1, plot = blank_plot, device = "tiff", width = 16, height = 18, units = "cm", dpi = 600)
    
  } else {
    # 'kegg_results' is not empty; proceed as normal
    my_plot <- clusterProfiler::dotplot(kegg_results, showCategory=10)
    ggplot2::ggsave(filename = filename1, plot = my_plot, device = "tiff", width = 16, height = 18, units = "cm", dpi = 600)
  }
  

  kegg_results <- as.data.frame(kegg_results)
  write.csv(kegg_results,paste0("P.",var.name, ".kegg.csv"))
}





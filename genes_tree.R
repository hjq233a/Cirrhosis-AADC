library(tidytree)
library(ggtree)
library(tidydr)
library(tidyr)
library(ggtreeExtra)
library(ggstar)
library(ggplot2)
library(ggnewscale)
tree <- ape::read.tree("./COG0080_1.416.nwk") 
tree$tip.label <- sapply(strsplit(tree$tip.label, "__"),"[", 2)
tax_anno <- readxl::read_xlsx("./tax.xlsx")

df <- readxl::read_xlsx("./genes_profile04018.xlsx")
df_new <- df %>% as.matrix %>% t %>% as.data.frame
colnames(df_new) <- df_new[1,]
df_new <- df_new[-1,]


leaf_label <- tree$tip.label
extra <- rownames(df_new)[!rownames(df_new) %in% leaf_label]
meta_da <- df_new[!rownames(df_new) %in% extra,]

meta_da$id <- meta_da %>% rownames(meta_da)
meta_da <- meta_da %>% gather(key = "genes", value = "count", -id)
meta_da$count <- meta_da$count %>% as.numeric() 
meta_da$log_count <- (meta_da$count + 1) %>% log10()

id_uniqe <- meta_da$id %>% unique()
col_data <- meta_da %>% filter(count > 0) %>% group_by(id) %>% count(id) %>% data.frame()
count0_id <- id_uniqe[!id_uniqe %in% col_data$id]
count0_id_df <- data.frame(id = count0_id, n = 0)
col_data <- col_data %>% rbind(count0_id_df)
col_data <- col_data %>% left_join(tax_anno, by = c("id" = "ID"))
names(col_data)[3] <- "Phylum_group"

colors <- c("#9ACD32", "#EE6A50", "#87CEFA", "#FFC125", "#D15FEE", "#8DEEEE", "#800000",
            "#006400", "#800080", "#808080", "#B0171F", "#191970", "#7B68EE",
            "#00CD00", "Black")

# colors <- c( "#800000","#800080",
#             "#006400",  "#B0171F", "#191970", "#7B68EE",
#             "#00CD00", "Black")


p <- ggtree(tree,size = 0.15, layout = "circular") #+ geom_text2(aes(x= branch, label = node), size = 1.5)
p 
p <- p %<+% tax_anno


p_hig <- p + 
    geom_hilight(node=630,  alpha=.2,extendto=1.8,fill = "#e57055") +
    geom_cladelab(node=630, label="Bacteroidetes",
                  barsize=NA, fontsize=3.3, angle="auto", offset.text = 0.2,
                  hjust=0.5, horizontal=FALSE, fontface="italic"
    )+
    geom_hilight(node=276, fill="#9ccc3c", alpha=.2,extendto=1.8) + 
    geom_hilight(node=565, fill="#9ccc3c", alpha=.2,extendto=1.8) +
    geom_cladelab(node=565, label="Actinobacteria",
                  barsize=NA, fontsize=3.3, angle="auto",
                  hjust=0.5, horizontal=FALSE, fontface="italic",offset.text = -0.1
    )+
    geom_hilight(node=626, fill="#9aeced", alpha=.2,extendto=1.8) + 
    geom_hilight(node=203, fill="#7a1208", alpha=.2,extendto=1.8) + 
    geom_hilight(node=208, fill="#176307", alpha=.2,extendto=1.8) + 
    geom_hilight(node=c(696,487,64,729), fill="#f8c337", alpha=.2,extendto=1.8) + 
    geom_cladelab(node=487, label="Firmicutes",
                  barsize=NA, fontsize=3.3, angle="auto",
                  hjust=0.5, horizontal=FALSE, fontface="italic",offset.text = 0.45
    ) + 
    geom_hilight(node=749, fill="#cd64ec", alpha=.2,extendto=1.8) + 
    geom_cladelab(node=749, label="Fusobacteria",
                  barsize=NA, fontsize=2.5, angle="auto",
                  hjust=0.5, horizontal=FALSE, fontface="italic",offset.text = 0.8
    ) + 
    geom_hilight(node=c(699,756), fill="#9aeced", alpha=.2,extendto=1.8) +
    geom_cladelab(node=699, label="Proteobacteria",
                  barsize=NA, fontsize=3, angle="auto",
                  hjust=0.5, horizontal=FALSE, fontface="italic",offset.text = 0.4
    ) + 
    geom_hilight(node=c(761,361,444,418,824,819,777,449,466), fill="#f8c337", alpha=.2,extendto=1.8)  +
    new_scale_fill() + 
    geom_tippoint(
        mapping=aes(
            fill=Phylum
        ),
        shape=21,
        size=1.5,
        stroke=0.5,
        position="identity",
        show.legend = T
        
    ) + scale_fill_manual(values=colors,
                          guide=guide_legend(keywidth = 0.5, keyheight = 0.5, order=4,
                                             override.aes = list(starstroke=0.3))) +
    new_scale_fill()

#Calling colors
library(paletteer)
load("./sysdata.rda")
color_filter <- function(names = d_names, min = 15, max = 50) {
    
    leng_ma <- lapply(names, function(i) {
        length(paletteer_d(i))
    }) %>% do.call(rbind, .)
    
    len_df <- as.data.frame(leng_ma) ##配色方案的长度值
    len_df$paletteer <- d_names #赋值命名
    paletterr_d_15 <- len_df[len_df[1] > min & len_df[1] < max,]##筛选
    
    color_15 <- lapply(paletterr_d_15$paletteer, function(i) paletteer_d(i) ) ##获取每个配色方案的具体值
    names(color_15) <- paletterr_d_15$paletteer ##命名
    return(color_15)
}
color_12 <- color_filter(min = 12)

p1 <- p_hig + geom_fruit(
    data=meta_da,
    geom=geom_tile,
    mapping=aes(y=id, x=genes, fill=genes,alpha=log_count),
    color = NA,
    offset = 0.04,
    size = 0,
    pwidth = 0.5
    
) +
    scale_fill_manual(
        values = color_12$`miscpalettes::brightPastel`,
        #values=c(RColorBrewer::brewer.pal(12,"Paired")),
        #values=c("#800080","#FFA500","#FF0000","#800000","#006400","#0000FF","#696969","#d1651a","#de3b8b","#72c0a5","#fcff46","#87b0d2"),
        guide=guide_legend(keywidth=0.7, keyheight=0.7),
        name = "Genes"
    ) + 
    scale_alpha_continuous(range=c(0, 1.5),guide=guide_legend(keywidth=0.7, keyheight=0.7),
                           name=bquote(paste(Log[10],"(",.("Gene Copy Number + 1"), ")"))) +
    ggnewscale::new_scale_fill()

#ggsave("p1.pdf", p1, width = 10, height = 10)


p2 <- p1 + geom_fruit(
    data=col_data,
    geom=geom_col,
    mapping=aes(y=id, x=n, fill = Phylum_group),
    pwidth=0.38,
    orientation="y",
    position=position_stackx(),
    color = "white",
    show.legend = FALSE
) + scale_fill_manual(values=colors, guide="none")

ggsave("p3.pdf", p2, width = 10, height = 10)


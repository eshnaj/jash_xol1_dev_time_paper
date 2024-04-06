library('tidyr')
library('ggplot2')
library('ggprism')
library('data.table')
library('biomaRt')

time_resolved <- read.table("time_resolved_embryo_RNAseq.txt", 
                            header = TRUE)

time_resolved_embryo <- time_resolved[,1:19]
time_resolved_embryo <- drop_na(time_resolved_embryo)
time_resolved_embryo <- time_resolved_embryo[rowSums(time_resolved_embryo[2:19]) > 0,]
time_resolved_direction_binned <- data.frame("WormbaseID" = time_resolved_embryo$WormbaseName,
                                           "direction" = rep(NA, nrow(time_resolved_embryo)))

for (i in 1:nrow(time_resolved_embryo)) {
  if ( rowSums(time_resolved_embryo[i,2:10]) > rowSums(time_resolved_embryo[i,11:19]) ) {
    time_resolved_direction_binned[i,2] <- "up" 
  } 
  else if ( rowSums(time_resolved_embryo[i,2:10]) < rowSums(time_resolved_embryo[i,11:19]) ) {
    time_resolved_direction_binned[i,2] <- "down"
  } else {
    time_resolved_direction_binned[i,2] <- "unchanged"
  }
}

time_resolved_direction_binned$direction <- as.factor(time_resolved_direction_binned$direction)
summary(time_resolved_direction_binned$direction)

him8_data_ee <- read.table("DEraw_him8.csv")

#male-biased filtered
him8_male_bias <- subset(him8_data_ee, baseMean > 1 & padj < 0.05 & log2FoldChange > 0)

mart<-useDataset("celegans_gene_ensembl", 
                 useMart("ENSEMBL_MART_ENSEMBL", host="https://www.ensembl.org"))
attributes <- as.data.frame(mart@attributes[["name"]])

genes<-rownames(him8_male_bias)
gene_list<- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "wormbase_cds"),
  filters = "ensembl_gene_id",
  values = genes,
  mart = mart, 
  useCache = FALSE)
head(gene_list)
gene_list<-data.frame(gene_list)

gene_list_cds <- gene_list[,"wormbase_cds"]
gene_list_cds <- gsub('[abcdefghijklmnopqrstuvwxyz]', '', gene_list_cds)
gene_list$wormbase_cds <- gene_list_cds
gene_list <- gene_list[!duplicated(gene_list[,1:2]),]

him8_male_bias <-merge.data.frame(him8_male_bias, 
                                  gene_list, 
                                  by.x="row.names", 
                                  by.y = "ensembl_gene_id")

him8_male_bias_filter_binned <- him8_male_bias[!him8_male_bias$wormbase_cds %in% time_resolved_direction_binned$WormbaseID[time_resolved_direction_binned$direction == "down"],] 
him8_male_bias_filter_binned <- subset(him8_male_bias_filter_binned, log2FoldChange > 1.5)

xol1_ee <- read.table("xol1_N2_ee_raw.csv")
xol1_ee <- subset(xol1_ee, log2FoldChange != "NA")

xol1_male_bias_binned <- xol1_ee[rownames(xol1_ee) %in% him8_male_bias_filter_binned$Row.names,]
xol1_male_bias_binned$geneset <- rep("Male-biased Genes", nrow(xol1_male_bias_binned))
xol1_male_bias_binned_rest <- xol1_ee[!rownames(xol1_ee) %in% him8_male_bias_filter_binned$Row.names,]
xol1_male_bias_binned_rest$geneset <- rep("All Other Genes", nrow(xol1_male_bias_binned_rest))
xol1_male_bias_binned_ggplot <- rbind(xol1_male_bias_binned, xol1_male_bias_binned_rest)
xol1_male_bias_binned_ggplot$geneset <- as.factor(xol1_male_bias_binned_ggplot$geneset)
xol1_male_bias_binned_ggplot$geneset <- factor(xol1_male_bias_binned_ggplot$geneset, 
                                               levels = c("Male-biased Genes", "All Other Genes"))

plot_malebiasbinned <- ggplot(xol1_male_bias_binned_ggplot, 
                              aes(x = geneset, y = log2FoldChange, fill = geneset)) + 
  stat_boxplot(geom = "errorbar", width = 0.25) +  
  geom_boxplot(outlier.shape = NA, show.legend = FALSE, width = 0.5) +
  scale_fill_manual(values = c("#5065A7","lightgrey")) +
  coord_cartesian(ylim = c(-5,2.5)) +
  geom_hline(linetype = 2, 
             alpha = 0.3, 
             yintercept = median(xol1_male_bias_binned_ggplot$log2FoldChange[xol1_male_bias_binned_ggplot$geneset == "All Other Genes"])) +
  theme_prism(base_line_size = 0.5, base_fontface = "bold", base_family = "Arial") +
  theme(axis.text.x = element_text(size = 16, color = "black"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 18), 
        title = element_text(size = 18)) +
  theme(axis.title.y = element_blank()) +
  scale_x_discrete(labels = c("Male-biased \n Genes", " All Other \n Genes"))

plot_malebiasbinned

ggsave(plot = plot_malebiasbinned, 
       filename = "xol1_malebias.tiff", 
       path = "xol-1 paper/Data/", 
       height = 4, 
       width = 5)

wilcox.test(xol1_male_bias_binned_ggplot$log2FoldChange[xol1_male_bias_binned_ggplot$geneset == "Male-biased Genes"],
            xol1_male_bias_binned_ggplot$log2FoldChange[xol1_male_bias_binned_ggplot$geneset == "All Other Genes"])

#hermaphrodite-biased filtered
him8_herm_bias <- subset(him8_data_ee, baseMean > 1 & padj < 0.05 & log2FoldChange < 0)

mart<-useDataset("celegans_gene_ensembl", 
                 useMart("ENSEMBL_MART_ENSEMBL", host="https://www.ensembl.org"))
attributes <- as.data.frame(mart@attributes[["name"]])
genes<-rownames(him8_herm_bias)
gene_list<- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "wormbase_cds"),
  filters = "ensembl_gene_id",
  values = genes,
  mart = mart, 
  useCache = FALSE)
head(gene_list)
gene_list<-data.frame(gene_list)

gene_list_cds <- gene_list[,"wormbase_cds"]
gene_list_cds <- gsub('[abcdefghijklmnopqrstuvwxyz]', '', gene_list_cds)
gene_list$wormbase_cds <- gene_list_cds
gene_list <- gene_list[!duplicated(gene_list[,1:2]),]

him8_herm_bias <-merge.data.frame(him8_herm_bias, 
                                  gene_list, 
                                  by.x="row.names" , 
                                  by.y = "ensembl_gene_id")

him8_herm_bias_filter_binned <- him8_herm_bias[!him8_herm_bias$wormbase_cds %in% time_resolved_direction_binned$WormbaseID[time_resolved_direction_binned$direction == "up"],]
him8_herm_bias_filter_binned <- subset(him8_herm_bias_filter_binned, log2FoldChange < -1.5)

xol1_herm_bias_binned <- xol1_ee[rownames(xol1_ee) %in% him8_herm_bias_filter_binned$Row.names,]
xol1_herm_bias_binned$geneset <- rep("Hermaphrodite-biased Genes", nrow(xol1_herm_bias_binned))
xol1_herm_bias_binned_rest <- xol1_ee[!rownames(xol1_ee) %in% him8_herm_bias_filter_binned$Row.names,]
xol1_herm_bias_binned_rest$geneset <- rep("All Other Genes", nrow(xol1_herm_bias_binned_rest))
xol1_herm_bias_binned_ggplot <- rbind(xol1_herm_bias_binned, xol1_herm_bias_binned_rest)
xol1_herm_bias_binned_ggplot$geneset <- as.factor(xol1_herm_bias_binned_ggplot$geneset)
xol1_herm_bias_binned_ggplot$geneset <- factor(xol1_herm_bias_binned_ggplot$geneset, 
                                             levels = c("Hermaphrodite-biased Genes", "All Other Genes"))

plot_hermbiasbinned <- ggplot(xol1_herm_bias_binned_ggplot, 
                            aes(x = geneset, y = log2FoldChange, fill = geneset)) + 
  stat_boxplot(geom = "errorbar", width = 0.25) +  
  geom_boxplot(outlier.shape = NA, show.legend = FALSE, width = 0.5) +
  scale_fill_manual(values = c("#A367B1","lightgrey")) +
  coord_cartesian(ylim = c(-2.5,2.5)) +
  geom_hline(linetype = 2, 
             alpha = 0.3, 
             yintercept = median(xol1_herm_bias_binned_ggplot$log2FoldChange[xol1_herm_bias_binned_ggplot$geneset == "All Other Genes"])) +
  theme_prism(base_line_size = 0.5, base_fontface = "bold", base_family = "Arial") +
  theme(axis.text.x = element_text(size = 16, color = "black"), axis.title.x = element_blank(), axis.title.y = element_text(size = 18), title = element_text(size = 18)) +
  theme(axis.title.y = element_blank()) +
  scale_x_discrete(labels = c("Hermaphrodite-biased \n Genes", " All Other \n Genes"))

plot_hermbiasbinned

ggsave(plot = plot_hermbiasbinned, 
       filename = "xol1_hermbias.tiff", 
       path = "xol-1 paper/Data/", 
       height = 4, 
       width = 5)

wilcox.test(xol1_herm_bias_binned_ggplot$log2FoldChange[xol1_herm_bias_binned_ggplot$geneset == "Hermaphrodite-biased Genes"],
            xol1_herm_bias_binned_ggplot$log2FoldChange[xol1_herm_bias_binned_ggplot$geneset == "All Other Genes"])

### gmx for GSEA 
#filtered male-biased
him8_ee_up_filtered <- him8_male_bias_filter_binned[order(him8_male_bias_filter_binned$log2FoldChange, 
                                                          decreasing = TRUE),]
geneset_him8_up_filtered <- as.data.frame(him8_ee_up_filtered[,"Row.names"])
geneset_him8_up_filtered[nrow(geneset_him8_up_filtered)+1,1] <- "NA" 
tail(geneset_him8_up_filtered)
geneset_him8_up_filtered <- as.data.frame(geneset_him8_up_filtered[c(nrow(geneset_him8_up_filtered), 
                                                                     1:c(nrow(geneset_him8_up_filtered)-1)),])
colnames(geneset_him8_up_filtered) <- "him-8_upregulated_genes_filtered"

write.table(geneset_him8_up_filtered, 
            "geneset_him8_malebias_filtered.gmx", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

#filtered herm-biased
him8_ee_down_filtered <- him8_herm_bias_filter_binned[order(him8_herm_bias_filter_binned$log2FoldChange,
                                                            decreasing = TRUE),]
geneset_him8_down_filtered <- as.data.frame(him8_ee_down_filtered[,"Row.names"])
geneset_him8_down_filtered[nrow(geneset_him8_down_filtered)+1,1] <- "NA" 
tail(geneset_him8_down_filtered)
geneset_him8_down_filtered <- as.data.frame(geneset_him8_down_filtered[c(nrow(geneset_him8_down_filtered), 1:c(nrow(geneset_him8_down_filtered)-1)),])
colnames(geneset_him8_down_filtered) <- "him-8_downregulated_genes_filtered"
geneset_him8_down_filtered <- as.data.frame(geneset_him8_down_filtered[1:499,])

write.table(geneset_him8_down_filtered, 
            "geneset_him8_hermbias_filtered.gmx", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

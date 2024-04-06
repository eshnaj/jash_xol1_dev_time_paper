library('tidyr')
library('ggplot2')
library('ggprism')
library('data.table')

DEgenes_him8_ee <- read.table("DEgenes_him8_ee.csv")

him8_ee_up <- subset(DEgenes_him8_ee, log2FoldChange > 0.5 & padj < 0.05)
him8_ee_down <- subset(DEgenes_him8_ee, log2FoldChange < -0.5 & padj < 0.05)

write.table(him8_ee_up, file = "him8_ee_up.csv")
write.table(him8_ee_down, file = "him8_ee_down.csv")

xol1_ee_df <- read.table("xol1_N2_ee_raw.csv")

#male-biased unfiltered
xol1_ee_df$EnsemblID <- row.names(xol1_ee_df)
him8_ee_up$EnsemblID <- him8_ee_up$Row.names
him8_ee_up_above2 <- subset(him8_ee_up, log2FoldChange > 2 & padj < 0.05)

xol1_ee_him8up_above2 <- xol1_ee_df[xol1_ee_df$EnsemblID %in% him8_ee_up_above2$EnsemblID,]
xol1_ee_him8up_above2 <- subset(xol1_ee_him8up_above2, log2FoldChange != "NA")

xol1_ee_him8up_rest_above2 <- xol1_ee_df[!xol1_ee_df$EnsemblID %in% him8_ee_up_above2$EnsemblID,]
xol1_ee_him8up_rest_above2 <- subset(xol1_ee_him8up_rest_above2, log2FoldChange != "NA")

median(xol1_ee_him8up_above2$log2FoldChange)
#-0.5109585
median(xol1_ee_him8up_rest_above2$log2FoldChange)
#-0.04116502

xol1_ee_him8up_above2$genotype <- c(rep("Male-Biased Genes", nrow(xol1_ee_him8up_above2)))
xol1_ee_him8up_rest_above2$genotype <- c(rep("All Other Genes", nrow(xol1_ee_him8up_rest_above2)))
xol1_him8_up_above2 <- rbind(xol1_ee_him8up_above2,
                             xol1_ee_him8up_rest_above2)

xol1_him8_up_above2$genotype <- factor(xol1_him8_up_above2$genotype, 
                                       levels = c("Male-Biased Genes", "All Other Genes"))

plot_5 <- ggplot(xol1_him8_up_above2, aes(x=genotype, y=log2FoldChange, fill = genotype)) +
  stat_boxplot(geom = "errorbar", width = 0.25) +  
  geom_boxplot(outlier.shape = NA, width = 0.5, show.legend = FALSE) +
  coord_cartesian(ylim=(c(-5.4,4.2))) +
  ggtitle("unfiltered male-biased genes in xol-1/WT, > 2 log2FC") +
  geom_hline(yintercept = median(xol1_him8_up_above2$log2FoldChange[xol1_him8_up_above2$genotype == "All Other Genes"]), 
             linetype = 2, 
             alpha = 0.3) +
  scale_fill_manual(values = c("#5065A7","lightgrey")) +
  theme_prism(base_line_size = 0.5, base_fontface = "bold", base_family = "Arial") +
  theme(axis.text.x = element_text(size = 16, color = "black"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 18), 
        title = element_text(size = 8)) +
  theme(axis.title.y = element_blank()) +
  scale_x_discrete(labels = c("Male-biased \n Genes", " All Other \n Genes"))
plot_5

ggsave(filename = "xol1_malebias_unfiltered.tiff", 
       plot = plot_5, 
       path = "xol-1 paper/Data/", 
       height = 4, 
       width = 4)

#wilcox test
wilcox.test(xol1_ee_him8up_above2$log2FoldChange, xol1_ee_him8up_rest_above2$log2FoldChange)

#hermaphrodite-biased unfiltered
him8_ee_down$EnsemblID <- him8_ee_down$Row.names
him8_ee_down_below2 <- subset(him8_ee_down, log2FoldChange < -2 & padj < 0.05)

xol1_ee_him8down_below2 <- xol1_ee_df[xol1_ee_df$EnsemblID %in% him8_ee_down_below2$EnsemblID,]
xol1_ee_him8down_below2 <- subset(xol1_ee_him8down_below2, log2FoldChange != "NA")
xol1_ee_him8down_rest_below2 <- xol1_ee_df[!xol1_ee_df$EnsemblID %in% him8_ee_down_below2$EnsemblID,]
xol1_ee_him8down_rest_below2 <- subset(xol1_ee_him8down_rest_below2, log2FoldChange != "NA")

median(xol1_ee_him8down_below2$log2FoldChange)
#0.8128171
median(xol1_ee_him8down_rest_below2$log2FoldChange)
#-0.05036803

xol1_ee_him8down_below2$genotype <- c(rep("Hermaphrodite-Biased Genes", nrow(xol1_ee_him8down_below2)))
xol1_ee_him8down_rest_below2$genotype <- c(rep("All Other Genes", nrow(xol1_ee_him8down_rest_below2)))
xol1_him8_down_below2 <- rbind(xol1_ee_him8down_below2, xol1_ee_him8down_rest_below2)
xol1_him8_down_below2$genotype <- factor(xol1_him8_down_below2$genotype, 
                                         levels = c("Hermaphrodite-Biased Genes", "All Other Genes"))

plot_7 <- ggplot(xol1_him8_down_below2, aes(x=genotype, y=log2FoldChange, fill = genotype)) +
  stat_boxplot(geom = "errorbar", width = 0.25) +  
  geom_boxplot(outlier.shape = NA, width = 0.5, show.legend = FALSE) +
  coord_cartesian(ylim=(c(-2.7,2.5))) +
  ggtitle("unfiltered hermaphrodite-biased genes in xol-1/WT, > 2 log2FC") +
  geom_hline(yintercept = median(xol1_him8_down_below2$log2FoldChange[xol1_him8_down_below2$genotype == "All Other Genes"]), linetype = 2, alpha = 0.3) +
  scale_fill_manual(values = c("#A367B1","lightgrey")) +
  theme_prism(base_line_size = 0.5, base_fontface = "bold", base_family = "Arial") +
  theme(axis.text.x = element_text(size = 14, color = "black"), axis.title.x = element_blank(), axis.title.y = element_text(size = 18), title = element_text(size = 8)) +
  theme(axis.title.y = element_blank()) +
  scale_x_discrete(labels = c("Hermaphrodite-biased \n Genes", " All Other \n Genes"))
plot_7

ggsave(filename = "xol1_hermbias_unfiltered.tiff", 
       plot = plot_7, 
       path = "xol-1 paper/Data/", 
       height = 4, 
       width = 4)

wilcox.test(xol1_ee_him8down_below2$log2FoldChange, xol1_ee_him8down_rest_below2$log2FoldChange)

### gmx for GSEA 
#unfiltered male-biased
him8_ee_up_unfiltered <- him8_ee_up_above2[order(him8_ee_up_above2$log2FoldChange, 
                                                 decreasing = TRUE),]
geneset_him8_up_unfiltered <- as.data.frame(him8_ee_up_unfiltered[,"EnsemblID"])
geneset_him8_up_unfiltered[nrow(geneset_him8_up_unfiltered)+1,1] <- "NA" 
tail(geneset_him8_up_unfiltered)
geneset_him8_up_unfiltered <- as.data.frame(geneset_him8_up_unfiltered[c(nrow(geneset_him8_up_unfiltered), 
                                                                         1:c(nrow(geneset_him8_up_unfiltered)-1)),])
colnames(geneset_him8_up_unfiltered) <- "him-8_upregulated_genes_unfiltered"

write.table(geneset_him8_up_unfiltered, 
            "geneset_him8_malebias_unfiltered.gmx", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

#unfiltered hermaphrodite-biased
him8_ee_down_unfiltered <- him8_ee_down_below2[order(him8_ee_down_below2$log2FoldChange, 
                                                     decreasing = TRUE),]
geneset_him8_down_unfiltered <- as.data.frame(him8_ee_down_unfiltered[,"EnsemblID"])
geneset_him8_down_unfiltered[nrow(geneset_him8_down_unfiltered)+1,1] <- "NA" 
tail(geneset_him8_down_unfiltered)
geneset_him8_down_unfiltered <- as.data.frame(geneset_him8_down_unfiltered[c(nrow(geneset_him8_down_unfiltered), 
                                                                             1:c(nrow(geneset_him8_down_unfiltered)-1)),])
colnames(geneset_him8_down_unfiltered) <- "him-8_downregulated_genes_unfiltered"

write.table(geneset_him8_down_unfiltered, 
            "geneset_him8_hermbias_unfiltered.gmx", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)
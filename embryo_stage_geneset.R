library('ggplot2')
library('data.table')
library('tidyr')
library('tidyverse')

xol1_ee_data <- read.table("xol1_N2_ee_raw.csv")
xol1_ee_data <- subset(xol1_ee_data, baseMean > 1)
xol1_ee_data$EnsemblID <- row.names(xol1_ee_data)

#importing gene expression datasets from https://www.vanderbilt.edu/wormdoc/wormmap/Welcome.html
EE_genes <- read.table("/Volumes/lsa-gyorgyi/Eshna/Datasets/EE_EGs.txt", header = TRUE)
LE_genes <- read.table("/Volumes/lsa-gyorgyi/Eshna/Datasets/LE_EGs.txt", header = TRUE)
L1_genes <- read.table("/Volumes/lsa-gyorgyi/Eshna/Datasets/L1_EGs.txt", header = TRUE)
L2_genes <- read.table("/Volumes/lsa-gyorgyi/Eshna/Datasets/L2_EGs.txt", header = TRUE)
L3_genes <- read.table("/Volumes/lsa-gyorgyi/Eshna/Datasets/L3_EGs.txt", header = TRUE)
L4_genes <- read.table("/Volumes/lsa-gyorgyi/Eshna/Datasets/L4_EGs.txt", header = TRUE)

common_embryo_genes <- EE_genes$gene_id[EE_genes$gene_id %in% LE_genes$gene_id] 
common_larval_genes <- L1_genes$gene_id[L1_genes$gene_id %in% L2_genes$gene_id &
                                          L1_genes$gene_id %in% L3_genes$gene_id &
                                          L1_genes$gene_id %in% L4_genes$gene_id]
common_embryo_larval_genes <- EE_genes$gene_id[EE_genes$gene_id %in% common_embryo_genes &
                                                 EE_genes$gene_id %in% common_larval_genes]
LE_genes_core <- LE_genes[!LE_genes$gene_id %in% common_embryo_genes,]
EE_genes_core <- EE_genes[!EE_genes$gene_id %in% common_embryo_genes,]
L1_genes_core <- L1_genes[!L1_genes$gene_id %in% common_embryo_larval_genes,]

# early embryonic genes
xol1_ee_eegenes <- xol1_ee_data[xol1_ee_data$EnsemblID %in% EE_genes_core$gene_id,]
xol1_ee_eegenes$geneset <- rep("Early Embryonic Genes", nrow(xol1_ee_eegenes))
xol1_ee_noteegenes <- xol1_ee_data[!xol1_ee_data$EnsemblID %in% EE_genes_core$gene_id,]
xol1_ee_noteegenes$geneset <- rep("All Other Genes", nrow(xol1_ee_noteegenes))
xol1_ee_eegenes_labeled <- rbind(xol1_ee_eegenes, xol1_ee_noteegenes)
xol1_ee_eegenes_labeled$geneset <- as.factor(xol1_ee_eegenes_labeled$geneset)
xol1_ee_eegenes_labeled$geneset <- factor(xol1_ee_eegenes_labeled$geneset, 
                                          levels = c("Early Embryonic Genes", "All Other Genes"))

plot_EE <- ggplot(xol1_ee_eegenes_labeled, 
                  aes(x = as.factor(geneset), y = log2FoldChange, fill = geneset)) + 
  stat_boxplot(geom = "errorbar", width = 0.25) +  
  geom_boxplot(outlier.shape = NA, show.legend = FALSE, width = 0.5) + 
  scale_fill_manual(values = c("#b7c1de", "lightgrey")) +
  coord_cartesian(ylim = c(-5,3.5)) +
  geom_hline(linetype = 2, 
             alpha = 0.5, 
             yintercept = median(xol1_ee_eegenes_labeled$log2FoldChange[xol1_ee_eegenes_labeled$geneset == "All Other Genes"])) +
  theme_prism(base_line_size = 0.5, 
              base_fontface = "bold", 
              base_family = "Arial") +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 16, color = "black"), 
        axis.title.y = element_text(size = 18)) +
  labs(y = substitute(paste('Log2 Fold Change (',italic('xol-1'), ' vs WT)'))) +
  scale_x_discrete(labels = c("Early Embryonic \n Genes", " All Other \n Genes"))

plot_EE

ggsave(plot = plot_EE, 
       filename = "ee_enriched_xol1.tiff", 
       path = "xol-1 paper/Data/", 
       height = 4, 
       width = 5)

wilcox.test(xol1_ee_eegenes$log2FoldChange, xol1_ee_noteegenes$log2FoldChange)

# late embryonic genes
xol1_ee_legenes <- xol1_ee_data[xol1_ee_data$EnsemblID %in% LE_genes_core$gene_id,]
xol1_ee_legenes$geneset <- rep("Late Embryonic Genes", nrow(xol1_ee_legenes))
xol1_ee_notlegenes <- xol1_ee_data[!xol1_ee_data$EnsemblID %in% LE_genes_core$gene_id,]
xol1_ee_notlegenes$geneset <- rep("All Other Genes", nrow(xol1_ee_notlegenes))
xol1_ee_legenes_labeled <- rbind(xol1_ee_legenes, xol1_ee_notlegenes)
xol1_ee_legenes_labeled$geneset <- as.factor(xol1_ee_legenes_labeled$geneset)
xol1_ee_legenes_labeled$geneset <- factor(xol1_ee_legenes_labeled$geneset, 
                                          levels = c("Late Embryonic Genes", "All Other Genes"))


plot_LE <- ggplot(xol1_ee_legenes_labeled, 
                  aes(x = geneset, y = log2FoldChange, fill = geneset)) + 
  stat_boxplot(geom = "errorbar", width = 0.25) +  
  geom_boxplot(outlier.shape = NA, show.legend = FALSE, width = 0.5) + 
  coord_cartesian(ylim = c(-2.7,2.7)) +
  scale_fill_manual(values = c("#b7c1de", "lightgrey")) +
  geom_hline(linetype = 2, 
             alpha = 0.5, 
             color = "black", 
             yintercept = median(xol1_ee_legenes_labeled$log2FoldChange[xol1_ee_legenes_labeled$geneset == "All Other Genes"])) +
  theme_prism(base_line_size = 0.5, base_fontface = "bold", base_family = "Arial") +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 16, color = "black"), 
        axis.title.y = element_text(size = 18)) +
  labs(y = substitute(paste('Log2 Fold Change (',italic('xol-1'), ' vs WT)'))) +
  scale_x_discrete(labels = c("Late Embryonic \n Genes", " All Other \n Genes"))

plot_LE

ggsave(plot = plot_LE, 
       filename = "le_enriched_xol1.tiff", 
       path = "xol-1 paper/Data/", 
       height = 5, 
       width = 5)

wilcox.test(xol1_ee_legenes$log2FoldChange, xol1_ee_notlegenes$log2FoldChange)

# L1 larval genes
xol1_ee_l1genes <- xol1_ee_data[xol1_ee_data$EnsemblID %in% L1_genes_core$gene_id,]
xol1_ee_l1genes$geneset <- rep("L1 Larval Genes", nrow(xol1_ee_l1genes))
xol1_ee_notl1genes <- xol1_ee_data[!xol1_ee_data$EnsemblID %in% L1_genes_core$gene_id,]
xol1_ee_notl1genes$geneset <- rep("All Other Genes", nrow(xol1_ee_notl1genes))
xol1_ee_l1genes_labeled <- rbind(xol1_ee_l1genes, xol1_ee_notl1genes)
xol1_ee_l1genes_labeled$geneset <- as.factor(xol1_ee_l1genes_labeled$geneset)
xol1_ee_l1genes_labeled$geneset <- factor(xol1_ee_l1genes_labeled$geneset, 
                                          levels = c("L1 Larval Genes", "All Other Genes"))

plot_L1 <- ggplot(xol1_ee_l1genes_labeled, 
                  aes(x = as.factor(geneset), y = log2FoldChange, fill = geneset)) + 
  stat_boxplot(geom = "errorbar", width = 0.25) +  
  geom_boxplot(outlier.shape = NA, show.legend = FALSE, width = 0.5) +
  scale_fill_manual(values = c("#b7c1de", "lightgrey")) +
  coord_cartesian(ylim = c(-4,4)) +
  geom_hline(linetype = 2, 
             alpha = 0.5, 
             yintercept = median(xol1_ee_l1genes_labeled$log2FoldChange[xol1_ee_l1genes_labeled$geneset == "All Other Genes"])) +
  theme_prism(base_line_size = 0.5, base_fontface = "bold", base_family = "Arial") +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 16, color = "black"), 
        axis.title.y = element_text(size = 18)) +
  labs(y = substitute(paste('Log2 Fold Change (',italic('xol-1'), ' vs WT)'))) +
  scale_x_discrete(labels = c("L1 Larval \n Genes", " All Other \n Genes"))

plot_L1

ggsave(plot = plot_L1, 
       filename = "l1_enriched_xol1.tiff", 
       path = "xol-1 paper/Data/", 
       height = 4.5, 
       width = 5)

wilcox.test(xol1_ee_l1genes$log2FoldChange, xol1_ee_notl1genes$log2FoldChange)
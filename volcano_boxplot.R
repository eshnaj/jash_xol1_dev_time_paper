library('biomaRt')
library('AnnotationHub')
library('ensembldb')
library('data.table')
library('tidyr')
library('ggplot2')
library('ggprism')
library('tidyverse')

xol1_ee_data <- read.table("xol1_N2_ee_raw.csv")

mart <- useDataset("celegans_gene_ensembl", 
                 useMart("ENSEMBL_MART_ENSEMBL", 
                 host="https://www.ensembl.org"))
genes <- row.names(xol1_ee_data)
gene_list <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "description"),
  filters = "ensembl_gene_id",
  values = genes,
  mart = mart, 
  useCache = FALSE)
head(gene_list)

gene_list <- data.frame(gene_list)
xol1_ee_genes <- merge.data.frame(xol1_ee_data, 
                                  gene_list, 
                                  by.x=0, 
                                  by.y = "ensembl_gene_id")

xol1_l1_data <- read.table(file = "xol1_N2_l1_raw.csv")

mart <- useDataset("celegans_gene_ensembl", 
                   useMart("ENSEMBL_MART_ENSEMBL", 
                           host="https://www.ensembl.org"))
genes <- rownames(xol1_l1_data)
gene_list <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "description"),
  filters = "ensembl_gene_id",
  values = genes,
  mart = mart, 
  useCache = FALSE)
head(gene_list)

gene_list <- data.frame(gene_list)
xol1_l1_genes <- merge.data.frame(xol1_l1_data, 
                                  gene_list,by.x="row.names", 
                                  by.y = "ensembl_gene_id")
head(xol1_l1_genes)

#early embryo volcano plot
volcano_ee <- ggplot(data = xol1_ee_genes, 
                     aes(x = log2FoldChange, y = -log10(padj), 
                         color = abs(log2FoldChange) > 1.5 & padj<0.05)) +
  geom_point(show.legend = FALSE, alpha = 0.3, size = 2) +
  scale_color_manual(values = c("lightgrey", "#092047")) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, alpha = 0.7) +
  geom_vline(xintercept = c(-1.5,1.5), linetype = 2, alpha = 0.7) +
  coord_cartesian(ylim = c(0,100), xlim = c(-12,12)) +
  theme_prism(base_line_size = 0.5, base_fontface = "bold", base_family = "Arial") +
  labs(x = "Log2 Fold Change (Early Embryos)", y = "-log10 (adjusted p-value)") +
  theme(axis.title.x = element_text(size = 17), axis.title.y = element_text(size = 17))

volcano_ee

ggsave(plot = volcano1_ee, 
       filename = "ee_volcano_xol1.tiff", 
       path = "/xol-1 paper/Data/", 
       height = 5, 
       width = 6)

##L1 volcano plot
volcano_l1 <- ggplot(xol1_l1_genes, 
                     aes(x = log2FoldChange, y = -log10(padj), 
                         color = abs(log2FoldChange) > 1.5 & padj<0.05)) + 
  geom_point(show.legend = FALSE, alpha = 0.4, size = 2) + 
  coord_cartesian(ylim=c(0,100), xlim = c(-15,15)) + 
  scale_color_manual(values = c("lightgrey", "#092047")) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, alpha = 0.7) +
  geom_vline(xintercept = c(-1.5,1.5), linetype = 2, alpha = 0.7) +
  theme_prism(base_line_size = 0.5) +
  labs(x = "Log2 Fold Change (L1 Larvae)", y = ("-log10 (adjusted p-value)")) +
  theme(axis.title.x = element_text(size = 17), axis.title.y = element_text(size = 17))
volcano_l1

ggsave(plot = volcano_l1, 
       filename = "l1_volcano_xol1.tiff", 
       path = "xol-1 paper/Data/" , height = 5, width = 6)

#X-chromosome expression
add_XorA <- function (inputdataframe) {
XorA <- data.frame(XorA = character(), stringsAsFactors = FALSE)
for (i in 1:nrow(inputdataframe)) {
  if (inputdataframe[i,"chromosome_name"] == "X"){
    XorA[i,1] <- "X"
  } else {
    XorA[i,1] <- "A"
  }
}
row.names(XorA) <- inputdataframe$Row.names
outputdataframe <- merge.data.frame(inputdataframe, 
                                  XorA, 
                                  by.x = "Row.names", 
                                  by.y = 0)
return(outputdataframe)
}

#early embryo
xol1_ee_genes <- add_XorA(xol1_ee_genes)
xol1_ee_genes$XorA <- as.factor(xol1_ee_genes$XorA)
xol1_ee_genes <- xol1_ee_genes %>% drop_na(log2FoldChange)
xol1_ee_genes <- subset(xol1_ee_genes, chromosome_name != "MtDNA")

ggplot(data = subset(xol1_ee_genes, baseMean > 1), 
       aes(x = XorA, y = log2FoldChange)) + 
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(-3,3)) +
  geom_hline(yintercept = median(xol1_ee_genes$log2FoldChange), 
             linetype = 2, 
             alpha = 0.5)

ggsave(filename = "X_A_boxplot_xol1_WT_ee.png", 
       path = "xol1/")

#L1
xol1_l1_genes <- add_XorA(xol1_l1_genes)
xol1_l1_genes$XorA <- as.factor(xol1_l1_genes$XorA)
xol1_l1_genes <- xol1_l1_genes %>% drop_na(log2FoldChange)
xol1_l1_genes <- subset(xol1_l1_genes, chromosome_name != "MtDNA")

ggplot(data = subset(xol1_l1_genes, baseMean > 1), 
       aes(x = XorA, y = log2FoldChange)) + 
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(-3,3)) +
  geom_hline(yintercept = median(xol1_l1_genes$log2FoldChange), 
             linetype = 2, 
             alpha = 0.5)

ggsave(filename = "X_A_boxplot_xol1_WT_l1.png", 
       path = "xol1/")
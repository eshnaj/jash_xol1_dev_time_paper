library('DESeq2')
library('tximport')
library('ensembldb')
library('AnnotationHub')
library('biomaRt')
library('tidyr')
library('ggplot2')
library('ggprism')

#diff expression for him-8/WT
data_him8 <- data.frame(sampleID=rep("PH", 6),
                                 condition=rep("PH", 6),
                                 genotype_rep=rep("PH", 6))
data_him8$genotype_rep <- c("N2_rep1",
                            "N2_rep2",
                            "N2_rep3",
                            "him8_rep1",
                            "him8_rep2",
                            "him8_rep3")
for (i in c(1:6)) {
  y=paste0("8125-EJ-",i)
  data_him8[i,1]= y
}
data_him8$condition <- c(rep("N2", 3), 
                         rep("him8", 3))

dir_him8 <- "salmon_reads/"
list.files(dir_him8)
files <- file.path(dir_him8, 
                   paste0("counts_salmon_", data_him8$sampleID), 
                   "quant.sf")
names(files) <- data_him8$sampleID
files
all(file.exists(files))

celegansdb_formalclassobject <- query(AnnotationHub(), 
                                      pattern = c("Caenorhabditis elegans", "EnsDb", "109"))
celegansdb <- celegansdb_formalclassobject[[1]]
genes <- genes(celegansdb)
tx2gene <- data.frame(TXNAME=genes$canonical_transcript,
                      GENEID=genes$gene_id)
head(tx2gene)

salmon_him8_N2 <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(salmon_him8_N2)
head(salmon_him8_N2$counts)

ddshim8 <- DESeqDataSetFromTximport(salmon_him8_N2,
                                    colData = data_him8,
                                    design = ~ condition)
ddshim8$condition <- relevel(ddshim8$condition, ref = "N2")
DErunhim8 <- DESeq(ddshim8)
DErun_him8_res <- results(DErunhim8)
head(DErun_him8_res)

#fpkm
fpkm_him8 <- fpkm(DErunhim8)
colnames(fpkm_him8) <- c("WT_rep1", 
                         "WT_rep2",
                         "WT_rep3", 
                         "him8_rep1", 
                         "him8_rep2", 
                         "him8_rep3")
write.table(fpkm_him8, 
            file = "xol1_paper/fpkm_him8_ee.txt",
            sep = "\t",
            quote = FALSE, 
            na = "NA", 
            row.names = TRUE,
            col.names = TRUE)

rld <- rlog(DErunhim8, blind = TRUE)
plotPCA(rld)

write.table(DErun_him8_res, 
            file = "DEraw_him8.csv")

#volcano plot
him8_ee <- read.table(file = "DEraw_him8.csv")

volcano <- ggplot(subset(him8_ee,c(baseMean > 1 & abs(log2FoldChange) < 1.5 | padj > 0.05)), 
                  aes(x = log2FoldChange,
                      y = -log10(padj))) +
  geom_point(show.legend = FALSE, alpha = 0.4, size = 2, color = "lightgrey") +
  geom_point(color = "#5065A7", alpha = 0.4, size = 2, data = subset(him8_ee, c(baseMean>1 & log2FoldChange > 1.5 & padj < 0.05))) +
  geom_point(color = "#A367B1", alpha = 0.4, size = 2, data = subset(him8_ee, c(baseMean>1 & log2FoldChange < -1.5 & padj < 0.05))) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = c(-1.5,1.5), linetype = 2, alpha = 0.5) +
  coord_cartesian(ylim=c(0,50), xlim = c(-12,12)) +
  labs(x = substitute(paste('Log2 Fold Change (', italic('him-8'), ' vs WT)')), y = "-log10 (adjusted p-value)") +
  theme_prism(base_line_size = 0.5, base_fontface = "bold", base_family = "Arial") +
  theme(axis.title.x = element_text(size = 17), axis.title.y = element_text(size = 17))

volcano

ggsave("him8_ee_volcano.png", plot = last_plot())

##chromosome and gene name
him8_ee <- read.table(file = "DEraw_him8.csv")
mart<-useDataset("celegans_gene_ensembl", 
                 useMart("ENSEMBL_MART_ENSEMBL", 
                         host="https://www.ensembl.org"))

genes<-rownames(him8_ee)
gene_list<- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "description"),
  filters = "ensembl_gene_id",
  values = genes,
  mart = mart, 
  useCache = FALSE)

gene_list<-data.frame(gene_list)
him8_ee_genes<-merge.data.frame(him8_ee, 
                                gene_list,
                                by.x="row.names" , 
                                by.y = "ensembl_gene_id")
head(him8_ee_genes)
write.table(him8_ee_genes, file="DEgenes_him8_ee.csv")

#chromosome-wide differential expression analysis
DEgenes_him8_ee <- read.table("DEgenes_him8_ee.csv")
row.names(DEgenes_him8_ee) <- DEgenes_him8_ee$Row.names
DEgenes_him8_ee <- DEgenes_him8_ee[,-1]
DEgenes_him8_ee <- DEgenes_him8_ee %>% drop_na(log2FoldChange)
DEgenes_him8_ee <- subset(DEgenes_him8_ee, chromosome_name != "MtDNA")

XorA <- data.frame(XorA = character(), stringsAsFactors = FALSE)
for (i in 1:nrow(DEgenes_him8_ee)) {
  if (DEgenes_him8_ee[i,8] == "X"){
    XorA[i,1] <- "X"
  } else {
    XorA[i,1] <- "A"
  }
}
row.names(XorA) <- row.names(DEgenes_him8_ee)

DEgenes_him8_ee <- merge.data.frame(DEgenes_him8_ee, XorA, by = c(0,0))
row.names(DEgenes_him8_ee) <- DEgenes_him8_ee$Row.names
DEgenes_him8_ee <- DEgenes_him8_ee[,-1]
DEgenes_him8_ee$XorA <- as.factor(DEgenes_him8_ee$XorA)
DEgenes_him8_ee$chromosome_name <- as.factor(DEgenes_him8_ee$chromosome_name)

ggplot(data=subset(DEgenes_him8_ee, baseMean > 1), 
       aes(x=XorA, y=log2FoldChange)) + 
  geom_boxplot(outlier.shape = NA) +
  ggtitle("him-8/N2 (early embryos)") +
  coord_cartesian(ylim=c(-1.9,1.9))

median_X <- subset(DEgenes_him8_ee, XorA=="X" & baseMean>1)
median(median_X$log2FoldChange)
#-0.2545882
median_A <- subset(DEgenes_him8_ee, XorA=="A" & baseMean>1)
median(median_A$log2FoldChange)
#-0.01457317
wilcox.test(median_A$log2FoldChange, median_X$log2FoldChange, 
            correct = TRUE, conf.int = TRUE)
#p-value = 2.2e-16
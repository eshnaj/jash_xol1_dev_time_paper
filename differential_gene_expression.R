library('biomaRt')
library('tidyr')
library('AnnotationHub')
library('ensembldb')
library('data.table')
library('DESeq2')

#early embryos
data_xol1_ee <- data.frame(sampleID=rep("PH", 6),
                           condition=rep("PH", 6),
                           genotype_rep=rep("PH", 6))
data_xol1_ee$genotype_rep <- c("N2_rep1",
                               "N2_rep2",
                               "N2_rep3",
                               "xol1_rep1",
                               "xol1_rep2",
                               "xol1_rep3")
#setting sampleIDs
for (i in c(1:6)) {
  y=paste0("3489-EJ-",i)
  data_xol1_ee[i,1]= y
}
data_xol1_ee$condition <- c(rep("N2", 3), rep("xol1", 3))

dir_xol1_ee <- "salmon_reads/"
list.files(dir_xol1_ee)
files <- file.path(dir_xol1_ee, 
                   paste0("counts_salmon_Sample_", data_xol1_ee$sampleID), 
                   "quant.sf")
names(files) <- data_xol1_ee$sampleID
files
all(file.exists(files))

celegansdb_formalclassobject <- query(AnnotationHub(), 
                                      pattern = c("Caenorhabditis elegans", "EnsDb", "109"))
celegansdb <- celegansdb_formalclassobject[[1]]
genes <- genes(celegansdb)
tx2gene <- data.frame(TXNAME=genes$canonical_transcript,
                      GENEID=genes$gene_id)
head(tx2gene)

salmon_xol1_ee <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(salmon_xol1_ee)
head(salmon_xol1_ee$counts)

ddsxol1ee <- DESeqDataSetFromTximport(salmon_xol1_ee,
                                      colData = data_xol1_ee,
                                      design = ~ condition)
ddsxol1ee$condition <- relevel(ddsxol1ee$condition, ref = "N2")
DErun_xol1_ee <- DESeq(ddsxol1ee)
DErun_xol1_ee_res <- results(DErun_xol1_ee)
head(DErun_xol1_ee_res)

#RDS
saveRDS(DErun_xol1_ee,
        file = "DE_run_xol1_ee.rds")

#fpkm
fpkm_xol1_ee <- fpkm(DErun_xol1_ee)
colnames(fpkm_xol1_ee) <- c("WT_rep1", 
                            "WT_rep2", 
                            "WT_rep3", 
                            "xol1_rep1", 
                            "xol1_rep2", 
                            "xol1_rep3")
write.table(fpkm_xol1_ee, 
            file = "xol1_paper/fpkm_xol1_ee.txt",
            sep = "\t",
            quote = FALSE, 
            na = "NA", 
            row.names = TRUE, 
            col.names = TRUE)

rld <- rlog(DErunxol1ee, blind = TRUE)
plotPCA(rld)

write.table(DErun_xol1_ee_res,
            "xol1_N2_ee_raw.csv")

#L1
data_xol1_L1 <- data.frame(sampleID=rep("PH", 12),
                           condition=rep("PH", 12),
                           genotype_rep=rep("PH", 12))
data_xol1_L1$genotype_rep <- c("N2_rep1",
                               "N2_rep2",
                               "N2_rep3",
                               "N2_rep4",
                               "xol1sex1_rep1",
                               "xol1sex1_rep2",
                               "xol1sex1_rep3",
                               "xol1sex1_rep4",
                               "xol1_rep1",
                               "xol1_rep2",
                               "xol1_rep3",
                               "xol1_rep4")
#setting sampleIDs
for (i in c(1:12)) {
  y=paste0("4818-EJ-",i)
  data_xol1_L1[i,1]= y
}
data_xol1_L1$condition <- c(rep("N2", 4), 
                            rep("xol1sex1", 4), 
                            rep("xol1", 4))
data_xol1_L1 <- data_xol1_L1[-11,]

dir_xol1_L1 <- "L1s_salmon_reads/"
list.files(dir_xol1_L1)
files <- file.path(dir_xol1_L1, 
                   paste0("counts_salmon_Sample_", data_xol1_L1$sampleID), 
                   "quant.sf")
names(files) <- data_xol1_L1$sampleID
files
all(file.exists(files))

celegansdb_formalclassobject <- query(AnnotationHub(), 
                                      pattern = c("Caenorhabditis elegans", "EnsDb", "109"))
celegansdb <- celegansdb_formalclassobject[[1]]
genes <- genes(celegansdb)
tx2gene <- data.frame(TXNAME=genes$canonical_transcript,
                      GENEID=genes$gene_id)
head(tx2gene)

salmon_xol1_L1 <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(salmon_xol1_L1)
head(salmon_xol1_L1$counts)

#piping into DESeq2
ddsL1nxxs <- DESeqDataSetFromTximport(salmon_xol1_L1,
                                      colData = data_xol1_L1,
                                      design = ~ condition)
ddsL1nxxs$condition <- relevel(ddsL1nxxs$condition, ref = "N2")
DErun_xol1_L1 <- DESeq(ddsL1nxxs)
DErun_xol1_L1_res <- results(DErun_xol1_L1)
head(DErun_xol1_L1_res)

#RDS
saveRDS(DErun_xol1_L1,
        file = "DE_run_xol1_l1.rds")
#fpkm
fpkm_xol1_l1 <- fpkm(DErun_xol1_L1)
fpkm_xol1_l1 <- fpkm_xol1_l1[,-5:-8]
colnames(fpkm_xol1_l1) <- c("WT_rep1", 
                            "WT_rep2", 
                            "WT_rep3", 
                            "WT_rep4", 
                            "xol1_rep1", 
                            "xol1_rep2", 
                            "xol1_rep3")
write.table(fpkm_xol1_l1, 
            file = "/Volumes/lsa-gyorgyi/Eshna/xol1_paper/fpkm_xol1_l1.txt",
            sep = "\t",
            quote = FALSE, 
            na = "NA", 
            row.names = TRUE, 
            col.names = TRUE)

rld <- rlog(DErun_xol1_L1, blind = TRUE)
plotPCA(rld)

write.table(DErun_xol1_L1_res, 
            file = "DEraw_xol1_L1.csv")
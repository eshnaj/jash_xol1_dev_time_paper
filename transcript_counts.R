library('biomaRt')
library('ggtext')
library('stringr')
library('ggprism')
library('ggplot2')

time_resolved <- read.table("time_resolved_embryo_RNAseq.txt", 
                            header = TRUE)

dim(time_resolved)
time_resolved_embryo <- time_resolved[,-20:-32]

mart<-useDataset("celegans_gene_ensembl", 
                   useMart("ENSEMBL_MART_ENSEMBL", 
                   host="https://www.ensembl.org"))
genes<-time_resolved_embryo$WormbaseName
gene_list<- getBM(
attributes = c("ensembl_gene_id", "external_gene_name", "wormbase_cds"),
filters = "wormbase_cds",
values = genes,
mart = mart, 
useCache = FALSE)

time_resolved_embryo <- merge.data.frame(time_resolved_embryo, 
                                         gene_list, 
                                         by.x = "WormbaseName", 
                                         by.y = "wormbase_cds", 
                                         all.x = TRUE)
  
time_resolved_embryo <- time_resolved_embryo[rowSums(time_resolved_embryo[,c(-1, -20:-21)]) > 0,]

#xol-1
xol1 <- time_resolved_embryo[time_resolved_embryo$WormbaseName == "C18A11.5",]
xol1 <- xol1[1,]
xol1 <- as.data.frame(t(xol1))
xol1$embryo_time <- rownames(xol1)
colnames(xol1) <- c("xol1", "embryo_time")
xol1 <- xol1[-1,]
xol1 <- xol1[-19:-20,]

xol1$embryo_time <- gsub("X", "", xol1$embryo_time)
xol1$embryo_time <- as.factor(xol1$embryo_time)
xol1$xol1 <- as.numeric(xol1$xol1)
xol1$xol1 <- round(xol1$xol1, digits = 2)
xol1$embryo_time <- gsub("min", " min", xol1$embryo_time)
xol1$embryo_time <- gsub("cell", " cell", xol1$embryo_time)
xol1$embryo_time = factor(xol1$embryo_time, levels = c("4 cell", 
                                                       "44 min", 
                                                       "83 min", 
                                                       "122 min", 
                                                       "161 min", 
                                                       "199 min", 
                                                       "238 min", 
                                                       "277 min", 
                                                       "316 min", 
                                                       "355 min", 
                                                       "393 min", 
                                                       "432 min", 
                                                       "471 min", 
                                                       "510 min", 
                                                       "548 min", 
                                                       "587 min", 
                                                       "626 min", 
                                                       "665 min"))

xol1_counts <- ggplot(xol1, aes(y = xol1, x = embryo_time, group = 1)) + 
  geom_line(color = "#5065A7", alpha = 0.5) + 
  geom_point(color = "#5065A7", size = 3, alpha = 0.9) +
  theme_prism(base_line_size = 0.5) +
  labs(y =  substitute(paste(italic('xol-1  '), 'transcripts','\n')), 
       x = " Time In Embryo Development") +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1, 
                                   size = 15),
        axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18, family = "Arial Bold")) 

xol1_counts

ggsave(plot = xol1_counts, 
       filename = "xol1_transcripts.tiff", 
       path = "xol-1 paper/Data/", 
       height = 5, 
       width = 7)

#met-2 and arle-14
met2_arle14_wormbaseID <- c("arle14" = "B0336.5", 
                            "met2" = "R05D3.11")
met2_arle14 <- time_resolved_embryo[time_resolved_embryo$WormbaseName %in% met2_arle14_wormbaseID,]
rownames(met2_arle14) <- names(met2_arle14_wormbaseID)
met2_arle14 <- as.data.frame(t(met2_arle14))
met2_arle14 <- met2_arle14[-1,]
met2_arle14 <- met2_arle14[-19:-20,]
met2_arle14$embryo_time <- rownames(met2_arle14)

arle14_ggplot <- data.frame("counts" = met2_arle14$arle14,
                                 "embryo_time" = met2_arle14$embryo_time,
                                 "gene_name" = rep("arle-14", nrow(met2_arle14)))
met2_ggplot <- data.frame("counts" = met2_arle14$met2,
                               "embryo_time" = met2_arle14$embryo_time,
                               "gene_name" = rep("met-2", nrow(met2_arle14)))
met2_arle14_ggplot <- rbind(arle14_ggplot, 
                            met2_ggplot)

met2_arle14_ggplot$embryo_time <- gsub("X", "", met2_arle14_ggplot$embryo_time)
met2_arle14_ggplot$embryo_time <- gsub("min", " min", met2_arle14_ggplot$embryo_time)
met2_arle14_ggplot$embryo_time <- gsub("cell", " cell", met2_arle14_ggplot$embryo_time)
met2_arle14_ggplot$embryo_time <- as.factor(met2_arle14_ggplot$embryo_time)
met2_arle14_ggplot$embryo_time = factor(met2_arle14_ggplot$embryo_time, levels = c("4 cell", 
                                                                                   "44 min", 
                                                                                   "83 min", 
                                                                                   "122 min", 
                                                                                   "161 min", 
                                                                                   "199 min", 
                                                                                   "238 min", 
                                                                                   "277 min", 
                                                                                   "316 min", 
                                                                                   "355 min", 
                                                                                   "393 min", 
                                                                                   "432 min", 
                                                                                   "471 min", 
                                                                                   "510 min", 
                                                                                   "548 min", 
                                                                                   "587 min", 
                                                                                   "626 min", 
                                                                                   "665 min"))
met2_arle14_ggplot$counts <- as.numeric(met2_arle14_ggplot$counts)
met2_arle14_ggplot$counts <- round(met2_arle14_ggplot$counts, digits = 2)
met2_arle14_ggplot$gene_name = factor(met2_arle14_ggplot$gene_name, levels = c("met-2", 
                                                                               "arle-14"))


met2_arle14_counts <- ggplot(met2_arle14_ggplot,
                              aes(y = counts, x = embryo_time, group = gene_name, color = gene_name)) + 
  geom_line(alpha = 0.5) + 
  geom_point(size = 3, alpha = 0.9) +
  theme_prism(base_line_size = 0.5) +
  scale_color_manual(values = c("#5065A7", "#A367B1"), 
                     labels = c(substitute(italic('met-2')), substitute(italic('arle-14')))) +
  labs(y =  'transcripts', x = " \nTime In Embryo Development") +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1, 
                                   size = 15, 
                                   color = "black"), 
        axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 15)) 
  
met2_arle14_counts

ggsave(plot = met2_arle14_counts, 
       filename = "met2_arle14_transcripts.tiff", 
       path = "xol-1 paper/Data/", 
       height = 5, 
       width = 7)

#H3K9me co-factors and regulators
H3K9_cofactors_wormbaseID <- c("lin65" = "Y71G12B.9", 
                               "lin61" = "R06C7.7", 
                               "nrde3" = "R04A9.2", 
                               "cec4" = "F32E10.2",
                               "set25" = "Y43F4B.3",
                               "set32" = "C41G7.4")
H3K9_cofactors <- time_resolved_embryo[time_resolved_embryo$WormbaseName %in% H3K9_cofactors_wormbaseID,]
rownames(H3K9_cofactors) <- names(H3K9_cofactors_wormbaseID)
H3K9_cofactors <- as.data.frame(t(H3K9_cofactors))
H3K9_cofactors <- H3K9_cofactors[-1,]
H3K9_cofactors <- H3K9_cofactors[-19:-20,]
H3K9_cofactors$embryo_time <- rownames(H3K9_cofactors)

lin65_ggplot <- data.frame("counts" = H3K9_cofactors$lin65,
                            "embryo_time" = H3K9_cofactors$embryo_time,
                            "gene_name" = rep("lin-65", nrow(H3K9_cofactors)))
lin61_ggplot <- data.frame("counts" = H3K9_cofactors$lin61,
                          "embryo_time" = H3K9_cofactors$embryo_time,
                          "gene_name" = rep("lin-61", nrow(H3K9_cofactors)))
nrde3_ggplot <- data.frame("counts" = H3K9_cofactors$nrde3,
                           "embryo_time" = H3K9_cofactors$embryo_time,
                           "gene_name" = rep("nrde-3", nrow(H3K9_cofactors)))
cec4_ggplot <- data.frame("counts" = H3K9_cofactors$cec4,
                           "embryo_time" = H3K9_cofactors$embryo_time,
                           "gene_name" = rep("cec-4", nrow(H3K9_cofactors)))
set25_ggplot <- data.frame("counts" = H3K9_cofactors$set25,
                           "embryo_time" = H3K9_cofactors$embryo_time,
                           "gene_name" = rep("set-25", nrow(H3K9_cofactors)))
set32_ggplot <- data.frame("counts" = H3K9_cofactors$set32,
                           "embryo_time" = H3K9_cofactors$embryo_time,
                           "gene_name" = rep("set-32", nrow(H3K9_cofactors)))
H3K9_cofactors_ggplot <- rbind(set25_ggplot,
                               set32_ggplot,
                               lin65_ggplot,
                               lin61_ggplot,
                               nrde3_ggplot,
                               cec4_ggplot)

H3K9_cofactors_ggplot$embryo_time <- gsub("X", "", H3K9_cofactors_ggplot$embryo_time)
H3K9_cofactors_ggplot$embryo_time <- gsub("min", " min", H3K9_cofactors_ggplot$embryo_time)
H3K9_cofactors_ggplot$embryo_time <- gsub("cell", " cell", H3K9_cofactors_ggplot$embryo_time)
H3K9_cofactors_ggplot$embryo_time <- as.factor(H3K9_cofactors_ggplot$embryo_time)
H3K9_cofactors_ggplot$embryo_time = factor(H3K9_cofactors_ggplot$embryo_time, levels = c("4 cell", 
                                                                                         "44 min", 
                                                                                         "83 min", 
                                                                                         "122 min", 
                                                                                         "161 min", 
                                                                                         "199 min", 
                                                                                         "238 min", 
                                                                                         "277 min", 
                                                                                         "316 min", 
                                                                                         "355 min", 
                                                                                         "393 min", 
                                                                                         "432 min", 
                                                                                         "471 min", 
                                                                                         "510 min", 
                                                                                         "548 min", 
                                                                                         "587 min", 
                                                                                         "626 min", 
                                                                                         "665 min"))
H3K9_cofactors_ggplot$gene_name = factor(H3K9_cofactors_ggplot$gene_name, levels = c("set-25", 
                                                                                     "set-32", 
                                                                                     "lin-65", 
                                                                                     "lin-61", 
                                                                                     "nrde-3", 
                                                                                     "cec-4"))
H3K9_cofactors_ggplot$counts <- as.numeric(H3K9_cofactors_ggplot$counts)
H3K9_cofactors_ggplot$counts <- round(H3K9_cofactors_ggplot$counts, digits = 2)

H3K9_cofactors_counts <- ggplot(H3K9_cofactors_ggplot, 
                                 aes(y = counts, x = embryo_time, group = gene_name, color = gene_name)) + 
  geom_line(alpha = 0.5) + 
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = c("#5065A7", "#A367B1", "#092047","#63345e", "#b7c1de","#3DB23D"), 
                     labels = c(substitute(italic('set-25')),
                                substitute(italic('set-32')),
                                substitute(italic('lin-65')),
                                substitute(italic('lin-61')),
                                substitute(italic('nrde-3')),
                                substitute(italic('cec-4')))) +
  labs(y =  'Transcripts\n', x = " \nTime In Embryo Development") +
  theme_prism(base_line_size = 0.5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 15, color = "black"), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), legend.title = element_blank(), legend.text = element_text(size = 15)) 

H3K9_cofactors_counts + theme(axis.title.x = element_text(family = "Arial Bold"))

ggsave(plot = H3K9_cofactors_counts, 
       filename = "H3K9_cofactors_transcripts.tiff", 
       path = "xol-1 paper/Data/", 
       height = 5, 
       width = 7)


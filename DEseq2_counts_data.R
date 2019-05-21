##############################
## R code to generate heatmap for RNAseq and RPPA data 
## using RSEM or HTseq counts for RNAseq data. 
## RNAseq HTSeq counts were normalized using vst normalization in
## DESeq2 package.
## Author: Sammed Mandape
## Contact: smandape@email.arizona.edu
## Date: Apr 16, 2019
##############################

library("DESeq2")
library(pheatmap)

directory<- "C:/Users/smandape/Box/Kraft_Sathish_Sammed/heatmap"

outputPrefix <- "TLL_DESeq2"

sampleFiles <- c("Azdr1-1.counts",
                 "Azdr1-2.counts",
                 "Azdr1-3.counts",
                 "Hsb2Naive-1.counts",
                 "HSB2Naive-2.counts",
                 "HSB2Naive-3.counts")

sampleNames <- c("AZDR1-1","AZDR1-2","AZDR1-3","HSB2Naive-1","HSB2Naive-2","HSB2Naive-3")
sampleCondition <- c("Resistant","Resistant","Resistant","Naive","Naive","Naive")
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)
treatments = c("Resistant","Naive")
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design = ~ condition)
dds <- DESeq(ddsHTSeq)
res <- results(dds)
res= subset(res, padj<0.05)
res <- res[order(res$padj),]

# save data results and normalized reads to csv
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata)[1] <- 'gene'
write.csv(resdata, file = paste0(outputPrefix, "-results-with-normalized.csv"))

# send normalized counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds),normalized=T), file = paste0(outputPrefix, "_normalized_counts.txt"), sep = '\t')

# produce DataFrame of results of statistical tests
mcols(res, use.names = T)
write.csv(as.data.frame(mcols(res, use.name = T)),file = paste0(outputPrefix, "-test-conditions.csv"))

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# genes with padj < 0.1 are colored Red
plotMA(dds, ylim=c(-8,8),main = "RNAseq experiment")
dev.copy(png, paste0(outputPrefix, "-MAplot_initial_analysis.png"))
dev.off()

# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
rld <- rlogTransformation(dds, blind=T)
vsd <- varianceStabilizingTransformation(dds, blind=T)

# save normalized values
write.table(as.data.frame(assay(rld)),file = paste0(outputPrefix, "-rlog-transformed-counts.txt"), sep = '\t')
write.table(as.data.frame(assay(vsd)),file = paste0(outputPrefix, "-vst-transformed-counts.txt"), sep = '\t')

# plot to show effect of transformation
# axis is square root of variance over the mean for all samples
par(mai = ifelse(1:4 <= 2, par('mai'),0))
px <- counts(dds)[,1] / sizeFactors(dds)[1]
ord <- order(px)
ord <- ord[px[ord] < 150]
ord <- ord[seq(1,length(ord),length=50)]
last <- ord[length(ord)]
vstcol <- c('blue','black')
matplot(px[ord], cbind(assay(vsd)[,1], log2(px))[ord, ],type='l', lty = 1, col=vstcol, xlab = 'n', ylab = 'f(n)')
legend('bottomright',legend=c(expression('variance stabilizing transformation'), expression(log[2](n/s[1]))), fill=vstcol)
dev.copy(png,paste0(outputPrefix, "-variance_stabilizing.png"))
dev.off()
            

# clustering analysis
# excerpts from http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/
library("RColorBrewer")
library("gplots")
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(condition, sampleNames, sep=" : "))
#Or if you want conditions use:
#rownames(mat) <- colnames(mat) <- with(colData(dds),condition)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
dev.copy(png, paste0(outputPrefix, "-clustering.png"))
heatmap.2(mat, trace = "none", col = rev(hmcol), margin = c(13,13))
dev.off()


#Principal components plot shows additional but rough clustering of samples
library("genefilter")
library("ggplot2")
library("grDevices")

rv <- rowVars(assay(rld))
select <- order(rv, decreasing=T)[seq_len(min(500,length(rv)))]
pc <- prcomp(t(assay(vsd)[select,]))

# set condition
condition <- treatments
scores <- data.frame(pc$x, condition)

(pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(condition))))
  + geom_point(size = 5)
  + ggtitle("Principal Components")
  + scale_colour_brewer(name = " ", palette = "Set1")
  + theme(
    plot.title = element_text(face = 'bold'),
    legend.position = c(.9,.2),
    legend.key = element_rect(fill = 'NA'),
    legend.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(colour = "Black"),
    axis.text.x = element_text(colour = "Black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = 'bold'),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(color = 'black',fill = NA)
  ))

ggsave(pcaplot,file=paste0(outputPrefix, "-ggplot2.pdf"))


# scatter plot of rlog transformations between Sample conditions
# nice way to compare control and experimental samples
head(assay(rld))
# plot(log2(1+counts(dds,normalized=T)[,1:2]),col='black',pch=20,cex=0.3, main='Log2 transformed')
plot(assay(rld)[,1:3],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,2:4],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,6:5],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")

# heatmap of data
library("RColorBrewer")
library("gplots")
# 1000 top expressed genes with heatmap.2
select <- order(rowMeans(counts(dds,normalized=T)),decreasing=T)[1:1000]
my_palette <- colorRampPalette(c("blue",'white','red'))(n=1000)
heatmap.2(assay(vsd)[select,], col=my_palette,
          scale="row", key=T, keysize=1, symkey=T,
          density.info="none", trace="none",
          cexCol=0.6, labRow=F,
          main="1000 Top Expressed Genes Heatmap")

heatmap.2(assay(vsd)[resdata$gene,], col=my_palette,
          scale="row", key=T, keysize=1, symkey=T,
          density.info="none", trace="none",
          cexCol=0.6, labRow=F,
          main="Differentially expressed genes (FDR < 0.05) Heatmap")

# write vsd transformed file with genes FDR<0.05 from resdata object
write.table(as.data.frame(assay(vsd)[resdata$gene,]), file = paste0(outputPrefix, "-vst-transformed-counts-FDR-0.05.txt"), sep = "\t" )

pheatmap(assay(vsd)[resdata$gene,], scale = "row", show_rownames = F, fontsize = 11, main = "Differentially expressed genes (FDR < 0.05) Heatmap", angle_col = 45,
         fontsize_col = 10)


##################################
## The following code is to plot heatmap of those genes
## which intersect with rppa data
##################################

# converting ensemble gene id to gene symbol. Here, resdata is considered that is FDR < 0.05.
# With resdata the final heatmap is going to be for DE (FDR<0.05) genes common in RPPA.

library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ensembl2gene<-getBM(filters = 'ensembl_gene_id', attributes = c('ensembl_gene_id','hgnc_symbol'),values = resdata$gene, mart=mart)
nrow(ensembl2gene) # 3300

# to see which rows are empty (ensemble id not matched to any gene symbol)
head(ensembl2gene[which(ensembl2gene$hgnc_symbol == ""),])
nrow(ensembl2gene[which(ensembl2gene$hgnc_symbol == ""),]) #81


# remove rows where hgnc_genesymbol is empty
ensembl2gene_nonempty <- ensembl2gene[!(ensembl2gene$hgnc_symbol == ""),]
nrow(ensembl2gene_nonempty) #3219

#write.table(ensembl2gene_nonempty, file = "TLL_DESeq2_FDR-0.05_ensembl2gene.csv", sep = ",", row.names = F)
# read rppa data file
rppagene <- read.delim("RPPA_gene_names_James.txt", header = F)

#ensembl2gene_nonempty_exp <- merge(ensembl2gene_nonempty, resdata[,-c(2:7)], by.x = "ensembl_gene_id", by.y = "gene")
#rm(ensembl2gene_nonempty_exp)

# merge rppa data with rnaseq data to get intersection of those two #45
rnaexp_rppa <- merge(ensembl2gene_nonempty, rppagene, by.x = "hgnc_symbol", by.y = "V1" )

# pull expression values for the intersected data
rnaexp_rppa_val <- cbind(rnaexp_rppa$hgnc_symbol, data.frame(assay(vsd)[rnaexp_rppa$ensembl_gene_id,], row.names = NULL))

# remove gene symbol column and add it back as rownames replacing ensemblids as rownames
rnaexp_rppa_val1<-as.matrix(rnaexp_rppa_val[,-c(1)])
rownames(rnaexp_rppa_val1) <-  rnaexp_rppa_val$`rnaexp_rppa$hgnc_symbol`

pheatmap(assay(vsd)[rnaexp_rppa$ensembl_gene_id,], scale = "row", show_rownames = T, fontsize = 11, main = "Differentially expressed genes (FDR < 0.05) common to RPPA data", angle_col = 45,
         fontsize_col = 10, fontsize_row = 10)

pheatmap(rnaexp_rppa_val1, scale = "row", show_rownames = T, fontsize = 11, angle_col = 45,
         fontsize_col = 10, fontsize_row = 10)

##################################

gsea_path <- read.delim("GSEA_pathways_of_interest.txt", header = T, stringsAsFactors = F)
library(biomaRt)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
gene2ensemb <- getBM(filters = 'hgnc_symbol', attributes = c('ensembl_gene_id','hgnc_symbol'),values = gsea_path$Genes, mart=mart)

assay(vsd)[gene2ensemb$ensembl_gene_id,]
write.table(gene2ensemb, file="Genesymbol_2_Ensembl_GSEA_Pathways.txt", row.names = F, sep = '\t')

grep('ENSG00000000003',assay(vsd),ignore.case = F)
head(assay(vsd))

fgsea_input <- read.delim("fgsea_AZDR1_rnaseq_input_file_heatmap.txt", header = T, row.names = "hgnc_symbol", sep = "\t", stringsAsFactors = F)
fgsea_input_sub <- fgsea_input[,-c(7)]

path_anno <- data.frame(fgsea_input[,7])
rownames(path_anno) <- rownames(fgsea_input)
colnames(path_anno) <- "Pathways"
pheatmap(fgsea_input_sub, scale = "row", cluster_rows = F, annotation_row = path_anno, angle_col = 90, fontsize = 11)


##################################
## Heatmap of pathways of interest picked by biologists
##################################

merged_uniq_genes <- read.delim("PathwaysOfInterest_uniq_geneList.txt", header = T)
library(biomaRt)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
gene2ensemb_mergeList <- getBM(filters = 'hgnc_symbol', attributes = c('ensembl_gene_id','hgnc_symbol'),values = merged_uniq_genes$Merged_Uniq_list, mart=mart)
#write.table(gene2ensemb_mergeList, file = "PathwaysOfInterest_uniq_geneList2Ensembl.csv", sep = ",", row.names = F)

genesFDR05 <- as.data.frame(assay(vsd)[resdata$gene,])

genesFDR05$ensemblId <- rownames(genesFDR05)
merged_uniq_genes_expval <- merge(gene2ensemb_mergeList, genesFDR05, by.x = "ensembl_gene_id", by.y= "ensemblId")
write.table(merged_uniq_genes_expval, file= "PathwaysOfInterest_uniq_genes_expval.csv", sep = ",", row.names = F)


nfkb <- read.delim("Biocarta_NFKB.txt", sep = '\t', header = T, row.names = 1)
pheatmap(nfkb, scale = "row", angle_col = 90)

hall_Pi3kakt <- read.delim("Hallmark_PI3K_AKT.txt", sep = '\t', header = T, row.names = 1)
pheatmap(hall_Pi3kakt, scale = "row", angle_col = 90)

jak_stat <- read.delim("KEGG_JAK_STAT.txt", sep = '\t', header = T, row.names = 1)
pheatmap(jak_stat, scale = "row", angle_col = 90)

notch <- read.delim("KEGG_Notch.txt", sep = '\t', header = T, row.names = 1)
pheatmap(notch, scale = "row", angle_col = 90)

tgf_beta <- read.delim("KEGG_TGF_Beta.txt", sep = '\t', header = T, row.names = 1)
pheatmap(tgf_beta, scale = "row", angle_col = 90)

pid_pi3kakt <- read.delim("PID_PI3K_AKT.txt", sep = '\t', header = T, row.names = 1)
pheatmap(pid_pi3kakt, scale = "row", angle_col = 90)








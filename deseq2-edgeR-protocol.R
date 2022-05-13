# R script to perform differential expression analysis on RNA-seq data using 
# Rsubread/featureCounts and DESeq2
# Author: Harry Pollitt
# Email: hap39@aber.ac.uk / hp508@exeter.ac.uk



# load required packages ------------------------------------------------------

library(Rsubread)
library(DESeq2)
library(ggplot2)
library(apeglm)



# Create table of metadata ----------------------------------------------------

# grab file names from working directory
sample_names <- list.files(pattern= "*.bam")

# get sample IDs from bam filenames
samples <- strsplit(sample_names, "_")
samples <- sapply(samples, '[', 2)

# factors - order needs to match sample_names order
condition <- factor(c("Ctrl","Ctrl","Nys","Nys","Ctrl","Ctrl","Nys","Nys","Ctrl","Nys"))
batch <- factor(c("8","8","8","8","5","5","5","5","6","6" ))

# Combine into table
metadata <- data.frame(row.names = samples, sample_names, condition, batch)



# Run featureCounts from Rsubread ---------------------------------------------

counts <- featureCounts(files=metadata$sample_names,  annot.ext = "El_Paco_V3_gene_annotations.gtf",
                        isGTFAnnotationFile = TRUE, GTF.featureType = "exon", GTF.attrType = "transcript_id",
                        isPairedEnd = TRUE, minFragLength = 50, maxFragLength = 600, strandSpecific = 2)
						
# Ready table for use by DESeq2 and remove features with no reads
raw_counts_df <- as.data.frame(counts$counts)
counts_df <- raw_counts_df[which(rowSums(raw_counts_df) > 0),]



# Running DESeq2 --------------------------------------------------------------

# create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData=counts_df, colData = metadata, design = ~ batch + condition)

# run DESeq2 
dds <- DESeq(dds) 

# obtaining results
resultsNames(dds)  # will tell you coef or contrast names if needed
res <- results(dds, contrast=c("condition", "Nys", "Ctrl")



# Exploring Quality of DESeq2 results -----------------------------------------

# Log fold change shrinkage
resLFC <- lfcShrink(dds, coef="condition_Nys_vs_Ctrl", type="apeglm")  # if apeglm doesn't work use "normal" or "ashr"

# View log fold change shrinkage with MAplot 
plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))

# Variance stabilisating transformation VST
vsd <- vst(dds, blind=FALSE)

# View vsd as PCA
pca_norm <- plotPCA(vsd, intgroup="condition")
pca_norm + geom_label(aes(label = metadata$sample)) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
   panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(expand = c(0.035, 0))
   
# Dispersion plot
plotDispEsts(dds)



# Exploring significantly differentially expressed genes ----------------------

# Only keep genes with padj below threshold
sigs <- na.omit(res)
deseq_sigs <- sigs[sigs$padj < 0.1]
deseq_ordered_sigs <- sigs[order(sigs$padj),]
   
# EdgeR -----------------------------------------------------------------------

d <- DGEList(counts=counts_df, group=factor(condition))
design <- model.matrix(~batch + condition)
rownames(design) <- colnames(counts_df)

d <- calcNormFactors(d)

# Estimate dispersion
d <- estimateGLMCommonDisp(d, design, verbose=TRUE)
d <- estimateGLMTrendedDisp(d, design)
d <- estimateGLMTagwiseDisp(d, design)

# Explore results
plotBCV(d)
plotMDS(d, gene.selection="common")
plotMDS(d, col=rep(1:2, each=4))

# Fit to glm 
fit <- glmFit(d, design)
lrt <- glmLRT(fit, coef=2:3)
topTags(lrt)

# Adjust p value for multiple testing
FDR <- p.adjust(lrt$table$PValue, method="BH")
sum(FDR < 0.05)

# Look at results
lrt <- glmLRT(fit)
topTags(lrt)

top <- rownames(topTags(lrt))
cpm(d)[top,]

summary(dt <- decideTestsDGE(lrt))

isDE <- as.logical(dt)
DEnames <- rownames(d)[isDE]

plotSmear(lrt, de.tags=DEnames)
abline(h=c(-1,1), col="blue")

edgeR-deg <- topTags(lrt, n = Inf, p = 0.05)$table
up <- row.names(deg[deg$logFC > 0,])
down <- row.names(deg[deg$logFC < 0,])

# look for matches in gene lists from DESeq2 and EdgeR ------------------------
matches <- deseq_ordered_sigs[match(rownames(edgeR-deg), deseq_ordered_sigs)]
matches <- na.omit(matches)

de <- as.data.frame(DEnames)




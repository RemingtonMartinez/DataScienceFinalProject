
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(EnsDb.Hsapiens.v86)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)

samplefiles <- list.files("/scratch/rmm2gq/RNASeqProject/samtoolsout/macs2", pattern= ".bed", full.names=T)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- c("DiseasedVsControl")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, 
                       tssRegion=c(-1000, 1000), verbose=TRUE)

plotAnnoBar(peakAnnoList)

plotDistToTSS(peakAnnoList, title="Distribution of transcription factor-binding loci \n relative to TSS")

diseased_annot <- data.frame(peakAnnoList[["DiseasedVsControl"]]@anno)

# Get the entrez IDs
entrez <- diseased_annot$geneId

# Return the gene symbol for the set of Entrez IDs
annotations_edb <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                         keys = entrez,
                                         columns = c("GENENAME"),
                                         keytype = "ENTREZID")

# Change IDs to character type to merge
annotations_edb$ENTREZID <- as.character(annotations_edb$ENTREZID)

# Write to file
diseased_annot %>% 
  left_join(annotations_edb, by=c("geneId"="ENTREZID")) %>% 
  write.table(file="results/peak_annotation.txt", sep="\t", quote=F, row.names=F)

# Run GO enrichment analysis 
ego <- enrichGO(gene = entrez, 
                    keyType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db, 
                    ont = "BP", 
                    pAdjustMethod = "BH", 
                    qvalueCutoff = 0.05, 
                    readable = TRUE)

# Output results from GO analysis to a table
cluster_summary <- data.frame(ego)
write.csv(cluster_summary, "results/clusterProfiler.csv")


# Dotplot visualization
dotplot(ego, showCategory=50)

ekegg <- enrichKEGG(gene = entrez,
                 organism = 'hsa',
                 pvalueCutoff = 0.05)

dotplot(ekegg)

# Create a list with genes from each sample
genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)

# Run KEGG analysis
compKEGG <- compareCluster(geneCluster = genes, 
                         fun = "enrichKEGG",
                         organism = "human",
                         pvalueCutoff  = 0.05, 
                         pAdjustMethod = "BH")
dotplot(compKEGG, showCategory = 20, title = "KEGG Pathway Enrichment Analysis")


               
               
# failed attempts at using motif enrichment through GADEM and memes was incompatible for the version of R that I was implementing               
```{r}
library(rGADEM)
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
BedFile<- "/scratch/rmm2gq/RNASeqProject/samtoolsout/macs2/pleasesendhelp_summits.bed"
Sequences<-import(BedFile)
gadem<-GADEM(Sequences,verbose=1,genome=Hsapiens)
```

```{r}
library(rGADEM)
library(ChIPpeakAnno)
library(BSgenome.Hsapiens.UCSC.hg38)
bed <- "/scratch/rmm2gq/RNASeqProject/samtoolsout/macs2/pleasesendhelp_summits.bed"
gr1 <- toGRanges(bed, format="BED", header=FALSE) 
peaks <- GRangesList(rep1=gr1)
gadem<-GADEM(gr1,verbose=1,genome=Hsapiens)
```


library(RColorBrewer)
library(openxlsx)
library(edgeR)


###########
### DEG ###
###########

###VOLCANO###

VolcanoFig<-function(data,file,main){
  library (ggrepel)
  data<-data
  data$Gene<-rownames(data)
  data$threshold<-"nDEG"
  data$threshold[data$P.Value <= 0.05 & (data$logFC >= 0.58 | data$logFC <= -0.58)] <- "DEG"
  data[data$P.Value==0,"P.Val"]<-1e-318
  
  
  ##Top 10 with the best FC
  real_DE<-data[data$P.Value<=0.05,]
  selected_FC<-head(real_DE[order(abs(real_DE$logFC),decreasing = T),], 10)
  pdf(file,paper="a4")
  g = ggplot(data=data, aes(x=logFC, y=-log10(P.Value), colour=threshold)) + coord_cartesian(xlim = c(-10, 10)) +
    geom_point(alpha=0.7, size=1.75) +
    xlab("log2 fold change") + ylab("-log10 P.Value") +
    geom_text_repel(data=selected_FC, aes(label=Gene),colour="black",size=3, max.overlaps = 30) + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5)) +#adding text for the top1 10 genes
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +  # P.Value threshold (0.05)
    geom_vline(xintercept = 0.58, linetype = "dashed", color = "black", size = 0.5) +  # logFC positive threshold
    geom_vline(xintercept = -0.58, linetype = "dashed", color = "black", size = 0.5) +  # logFC negative threshold
    
    scale_color_manual(values = c("nDEG" = "lightgray", "DEG" = "skyblue"))
  
  plot(g)
  dev.off()
}

                  ##########
                  ###DATA###
                  ##########

# Comparisons #
label<-("Samples")

#reading counts
rawdata<-read.delim(file="str-ReadCount.tab",header=TRUE,row.names=1,check.names = F)
colnames(rawdata)<-gsub("S\\d\\d.*","",colnames(rawdata))
colnames(rawdata)<-gsub("S\\d.*","",colnames(rawdata))
colnames(rawdata)<-gsub("-","_",colnames(rawdata))
colnames(rawdata)<-gsub("_nat","",colnames(rawdata))
#Delete sample c48
rawdata <- subset(rawdata, select = -C48)
rawdata <- subset(rawdata, select = -A11)
rawdata <- subset(rawdata, select = -A12)
rawdata <- subset(rawdata, select = -A14)
rawdata <- subset(rawdata, select = -A15)
rawdata <- subset(rawdata, select = -A17)
rawdata <- subset(rawdata, select = -A18)
rawdata <- subset(rawdata, select = -A20)
rawdata <- subset(rawdata, select = -A21)
rawdata <- subset(rawdata, select = -A23)
rawdata <- subset(rawdata, select = -A25)
rawdata <- subset(rawdata, select = -A26)
rawdata <- subset(rawdata, select = -A28)
rawdata <- subset(rawdata, select = -A3)
rawdata <- subset(rawdata, select = -A31)
rawdata <- subset(rawdata, select = -A32)
rawdata <- subset(rawdata, select = -A33)
rawdata <- subset(rawdata, select = -A34)
rawdata <- subset(rawdata, select = -A38)
rawdata <- subset(rawdata, select = -A39)
rawdata <- subset(rawdata, select = -A4)
rawdata <- subset(rawdata, select = -A40)
rawdata <- subset(rawdata, select = -A41)
rawdata <- subset(rawdata, select = -A42)
rawdata <- subset(rawdata, select = -A43)
rawdata <- subset(rawdata, select = -A44)
rawdata <- subset(rawdata, select = -A47)
rawdata <- subset(rawdata, select = -A49)
rawdata <- subset(rawdata, select = -A50)
rawdata <- subset(rawdata, select = -A56)
rawdata <- subset(rawdata, select = -A59)
rawdata <- subset(rawdata, select = -A60)
rawdata <- subset(rawdata, select = -A61)
rawdata <- subset(rawdata, select = -A9)
rawdata <- subset(rawdata, select = -C100)
rawdata <- subset(rawdata, select = -C105)
rawdata <- subset(rawdata, select = -C13)
rawdata <- subset(rawdata, select = -C19)
rawdata <- subset(rawdata, select = -C2)
rawdata <- subset(rawdata, select = -C20)
rawdata <- subset(rawdata, select = -C31)
rawdata <- subset(rawdata, select = -C32)
rawdata <- subset(rawdata, select = -C39)
rawdata <- subset(rawdata, select = -C41)
rawdata <- subset(rawdata, select = -C42)
rawdata <- subset(rawdata, select = -C44)
rawdata <- subset(rawdata, select = -C49)
rawdata <- subset(rawdata, select = -C5)
rawdata <- subset(rawdata, select = -C52)
rawdata <- subset(rawdata, select = -C54)
rawdata <- subset(rawdata, select = -C56)
rawdata <- subset(rawdata, select = -C6)
rawdata <- subset(rawdata, select = -C65)
rawdata <- subset(rawdata, select = -C67)
rawdata <- subset(rawdata, select = -C68)
rawdata <- subset(rawdata, select = -C72)
rawdata <- subset(rawdata, select = -C79)
rawdata <- subset(rawdata, select = -C84)
rawdata <- subset(rawdata, select = -C90)
rawdata <- subset(rawdata, select = -C97)
rawdata<-rawdata[,sort(colnames(rawdata))]


#reading pheno
targets <- read.table(file = "targets_comp.txt", header = T, stringsAsFactors = F)
filenames_to_remove <- c("A11", "A12", "A14", "A15", "A17", "A18", "A20", "A21", "A23", "A25", "A26", "A28", "A3", 
                         "A31", "A32", "A33", "A34", "A38", "A39", "A4", "A40", "A41", "A42", "A43", "A44", "A47",
                         "A49", "A50", "A56", "A59", "A60", "A61", "A9", "C100", "C105", "C13", "C19", "C2", "C20", 
                         "C31", "C32", "C39", "C41", "C42", "C44", "C49", "C5", "C52", "C54", "C56", "C6", "C65", 
                         "C67", "C68", "C72", "C79", "C84", "C90", "C97")
targets <- targets[!(targets$Filename %in% filenames_to_remove), ]
print(targets)

colnames(rawdata)
targets$Filename

targets$Filename<-as.vector(targets$Filename)
rawdata<-rawdata[,targets$Filename]

#group=targets$Type

#reading gene length
gene_lengths<-read.table(file="str-Size.tab",header=T)
idx<-match(rownames(rawdata),gene_lengths$Gene) #matches the row names of the data frame read_count to the Gene column of the data frame gene_lengths.
results_counts<-gene_lengths[idx,] #reorder gene_lengths to match the order of the read_count row names.
results_counts[is.na(results_counts$Length),"Length"]<-0 #looks for NA values in the Length column of results_counts and replaces them with 0.
nrow(results_counts)

#creating DGE object
x<-DGEList(counts=rawdata, group=targets$diagnosis_code, genes=results_counts) 
x$sample

table(targets$diagnosis_code)

#DEG for groups:
Group<-factor(targets$Group)
Weight<-factor(targets$weight_kg)
IBM<-factor(targets$IBM_code)
Origin<-factor(targets$origin)
Smoker<-factor(targets$smoker)
Alcohol<-factor(targets$alcohol)
Progressive<-factor(targets$progressive)
design<-model.matrix(~0+Group+IBM+Origin+Smoker+Alcohol+Progressive)

#colnames(design)<-levels(Group)
rownames(design)<-colnames(x)

#filtering low expressed genes
keep<-filterByExpr(x,design)
x<-x[keep,,keep.lib.size=FALSE]
nrow(x)
x <- calcNormFactors(x)

#normalizing
# Differential expresion
v <- voom(x, design, plot = TRUE)


###############
### QUALITY ###
###############


label <- "tfm_SIvsCI"
group <- targets$Group
x<-DGEList(counts=rawdata, genes=results_counts) #Creates an object of type DGEList (Differential Gene Expression List) with the count data (read_count) and the gene information (results_counts).
n <- 3

cpm <- cpm(x) #Calculates the Counts Per Million (CPM) values for object x, which contains the count data and gene information.
lcpm <- cpm(x, log=TRUE) #Calculates the logarithmic CPM values for object x

L <- mean(x$samples$lib.size) * 1e-6 #Calculates the average library size of the samples in object x, and converts it to millions by dividing by 10^6 .
M <- median(x$samples$lib.size) * 1e-6 #Calculates the median library size of the samples in object x, and converts it to millions by dividing by 10^6 .
c(L, M) #Combines L and M into a vector and returns it.

summary(lcpm)
table(rowSums(x$counts==0)==8) #Create a table that counts how many genes have counts of zero in all samples (assuming there are 8 samples).

#This code block is part of a quality control (QC) process to visualise and evaluate gene expression data.
pdf(paste(label,"_QC.pdf",sep=""),paper="A4")
#tiff(paste(label,"_QC.tiff",sep=""),width = 62, height = 42, units = "cm", res = 300)
lcpm.cutoff <- log2(10/M + 2/L) #Calculates a cutoff point (cutoff) at the logarithmic values of CPM.
nsamples <- ncol(x) #Calculates the number of samples (nsamples) in object x by counting the number of columns in x.
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(rawdata[,1]), col=col[1], lwd=2, ylim=c(0,0.4), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright",legend=targets$diagnosis_code, text.col=col, bty = "n", cex = 0.5)

lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.4), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", legend=targets$diagnosis_code, text.col=col, bty="n", cex = 0.5)

x <- calcNormFactors(x) #Calculates normalisation factors for object x using the calcNormFactors function.
x$samples
x <- estimateCommonDisp(x, robust=TRUE) #Estimates the common dispersion for object x using the estimateCommonDisp function. The argument robust=TRUE indicates that a robust method should be used for estimation, which can help reduce the impact of outliers.
x <- estimateTagwiseDisp(x)

x2 <- x #Creates a copy of object x and stores it in x2
x2$samples$norm.factors <- 1 #Adjusts the normalisation factors of all samples in x2 to 1

#Box plots of log CPM values before and after normalisation are created, 
#allowing visual comparison of the distribution of gene expression data.
col.group <- as.factor(group)
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)

par(mfrow=c(1,2))
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col.group, main="", names=targets$diagnosis_code, cex.axis=0.4)
title(main="Unnormalised data",ylab="Log-cpm")

lcpm <- cpm(x, log=TRUE)
boxplot(lcpm, las=2, col=col.group, main="", names=targets$diagnosis_code, cex.axis=0.4)
title(main="Normalised data",ylab="Log-cpm")

#A multidimensional scaling (MDS) plot is created that shows how samples are grouped according to their gene expression profiles. 
#This helps to identify patterns in data and to verify the consistency of groups of samples.
lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,1))
col.group <- as.factor(group)
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
levels(col.group)
plotMDS(lcpm, labels=targets$diagnosis_code, col=col.group)
title(main="A. Sample groups")

lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,1))
col.group <- as.factor(group)
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
levels(col.group)
plotMDS(lcpm, labels=targets$diagnosis_code, col=col.group)
title(main="A. Sample groups")

lcpm <- cpm(x, log= FALSE)
par(mfrow=c(1,1))
col.group <- as.factor(group)
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
levels(col.group)
plotMDS(lcpm, labels=targets$diagnosis_code, col=col.group)
title(main="A. Sample groups")

library(ggfortify)

#PCA by type
data_pca<-as.matrix(x) #Converts object x to an array.
data_pca<-as.data.frame(t(data_pca)) #Transposes the matrix (swapping rows and columns) and converts it into a data frame.
rownames(data_pca)<-targets$diagnosis_code
data_pca.PC = prcomp(data_pca) #Perform a Principal Component Analysis (PCA) on gene expression data using the prcomp function.
data_pca$Type<-targets$diagnosis_code #Adds a new Type column to the data frame data_pca, which contains the diagnosis codes (diagnosis$diagnosis_code) for each sample.

autoplot(data_pca.PC,label=T,data=data_pca,colour='Type')

library(sva)
library(ggrepel)
#visualisation of a heat map for the most variable genes
rsd <- rowSds(as.matrix(x))
sel <- order(rsd, decreasing=TRUE)[1:250]
samplenames<-as.character(targets$diagnosis_code)
heatmap(na.omit(as.matrix(x[sel,])),margins=c(10,8),main="Heatmap 250 most DE entities",cexRow=0.01,cexCol=0.5,labCol=samplenames)

#CLUSTER
#performs a hierarchical clustering analysis on the gene expression data 
#to compare how samples are grouped before and after normalisation.
library(cluster)
library(factoextra)
library(ggplot2)


par(mfrow=c(1,1), col.main="royalblue4", col.lab="royalblue4", col.axis="royalblue4", bg="white", fg="royalblue4", font=1, cex.axis=0.6, cex.main=0.8)
pr.hc.c<- hclust(na.omit(dist(t(rawdata))))
plot(pr.hc.c, xlab="Sample Distance",main=paste("Hierarchical Clustering of ", label, sep=""), labels=targets$diagnosis_code)
pr.hc.c<- hclust(na.omit(dist(t(cpm(x$counts,log=T)),method = "euclidean")))
plot(pr.hc.c, xlab="Sample Distance",main=paste("Hierarchical Clustering of Normalized samples of ", label, sep=""), labels=targets$diagnosis_code)


#tSNE
#performs a t-SNE analysis to reduce the dimensionality of gene expression data and visualise 
#the relationships between samples in a two-dimensional space.
library(M3C)

a<-tsne(x$counts,seed=100,labels=as.factor(targets$diagnosis_code), perplex=1, legendtitle="Types",text=targets$diagnosis_code ,dotsize=3, legendtextsize = 6) + ggtitle("Tsne") + theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))

plot(a)

dev.off() #Closes the current graphics device (pdf or tiff)

#######################
###DEG VISUALIZATION###
#######################

barplot(x$samples$lib.size,names.arg = targets$diagnosis_code,las=2)
abline(h=5.6e6, col="red")

#saving RPKM
rpkm<-rpkm(x,normalized.lib.sizes=TRUE)
colnames(rpkm)<-targets$Group
resultsfile<-file.path(paste(label,"_RPKM_FvsSI.xls", sep=""))
write.table(rpkm, file=resultsfile, sep = "\t", col.names = NA , qmethod = "double") 
resultsfile

fit <- lmFit(v, design) #
contrast<-makeContrasts(
  b=(GroupCOMBINED_INFERTILITY) - (GroupSIMPLE_INFERTILITY),
  levels=design
)

# DEG b 
contrast
tmp <- contrasts.fit(fit, contrast[,"b"])
tmp <- eBayes(tmp)
tmp
top1.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top1.table, 50)

#Up
top1.table[top1.table$Gene=="RN7SL2","logFC"]
control<-mean(as.numeric(rpkm["RN7SL2",targets[targets$Group=="SIMPLE_INFERTILITY","Group"]]), na.rm = TRUE)
sample<-mean(as.numeric(rpkm["RN7SL2",targets[targets$Group=="COMBINED_INFERTILITY","Group"]]), na.rm = TRUE)
#logFC
log2(sample/control)


write.table(top1.table,file="SIMPLE_INFERTILITY_vs_COMBINED_INFERTILITY.xls",quote = F,sep="\t",col.names = NA)

VolcanoFig(top1.table, "SIMPLE_INFERTILITY_vs_COMBINED_INFERTILITY.pdf", "SIMPLE_INFERTILITY_vs_COMBINED_INFERTILITY")

###ENRICHMENT###

format_enrichment<-function(ego,de_data){
  data<-data.frame(ego)
  colnames(data)<-gsub("Count","Total",colnames(data))
  genes_down<-de_data[de_data$P.Value<=0.05 & de_data$logFC<=-0.58,"Gene"]
  genes_up<-de_data[de_data$P.Value<=0.05 & de_data$logFC>=0.58,"Gene"]
  data$GenesID_up<-NA
  data$GenesID_dn<-NA
  data$Total_up<-NA
  data$Total_dn<-NA
  
  for(n in 1:nrow(data)){
    up<-genes_up[genes_up %in% unlist(strsplit(data[n,"geneID"],"/"))]
    dn<-genes_down[genes_down %in% unlist(strsplit(data[n,"geneID"],"/"))]
    data[n,"GenesID_up"]<-paste(up,collapse = "/")
    data[n,"GenesID_dn"]<-paste(dn,collapse = "/")
    data[n,"Total_up"]<-length(up)
    data[n,"Total_dn"]<-length(dn)
  }
  return(data[,c(1:8,10,11,9,12,13)])
}

library(clusterProfiler)

my_enrichment <- function(de_data, P.Value = 0.05, FA_label, cutoff = 0.05, organism = "human", FC = 1.5) {
  
  gene <- de_data[de_data$P.Value <= P.Value & abs(de_data$logFC) >= 1.5, "Gene"]
  
  dbName <- NULL
  keggDB <- NULL
  if (organism == "human") {
    library("org.Hs.eg.db")
    dbName <- org.Hs.eg.db
    keggDB <- "hsa"
  } else {
    stop("At the moment only human and mouse are supported")
  }
  
  # Conversion of gene symbols to ENTREZ IDs
  eg <- bitr(gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = dbName)
  all_eg <- bitr(de_data$Gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = dbName)
  
  # Perform GO enrichment for BP, MF, CC
  ego_BP <- enrichGO(gene = as.vector(eg$ENTREZID),
                     universe = as.vector(all_eg$ENTREZID),
                     OrgDb = dbName,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = cutoff,
                     qvalueCutoff = cutoff,
                     readable = TRUE)
  
  ego_MF <- enrichGO(gene = as.vector(eg$ENTREZID),
                     universe = as.vector(all_eg$ENTREZID),
                     OrgDb = dbName,
                     ont = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff = cutoff,
                     qvalueCutoff = cutoff,
                     readable = TRUE)
  
  ego_CC <- enrichGO(gene = as.vector(eg$ENTREZID),
                     universe = as.vector(all_eg$ENTREZID),
                     OrgDb = dbName,
                     ont = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff = cutoff,
                     qvalueCutoff = cutoff,
                     readable = TRUE)
  
  # Save results and generate graphs for BP, MF, and CC
  if (nrow(as.data.frame(ego_BP)) > 0) {
    write.table(format_enrichment(ego_BP, de_data),
                file = paste("Enrichment_GO_BP_", FA_label, ".xls", sep = ""),
                sep = "\t", row.names = F)
    
    pdf(file = paste("Enrichment_GO_BP_", FA_label, ".pdf", sep = ""), paper = "a4")
    barplot(ego_BP, main = "GO terms from BP")
    plot(dotplot(ego_BP) + ggtitle("GO terms from BP"))
    plotGOgraph(ego_BP, firstSigNodes = 5, sigForAll = TRUE, useFullNames = TRUE)
    dev.off()
    
    png(file = paste("GO_BP_dotplot_", FA_label, ".png", sep = ""))
    plot(dotplot(ego_BP) + ggtitle("GO terms from BP"))
    dev.off()
  }
  
  if (nrow(as.data.frame(ego_MF)) > 0) {
    write.table(format_enrichment(ego_MF, de_data),
                file = paste("Enrichment_GO_MF_", FA_label, ".xls", sep = ""),
                sep = "\t", row.names = F)
    
    pdf(file = paste("Enrichment_GO_MF_", FA_label, ".pdf", sep = ""), paper = "a4")
    barplot(ego_MF, main = "GO terms from MF")
    plot(dotplot(ego_MF) + ggtitle("GO terms from MF"))
    plotGOgraph(ego_MF, firstSigNodes = 5, sigForAll = TRUE, useFullNames = TRUE)
    dev.off()
    
    png(file = paste("GO_MF_dotplot_", FA_label, ".png", sep = ""))
    plot(dotplot(ego_MF) + ggtitle("GO terms from MF"))
    dev.off()
  }
  
  if (nrow(as.data.frame(ego_CC)) > 0) {
    write.table(format_enrichment(ego_CC, de_data),
                file = paste("Enrichment_GO_CC_", FA_label, ".xls", sep = ""),
                sep = "\t", row.names = F)
    
    pdf(file = paste("Enrichment_GO_CC_", FA_label, ".pdf", sep = ""), paper = "a4")
    barplot(ego_CC, main = "GO terms from CC")
    plot(dotplot(ego_CC) + ggtitle("GO terms from CC"))
    plotGOgraph(ego_CC, firstSigNodes = 5, sigForAll = TRUE, useFullNames = TRUE)
    dev.off()
    
    png(file = paste("GO_CC_dotplot_", FA_label, ".png", sep = ""))
    plot(dotplot(ego_CC) + ggtitle("GO terms from CC"))
    dev.off()
  }
}

my_enrichment(top1.table, P.Value=0.05, FA_label="enrichment_SIMPLE_INFERTILITY_vs_COMBINED_INFERTILITY", cutoff=1,organism="human",FC=1.5)


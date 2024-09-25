### Loading the R packages
library(deco)
library(BiocParallel) # for parallel computation
library(Homo.sapiens)
library(edgeR)


rawdata<-read.delim(file="str-ReadCount.tab",header=TRUE,row.names=1,check.names = F)

colnames(rawdata)<-gsub("S\\d\\d.*","",colnames(rawdata))
colnames(rawdata)<-gsub("S\\d.*","",colnames(rawdata))
colnames(rawdata)<-gsub("-","_",colnames(rawdata))
colnames(rawdata)<-gsub("_nat","",colnames(rawdata))
rawdata<-rawdata[,sort(colnames(rawdata))]


# Create subset for FERTILITY vs COMBINED INFERTILITY of rawdata
selected_samples_FvsCI <- c(
  "C56", "C6", "C100", "C84", "C39", "C41", "C54", "C65", "C2", "C44",
  "C31", "C105", "C19", "C32", "C68", "C97", "C72", "C13", "C42", "C20", 
  "C49", "C52", "C79", "C90", "C67", "C5", "A2", "A3", "A4", "A8", "A9",
  "A11", "A12", "A13", "A14", "A15", "A17", "A18", "A19", "A20", "A21", 
  "A23", "A25", "A27", "A26", "A28", "A29", "A30", "A31", "A32", "A33", 
  "A34", "A35", "A38", "A39", "A40", "A41", "A42", "A43", "A44", "A47", 
  "A49", "A50", "A56", "A59", "A60", "A61"
)

#reading pheno
targets<-read.table(file = "targets.txt",header = T,stringsAsFactors=F)
targets
colnames(rawdata)
targets$Filename

targets$Filename<-as.vector(targets$Filename)
rawdata<-rawdata[,targets$Filename]

#reading gene length
gene.length<-read.table(file="str-Size.tab",header=T)
idx<-match(rownames(rawdata),gene.length$Gene)
results_counts<-gene.length[idx,]
results_counts[is.na(results_counts$Length),"Length"]<-0
nrow(results_counts)

# Differential expresion

v <- limma::voom(rawdata)$E

#########################################
### FERTILITY vs COMBINED INFERTILITY ###
#########################################

v <- v[, selected_samples_FvsCI, drop=FALSE]

#reading pheno
targets_FvsCI<-read.table(file = "targets_ferti_combi.txt",header = T,stringsAsFactors=F)
targets_FvsCI
colnames(rawdata)
targets_FvsCI$Filename
Group<-factor(targets_FvsCI$Group)
design<-model.matrix(~0+Group)

gr <- GRanges(seqnames = Rle("chr1", nrow(v)), # Assuming that all genes are in "chr1".
              ranges = IRanges(start = 1:nrow(v), width = 1),  # Fictitious positions
              strand = Rle("*", nrow(v))) # No information on the chain
names(gr) <- rownames(v)

# Ensure that the names of the rows in v$E match the names of the rows in gr
rownames(v) <- rownames(gr)

# Ensure that the names of the columns of v$E match
design_df <- DataFrame(design)
rownames(design_df) <- colnames(v)

# Create a SummarisedExperiment object from voom_data
se <- SummarizedExperiment(
  assays = list(counts = v),  # E contains the standardised expression data
  rowRanges = gr, # Gene names
  colData = design_df  # Experimental design
)

# Check the SummarizedExperiment object
print(se)

## UNSUPERVISED analysis
sub.ma.3r.1K.uns <- decoRDA(data = se, q.val = 0.05,
                            rm.xy = TRUE, r = 10, annot = TRUE,
                            id.type = "SYMBOL", iterations = 1000,
                            pack.db = "Homo.sapiens")

## Group-vs-group comparison
classes.DIAGNOSIS <- targets_FvsCI$Group
names(classes.DIAGNOSIS) <- targets_FvsCI$Filename

sub_binary <- decoRDA(
  data = se, classes = classes.DIAGNOSIS, q.val = 0.05,
  iterations = 1000, rm.xy = TRUE, r = NULL,
  control = "FERTILE", annot = TRUE, bpparam = SerialParam(),
  id.type = "SYMBOL", pack.db = "Homo.sapiens"
)

deco.results.ma <- decoNSCA(sub = sub_binary, v = 80, method = "ward.D", bpparam = SerialParam(),
                            k.control = NULL, k.case = NULL, samp.perc = 0.05, rep.thr = 1)

summary(deco.results.ma)

featTable <- featureTable(deco.results.ma)

head(featTable)

library(writexl)
write_xlsx(featTable, "FERTILEvsCOMBINED_INFERTILITY_featTable.xlsx")

decoReport(deco.results.ma, sub_binary , pdf.file = "Report_FERTILITYvsCOMBINED_INFERTILITY.pdf", cex.names = 0.3, print.annot=T)
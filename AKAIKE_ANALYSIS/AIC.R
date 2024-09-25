library(readr)
library(dplyr)
library(caret)
library(edgeR)

# Reading patients data
data <- read.delim("data.txt", ,header=TRUE,row.names=1,check.names = F)

# Reading count data
rawdata<-read.delim(file="str-ReadCount.tab",header=TRUE,row.names=1,check.names = F)

colnames(rawdata)<-gsub("S\\d\\d.*","",colnames(rawdata))
colnames(rawdata)<-gsub("S\\d.*","",colnames(rawdata))
colnames(rawdata)<-gsub("-","_",colnames(rawdata))
colnames(rawdata)<-gsub("_nat","",colnames(rawdata))
rawdata<-rawdata[,sort(colnames(rawdata))]


# Create subset for FERTILITY vs SIMPLE INFERTILITY of rawdata
selected_samples <- c(
  "C56", "C6", "C100", "C84", "C39", "C66", "C45", "C8", "C41", "C54",
  "C65", "C2", "C44", "C31", "C105", "C19", "C32", "C68", "C97", "C72",
  "C13", "C42", "C20", "C49", "C4", "C52", "C79", "C90", "C67", "C5",
  "A1", "A2", "A3", "A4", "A7", "A8", "A9", "A11", "A12", "A13","A14", "A15", "A16", "A17",
  "A18", "A19","A20", "A21", "A23", "A25", "A26", "A27","A28", "A29","A30","A31", "A32", "A33",
  "A34", "A35", "A38", "A39", "A40", "A41", "A42", "A43", "A44", "A45", "A46",
  "A47", "A48", "A49", "A50", "A56", "A59", "A60", "A61"
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

#creating DGE object
x<-DGEList(counts=rawdata, group=targets$Group, genes=results_counts)
x$sample

table(targets$Group)

Group<-factor(targets$Group)
design<-model.matrix(~0+Group)

rownames(design)<-colnames(x)
#filtering low expressed genes
keep<-filterByExpr(x,design)
x<-x[keep,,keep.lib.size=FALSE]
nrow(x)
x <- calcNormFactors(x)

#normalizing

# Differential expresion
v <- voom(x, design)

norm_rawdata <- data.frame(v$E)
norm_rawdata<-norm_rawdata[,sort(colnames(norm_rawdata))]


# Combining data
combined <- rbind(data, norm_rawdata)
transposed_combined <- t(combined)
transposed_combined <- as.data.frame(transposed_combined)

###AKAIKE###

# Function to impute missing values with the mean for numeric columns
impute_mean <- function(x) {
  if (is.numeric(x)) {
    replace(x, is.na(x), mean(x, na.rm = TRUE))
  } else {
    x
  }
}

# Function to impute missing values with the mode for categorical columns
impute_mode <- function(x) {
  if (is.factor(x) || is.character(x)) {
    mode_value <- as.character(names(sort(table(x), decreasing = TRUE)[1]))
    replace(x, is.na(x), mode_value)
  } else {
    x
  }
}

# Apply imputation to numeric columns
transposed_combined <- transposed_combined %>%
  mutate(across(where(is.numeric), impute_mean))

# Convert appropriate columns to factors (categorical)
transposed_combined <- transposed_combined %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.factor), impute_mode))

# Check for missing values
cat("Number of missing values:", sum(is.na(transposed_combined)), "\n")

# Convert the classification variable to factor
transposed_combined$Group <- as.factor(transposed_combined$diagnosis)

# Create a model with all variables and perform stepwise AIC-based variable selection
full_model <- glm(Group ~ ., data = transposed_combined, family = binomial)

# Perform AIC-based backward variable selection
stepwise_model <- step(full_model, direction = "backward")

# Summary of the final model
summary(stepwise_model)

# Obtain the AIC of the final model
aic_stepwise <- AIC(stepwise_model)
print(paste("AIC:", aic_stepwise))


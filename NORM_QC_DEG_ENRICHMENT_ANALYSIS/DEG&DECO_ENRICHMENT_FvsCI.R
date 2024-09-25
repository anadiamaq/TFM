library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(readxl)


# Loading and filtering DECO data
deco_data <- readxl::read_excel("FERTILEvsCOMBINED_INFERTILITY_featTable.xlsx")
filtered_genes <- deco_data %>% filter(Profile %in% c("Majority", "Complete"))
genes_of_interest <- filtered_genes$SYMBOL
print(genes_of_interest)

# Loading DEG data
deg_data <-  readxl::read_excel("FERTILE_vs_COMBINED_INFERTILITY.xlsx")

# Filter deg_data to contain only the genes of interest
filtered_deg_data <- deg_data %>% filter(Gene %in% genes_of_interest)
print(filtered_deg_data)

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
  gene <- gene$Gene
  
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

my_enrichment(filtered_deg_data, P.Value=0.05, FA_label="enrichment_FERTILE_vs_COMBINED_INFERTILITY_DEG&DECO", cutoff=1,organism="human",FC=1.5)


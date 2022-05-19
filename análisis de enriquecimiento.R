if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
install.packages("knitr")
install.packages("colorspace")
install.packages("gplots")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("htmlTable")
install.packages("prettydoc")
install.packages("devtools")
install.packages("BiocManager")
BiocManager::install("oligo")
BiocManager::install("pd.clariom.s.human")
BiocManager::install("arrayQualityMetrics")
BiocManager::install("pvca")
# NOT NEEDED UNTIL ANALYSES ARE PERFORMED
BiocManager::install("limma")
BiocManager::install("genefilter")
BiocManager::install("clariomshumantranscriptcluster.db")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("ReactomePA")

dato<-read.csv("C:/Users/Esteban Vargas Parra/Downloads/ExpressAndTopRealData.csv", sep = ";")
gene4<-dato$X
gene4<-as.character(gene4)
str(gene4)
library(clusterProfiler)
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
str(gene)
summary(gene)
head(gene)
library(org.Hs.eg.db)

ejemplo<-bitr(gene4, fromType = "ENSEMBL",
              toType = "ENTREZID",
              OrgDb = "org.Hs.eg.db")
write.csv(ejemplo,"C:/Users/Esteban Vargas Parra/Downloads/ejemplo.csv")

ejemplo2 <- enrichGO(gene         = ejemplo$ENSEMBL,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

goplot(ejemplo2)

paste0(substring(gene4, 1, 1), tolower(substring(gene4, 2))) -> gene42
gene42[gene42=="Pvr"] = "PVR"
conversion<-bitr(gene42, fromType = 'ENSEMBL', toType = c("ENTREZID"),  OrgDb = org.Hs.eg.db)


#Tesis

gseGO()

go_gse <- gseGO(geneList=gene_list, 
                ont = input$ontology, 
                keyType = input$keytype, 
                nPerm = input$nPerm, 
                minGSSize = input$minGSSize, 
                maxGSSize = input$maxGSSize, 
                pvalueCutoff = input$pvalCuttoff, 
                verbose = T, 
                OrgDb = orgDb.obj, 
                pAdjustMethod = input$pAdjustMethod)

kegg_gse <- gseKEGG(geneList=kegg_gene_list, 
                    organism=organismsDbKegg[input$organismDb],
                    nPerm = input$nPerm,
                    minGSSize = input$minGSSize, 
                    maxGSSize = input$maxGSSize, 
                    pvalueCutoff = input$pvalCuttoff,
                    pAdjustMethod = input$pAdjustMethod,
                    keyType = input$keggKeyType)
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
gene4<-dato$Entrez


library(clusterProfiler)
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
View(gene)
str(gene)
summary(gene)
head(gene)
library(org.Hs.eg.db)
ggo <- groupGO(gene     = gene,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)
head(ggo)
ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)
gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
ego2 <- enrichGO(gene         = dato$Entrez,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
head(ego2, 3)  
ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "CC",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

ego4 <- enrichGO(gene         = dato$Gene.Symbol,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
goplot(ego2)

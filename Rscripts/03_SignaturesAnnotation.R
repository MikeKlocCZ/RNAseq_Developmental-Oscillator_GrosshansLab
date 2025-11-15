

# upload of time-resolved RNA seq on c. elegans capturing developmental clock
# 2020 Mollecular systems biology, MWM Meeuse et al

# Here, using clusterProfiler we annotate the extracted signature using NMF



suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(limma)
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(ComplexHeatmap)
  library(viridis)
  library(clusterProfiler)
  library(msigdbr)
  library(org.Ce.eg.db)
})



#Set wd to the directory of the script (in R studio)
setwd(sub("03_SignaturesAnnotation\\.R", "", getSourceEditorContext()$path) )


Sign1.genes <- read.table("tables/Signature1.csv",  header =FALSE) %>% pull
Sign2.genes <- read.table("tables/Signature2.csv",  header =FALSE) %>% pull
Sign3.genes <- read.table("tables/Signature3.csv",  header =FALSE) %>% pull
Sign4.genes <- read.table("tables/Signature4.csv",  header =FALSE) %>% pull
Sign5.genes <- read.table("tables/Signature5.csv",  header =FALSE) %>% pull
Sign6.genes <- read.table("tables/Signature6.csv",  header =FALSE) %>% pull
Sign7.genes <- read.table("tables/Signature7.csv",  header =FALSE) %>% pull
Sign8.genes <- read.table("tables/Signature8.csv",  header =FALSE) %>% pull

#SIGNATURES
#maturation Together: Sign7 : final maturation, Sign1 meand maturation , Sign 6: initial form of maturation, 
#Sign3. Transition to maturation
# exit: Signature8
# phase1: Sign5
# phase2: Sign2
# between phases: Sign.4

Sign.Phase1 <- Sign5.genes
Sign.Phase2 <- Sign2.genes
Sign.Phase1to2 <- Sign4.genes
Sign.CycleExit <- Sign8.genes
Sign.InitialMaturation <- Sign3.genes
Sign.Maturation <- c(Sign7.genes,Sign1.genes,Sign6.genes)

#database
#msigdbr_species()


#C elegans symbol
ce <- org.Ce.eg.db

dictionary.Phase1 <- select(ce, 
                            keys = Sign.Phase1,
                            columns = c("ENTREZID", "ENSEMBL"),
                            keytype = "ENSEMBL")
dictionary.Phase1 <- na.omit(dictionary.Phase1)


dictionary.Phase2 <- select(ce, 
                            keys = Sign.Phase2,
                            columns = c("ENTREZID", "ENSEMBL"),
                            keytype = "ENSEMBL")
dictionary.Phase2 <- na.omit(dictionary.Phase2)


dictionary.Phase1to2 <- select(ce, 
                               keys = Sign.Phase1to2 ,
                               columns = c("ENTREZID", "ENSEMBL"),
                               keytype = "ENSEMBL")
dictionary.Phase1to2 <- na.omit(dictionary.Phase1to2)


dictionary.CycleExit <- select(ce, 
                               keys = Sign.CycleExit ,
                               columns = c("ENTREZID", "ENSEMBL"),
                               keytype = "ENSEMBL")
dictionary.CycleExit <- na.omit(dictionary.CycleExit)

dictionary.InitialMaturation <- select(ce, 
                               keys = Sign.InitialMaturation ,
                               columns = c("ENTREZID", "ENSEMBL"),
                               keytype = "ENSEMBL")
dictionary.InitialMaturation <- na.omit(dictionary.InitialMaturation)


dictionary.Maturation <- select(ce, 
                               keys = Sign.Maturation ,
                               columns = c("ENTREZID", "ENSEMBL"),
                               keytype = "ENSEMBL")
dictionary.Maturation <- na.omit(dictionary.Maturation)



input.overrepr <- list(dictionary.Phase1$ENTREZID,
                       dictionary.Phase2$ENTREZID,
                       dictionary.Phase1to2$ENTREZID,
                       dictionary.CycleExit$ENTREZID,
                       dictionary.InitialMaturation$ENTREZID,
                       dictionary.Maturation$ENTREZID)

names(input.overrepr) <- c("Sig.Phase1","Sig.Phase2","Sig.Phase1to2","Sig.CycleExit","Sig.IniMat","Sig.LateMat")

###### Enrichment analysis
m_t2g <- msigdbr::msigdbr(species = "Caenorhabditis elegans", category = c("H")) %>% 
  dplyr::select(gs_name, entrez_gene) %>% dplyr::distinct(gs_name, entrez_gene)
xx <- compareCluster(input.overrepr, enricher, TERM2GENE=m_t2g, pvalueCutoff=0.1, pAdjustMethod="BH") #p-value cutoff 0.1
xx <- as.data.frame(xx)

write.csv(as.matrix(xx), paste0("tables/H_Enrichment_largerCutOff.csv"),row.names = T, quote = F)
head(xx, n = 3)
# Cluster                               ID
# 1    Sig.Phase1              HALLMARK_GLYCOLYSIS
# 2    Sig.Phase1 HALLMARK_CHOLESTEROL_HOMEOSTASIS
# 3 Sig.Phase1to2      HALLMARK_HEDGEHOG_SIGNALING
# Description GeneRatio  BgRatio       pvalue   p.adjust
# 1              HALLMARK_GLYCOLYSIS      6/19 126/2093 0.0005992923 0.01018797
# 2 HALLMARK_CHOLESTEROL_HOMEOSTASIS      3/19  53/2093 0.0111559539 0.09482561
# 3      HALLMARK_HEDGEHOG_SIGNALING       1/3  19/2093 0.0269999471 0.07084829
# qvalue                                    geneID Count
# 1 0.008200842 260026/173326/266833/173171/178170/183790     6
# 2 0.076330211                      260026/172715/173326     3
# 3 0.037288573                                    186880     1

#bioprocesses
m_t2g <- msigdbr::msigdbr(species = "Caenorhabditis elegans", category = "C5", subcategory = "BP") %>% 
  dplyr::select(gs_name, entrez_gene) %>% dplyr::distinct(gs_name, entrez_gene)
xx <- compareCluster(input.overrepr, enricher, TERM2GENE=m_t2g, pvalueCutoff=0.1, pAdjustMethod="BH")
xx <- as.data.frame(xx)


write.csv(as.matrix(xx), paste0("tables/GO_BP_Enrichment.csv"),row.names = T, quote = F)

#celullar components
m_t2g <- msigdbr::msigdbr(species = "Caenorhabditis elegans", category = "C5", subcategory = "CC") %>% 
  dplyr::select(gs_name, entrez_gene) %>% dplyr::distinct(gs_name, entrez_gene)
xx <- compareCluster(input.overrepr, enricher, TERM2GENE=m_t2g, pvalueCutoff=0.1, pAdjustMethod="BH")
xx <- as.data.frame(xx)

write.csv(as.matrix(xx), paste0("tables/GO_CC_Enrichment.csv"),row.names = T, quote = F)


#molecular functional ontology
m_t2g <- msigdbr::msigdbr(species = "Caenorhabditis elegans", category = "C5", subcategory = "MF") %>% 
  dplyr::select(gs_name, entrez_gene) %>% dplyr::distinct(gs_name, entrez_gene)
xx <- compareCluster(input.overrepr, enricher, TERM2GENE=m_t2g, pvalueCutoff=0.1, pAdjustMethod="BH")
xx <- as.data.frame(xx)

write.csv(as.matrix(xx), paste0("tables/GO_MF_Enrichment.csv"),row.names = T, quote = F)

#check oscillations of the signatures (to understand better how biological processes alter through the oscilaltions)
library(ButchR)
load(file = "Rdata/NMF_mydata.rda")
Hk8 <- HMatrix(procr_nmf_exp, k = 8)
rownames(Hk8) <- c("Sign.1", "Sign.2", "Sign.3", "Sign.4", "Sign.5",
                   "Sign.6", "Sign.7", "Sign.8")

Hk8 <- as.matrix(Hk8)


df.anno <- read.table("tables/annotation.csv", row.names = 1, sep = ",", header = TRUE) %>% as.data.frame

#oscillating signatures 2,3 and 5
df.plot.sigs <- as.data.frame(t(Hk8) )
df.plot.sigs <- cbind(df.plot.sigs,df.anno$timepoint)
colnames(df.plot.sigs)[length(colnames(df.plot.sigs))] <- "timepoint"

df.osc <- df.plot.sigs %>% filter(timepoint >= 8 & timepoint <=35) 
  
#2 5 4
ggplot(df.osc,aes(x=timepoint,y = Sign.4)) + geom_point()

#regression with 8 hours period
omega = 2*pi/7.5
Time = df.osc$timepoint
xc <- cos(omega*Time)
xs <- - sin(omega*Time)




pdf("plots/cyclingSignatures.pdf", width = 17, height = 5.5)
for (i in c(2,4,5,8) ){
  fit.lm <- lm(df.osc[,i] ~  xc+xs )
  pred <- stats::predict(fit.lm, newdata=data.frame(Time=Time))   
  plot(df.plot.sigs[,i] ~ df.plot.sigs[,"timepoint"],  main=paste("row id", i, colnames(df.osc)[i]))
  lines( pred  ~ Time, col="blue")
}
dev.off()

pdf("plots/TimeDepSignatures.pdf", width = 17, height = 5.5)
for (i in seq(8) ){
  plot(df.plot.sigs[,i] ~ df.plot.sigs[,"timepoint"],  main=paste("row id", i, colnames(df.osc)[i]))
}
dev.off()

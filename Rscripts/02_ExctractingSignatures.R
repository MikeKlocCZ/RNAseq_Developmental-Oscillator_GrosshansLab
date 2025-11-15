

# upload of time-resolved RNA seq on c. elegans capturing developmental clock
# 2020 Mollecular systems biology, MWM Meeuse et al

# Here, we extract molecular signature using NMF

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(limma)
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(ComplexHeatmap)
  library(viridis)
})

#Set wd to the directory of the script (in R studio)
setwd(sub("02_ExctractingSignatures\\.R", "", getSourceEditorContext()$path) )

#read-in data, phase-ordered expression
dge <- readRDS("Rdata/ordered_expression.rds")


#observation about the developmental stages
# time 0 - 5 : initial phase
#time 35 and more : final phase

dge$samples$stages <- NA
dge$samples$stages[dge$samples$timepoint <= 5] <- "initial"
dge$samples$stages[dge$samples$timepoint %in% c(6,7,8,9,13,14,15,19,20,21,22,27,28,29,30)] <- "phase_1"
dge$samples$stages[dge$samples$timepoint %in% c(10,11,12,16,17,18,23,24,25,26,31,32,33,34)] <- "phase_2"
dge$samples$stages[dge$samples$timepoint >=35 ] <- "mature"

df.anno <- as.data.frame(dge$samples)
df.anno$batch <- as.factor(df.anno$batch )
df.anno <- df.anno[,-1]

write.csv(df.anno, "tables/annotation.csv", quote = TRUE)

set.seed(123)
ha = HeatmapAnnotation(df = df.anno,  col = list(batch = c("1" = "#4C028BFF",
                                                                "2" = "#DE6106FF") ,
                                                 stages = c("initial" = "blue3",
                                                               "phase_1" = "cyan2",
                                                               "phase_2" = "#AE04D5FF",
                                                               "mature" = "darkgreen")))

#get_color_mapping_list(ha)


pdf("plots/phased_genes_HM_annot.pdf", width = 10, height = 8)
Heatmap(edgeR::cpm(dge, log = TRUE), 
        name = "log2cpm",
        row_names_gp = grid::gpar(fontsize = 5),
        column_names_gp = grid::gpar(fontsize = 5),
        cluster_rows=FALSE,
        row_title = NULL,
         top_annotation = ha,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = TRUE)
dev.off()


pdf("plots/phased_genes_HM_annot_loLabs.pdf", width = 7, height = 6.5)
Heatmap(edgeR::cpm(dge, log = TRUE), 
        name = "log2cpm",
        row_names_gp = grid::gpar(fontsize = 5),
        column_names_gp = grid::gpar(fontsize = 5),
        cluster_rows=FALSE,
        row_title = NULL,
        top_annotation = ha,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE)
dev.off()


#try to extract features
#STEP: try to do the NMF
# ButchR: https://github.com/wurst-theke/ButchR
library(reticulate)
reticulate::use_condaenv("tf_4.2", required = TRUE) #tf_4.2 is a conda environment containing tensorflow
#py_config()

#conda env needs to be activated before ButchR
library(ButchR)

#get counts 
matrix.cts <- edgeR::cpm(dge, log = FALSE) #normalized counts, but not log


#latent dimensions - span
k_min <- 3
k_max <- 10

seed <- 102
set.seed(seed)

procr_nmf_exp <- run_NMF_tensor(X = matrix.cts,
                                ranks = k_min:k_max,
                                method = "NMF",
                                n_initializations = 10,  #to make it run faster
                                extract_features = TRUE,
                                seed = seed)

# [1] "2025-11-15 19:44:25 CET"
# Factorization rank:  3 
# [1] "NMF converged after  117,77,163,115,198,88,93,73,80,135 iterations"
# [1] "2025-11-15 19:44:27 CET"
# Factorization rank:  4 
# [1] "NMF converged after  116,135,61,102,88,211,168,316,140,85 iterations"
# [1] "2025-11-15 19:44:29 CET"
# Factorization rank:  5 
# [1] "NMF converged after  247,322,102,156,231,188,153,74,250,135 iterations"
# [1] "2025-11-15 19:44:32 CET"
# Factorization rank:  6 
# [1] "NMF converged after  111,130,210,292,229,182,123,92,121,85 iterations"
# [1] "2025-11-15 19:44:35 CET"
# Factorization rank:  7 
# [1] "NMF converged after  300,301,222,373,522,115,163,279,128,263 iterations"
# [1] "2025-11-15 19:44:39 CET"
# Factorization rank:  8 
# [1] "NMF converged after  261,223,296,273,156,319,237,232,187,264 iterations"
# [1] "2025-11-15 19:44:43 CET"
# Factorization rank:  9 
# [1] "NMF converged after  242,429,120,174,275,157,122,238,251,195 iterations"
# [1] "2025-11-15 19:44:47 CET"
# Factorization rank:  10 
# [1] "NMF converged after  125,158,213,266,319,112,193,89,325,198 iterations""

## SAVING THE DECOMPOSITION
#save(procr_nmf_exp , file = "Rdata/NMF_mydata.rda")


#stability of the signatures
plt <- generateRiverplot(procr_nmf_exp, edges.cutoff = 0.2)
pdf("plots/NMF_riverplot.pdf", width = 8, height = 10)
plot(plt, plot_area = 0.8, yscale = 0.5, nodewidth = 0.5) + theme(axis.text = element_text(size = 1)) + geom_text_repel()
dev.off()

#signature 8 most favourable/stable
pdf("plots/NMF_k_selection.pdf", width = 5, height = 5)
gg_plotKStats(procr_nmf_exp)
dev.off()


# Follow with the latent dimension 8
Hk8 <- HMatrix(procr_nmf_exp, k = 8)
rownames(Hk8) <- c("Sign.1", "Sign.2", "Sign.3", "Sign.4", "Sign.5",
                   "Sign.6", "Sign.7", "Sign.8")

Hk8 <- as.matrix(Hk8)
h_heatmap_8 <- Heatmap(Hk8,
                         col = viridis(100),
                         name = "Exposure",
                         clustering_distance_columns = 'pearson',
                         show_column_dend = TRUE,
                         top_annotation = ha,
                         show_column_names = FALSE,
                         show_row_names = TRUE,
                         cluster_columns = TRUE,
                         cluster_rows = FALSE)


h_heatmap_8.b <- Heatmap(Hk8,
                       col = viridis(100),
                       name = "Exposure",
                       clustering_distance_columns = 'pearson',
                       show_column_dend = TRUE,
                       top_annotation = ha,
                       show_column_names = FALSE,
                       show_row_names = TRUE,
                       cluster_columns = FALSE,
                       cluster_rows = FALSE)


pdf("plots/NMF_H.pdf", width = 8.5, height = 5)
print(h_heatmap_8)
print(h_heatmap_8.b)
dev.off()

# pca and umap clustering?
library(umap)
Hk8.umap <- umap(t(Hk8))

umap_df <- Hk8.umap$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2")
umap_df$anno.stages <- df.anno$stages 
umap_df$timepoint <- df.anno$timepoint

#the phases show specific embeddings
pdf(paste0("plots/UMAP_Hk8.pdf"), width = 6,height = 5)
umap_df %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2, 
             color = anno.stages))+
  geom_point()+ theme_bw(18) +
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "UMAP embedding of the exposure matrix H")

umap_df %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2, 
             color = timepoint))+
  geom_point()+ theme_bw(18) +
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "UMAP embedding of the exposure matrix H")
dev.off()


# PCA
#
pca <- prcomp(t(Hk8), center = TRUE, scale. = FALSE)
fov <- pca$sdev^2/sum(pca$sdev^2)

pcaData <- pca$x |>
  as_tibble() |>
  mutate(Sample = rownames(pca$x),
         Evolution = "evo",
         batch = as.factor(dge$samples$batch)) 


# Principal components
pc1 <- ggplot(data = pcaData,
              mapping = aes(x = PC1, y = PC2, color=Sample)) +
  geom_point() + 
  geom_path(aes(group = Evolution)) +
  labs(x = sprintf('PC1 (%.0f%%)', 100*fov[1]),
       y = sprintf('PC2 (%.0f%%)', 100*fov[2])) +
  geom_text_repel(aes(label = Sample, colour = Sample, segment.size = 0.3)) +
  guides(colour = guide_legend(override.aes = list(label=""))) +
  theme_bw()

pdf(paste0("plots/PCA_Hk8.pdf"), width = 10,height = 5)
print(pc1)
dev.off()




######################
#Signatures
#maturation: Sign7 : final maturation, Sign1 meand maturation , Sign 6: initial form of maturation, Sign3. Transition to maturation
# exit: Signature8
# phase1: Sign5
# phase2: Sign2
# between phases: Sign.4


#Association plots
pdf(paste0("plots/AssociationWithGroups_8.pdf"), width = 8,height = 5)
recovery_plot(Hk8,
              df.anno$stages )
dev.off()


# what are the signatures
#W matrix
library(viridis)
library(knitr)

#load( "Rdata/NMF_mydata.rda")
Wk8 <- WMatrix(procr_nmf_exp, k = 8)
dim(Wk8) #4407    8
colnames(Wk8) <- paste0("Sign.", 1:8)

#SIGNATURES
features <- SignatureSpecificFeatures(procr_nmf_exp,
                                      k = 8, 
                                      return_all_features = TRUE)
colnames(features) <- paste0("Sign.", 1:8)
head(features)
#show

#which unique genes constitute the signature
names_specific <- rownames(features)[rowSums(features) == 1]
length(names_specific) #  992

#reduce to signature-specific genes
Wk8.specific <- Wk8[names_specific, ]
norm_Wspecific <- Wk8.specific/matrixStats::rowMaxs(Wk8.specific)

w_heatmap <- Heatmap(norm_Wspecific,
                     col = inferno(100),
                     name = "W matrix",
                     clustering_distance_columns = 'pearson',
                     show_column_dend = TRUE,
                     show_column_names = TRUE,
                     show_row_names = FALSE,
                     cluster_rows = TRUE,
                     cluster_columns = FALSE)

pdf("plots/NMF_W8.pdf", width = 5, height = 6)
print(w_heatmap)
dev.off()

#save the signatures
Sign1.genes <- rownames(norm_Wspecific)[norm_Wspecific[,1] == 1 ] # 1
Sign2.genes <- rownames(norm_Wspecific)[norm_Wspecific[,2] == 1 ] #530
Sign3.genes <- rownames(norm_Wspecific)[norm_Wspecific[,3] == 1 ] #29
Sign4.genes <- rownames(norm_Wspecific)[norm_Wspecific[,4] == 1 ]  #115
Sign5.genes <- rownames(norm_Wspecific)[norm_Wspecific[,5] == 1 ]  #229
Sign6.genes <- rownames(norm_Wspecific)[norm_Wspecific[,6] == 1 ]  # 11
Sign7.genes <- rownames(norm_Wspecific)[norm_Wspecific[,7] == 1 ] #1
Sign8.genes <- rownames(norm_Wspecific)[norm_Wspecific[,8] == 1 ]  #76

write.table(Sign1.genes, "tables/Signature1.csv", quote = TRUE, row.names = FALSE, col.names = FALSE)
write.table(Sign2.genes, "tables/Signature2.csv", quote = TRUE, row.names = FALSE, col.names = FALSE)
write.table(Sign3.genes, "tables/Signature3.csv", quote = TRUE, row.names = FALSE, col.names = FALSE)
write.table(Sign4.genes, "tables/Signature4.csv", quote = TRUE, row.names = FALSE, col.names = FALSE)
write.table(Sign5.genes, "tables/Signature5.csv", quote = TRUE, row.names = FALSE, col.names = FALSE)
write.table(Sign6.genes, "tables/Signature6.csv", quote = TRUE, row.names = FALSE, col.names = FALSE)
write.table(Sign7.genes, "tables/Signature7.csv", quote = TRUE, row.names = FALSE, col.names = FALSE)
write.table(Sign8.genes, "tables/Signature8.csv", quote = TRUE, row.names = FALSE, col.names = FALSE)

#notes: the single gene signature are not really meaningful

sessionInfo()
# R version 4.2.2 (2022-10-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 22.04.5 LTS
# 
# Matrix products: default
# BLAS/LAPACK: /scicore/soft/easybuild/apps/FlexiBLAS/3.2.1-GCC-12.2.0/lib/libflexiblas.so.3.2
# 
# locale:
#   [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8        LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
# [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C           LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
# 
# attached base packages:
#   [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] knitr_1.42                  umap_0.2.10.0               ButchR_1.0                  reticulate_1.28             viridis_0.6.2              
# [6] viridisLite_0.4.1           patchwork_1.1.2             ComplexHeatmap_2.14.0       vroom_1.6.1                 readr_2.1.4                
# [11] tidyr_1.3.0                 tibble_3.2.1                tidyverse_2.0.0             dplyr_1.1.1                 ggrepel_0.9.3              
# [16] ggplot2_3.4.2               edgeR_3.40.2                limma_3.54.2                SummarizedExperiment_1.28.0 Biobase_2.58.0             
# [21] GenomicRanges_1.50.2        GenomeInfoDb_1.34.9         IRanges_2.32.0              S4Vectors_0.36.2            BiocGenerics_0.44.0        
# [26] MatrixGenerics_1.10.0       matrixStats_0.63.0          rstudioapi_0.14            
# 
# loaded via a namespace (and not attached):
#   [1] bitops_1.0-7           bit64_4.0.5            doParallel_1.0.17      RColorBrewer_1.1-3     rprojroot_2.0.3        tools_4.2.2           
# [7] utf8_1.2.3             R6_2.5.1               colorspace_2.1-0       GetoptLong_1.0.5       withr_2.5.0            tidyselect_1.2.0      
# [13] gridExtra_2.3          bit_4.0.5              compiler_4.2.2         cli_3.6.1              csaw_1.32.0            DelayedArray_0.24.0   
# [19] labeling_0.4.2         scales_1.2.1           nnls_1.4               askpass_1.1            rappdirs_0.3.3         digest_0.6.31         
# [25] Rsamtools_2.14.0       XVector_0.38.0         pkgconfig_2.0.3        rlang_1.1.0            GlobalOptions_0.1.2    shape_1.4.6           
# [31] generics_0.1.3         farver_2.1.1           riverplot_0.10         jsonlite_1.8.4         BiocParallel_1.32.6    RCurl_1.98-1.12       
# [37] magrittr_2.0.3         GenomeInfoDbData_1.2.9 Matrix_1.5-4           Rcpp_1.0.10            munsell_0.5.0          fansi_1.0.4           
# [43] lifecycle_1.0.3        zlibbioc_1.44.0        parallel_4.2.2         crayon_1.5.2           lattice_0.21-8         Biostrings_2.66.0     
# [49] cowplot_1.1.1          circlize_0.4.15        hms_1.1.3              locfit_1.5-9.7         magick_2.7.4           metapod_1.6.0         
# [55] pillar_1.9.0           rjson_0.2.21           codetools_0.2-19       glue_1.6.2             png_0.1-8              vctrs_0.6.1           
# [61] tzdb_0.3.0             foreach_1.5.2          openssl_2.0.6          gtable_0.3.3           purrr_1.0.1            clue_0.3-64           
# [67] xfun_0.38              RSpectra_0.16-1        iterators_1.0.14       cluster_2.1.4          here_1.0.1            
# > 
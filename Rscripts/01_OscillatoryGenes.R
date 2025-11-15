

# upload of time-resolved RNA seq on c. elegans capturing developmental clock
# 2020 Mollecular systems biology, MWM Meeuse et al


suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(limma)
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(ComplexHeatmap)
  library(ggrepel)
  library(patchwork)
  
  library(rstudioapi)
})

#Set wd to the directory of the script (in R studio)
setwd(sub("01_OscillatoryGenes\\.R", "", getSourceEditorContext()$path) )

#Generate folder structure to save intermediate results
directories <- c("plots","tables","Rdata")

for (directory in directories){
  if (!dir.exists(directory)){
    dir.create(directory)
    print(paste("dir", directory, "created"))
  } else {
    print(paste("dir", directory, "exists"))
  }
}

#read-in data
RNAseq.1 <- read.table("../input_data/GSE130782_expr_mRNA.tab")
RNAseq.2 <-  read.table("../input_data/GSE130811_expr_mRNA_CE10_coding.tab")

#data come in 2 batches, in the batch 2, some data points are overlapping
#clean data, first row reflects the exon width, rownames are gene names
counts.1 <- RNAseq.1[,2:dim(RNAseq.1)[2]]
exon.widths.1  <- RNAseq.1[,1]

counts.2 <- RNAseq.2[,2:dim(RNAseq.2)[2]]
exon.widths.2  <- RNAseq.2[,1]

#generate DGELists
dge.1 <- DGEList(counts = counts.1)
dge.2 <- DGEList(counts = counts.2)

#filter low expressed genes and compute  normalization factors
# filtering based on .2 batch,  many genes start later expression
keep.exprs.2 <- filterByExpr(dge.2, group = colnames(dge.2)) #each sample represents one timepoint (group)
length(keep.exprs.2); sum(keep.exprs.2)
#1] 20392
#[1] 17094

#use shared genes
names.to.keep <- rownames( dge.2[keep.exprs.2,])[rownames( dge.2[keep.exprs.2,]) %in% rownames(dge.1 )]
length(names.to.keep ) #16772

dge.2 <- dge.2[names.to.keep ,, keep.lib.sizes=FALSE]
dge.2 <- calcNormFactors(dge.2)
dge.2$samples$norm.factors

dge.1 <- dge.1[names.to.keep ,, keep.lib.sizes=FALSE]
dge.1 <- calcNormFactors(dge.1)
dge.1$samples$norm.factors

#unite the genes
logCPM <- edgeR::cpm(dge.2, log = TRUE)

#oscillations 10 to 25 hours
dge.cycle <- dge.2[,c(which(colnames(dge.2) == "X10hr"):which(colnames(dge.2) == "X25hr"))]

pdf("plots/corelations_data2.pdf", width = 6, height = 5)
Heatmap(cor(edgeR::cpm(dge.2,log = TRUE)), 
        name = "correlations",
        row_names_gp = grid::gpar(fontsize = 5),
        column_names_gp = grid::gpar(fontsize = 5),
        cluster_rows=FALSE,
        row_title = NULL,
        #   top_annotation = ha,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE)

Heatmap(cor(edgeR::cpm(dge.cycle,log = TRUE)), 
        name = "correlations",
        row_names_gp = grid::gpar(fontsize = 10),
        column_names_gp = grid::gpar(fontsize = 10),
        cluster_rows=FALSE,
        row_title = NULL,
        #   top_annotation = ha,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE)
dev.off()


#ok, identify oscillatory genes
EE <- as.matrix(edgeR::cpm(dge.2,log = TRUE)[,1:22])
            
pca <- prcomp(t(EE), center = TRUE, scale. = FALSE)
fov <- pca$sdev^2/sum(pca$sdev^2)

pcaData <- pca$x |>
  as_tibble() |>
  mutate(Sample = rownames(pca$x),
         Evolution = "evo") 



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

pc2 <- ggplot(data = pcaData,
              mapping = aes(x = PC2, y = PC3, color=Sample)) +
  geom_point() + 
  geom_path(aes(group = Evolution)) +
  labs(x = sprintf('PC2 (%.0f%%)', 100*fov[2]),
       y = sprintf('PC3 (%.0f%%)', 100*fov[3])) +
  geom_text_repel(aes(label = Sample, colour = Sample, segment.size = 0.3)) +
  guides(colour = guide_legend(override.aes = list(label=""))) +
  theme_bw()


pdf("plots/cycles_pca_plots.pdf", width = 17, height = 5.5)
pc1  + pc2 
dev.off()

#######################
#check individual ocillating genes
#regression with 7.5 hours period (experimentally known)
omega = 2*pi/7.5
Time = seq( 8,34)
xc <- cos(omega*Time)
xs <- - sin(omega*Time)


y <- edgeR::cpm(dge.2,log = TRUE)[,4:30]

pdf("plots/cyclingGenes.pdf", width = 17, height = 5.5)
for (i in seq(1500) ){
  fit.lm <- lm(y[i,] ~  xc+xs + Time)
  pred <- stats::predict(fit.lm, newdata=data.frame(Time=Time))   
  plot(y[i,] ~ Time,  main=paste("row id", i, rownames(y)[i]))
  lines( pred  ~ Time, col="blue")
}
dev.off()


######### Regression model
dsg <- model.matrix(~  xc + xs + Time)
limma.fit <- lmFit(y,design = dsg)

limma.fit$coefficients %>% head(5)
#                (Intercept)          xc          xs         Time
#WBGene00007063    2.770971  0.07186378 -0.31439996 -0.011419590
#WBGene00007064    5.588660 -0.12071314 -0.14268798 -0.031230437
#WBGene00044165    1.227348  0.19772444 -0.24557401 -0.047292633
#WBGene00007065    1.398517  0.03396921 -0.04969889  0.042402780
#WBGene00003525    1.928584  0.94907147 -2.50324600  0.006059336


phases <- atan2(limma.fit$coefficients[,"xc"],limma.fit$coefficients[,"xs"] )
amplitudes <- sqrt(limma.fit$coefficients[,"xc"]**2 + limma.fit$coefficients[,"xs"]**2)
trend <- limma.fit$coefficients[,"Time"]

summary(amplitudes)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0000  0.1092  0.2314  0.4680  0.4894  4.0242

#plot sorted genes that show some oscillatory trend (via amplitude)
oscillatory.genes <- names(amplitudes[amplitudes > mean(amplitudes)])

phases <-  sort(phases[oscillatory.genes])
amplitudes <- amplitudes[names(phases)]
trend <- trend[names(phases)]

length(phases ) 
#[1] 4407 ... 4 and half thousands of them

y.osc <- y[names(phases),]
pdf("plots/cyclingGenes_selected.pdf", width = 17, height = 5.5)
for (i in seq(dim(y.osc)[1]) ){
  fit.lm <- lm(y.osc[i,] ~  xc+xs + Time)
  pred <- stats::predict(fit.lm, newdata=data.frame(Time=Time))   
  plot(y.osc[i,] ~ Time,  main=paste("row id", i, rownames(y)[i]))
  lines( pred  ~ Time, col="magenta")
}
dev.off()


#plot the heatmap
#AMAZING!s
heat.osc <- Heatmap(y.osc, 
        name = "log2cpm",
        row_names_gp = grid::gpar(fontsize = 10),
        column_names_gp = grid::gpar(fontsize = 10),
        cluster_rows=FALSE,
        row_title = NULL,
        #   top_annotation = ha,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = TRUE)

pdf("plots/HM_oscillatoryGenes.pdf", width = 8, height = 7)
print(heat.osc)
dev.off()

#select these names for the whole data set
dge.1 <- dge.1[names(phases),]
dge.2 <- dge.2[names(phases),]
dge <- cbind(dge.1,dge.2)

dge$samples$batch <- 1
dge$samples$batch[17:length(dge$samples$batch)] <- 2

dge$samples$timepoint <- c(c(0:15),c(5:48),c(37,38,39,41,45,47,48))

#samples ordered by time, genes ordered by phase
dge <- dge[,order(dge$samples$timepoint)]

heat.whole <- Heatmap(edgeR::cpm(dge, log = TRUE), 
        name = "log2cpm",
        row_names_gp = grid::gpar(fontsize = 10),
        column_names_gp = grid::gpar(fontsize = 10),
        cluster_rows=FALSE,
        row_title = NULL,
        #   top_annotation = ha,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = TRUE)

pdf("plots/HM_allGenes-aligned.pdf", width = 12, height = 7)
print(heat.whole )
dev.off()

#pca on all
EE <- as.matrix(edgeR::cpm(dge, log = TRUE))

pca <- prcomp(t(EE), center = TRUE, scale. = FALSE)
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

pc2 <- ggplot(data = pcaData,
              mapping = aes(x = PC1, y = PC2, color=batch)) +
  geom_point() + 
  geom_path(aes(group = Evolution)) +
  labs(x = sprintf('PC1 (%.0f%%)', 100*fov[1]),
       y = sprintf('PC2 (%.0f%%)', 100*fov[2])) +
  geom_text_repel(aes(label = Sample, colour = batch, segment.size = 0.3)) +
  guides(colour = guide_legend(override.aes = list(label=""))) +
  theme_bw()

pdf("plots/cycles_pca_plots_allDatapoints_cycling_genes.pdf", width = 17, height = 5.5)
pc1  + pc2 
dev.off()

saveRDS(dge,"Rdata/ordered_expression.rds")

sessionInfo()
# R version 4.2.2 (2022-10-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 22.04.5 LTS
# 
# Matrix products: default
# BLAS/LAPACK: /scicore/soft/easybuild/apps/FlexiBLAS/3.2.1-GCC-12.2.0/lib/libflexiblas.so.3.2
# 
# locale:
#   [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8        LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8    LC_PAPER=C.UTF-8      
# [8] LC_NAME=C              LC_ADDRESS=C           LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
# 
# attached base packages:
#   [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] patchwork_1.1.2             ComplexHeatmap_2.14.0       vroom_1.6.1                 readr_2.1.4                 tidyr_1.3.0                
# [6] tibble_3.2.1                tidyverse_2.0.0             dplyr_1.1.1                 ggrepel_0.9.3               ggplot2_3.4.2              
# [11] edgeR_3.40.2                limma_3.54.2                SummarizedExperiment_1.28.0 Biobase_2.58.0              GenomicRanges_1.50.2       
# [16] GenomeInfoDb_1.34.9         IRanges_2.32.0              S4Vectors_0.36.2            BiocGenerics_0.44.0         MatrixGenerics_1.10.0      
# [21] matrixStats_0.63.0          rstudioapi_0.14            
# 
# loaded via a namespace (and not attached):
#   [1] bit64_4.0.5            foreach_1.5.2          GenomeInfoDbData_1.2.9 Rsamtools_2.14.0       pillar_1.9.0           lattice_0.21-8         glue_1.6.2            
# [8] digest_0.6.31          RColorBrewer_1.1-3     XVector_0.38.0         colorspace_2.1-0       Matrix_1.5-4           pkgconfig_2.0.3        GetoptLong_1.0.5      
# [15] csaw_1.32.0            magick_2.7.4           zlibbioc_1.44.0        purrr_1.0.1            scales_1.2.1           tzdb_0.3.0             BiocParallel_1.32.6   
# [22] farver_2.1.1           generics_0.1.3         withr_2.5.0            cli_3.6.1              magrittr_2.0.3         crayon_1.5.2           fansi_1.0.4           
# [29] doParallel_1.0.17      tools_4.2.2            hms_1.1.3              GlobalOptions_0.1.2    lifecycle_1.0.3        munsell_0.5.0          locfit_1.5-9.7        
# [36] cluster_2.1.4          DelayedArray_0.24.0    Biostrings_2.66.0      compiler_4.2.2         rlang_1.1.0            RCurl_1.98-1.12        iterators_1.0.14      
# [43] rjson_0.2.21           circlize_0.4.15        labeling_0.4.2         bitops_1.0-7           gtable_0.3.3           codetools_0.2-19       R6_2.5.1              
# [50] gridExtra_2.3          bit_4.0.5              utf8_1.2.3             clue_0.3-64            metapod_1.6.0          shape_1.4.6            parallel_4.2.2        
# [57] Rcpp_1.0.10            vctrs_0.6.1            png_0.1-8              tidyselect_1.2.0  
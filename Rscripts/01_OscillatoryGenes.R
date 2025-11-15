

# upload of time-resolved RNA seq on c. elegans capturing developmental clock
# 2020 Mollecular systems biology, MWM Meeuse et al

setwd("/scicore/home/bentires/GROUP/michal/Developmental__oscillator/Grosshans_Data")

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(limma)
  library(edgeR)
  library(stringr)
  library(ggplot2)
  library(ggfortify)
  library(ggrepel)
  library(dplyr)
  library(tidyverse)
  library(vroom)
  library(ComplexHeatmap)
})


#read-in data
RNAseq.1 <- read.table("input_data/GSE130782_expr_mRNA.tab")
RNAseq.2 <-  read.table("input_data/GSE130811_expr_mRNA_CE10_coding.tab")

counts.1 <- RNAseq.1[,2:dim(RNAseq.1)[2]]
exon.widths.1  <- RNAseq.1[,1]

counts.2 <- RNAseq.2[,2:dim(RNAseq.2)[2]]
exon.widths.2  <- RNAseq.2[,1]

dge.1 <- DGEList(counts = counts.1)
dge.2 <- DGEList(counts = counts.2)

#filter low expressed genes and compute  normalization factors
# filtering based on .2 batch,  many genes start later expression

keep.exprs.2 <- filterByExpr(dge.2, group = colnames(dge.2))
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

Heatmap(logCPM[1:3000,], 
        name = "log2cpm",
        row_names_gp = grid::gpar(fontsize = 10),
        column_names_gp = grid::gpar(fontsize = 10),
        cluster_rows=TRUE,
        row_title = NULL,
     #   top_annotation = ha,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = TRUE)


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

library(ggrepel)
library(patchwork)

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

#check inidvidual genes

#regression with 8 hours period
omega = 2*pi/7.5
#Time = seq(5,26)
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
Heatmap(y.osc, 
        name = "log2cpm",
        row_names_gp = grid::gpar(fontsize = 10),
        column_names_gp = grid::gpar(fontsize = 10),
        cluster_rows=FALSE,
        row_title = NULL,
        #   top_annotation = ha,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = TRUE)

#select these names for the whole data set

dge.1 <- dge.1[names(phases),]
dge.2 <- dge.2[names(phases),]
dge <- cbind(dge.1,dge.2)

dge$samples$batch <- 1
dge$samples$batch[17:length(dge$samples$batch)] <- 2

dge$samples$timepoint <- c(c(0:15),c(5:48),c(37,38,39,41,45,47,48))

#samples ordered by time, genes ordered by phase
dge <- dge[,order(dge$samples$timepoint)]

Heatmap(edgeR::cpm(dge, log = TRUE), 
        name = "log2cpm",
        row_names_gp = grid::gpar(fontsize = 10),
        column_names_gp = grid::gpar(fontsize = 10),
        cluster_rows=FALSE,
        row_title = NULL,
        #   top_annotation = ha,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = TRUE)

#pca on all
EE <- as.matrix(edgeR::cpm(dge, log = TRUE))

pca <- prcomp(t(EE), center = TRUE, scale. = FALSE)
fov <- pca$sdev^2/sum(pca$sdev^2)

pcaData <- pca$x |>
  as_tibble() |>
  mutate(Sample = rownames(pca$x),
         Evolution = "evo",
         batch = as.factor(dge$samples$batch)) 

library(ggrepel)
library(patchwork)

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

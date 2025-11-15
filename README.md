# Oscillatory Patterns in Expression During the Larval Development of C. Elegans.
The repo contains code that reproduces and extends results from the paper [Developmental function and state transitions of a gene expression oscillator in Caenorhabditis elegans](https://www.embopress.org/doi/full/10.15252/msb.20209498) published by the group of Helge Grosshans, Friedrich Miescher Institut, Novartis. 

Here is the description, what the code does:

## Cosine Fitting
The gene expression shows a periodic signal with a period around 7.5 hours.
![pca plot](https://github.com/MikeKlocCZ/RNAseq_Developmental-Oscillator_GrosshansLab/blob/main/figures/cycles_pca_plots.png "pca plot")

With this periodicity in mind, we can perform a cosine fitting that establishes a phase shift between the genes. Based on this phase shift, the genes can be sorted, which reveals the oscillatory pattern.
![HM oscillation](https://github.com/MikeKlocCZ/RNAseq_Developmental-Oscillator_GrosshansLab/blob/main/figures/phased_genes_HM_annot.png "oscillation - HM")


## Dimensionality Reduction: Non-Negative Matrix Factorization (NMF)
There are many oscillatory genes. In order to obtain a more "functional" description, the idea was to reduce the dimensionality with an NMF factorization, namely using the [ButchR package](https://github.com/wurst-theke/ButchR). It uses Tensorflow to perform the factorization, which is "called" using the `reticulate` package in R.

The `ButchR` package offers several measures to establish a suitable latent space dimensionality $k$. Here, the optimal option was $k = 8$. This is the exposure matrix.
![Exposure](https://github.com/MikeKlocCZ/RNAseq_Developmental-Oscillator_GrosshansLab/blob/main/figures/NMF_H-0.png "Exposure")

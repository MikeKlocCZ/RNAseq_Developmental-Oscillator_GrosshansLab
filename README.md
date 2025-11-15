# Oscillatory patterns in expression during the larval development of C. elegans.
The repo contains code that reproduces and extends results from the paper [Developmental function and state transitions of a gene expression oscillator in Caenorhabditis elegans](https://www.embopress.org/doi/full/10.15252/msb.20209498) published by the group of Helge Grosshans, Friedrich Miescher Institut, Novartis.

## Cosine fitting
The gene expression shows a periodic signal with a period around 7.5 hours.
![pca plot](https://github.com/MikeKlocCZ/RNAseq_Developmental-Oscillator_GrosshansLab/blob/main/figures/cycles_pca_plots.png "pca plot")

With this periodicity in mind, we can perform a cosine fitting that establishes a phase shift between the genes. Based on this phase shift, the genes can be sorted, which reveals the oscillatory pattern.
![HM oscillation](https://github.com/MikeKlocCZ/RNAseq_Developmental-Oscillator_GrosshansLab/blob/main/figures/phased_genes_HM_annot.png "oscillation - HM")

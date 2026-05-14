# Neuron_Seurat_code
Code for single-cell RNA-seq analysis in the Neuron Seurat project.
Title:
Analysis code for the single-cell RNA-seq study

Description:
This folder contains the R scripts used for the single-cell RNA-seq analyses reported in the manuscript.

Files:
1. 01_QC_merge_batch_correction.R
   - Data loading
   - Quality control
   - Sample merging
   - Batch correction / data integration

2. 02_clustering_annotation.R
   - Normalization
   - Variable feature selection
   - Scaling
   - PCA
   - Neighbor graph construction
   - Clustering
   - UMAP visualization
   - Marker gene identification
   - Cell-type annotation

3. 03_cell_proportion_analysis.R
   - Calculation of cell-type proportions
   - Group comparison
   - Visualization

4. 04_monocle2_pseudotime_analysis.R
   - Pseudotime analysis using Monocle 2
   - Cell trajectory reconstruction
   - Related visualizations

5. 05_mitochondrial_gene_score.R
   - Mitochondrial gene score calculation
   - Downstream comparison and visualization

Notes:
- The scripts contain the actual parameters used to generate the results in the manuscript.
- Local file paths were replaced with generic relative paths.
- Input data are described in the manuscript and associated data availability statement.

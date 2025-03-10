# WGCNA-Longevidad-Ratones

This project performs Weighted Gene Co-expression Network Analysis (WGCNA) using RNA sequencing count data. The objective is to identify gene modules associated with phenotypic traits such as condition, sex, and age. The analysis includes data preprocessing, normalization, variance filtering, clustering, module detection, and functional enrichment analysis.

Data Preparation

Count Data: Read and filter RNA-seq count data from SM539_countTable.txt.

Sample Information: Read sample metadata from SM539_sampleinfo.txt.

Normalization: Using DESeq2 to normalize read counts and remove lowly expressed genes.

Quality Control and Clustering

Variance Filtering: Remove genes with low variance to enhance clustering accuracy.

Sample Clustering: Generate dendrograms to visualize relationships between samples.

Trait Heatmap: Map phenotypic traits to sample clustering using heatmaps.

Network Construction (WGCNA)

Soft Thresholding: Determine the optimal power for network construction.

Adjacency Matrix & TOM Calculation: Convert expression data into a co-expression network.

Module Detection: Identify gene modules using hierarchical clustering.

Module-Trait Associations

Eigengene Correlation: Identify modules associated with phenotypic traits.

Heatmaps: Visualize correlations between gene modules and traits.

Functional Enrichment Analysis

Gene Ontology (GO) Analysis: Identify biological processes enriched in each module.

Data Conversion: Convert gene IDs to ENTREZ for enrichment analysis.

Visualization: Generate dot plots for top GO terms.

Output Files

GOEnrichmentTable.xlsx: Full GO enrichment analysis results.

Best_GO_Terms.xlsx: Top 5 GO terms per module.

Various plots and heatmaps for quality control, clustering, and module relationships.


MIT License

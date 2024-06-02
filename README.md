# Protein expression quantitative trait loci (pQTLs): software and analytic codes

## Overview
**Study Title**: Mutation Impact on mRNA Versus Protein Expression across Human Cancers

**Authors**: 
- Yuqi Liu¹
- Abdulkadir Elmas¹
- Kuan-lin Huang¹#

### Abstract
Cancer mutations are often assumed to alter proteins, thus promoting tumorigenesis. However, how mutations affect protein expression has rarely been systematically investigated. We conduct a comprehensive analysis of mutation impacts on mRNA- and protein-level expressions of 953 cancer cases with paired genomics and global proteomic profiling across six cancer types. Protein-level impacts are validated for 47.2% of the somatic expression quantitative trait loci (seQTLs), including mutations from likely “long-tail” driver genes. Devising a statistical pipeline for identifying somatic protein-specific QTLs (spsQTLs), we reveal several gene mutations, including NF1 and MAP2K4 truncations and TP53 missenses showing disproportional influence on protein abundance not readily explained by transcriptomics. Cross-validating with data from massively parallel assays of variant effects (MAVE), TP53 missenses associated with high tumor TP53 proteins were experimentally confirmed as functional. Our study demonstrates the importance of considering protein-level expression to validate mutation impacts and identify functional genes and mutations.

## File Descriptions
All files are located in the `analysis` directory.

- **analysis.ipynb**: Simple analyses.
- **DNA_Pro_Regression.R**: Find pQTLs in prospective data with multivariate linear regression.
- **DNA_RNA_Regression.R**: Find eQTLs in prospective data with multivariate linear regression.
- **Figs.R**: Generate the figures in the manuscript.
- **Fisher_Exact.R**: Conduct Fisher exact test to protein data to see if there is significant pQTL enrichment in truncation.
- **heatmap_data.ipynb**: Extract the data used to make heatmap.
- **HGVSp.R**: Extract HGVSp short.
- **Lolliplot_Data.R**: Prepare the data needed to generate lolliplots.
- **lrt.R**: Likelihood ratio test.
- **Pro_regression_retro_cptac.R**: Find pQTLs in retrospective data with multivariate linear regression.
- **pros_retro_overlap.ipynb**: Overlap retrospective data with prospective data.
- **retro_cptac_mutation_matrix.R**: Prepare mutation matrix of retrospective data.
- **RNA_regression_retro_cptac.R**: Find eQTLs in retrospective data with multivariate linear regression.
- **sig_pro_rna.ipynb**: Generate a table that has information of QTLs that are significant at no less than one level.
- **violin_plot_data.ipynb**: Prepare the data needed to generate violin plots.

## Contact
For any questions or further information, please open an issue. 

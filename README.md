# SCHoCO-seq-Paper
This repository contains all scripts to fully reproduce plots and statistical analysis.

The analysis of single-cell microbiome data is based on the phyloseq data system:

  phyloseq: An R package for reproducible interactive analysis and
  graphics of microbiome census data. Paul J. McMurdie and Susan
  Holmes (2013) PLoS ONE 8(4):e61217.

To run the provided scripts, necessary input data comprise an ASV count table resulting from the DADA2 pipeline (Callahan et al. 2016), a corresponding ASV taxonomy file that links any ASV to a respective taxonomic annotation, and a sample info/metadata file that links the sample IDs of the ASV count table to respective sample variables. Moreover, in this study a file containing potential contaminant genera is provided which is used to subset the taxonomy and thereby remove any potential contaminations.

This study is devided into two sections which are contained in the two folders of this repository: 
- SC-seq, which is representing diatom single-cell microbiome analysis without Cas9-optimized sequenzing
- SCHoCO-seq, containing microbiome data analysis for diatom microbiomes which underwent Cas9-optimized sequenzing

# Analyses and data of:

_Padfield et al. (2020) Evolution of diversity explains the impact of pre-adaptation of a focal species on the structure of a natural microbial community. ISME_ 

**If you would like access to the paper and do not have it, please email me on d.padfield@exeter.ac.uk**

Link to read paper: https://rdcu.be/b9q7X

DOI of paper: https://doi.org/10.1038/s41396-020-00755-3

DOI of analyses and dataset: 

Link to raw sequencing files: https://www.ebi.ac.uk/ena/browser/view/PRJEB40882

### Outline

This repository contains the final datasets, analyses and figures of the above-mentioned paper. It can recreate the all of the analyses and figures in the main text (Figure 2 to Figure 7) and all of the Supplementary Information.

### Feedback

- Please report any problems or bugs in the code in the [Issues](https://github.com/padpadpadpad/Padfield_2019_ISME_bact_phage_temperature) tab of the GitHub repository. Alternatively, please email _d.padfield@exeter.ac.uk_.

### Licensing

This code is licensed under GPL-3.

### Notes

- The bioinformatics pipeline for calling genetic variants from the raw whole genome sequencing data is provided in `scripts/clone_sequencing_pipeline.sh`. The vcf output was then imported to R and the final dataset `data/genetic_changes.csv` was produced. This R code is not present in the repository. Please contact me as above if you would like to see this code and the intermediate vcf files outputted by [freebayes](https://github.com/ekg/freebayes).
- The pipeline for processing amplicon 16S sequences is available [here](https://github.com/padpadpadpad/AB_dada2_pipeline_R).

### Running the scripts and analyses

- The project can be `cloned` or for those not familiar with GitHub, a zip file of this project can be downloaded using the "Clone or download" button at the top right of this page.
- Open the R project file in the downloaded folder. [R projects](https://support.rstudio.com/hc/en-us/articles/200526207-Using-Projects) automatically assigns the root directory to the directory in which the project resides. Consequently all of the analyses should be runnable without altering paths. These are very easy to open using RStudio. All of the scripts for the analyses can be found in `scripts/`.
- `resource_use_analysis.R` does the analyses of the biolog plates, recreating Figure 2 and Figure 3 as well as elements of the Supplementary Information.
- `phyloseq_trees.R` cleans the raw phyloseq object. It creates `figures/phyloseq_trees.pdf` which does not appear in the manuscript. Saves out the phyloseq object which has removed any instance where phylum is `NA`.
- `amplicon_16S_analysis.R` runs the analyses of the community 16S sequencing. Recreates Figure 4, Figure 5, and Figure 6 as well as elements of the Supplementary Information.
- `alpha_diversity_16S_analysis.R` analyses the alpha diversity of the community 16S sequencing. Recreates Figure S5.
- `final_abundance_analysis.R` runs the code to analyse final abundance of the focal species, _Pseudomonas fluorescens_. Recreates Figure 7.
- `extra_functions.R` contains functions for grabbing data out of the output of `vegan::betadisper` for plotting in `ggplot2`.
- All of the data needed to run the analyses are stored in `data`.
- All figures produced by the analysis are saved in `figures/` and are labelled as they are in the manuscript and Supplementary Information.

__All analyses are done in R version 4.0.2, on macOS High Sierra 10.15.6. I am unsure whether some older version of R will support all of the packages and whether the analyses will run exactly the same.__

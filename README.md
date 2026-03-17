# ds4EEB_RLP_hydra_gene

# Genomic Predictors of Duplication Mode in *Hydra vulgaris* AEP

This repository contains all data, scripts, and Quarto documents for the
analysis of genomic features predicting tandem versus dispersed gene
duplication in *Hydra vulgaris* AEP for the completion of final project of Data Science in Ecology and Evolutionary biology in the Winter quarter of 2026.

**Motivating question:** Do gene length, local gene density, distance to
the nearest gene, and scaffold size predict whether a duplicate gene
pair arose through tandem versus dispersed duplication? Are opsin genes
— mediators of extraocular phototransduction in *Hydra* — preferentially
tandemly duplicated?

**Key finding:** Gene length and genomic isolation significantly predict
tandem duplication 

## **Target Audience**
This analysis is is written for evolutionary and computational biologists familiar with logistic regression, Bayesian inference, and genomic feature analysis.


## Data Sources

Genome files must be downloaded separately due to size. Scripts will
prompt with download instructions if files are missing.

`HVAEP.GeneModels.gff3.gz` | <https://research.nhgri.nih.gov/HydraAEP/download/coordinates/hv_aep/> |
| `HVAEP.genome.fa.gz` | <https://research.nhgri.nih.gov/HydraAEP/download/sequences/hv_aep/> |
| `AEP.txt` | Lab annotation (Macias-Muñoz lab, UCSC) |
 `dup_pairsAEP.tsv` | Lab annotation (Macias-Muñoz lab, UCSC) |
 
 
## Reproducing the Analysis

### Requirements

``` r
# R >= 4.3.0
# Required packages:
pkgs <- c("tidyverse", "rtracklayer", "Biostrings", "GenomicRanges",
          "rstanarm", "bayesplot", "tidybayes", "pROC", "caret",
          "ggeffects", "ggridges", "patchwork", "broom",
          "scales", "viridis", "knitr", "kableExtra")
install.packages(pkgs)
# Bioconductor packages:
BiocManager::install(c("rtracklayer","Biostrings","GenomicRanges"))
```

## Future Directions

This analysis is limited to *H. vulgaris* AEP (n = 172 opsin pairs).
Planned extensions:

-   **H. vulgaris strain 105** — independent replication within species
-   **H. oligactis** — cross-species conservation test
-   Mixed-effects model: `dup_type ~ predictors + (1|species)`


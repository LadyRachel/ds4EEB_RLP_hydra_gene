# ds4EEB_RLP_hydra_gene

Which (any?) genomic features predict tandem vs dispersed duplication?

Predictor variables:
1-distance to nearest gene
2-chromosome/scaffold size
3-gene density region
4-gene identity (opsin vs. non-opsin)
  4.5-gene_length

NEED: Gene-level feature table

THEN: Run model dup_type~distance_to_gene+ gene_density+ gene_length+ chrom_size+ gene_ID

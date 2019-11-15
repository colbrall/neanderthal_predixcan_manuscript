## Code and data for "Inferred divergent gene regulation in archaic hominins reveals potential phenotypic differences"
https://doi.org/10.1038/s41559-019-0996-x

We anticipate that the most common quest that would bring people here is searching for the list of **DR genes** dicussed 
in the majority of the paper. To get the full list of DR genes in the Altai Neanderthal discussed in figures 2-4, look in 
`./data/altai_original_pvalues_2sided.txt`. This is a gene x tissue matrix containing empirical p-values calculated vs. 
1000 Genomes. A gene is DR in a particular tissue if $p=0$, and "NA" indicates that the gene was not succesfully modeled in 
that tissue. If you want specifically the **DR GWARRs** do an antijoin on gene id between this file and 
`./data/altai_intro_genes.txt`, then another antijoin with `./data/altai_original_missingModels.txt`; GWARRs are all 
genes that are not introgressed and are not missing all SNPs.

If you want DR genes to compare between archaic hominins, look in `./data/altai_update_pvalues_2sided.txt`, 
`./data/denisovan_update_pvalues_2sided.txt`, and `./data/vindija_pvalues_2sided.txt` and again identify genes for 
which $p=0$ in a given tissue.

Further description of files in this repo:
### scripts/
`comp_pop.jl` \
	calculates empirical p-values, runs PCA and clustering analyses
	requires population predictions and archaic predictions or pre-computed p-values

`comp_inds.jl` \
	comparisons of predicted expression between archaic individuals
	requires archaic predictions

`gene_overlap.jl`\
	intersects gene models and SNPs or regions of interest
	requires lists of regions of interest

`human_dr.jl`\
	for a set of human IDs, calculates empirical P-value and call DR genes

`intro_individuals.jl`\
	for each gene, identifies individuals with certain alleles within 1Mb

`metric_summary_by_gene.jl`\
	generalized version of recombination_rate_map.jl; compares a metric (i.e. gene density) across sets of genes

`model_weights.jl`\
	stats for comparing model weights between sets of SNPs
	requires model databases, SNP sets

`neanderthal_predixcan_plots.ipynb`\
	Jupyter notebook for plotting

`PrediXcan.py`\
	code for performing PrediXcan predictions on genetic information. 
	adapted from script downloaded [here](https://github.com/hakyimlab/PrediXcan/tree/master/Software) on May 19, 2017

`recombination_rate_map.jl`\
	calculates recombination rate and summarizes for gene regions, then compares R.R. across gene sets.

### data/
`*predExp.txt`\
	archaic predicted expression files-- genes x tissue matrix. altai_original indicated original genome called 
	for Altai neanderthal, used for most analyses. altai_update is the updated one used for comparisons to 
	other hominins.

`*pvalues_2sided.txt`\
	archaic 1kG empirical p-values-- genes x tissue matrix. matched to predExp file. A gene is DR in a tissue
	if this p-value equals zero. NA indicates that the gene was not successfully modelled in that tissue.

`altai_disgenet.txt`\
	Disgenet enrichments for DR GWARRS. corresponds to fig. 2d.

`altai_dr_phewas.txt`\
	phenotype associations with differential regulation of DR GWARRs. Corresponds to figure 3 and Table 1.

`altai_desert_genes.txt`\
	genes with at least 1 SNP in any of their tissue models that overlaps an introgression desert. Ensembl ids.

`altai_intro_genes.txt`\
	genes with at least 1 SNP in any of their tissue models that is on an introgressed haplotype (ie are non-GWARRs). 
	Anti-join this list and the genes in predixcan_model_summary.txt to get all GWARRs.

`altai_intro_genes_names.txt`\
	human-readable names of the genes in the file above.

`altai_original_nonIntro_DR_genes.txt`\
	ensembl ids for all DR GWARRs. locations plotted in figure 2b.

`altai_original_missingModels.txt`\
	all genes for which the Altai genome was missing all SNPs in at least one tissue model. These were excluded 
	from analyses.

`gene-dersnp_intersection.txt`\
	genes by ensembl ID and counts of archaic-derived SNPs within 1Mb of them

`gene_region_neanSNPs.txt.gz`\
	for each gene, list of introgressed SNPs (by rsID) within 1 Mb (ie the region considered by PrediXcan). 

`DR_crossSpeces_GO.tgz`\
	cross-species GO analyses- overlap of DR genes. Corresponds to figure 5b.

`diff_crossSpeces_GO.tgz`\
	cross-species GO analyses- absolute difference in regulation. corresponds to figure 5c.

`gencode.v12.annotation.bed.gz`\
	hg19 gene locations. A gene's location is considered to be all bases covered by at least 1 version of the gene.

`predixcan_model_summary.txt`\
	summary of all predixcan models in all tissues. models described by ensembl id, and includes model performance 
	(R2 for observed/predicted in training), number of snps in the model, number of those SNPs that were 
	missing in Altai_original, number of introgressed SNPs in the model, number of SNPs shared between AMH and 
	Altai, and the empirical p-value for divergent regulation (same as in altai_original_pvalues_2sided.txt)

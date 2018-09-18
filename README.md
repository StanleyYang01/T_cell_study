# T cell study
To study age, sex and region specific profile of T cell by RNA-seq

## Experimental design
1. Samples
 	- Aged male and female T cells with or without PMA stimulation were isolated from both brain and spleen
	- Young male T cells without PMA stimulation was isolated from brain
1. RNA from the T cells were extracted and exome enriched 


## Data pre-processing: 
1. fastq files of RNA sequencing was processed using standard JAX in-house pipeline (mrpe) and raw transcript count was generated at the end of the pipeline. 

## Data Analysis workflow
1. Raw transcript counts were converted into counts per million (CPM)
2. Used standard *edgeR*'s quasi-likelihood pipeline to determine differentially expressed (DE) (q<0.05) genes in different comparisons, dissecting sex, age, stimulation and region effect on T cell profiles. Comparisons documented in Rmd files

## Guidance for reviewing analysis
1. Quality control and generating DE genes in all group comparisons using *edgeR* check:
	- "scripts/T_cell_all.Rmd" file
	- "scripts/T_cell_all.pdf" file
	- "scripts/T_cell_all.html" file
1. The resulting of gene list and brief KEGG pathway overview  in each comparison from Rmd file, check: 
	- "analysis" folder for lists of genes with q<0.05
	- "analysis_unfiltered" folder for all the genes with q-value calculated but not filtered out 
1. RNA seq QC figures, check:
	- "figure" folder
1. Export files from IPA analysis, check: 
	- "IPA" folder

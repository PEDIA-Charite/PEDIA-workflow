# Project PEDIA

This folder contains the project data of the PEDIA study

## Description of the subfolders


* Annotate_PhenomizeBoca
	? Scripts for annotation of phenomizer/Boqua ?
* cases
	All available cases with pathogenic mutations in a multi VCF (combined cases: with variant file and only pathigenic mutations) 
* data
	Databases needed for annotation, like ExAC etc..

* exome_cases
	Project folder of cases with real exomes from collaborators. See Snakefile for the pipeline
	** annotated_vcfs
		VCFs annotated with ExAC, CADD etc...
	** vcfs
		Original VCFs from collaborators (but bgzipped and indexed)

* jannovar
	jannovar directory with binaries and databases. also avaiable (some version ) with the module system
	** data
	Transcript DBs (jannovar 0.19, 0.20 and 0.21-SNAPSHOT)

* json_cases
	JSON files of original and phenomized jsons

* simulation
	Prject folder of simulation with 1KG data. See Snakefile for the pipeline and config.json file for the used samples
	** 1kG_background
		1kG genome data annotated
	** json_simulation
		Final jsons where cadd scores of the simulations are appended. (folder json_cases, using phenomized_*.json together with vcf_simulation vcfs).
	** mutations
		multivcf of mutations (duplicated file to cases folder). And also annotated mVCF file.
	** OMIM
		OMIM genes
	** simulator
		The simulator program
	** testsamples
		just for testing
	** vcf_simulation
		Simulation of VCFs using one sample of the mVCF and spike it in one sample of the 1KG.
		*** IranianTrios: 100 Irian trios, the three sample columns are always index, mother, father. 
		    The index has an ID, but the parents are always unaffected. 
		    The ID in approximately half of the patients is due to AR mutations.
		    The other half is probably due to de novo mutations.
   		    This means that there will be in average more pathogenic alleles in this cohort than in the 1KGP data.
		    
 


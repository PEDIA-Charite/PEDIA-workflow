Project folder of simulation with 1KG data. See Snakefile for the pipeline and config.json file for the used samples
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

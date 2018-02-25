# Old jSON format structure

This file describes the old json format utilized by FDNA. This format should
work with the previous iteration of the pipeline.

## case_id

<numeric string>

Unique identifier in face2gene.

## submitter

[<Team stirng>, <user email>]

Information on the submitter. This is a list structure, not a dictionary.

## genomicData

Information on the described mutations.

### test information

{Mutation Type, Molecular Test, Gene Name, Genotype, Notation}

### Mutations

Description of the specific mutation. Contains two sub-fields.

#### Mutation 1

{HGVS-code, Inheritance Mode, result, Additional info}

#### Mutation 2

{HGVS-code, Inheritance Mode, result, Additional info}

## features

[<HPO term strings>]

List of HPO terms.

## geneList

{gene_id, gene_symbol, each score, syndrome_id}

List of genes with annotated scores.

## vcf

{original_filename, vcf_filename, ...}

VCF information including urls for accessing these.

## ranks

{OMIM_phenotype_ID, combined_rank, phenix_rank, PEDIA_rank}

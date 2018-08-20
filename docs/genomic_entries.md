# File structure of genomic entries

Genomic entries are references by the cases on the basis of the entry id.

## entry_id

Entry identifier, should be unique.

## gene

{'gene_id', 'gene_symbol', 'gene_omim_id'}

gene_id - Entrez id

gene_symbol - Alphanumeric identifier of the gene.

gene_omim_id - Gene omim id. Do not mistake this with phenotypic omim ids. (They have a n:n mapping.)

## result

Result of genetic diagnosis. Might be one of the following states.

```
'ABNORMAL'
'ABNORMAL_DIAGNOSTIC'
'DELETION_DUPLICATION'
'NORMAL'
'NORMAL_FEMALE'
'NORMAL_MALE'
'NO_DELETION_DUPLICATION'
'NO_SIGNIFICANT_VARIANTS'
'VARIANTS_DETECTED'
```

## test_type

Type of test performed, might be one of the following.

```
'CHROMOSOMAL_MICROARRAY'
'EXOME_SEQUENCING'
'FISH'
'KARYOTYPE'
'METHYLATION_TESTING'
'MULTIGENE_PANEL'
'OTHER'
'SINGLE_GENE_DELETION_DUPLICATION_TESTING'
'SINGLE_GENE_SEQUENCING'
'TARGETED_TESTING'
'WHOLE_GENE_SEQUENCING'
```


## variant_type

```
'EXON_DELETION'
'EXON_DUPLICATION'
'OTHER'
'SEQUENCE_CHANGE'
```

## variants

```
'chromosome'
'coordinates'
'finding'
	SECONDARY, PRIMARY
'gene'
'hgvs_variant_description' - should contain HGVS string but will almost never be correctly formatted
'interpretation'
	ABNORMAL, LIKELY_PATHOGENIC, PATHOGENIC, UNCERTAIN_SIGNIFICANCE
'mosaic'
	False | True
'mutation'
'mutation1'
'mutation2'
'notes', - some hgvs might even be found here
'rearrangement'
'variant_information', - type of variant
	CDNA_LEVEL, GENOMIC_DNA_LEVEL, PROTEIN_LEVEL, RS_NUMBER
'variant_type'
	TRANSLOCATION, DELETION, DUPLICATION, OTHER
'zygosity'
	COMPOUND_HETEROZYGOUS, HEMIZYGOUS, HETEROZYGOUS, HOMOZYGOUS
```

Beware that this field might also be ''. Always test type or check for empty string.

Do not take any of these fields for granted. They might or might not be there.

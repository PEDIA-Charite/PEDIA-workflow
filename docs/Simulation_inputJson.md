
# Simulation input Json format
\
This is the desired Json format for the Simulation step input files
\
**case_id**\
\
Unique id identifying the case on Face2Gene.\

**submitter\{\}**\
\
\{user_email, user_name, user_team}


**geneList[]**\
\
     [has_mask, feature_score, gene_symbol, combined_score, syndrome_name,  gestalt_score, gene_id, gene_omim_id]\
\
     &gt;&gt; gene information added after mapping gene to syndromes\

**features[]**\
\
&gt;&gt; List of phenotypic features specified using HPO terminology.



**vcf []**\
\
[file_url, vcf_filename]
\
**processing []**\
\
List of the scripts used to produce output file and background_sample\
\
**genomicData[]**\
\
[ Mutations \{result, Build, HGVS-code, additional info, Inheritance Mode\}\
\
, Test Information \{Gene Name, Genotype, Notation, Mutation Type, Molecular Test\}\
\
&gt;&gt; information extracted from genomic_entries in preprocessing step

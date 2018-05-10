
# Final results jSON format 
\
This is the desired Json format for the final results.
\
**submitter\{\}**\
\
\{user_email, user_name, user_team}\
\
**documents[]**\
\
[{document_name,is_vcf}]\
\
&gt;&gt; Additional information on the case, which might be provided by the user.\
\
&gt;&gt; Can be a provided vcf  file.\
\
**geneList[]**\
\
     [has_mask, feature_score, gene_symbol, combined_score, syndrome_name,  gestalt_score, gene_id, gene_omim_id]\
\
     &gt;&gt; gene information added after mapping gene to syndromes\
\
     &gt;&gt; Can contain extra information from  Simulation step : [cadd_phred_score, pheno_score, boqa_score, cadd_raw_score]\
\
**selected_syndromes[]**\
\
[has_mask, omim_id, syndrome_name]\
\
&gt;&gt;Syndrome as selected by the user in the Face2Gene interface.\
\
**detected_syndromes[]**\
\
[combined_score, feature_score, gestalt_score, has_mask, omim_id:|[], syndrome_name]\
\
&gt;&gt; Syndromes detected by the Face2Gene algorithm based on the provided phenotypic information as specified in features and the provided image.\
\
**features[]**\
\
&gt;&gt; List of phenotypic features specified using HPO terminology.\
\
**ranks []**\
\
[\{feature_score, value_pheno, disease-name_pheno,omim_id, combined_score , value_boqa, disease-name_boqa, syndrome_name , gene-id , gestalt_score , confirmed ,  gene-symbol}]\
\
&gt;&gt; Ranks generated after Phenomization\
\
**algo_deploy_version**\
\
     Algorithm version used in the generation of the jSON File. This might affect the detected syndromes and masks. This might change if the dumps are updated for single cases.\
\
**case_id**\
\
Unique id identifying the case on Face2Gene.\
\
**vcf []**\
\
Name of Vcf file when available\
\
**processing []**\
\
List of the scripts used to produce output file and background_sample\
\
**genomic_entries []**\
\
[&lt;entry_filename&gt;]\
\
References the filename of the genomic entry for the case. There might be more than one genomic entry per case.\
\
**genomicData[]**\
\
[ Mutations \{result, Build, HGVS-code, additional info, Inheritance Mode\}\
\
, Test Information \{Gene Name, Genotype, Notation, Mutation Type, Molecular Test\}\
\
&gt;&gt; information extracted from genomic_entries in preprocessing step

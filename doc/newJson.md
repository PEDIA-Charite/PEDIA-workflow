#JSON Data format

Case information have been provided by Face2Gene in several incompatible jSON formats.

This file will describe the most recent revision of the jSON provided by Face2Gene as of February 2018.

Beware that the entire information of a single case is now split in the original case information and separately stored mutation information (only filename reference to this data type).

## Main file fields

###algo_deploy_version

Algorithm version used in the generation of the jSON File. This might affect the detected syndromes and masks. This might change if the dumps are updated for single cases.

###case_id

Unique id identifying the case on Face2Gene.

###detected_syndromes

[combined_score, feature_score, gestalt_score, has_mask, omim_id:<ID>|[<ID>], syndrome_name]

Syndromes detected by the Face2Gene algorithm based on the provided phenotypic information as specified in features and the provided image.

###documents

[{document_name,is_vcf}]

Additional information on the case, which might be provided by the user.

Can be a provided vcs file.

###features

[<HPO features>]

List of phenotypic features specified using HPO terminology.

###genomic_entries

[<entry_filename>]

References the filename of the genomic entry for the case. There might be more than one genomic entry per case.

###selected_syndromes

[has_mask, omim_id, syndrome_name]

Syndrome as selected by the user in the Face2Gene interface.

###submitter

user_email, user_name, user_team

configfile: "config_gestalt.yml"

vcf_samples = config['VCF_SAMPLES']
single_samples = config['SINGLE_SAMPLES']
samples = []
samples.extend(single_samples)
samples.extend(vcf_samples)

preproc_dir = "preproc/"
merged_cases_dir = preproc_dir + "merged_dir/cases/"
pedia_dir = "data/PEDIA/"
converted_json_dir = pedia_dir + "jsons/phenomized/"
mutation_dir = pedia_dir + "mutations/"
vcf_dir = pedia_dir + "vcfs/original/"


rule all:
	input:
		cases = expand(converted_json_dir + "{case}.json", case = samples)

# orginal_json/ contains meta data
# final_json_per_case_280618 contains the gestalt and FM scores from data freeze
# We merge them and save to merged_dir/cases/ dir
rule merge:
	input:
		cases = expand(preproc_dir + "original_json/{case}.json", case = samples),
		syn = preproc_dir + "216_gestalt_syn_to_omim_final.csv",
		f2g = expand(preproc_dir + "final_json_per_case_280618/{case}.json", case = samples),
		list = preproc_dir + "config_gestalt.csv"
	output:
		cases = expand(merged_cases_dir + "{case}.json", case = samples)
	params:
		output = merged_cases_dir,
		input_dir = preproc_dir + "original_json/",
		f2g = preproc_dir + "final_json_per_case_280618"
	shell:
		"""
		python merge_pedia_dgfm.py -l {input.list} -o {params.output} -s {input.syn} -f {params.f2g} -i {params.input_dir}
		"""

rule convert:
	input:
		cases = merged_cases_dir + "{case}.json"
	output:
		cases = converted_json_dir + "{case}.json"
	params:
		output = converted_json_dir
	shell:
		"""
		python preprocessing.py -s {input.cases} -o {params.output} 
		"""


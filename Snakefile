classify_file = 'classifier/pedia.py'
mapping_file = 'classifier/lib/mapping.py'
mapping_vcf_file = 'classifier/lib/get_variant.py'

if 'sample_index' in config:
	sample_index = config['sample_index']
else:
	sample_index = '0'

if 'data_path' in config:
	data_path = config['data_path']
else:
	data_path = 'data'

if 'train_pickle' in config:
	train_pickle = config['train_pickle']
else:
	train_pickle = 'data/train/train_v1.2.p'

if 'param_c' in config:
	param_c = '--param-c ' + str(config['param_c'])
else:
	param_c = ''

# Check if we use phenomizer
if 'use_pheno' in config:
	if config['use_pheno']:
		exclude_pheno = ''
	else:
		exclude_pheno = '--exclude 3_4'
else:
	exclude_pheno = '--exclude 3_4'

rule decompress:
    input:
        "{output}/vcfs/original/{sample}.vcf.gz"
    output:
        temp("{output}/vcfs/original/{sample}.vcf")
    shell:
        """
        bgzip -d -c {input} | grep -v "##sgmutationstatistics=" | awk '{{gsub(/chr/,""); print}}' | awk '{{if($1!="M" && $5!=".") print $0}}' > {output}
        """

rule sort:
    input:
        "{output}/vcfs/original/{sample}.vcf"
    output:
        temp("{output}/vcfs/sorted/{sample}.vcf.gz")
    log: "{output}/logs/{sample}/sort.log"
    shell:
        """
        cat {input} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | "sort -k1,1V -k2,2n"}}' | bgzip -c > '{output}'
        """

rule index_sorted:
    input:
        "{output}/vcfs/sorted/{sample}.vcf.gz"
    output:
        temp("{output}/vcfs/sorted/{sample}.vcf.gz.tbi")
    shell:
        "tabix '{input}'"

rule filter:
    input:
        "{output}/vcfs/sorted/{sample}.vcf.gz"
    output:
        "{output}/vcfs/filtered_vcfs/{sample}.vcf.gz"
    log: "{output}/logs/{sample}/filter.log"
    shell:
        """
        zcat {input} | sed -e 's/nan/NaN/g' | vcftools --vcf - --bed {data_path}/referenceGenome/data/ncbi_refseq_exon_extend_100bp.bed --stdout --recode | bcftools view -i 'GT!~"\."' - | bcftools view -e 'QUAL<100' - | bgzip -c > {output} 2>&1 | tee {log}
       	"""

rule index_filter:
    input:
        "{output}/vcfs/filtered_vcfs/{sample}.vcf.gz"
    output:
        "{output}/vcfs/filtered_vcfs/{sample}.vcf.gz.tbi"
    shell:
        """
        tabix {input}
        """

rule annotate:
    input:
        vcf="{output}/vcfs/filtered_vcfs/{sample}.vcf.gz",
        vcf_index="{output}/vcfs/filtered_vcfs/{sample}.vcf.gz.tbi",
        db="{}/jannovar/data/hg19_refseq.ser".format(data_path),
        exac="{}/populationDBs/ExAC.r1.sites.vep.vcf.gz".format(data_path),
        kg="{}/populationDBs/1KG_ncbi_exon.vcf.gz".format(data_path),
        uk="{}/populationDBs/UK10K_ncbi_exon.vcf.gz".format(data_path),
        caddsnv="{}/pathogenicityScores/cadd_exon_snv.v1.4.tsv.gz".format(data_path),
        caddindel="{}/pathogenicityScores/cadd_exon_indel.v1.4.tsv.gz".format(data_path),
        ref="{}/referenceGenome/data/human_g1k_v37.fasta".format(data_path)
    output:
        "{output}/vcfs/annotated_vcfs/{sample}_annotated.vcf.gz"
    log: "{output}/logs/{sample}/annotation.log"
    shell:
        "java -jar -Xmx3g {data_path}/jannovar/jannovar-cli-0.21-SNAPSHOT.jar annotate-vcf -d {input.db} --exac-vcf {input.exac} --uk10k-vcf {input.uk} --1kg-vcf {input.kg} --tabix {input.caddsnv} {input.caddindel} --tabix-prefix CADD_SNV_ CADD_INDEL_ --ref-fasta {input.ref} -o '{output}' -i '{input.vcf}' 2>&1 | tee {log}"

rule index_annotated:
    input:
        "{output}/vcfs/annotated_vcfs/{sample}_annotated.vcf.gz"
    output:
        "{output}/vcfs/annotated_vcfs/{sample}_annotated.vcf.gz.tbi"
    shell:
        "tabix '{input}'"

rule json:
    input:
        vcf="{output}/vcfs/annotated_vcfs/{sample}_annotated.vcf.gz",
        vcf_index="{output}/vcfs/annotated_vcfs/{sample}_annotated.vcf.gz.tbi",
        omim="{}/omim/genemap2.txt".format(data_path),
        json="{output}/jsons/phenomized/{sample}.json",
        simulator="3_simulation/simulator/pedia-simulator-0.0.3-SNAPSHOT-jar-with-dependencies.jar"
    output:
        "{output}/jsons/test/{sample}.json"
    log: "{output}/logs/{sample}/extend_json.log"
    shell:
        """
        java -jar -Xmx20g {input.simulator} extendjson \
        -j {input.json} -v {input.vcf} -o {input.omim} -out {output} -s {sample_index} 2>&1 | tee {log}
        """

rule test:
    input:
        json = "{output}/jsons/test/{sample}.json"
    output:
        csv = "{output}/results/{sample}/{sample}.csv"
    params:
        label = "1KG",
        dir = "{output}/results/{sample}/",
        train = "3_simulation/jsons/1KG/CV"
    log: "{output}/logs/{sample}/classification.log"
    shell:
        """
        python {classify_file} '{params.train}' '{params.label}' -t {input.json} -o '{params.dir}' {param_c} --train-pickle {train_pickle} {exclude_pheno} 2>&1 | tee {log}
        """

rule map_pedia:
    input:
        csv = "{output}/results/{sample}/{sample}.csv",
        json = "{output}/jsons/test/{sample}.json"
    output:
        json = "{output}/results/{sample}/{sample}_pedia.json",
    params:
        dir = "{output}/results/{sample}/",
    log: "{output}/logs/{sample}/map_pedia.log"
    shell:
        """
        python {mapping_file} --input '{input.json}' --pedia '{input.csv}' --output '{output.json}' 2>&1 | tee {log}
        """

rule map_vcf:
    input:
        csv = "{output}/results/{sample}/{sample}.csv",
        vcf = "{output}/vcfs/annotated_vcfs/{sample}_annotated.vcf.gz"
    output:
        vcf = "{output}/results/{sample}/{sample}.vcf.gz",
    params:
        dir = "{output}/results/{sample}/",
        sample_index=sample_index
    log: "{output}/logs/{sample}/map_vcf.log"
    shell:
        """
        python {mapping_vcf_file} --input '{input.vcf}' --pedia '{input.csv}' --output '{output.vcf}' --sample-index {params.sample_index} 2>&1 | tee {log}
        """

rule map:
    input:
        vcf = "{output}/results/{sample}/{sample}.vcf.gz",
        json = "{output}/results/{sample}/{sample}_pedia.json"
    output:
        out = touch("{output}/results/{sample}/run.out")

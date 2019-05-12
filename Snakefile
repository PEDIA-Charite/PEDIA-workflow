classify_file = 'classifier/pedia.py'
mapping_file = 'classifier/lib/mapping.py'
mapping_vcf_file = 'classifier/lib/get_variant.py'

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
        (cat '{input}' | egrep "^#"; cat '{input}' | egrep -v "^#" | sort -k1,1 -k2,2n) | bgzip -c > '{output}'
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
    params:
        exon = "data/referenceGenome/data/ncbi_refseq_exon_extend_100bp.bed"
    log: "{output}/logs/{sample}/filter.log"
    shell:
        """
        bcftools view -e 'QUAL<100||GT="./."||GT="0/0"||GT=".|."||GT="."' {input} | sed -e 's/nan/NaN/g' | vcftools --vcf - --bed {params.exon} --stdout --recode --recode-INFO-all | bgzip -c > {output} 2> {log}
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
        db="data/jannovar/data/hg19_refseq.ser",
        exac="data/populationDBs/ExAC.r1.sites.vep.vcf.gz",
        uk="data/populationDBs/UK10K_COHORT.20160215.sites.vcf.gz",
        #dbsnp="data/dbSNPs/b147/All_20160601.vcf.gz",
        caddsnv="data/pathogenicityScores/cadd_exon_snv.v1.3.tsv.gz",
        caddindel="data/pathogenicityScores/cadd_exon_indel.v1.3.tsv.gz",
        ref="data/referenceGenome/data/human_g1k_v37.fasta"
    output:
        "{output}/vcfs/annotated_vcfs/{sample}_annotated.vcf.gz"
    log: "{output}/logs/{sample}/annotation.log"
    shell:
        "java -jar -Xmx3g data/jannovar/jannovar-cli-0.21-SNAPSHOT.jar annotate-vcf -d {input.db} --exac-vcf {input.exac} --uk10k-vcf {input.uk} --tabix {input.caddsnv} {input.caddindel} --tabix-prefix CADD_SNV_ CADD_INDEL_ --ref-fasta {input.ref} -o '{output}' -i '{input.vcf}' 2> {log}"

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
        omim="data/omim/genemap2.txt",
        json="{output}/jsons/phenomized/{sample}.json",
        simulator="3_simulation/simulator/pedia-simulator-0.0.2-SNAPSHOT-jar-with-dependencies.jar"
    output:
        "{output}/jsons/test/{sample}.json"
    log: "{output}/logs/{sample}/extend_json.log"
    shell:
        """
        java -jar -Xmx20g {input.simulator} extendjson \
        -j {input.json} -v {input.vcf} -o {input.omim} -out {output} 2> {log}
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
        python {classify_file} '{params.train}' '{params.label}' -t {input.json} -o '{params.dir}' 2> {log}
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
        python {mapping_file} --input '{input.json}' --pedia '{input.csv}' --output '{output.json}' 2> {log}
        """

rule map_vcf:
    input:
        csv = "{output}/results/{sample}/{sample}.csv",
        vcf = "{output}/vcfs/annotated_vcfs/{sample}_annotated.vcf.gz"
    output:
        vcf = "{output}/results/{sample}/{sample}.vcf.gz",
    params:
        dir = "{output}/results/{sample}/"
    log: "{output}/logs/{sample}/map_vcf.log"
    shell:
        """
        python {mapping_vcf_file} --input '{input.vcf}' --pedia '{input.csv}' --output '{output.vcf}' 2> {log}
        """

rule map:
    input:
        vcf = "{output}/results/{sample}/{sample}.vcf.gz",
        json = "{output}/results/{sample}/{sample}_pedia.json"
    output:
        out = touch("{output}/results/{sample}/run.out")

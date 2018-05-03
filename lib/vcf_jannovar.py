'''
Create VCF files from HGVS strings using jannovar.
'''
import os
import subprocess
import csv
import typing

import tempfile
import pandas

from lib.vcf_operations import handle_uncompressed


JANNOVAR_BINARY = "data/jannovar/jannovar_0.25/jannovar-cli-0.25-SNAPSHOT.jar"

REFSEQ_SER = "data/jannovar/jannovar_0.25/data/hg19_refseq.ser"
REF_FASTA = "data/referenceGenome/data/human_g1k_v37.fasta"


VCF_HEADER = (
    '##fileformat=VCFv4.1\n'
    '##INFO=<ID=HGVS,Number=1,Type=String,Description="HGVS-Code">\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    '##contig=<ID=1,assembly=b37,length=249250621>\n'
    '##contig=<ID=2,assembly=b37,length=243199373>\n'
    '##contig=<ID=3,assembly=b37,length=198022430>\n'
    '##contig=<ID=4,assembly=b37,length=191154276>\n'
    '##contig=<ID=5,assembly=b37,length=180915260>\n'
    '##contig=<ID=6,assembly=b37,length=171115067>\n'
    '##contig=<ID=7,assembly=b37,length=159138663>\n'
    '##contig=<ID=8,assembly=b37,length=146364022>\n'
    '##contig=<ID=9,assembly=b37,length=141213431>\n'
    '##contig=<ID=10,assembly=b37,length=135534747>\n'
    '##contig=<ID=11,assembly=b37,length=135006516>\n'
    '##contig=<ID=12,assembly=b37,length=133851895>\n'
    '##contig=<ID=13,assembly=b37,length=115169878>\n'
    '##contig=<ID=14,assembly=b37,length=107349540>\n'
    '##contig=<ID=15,assembly=b37,length=102531392>\n'
    '##contig=<ID=16,assembly=b37,length=90354753>\n'
    '##contig=<ID=17,assembly=b37,length=81195210>\n'
    '##contig=<ID=18,assembly=b37,length=78077248>\n'
    '##contig=<ID=19,assembly=b37,length=59128983>\n'
    '##contig=<ID=20,assembly=b37,length=63025520>\n'
    '##contig=<ID=21,assembly=b37,length=48129895>\n'
    '##contig=<ID=22,assembly=b37,length=51304566>\n'
    '##contig=<ID=X,assembly=b37,length=155270560>\n'
    '##contig=<ID=Y,assembly=b37,length=59373566>\n'
)


def create_vcf(
        variants: [str],
        zygosity: str,
        case_id: str,
        path: str,
) -> typing.Union[str, pandas.DataFrame]:
    '''Generates vcf dataframe. If an error occurs the error message is
    returned.
    '''
    with tempfile.NamedTemporaryFile(mode="w+", dir=path) as hgvsfile:
        hgvsfile.write("\n".join(variants))
        hgvsfile.seek(0)
        with tempfile.NamedTemporaryFile(
            mode="w+", dir=path, suffix=".vcf"
        ) as vcffile:
            try:
                process = subprocess.run(
                    [
                        "java", "-jar", JANNOVAR_BINARY, "hgvs-to-vcf",
                        "-d", REFSEQ_SER,
                        "-r", REF_FASTA,
                        "-i", hgvsfile.name,
                        "-o", vcffile.name,
                    ],
                    check=True,
                    universal_newlines=True,
                    stderr=subprocess.PIPE
                )
            except subprocess.CalledProcessError as error:
                return str(error)
            columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT',
                       'QUAL', 'FILTER', 'INFO', 'FORMAT', case_id]
            hgvs_data = pandas.read_table(
                vcffile.name, sep='\t', comment='#', names=columns
            )
            hgvs_data.ALT.fillna("NA", inplace=True)
            if any(hgvs_data.ALT == '<ERROR>'):
                return process.stderr
            if zygosity.lower() == 'hemizygous':
                genotype = '1'
            elif zygosity.lower() == 'homozygous':
                genotype = '1/1'
            elif (zygosity.lower() == 'heterozygous'
                  or zygosity.lower() == 'compound heterozygous'):
                genotype = '0/1'
            else:
                genotype = '0/1'
            hgvs_data[case_id] = genotype
            hgvs_data['FORMAT'] = 'GT'
            hgvs_data['INFO'] = ['HGVS="{}"'.format(v) for v in variants]
            hgvs_data = hgvs_data.sort_values(by=['#CHROM', "POS"])
            hgvs_data = hgvs_data.drop_duplicates()
            return hgvs_data


def write_vcf(data: pandas.DataFrame, path: str) -> None:
    '''Writes pandas dataframe to vcf file.'''
    # add header to vcf
    with open(path, 'w') as outfile:
        outfile.write(VCF_HEADER)

    # append vcf data
    data.to_csv(
        path, mode='a', sep='\t', index=False,
        header=True, quoting=csv.QUOTE_NONE
    )
    # compress to vcf.gz file
    handle_uncompressed(path, path+'.gz')
    os.remove(path)

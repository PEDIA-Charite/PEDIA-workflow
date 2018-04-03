'''
Case model created from json files.
'''
import logging
from typing import Union, Dict


import pandas
import csv
import subprocess
import tempfile
import os

from lib.model.json import OldJson, NewJson
from lib.vcf_operations import move_vcf
from lib.utils import explode_df_column


LOGGER = logging.getLogger(__name__)


def genes_to_single_cols(rowdata: pandas.Series) -> pandas.Series:
    '''Separate gene dict in genes column into separate columns.
    '''
    gene_dict = rowdata['genes']
    rowdata.drop('genes', inplace=True)
    genes = pandas.Series(gene_dict)
    rowdata = rowdata.append(genes)
    return rowdata


def create_gene_table(rowdata: pandas.Series, omim: 'Omim') -> pandas.Series:
    '''Get gene information from row information.
    This includes: all scores, gene_symbol, gene_id, gene_omim_id, syndrome_id
    '''
    disease_id = rowdata['omim_id']

    syndrome_name = rowdata["syndrome_name"] \
        or rowdata["disease-name_pheno"] \
        or rowdata["disease-name_boqa"]

    # get dict containing gene_id, gene_symbol, gene_omim_id
    genes = list(omim.mim_pheno_to_gene(disease_id).values())
    # get all three scores provided by face2gene
    gestalt_score = rowdata["gestalt_score"]
    feature_score = rowdata['feature_score']
    combined_score = rowdata['combined_score']
    pheno_score = rowdata['value_pheno']
    boqa_score = rowdata['value_boqa']

    resp = pandas.Series({
        "disease_id": disease_id,
        "syndrome_name": syndrome_name,
        "genes": genes,
        "gestalt_score": gestalt_score,
        "feature_score": feature_score,
        "combined_score": combined_score,
        "pheno_score": pheno_score,
        "boqa_score": boqa_score
    })
    return resp


class Case:
    '''
    Exposes the following properties:
    case_id - Unique identifier for the case in Face2Gene
    variants - list of hgvs objects describing valid hgvs variants
    analyzed_syndromes - dictionary of syndromes using omim id as key
    features - list of hpo terms
    diagnosis - list of syndromes selected as diagnosis
    submitter - submitter information containing fields for email, name, team
    realvcf - list of vcf filenames
    '''

    def __init__(self, data: Union[OldJson, NewJson],
                 error_fixer: "ErrorFixer"):
        self.algo_version = data.get_algo_version()
        self.case_id = data.get_case_id()
        # get both the list of hgvs variants and the hgvs models used in the
        # parsing
        self.variants, self.hgvs_models = data.get_variants(error_fixer)
        self.syndromes = data.get_syndrome_suggestions_and_diagnosis()
        self.features = data.get_features()
        self.submitter = data.get_submitter()
        self.realvcf = data.get_vcf()
        self.gene_scores = None
        # also save the json object to easier extract information from the
        # new format
        self.data = data
        LOGGER.debug("Creating case %s", self.case_id)

    def phenomize(self, pheno: 'PhenomizerService') -> bool:
        '''Add phenomization information to genes from boqa and phenomizer.
        Args:
            omim: Omim object to handle id translation.
            pheno: PhenomizerService to handle API calls for phenomizer and
                   boqa.
        '''
        pheno_boqa = pheno.disease_boqa_phenomize(self.features)

        pheno_boqa.index = pheno_boqa.index.astype(int)
        # merge pheno and boqa scores dataframe with our current syndromes
        # dataframe which contains face2gene scores
        self.syndromes = self.syndromes.merge(
            pheno_boqa, left_on='omim_id', how='outer', right_index=True)
        self.syndromes.reset_index(drop=True, inplace=True)

        # fill nans created by merge
        self.syndromes.fillna(
            {
                'combined_score': 0.0,
                'feature_score': 0.0,
                'gestalt_score': 0.0,
                'value_pheno': 0.0,
                'value_boqa': 0.0,
                'syndrome_name': '',
                'confirmed': False,
                'has_mask': 0,
                'gene-symbol': '',
                'gene-id': '',
                'disease-name_boqa': '',
                'disease-name_pheno': '',
                'disease-id_boqa': '',
                'disease-id_pheno': '',
            },
            inplace=True
        )

        LOGGER.debug("Phenomization case %s success", self.case_id)

        return True

    def get_gene_list(self, omim: 'Omim', recreate: bool = False,
                      filter_entrez_id: bool = True) \
            -> Dict[str, str]:
        '''Get list of genes from syndromes. Save them back to self.genes
        for faster lookup afterwards.
        Args:
            filter_entrez_id: Only include genes with existing entrez gene ids.
        '''
        # return existing gene list, if it already exists
        if self.gene_scores is not None and not recreate:
            return self.gene_scores

        LOGGER.debug("Generating geneList for case %s", self.case_id)
        gene_table = self.syndromes.apply(
            lambda x: create_gene_table(x, omim), axis=1)
        # explode the gene table on genes to separate the genetic entries
        gene_table = explode_df_column(gene_table, 'genes')
        gene_table = gene_table.apply(genes_to_single_cols, axis=1)
        # only select entries with non-empty gene ids
        if filter_entrez_id:
            gene_table = gene_table.loc[gene_table["gene_id"] != ""]
        gene_scores = gene_table.to_dict('records')
        self.gene_scores = gene_scores
        return gene_scores

    def check(self) -> bool:
        '''Check whether Case fulfills all provided criteria.

        The criteria are:
            picture has been provided - gestalt_score in detected_syndromes
            should be greater than 0
            clinical diagnosis - selected_syndromes should not be empty
            single monogenetic disease - not multiple syndromes selected and
            not multiple pathogenic mutations in different genes
            SNP mutations - no microdel/dup or other large scale aberrations

        Note that this function is more definitive than the json level check,
        as the validity of hgvs parsing has already been established.
        '''
        valid = True
        issues = []
        # check maximum gestalt score
        max_gestalt_score = max(self.syndromes['gestalt_score'])
        if max_gestalt_score <= 0:
            issues.append('Maximum gestalt score is 0. \
                          Probably no image has been provided.')
            valid = False

        # check that only one syndrome has been selected
        diagnosis = self.syndromes.loc[self.syndromes['confirmed']]
        if len(diagnosis) != 1:
            issues.append(
                '{} syndromes have been selected. Only 1 syndrome should be \
                selected for PEDIA inclusion.'.format(len(diagnosis)))
            valid = False

        # check that molecular information is available at all
        if len(self.variants) == 0:
            issues.append('No valid genomic entries available.')
            valid = False

        return valid, issues

    def eligible_training(self) -> bool:
        '''Eligibility of case for training.
        Exclusion criteria are:
            available real vcf
        '''
        return not self.realvcf

    def create_vcf(self, path: str) -> pandas.DataFrame:
        '''Generates vcf dataframe. If an error occurs the error message is returned.
        '''
        with tempfile.NamedTemporaryFile(mode="w+", dir=path) as hgvsfile:
            for v in self.variants:
                hgvsfile.write(str(v) + "\n")
            hgvsfile.seek(0)
            with tempfile.NamedTemporaryFile(mode="w+", dir=path, suffix=".vcf") as vcffile:
                try:
                    process=subprocess.run(["java", "-jar", 'data/jannovar/jannovar-cli-0.25-SNAPSHOT.jar', "hgvs-to-vcf", "-d",
                                                  'data/jannovar/data/hg19_refseq.ser', "-i", hgvsfile.name, "-o", vcffile.name, "-r", "data/jannovar/data/hg19/hg19.fa"], check=True, universal_newlines=True, stderr=subprocess.PIPE)
                except subprocess.CalledProcessError as e:
                    return str(e)
                columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT',
                           'QUAL', 'FILTER', 'INFO', 'FORMAT', self.case_id]
                df = pandas.read_table(
                    vcffile.name, sep='\t', comment='#', names=columns)
                if any(df.ALT == '<ERROR>'):
                    return(process.stderr)
                if self.hgvs_models[0].zygosity.lower() == 'hemizygous':
                    genotype = '1'
                elif self.hgvs_models[0].zygosity.lower() == 'homozygous':
                    genotype = '1/1'
                elif self.hgvs_models[0].zygosity.lower() == 'heterozygous' or self.hgvs_models[0].zygosity.lower() == 'compound heterozygous':
                    genotype = '0/1'
                else:
                    genotype = '0/1'
                df[self.case_id] = genotype
                df['FORMAT'] = 'GT'
                df['INFO'] = ['HGVS="' + str(v) + '"' for v in self.variants]
                df = df.sort_values(by=['#CHROM', "POS"])
                return df

    def dump_vcf(self, path: str, recreate: bool = False) -> None:
        '''Dumps vcf file to given path. Initializes vcf generation if none has yet been created.
        Created vcf is saved to self.vcf.
        '''
        if hasattr(self, 'vcf') and not recreate:
            if isinstance(self.vcf, str):
                LOGGER.error("VCF generation for case %s failed. Error message:%s",self.case_id,self.vcf)
            else:
                outputpath = os.path.join(path, self.case_id + '.vcf')
                # add header to vcf
                with open(outputpath, 'w') as outfile:
                    outfile.write(
                        '##fileformat=VCFv4.1\n##INFO=<ID=HGVS,Number=1,Type=String,Description="HGVS-Code">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                    outfile.write(
                        '##contig=<ID=1,assembly=b37,length=249250621>\n')
                    outfile.write(
                        '##contig=<ID=2,assembly=b37,length=243199373>\n')
                    outfile.write(
                        '##contig=<ID=3,assembly=b37,length=198022430>\n')
                    outfile.write(
                        '##contig=<ID=4,assembly=b37,length=191154276>\n')
                    outfile.write(
                        '##contig=<ID=5,assembly=b37,length=180915260>\n')
                    outfile.write(
                        '##contig=<ID=6,assembly=b37,length=171115067>\n')
                    outfile.write(
                        '##contig=<ID=7,assembly=b37,length=159138663>\n')
                    outfile.write(
                        '##contig=<ID=8,assembly=b37,length=146364022>\n')
                    outfile.write(
                        '##contig=<ID=9,assembly=b37,length=141213431>\n')
                    outfile.write(
                        '##contig=<ID=10,assembly=b37,length=135534747>\n')
                    outfile.write(
                        '##contig=<ID=11,assembly=b37,length=135006516>\n')
                    outfile.write(
                        '##contig=<ID=12,assembly=b37,length=133851895>\n')
                    outfile.write(
                        '##contig=<ID=13,assembly=b37,length=115169878>\n')
                    outfile.write(
                        '##contig=<ID=14,assembly=b37,length=107349540>\n')
                    outfile.write(
                        '##contig=<ID=15,assembly=b37,length=102531392>\n')
                    outfile.write(
                        '##contig=<ID=16,assembly=b37,length=90354753>\n')
                    outfile.write(
                        '##contig=<ID=17,assembly=b37,length=81195210>\n')
                    outfile.write(
                        '##contig=<ID=18,assembly=b37,length=78077248>\n')
                    outfile.write(
                        '##contig=<ID=19,assembly=b37,length=59128983>\n')
                    outfile.write(
                        '##contig=<ID=20,assembly=b37,length=63025520>\n')
                    outfile.write(
                        '##contig=<ID=21,assembly=b37,length=48129895>\n')
                    outfile.write(
                        '##contig=<ID=22,assembly=b37,length=51304566>\n')
                    outfile.write(
                        '##contig=<ID=X,assembly=b37,length=155270560>\n')
                    outfile.write(
                        '##contig=<ID=Y,assembly=b37,length=59373566>\n')

                self.vcf.to_csv(outputpath, mode='a', sep='\t', index=False,
                                header=True, quoting=csv.QUOTE_NONE)
                move_vcf(outputpath, outputpath + '.gz', 'text')
                os.remove(outputpath)
        else:
            LOGGER.debug("Generating VCF for case %s", self.case_id)
            self.vcf = self.create_vcf(path)
            self.dump_vcf(path)

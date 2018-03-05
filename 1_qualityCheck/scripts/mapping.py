#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@authors: hertzberg,leitheim
"""

import json
import os
import csv  # necessary for creating genedict
import requests
import getopt
import sys
import argparse


# METHODS


def makegenedict():  # creates a dictionary of genemap2.txt containing the omim id of a gene as key and genesymbol,entrez gene id as values
    genedict = {}
    with open("genemap2.txt") as genemap:
        genemap = genemap.readlines()[3:]
        genemap = csv.DictReader(genemap, delimiter="\t")
        for line in genemap:
            genesymbol = str(line["Gene Symbols"]).split(",")[0]
            geneid = line["Entrez Gene ID"]
            genedict[line["Mim Number"]] = genesymbol, geneid
    return genedict


def createentry(features, omimid, genedict):  # creates an entry for the genelist
    entry = features
    try:
        genesymbol = genedict[str(omimid)][0]
        geneid = genedict[str(omimid)][1]
        entry["gene_id"] = geneid
        entry["gene_symbol"] = genesymbol
        entry["gene_omim_id"] = omimid
    except KeyError:
        with open("Exceptions_by_missing_omim_in_genemap2.txt", "a") as exceptions:
            exceptions.write(features["syndrome_name"] + str(omimid) + "\n")
    return(entry)


def getTestInformation(genomicEntry):
    testInformation = {}

    testInformation['Molecular Test'] = getTestInformation_MolecularTest(
        genomicEntry)
    testInformation['Notation'] = getTestInformation_Notation(genomicEntry)
    testInformation['Genotype'] = getTestInformation_Genotype(genomicEntry)
    testInformation['Mutation Type'] = getTestInformation_MutationType(
        genomicEntry)
    testInformation['Gene Name'] = getTestInformation_GeneName(genomicEntry)

    return testInformation


def getMutations(genomicEntry):
    mutations = {}

    mutations['additional info'] = getMutations_AddInfo(genomicEntry)
    mutations['Build'] = getMutations_Build(genomicEntry)
    mutations['result'] = getMutations_Result(genomicEntry)
    mutations['Inheritance Mode'] = getMutations_InheritanceMode(genomicEntry)

    # deal with the possibility of two hgvs codes of differing quality
    primaryHGVSCode, alternativeHGVSCode = getMutations_HGVSCode(genomicEntry)
    mutations['HGVS-code'] = primaryHGVSCode
    if alternativeHGVSCode:
        mutations['alternative HGVS-code'] = alternativeHGVSCode

    return mutations


def getTestInformation_MolecularTest(genomicEntry):
    if 'test_type' in genomicEntry:
        return genomicEntry['test_type']
    return ''


def getTestInformation_Notation(genomicEntry):
    if 'variants' in genomicEntry:
        if 'variant_information' in genomicEntry['variants']:
            return genomicEntry['variants']['variant_information']
    return ''


def getTestInformation_Genotype(genomicEntry):
    if 'variants' in genomicEntry:
        if 'zygosity' in genomicEntry['variants']:
            return genomicEntry['variants']['zygosity']
    return ''


def getTestInformation_MutationType(genomicEntry):
    if 'variants' in genomicEntry:
        if 'notes' in genomicEntry['variants']:
            return findMutationTypeInString(genomicEntry['variants']['notes'])
    return ''


def getTestInformation_GeneName(genomicEntry):
    """ Get the name of the tested gene from a genomic entry.
        The gene symbol (e.g. 'RAB3GAP1') is preferred but if that is not
        given the OMIM-ID is used if present.
    """

    if 'gene' in genomicEntry:
        if 'gene_symbol' in genomicEntry['gene']:
            return genomicEntry['gene']['gene_symbol']
        elif 'gene_omim_id' in genomicEntry['gene']:
            return 'OMIM-ID:' + genomicEntry['gene']['gene_symbol']
        elif 'gene' in genomicEntry['variants']:
            return genomicEntry['variants']['gene']['gene_symbol']

    return ''


def getMutations_AddInfo(genomicEntry):
    infoSnippets = []
    addInfo = ''

    # add notes item
    if 'variants' in genomicEntry:
        if 'notes' in genomicEntry['variants']:
            addInfo = addInfo + genomicEntry['variants']['notes']

    # TODO add other relevant things to info snippets

    # combine info snippets
    for snippet in infoSnippets:
        addInfo = addInfo + ' || ' + snippet

    return addInfo


def getMutations_Build(genomicEntry):
    # This seems to apply only to variants described on the genomic level

    if 'variants' in genomicEntry:
        if 'mutation' in genomicEntry['variants']:
            if 'chromosome' in genomicEntry['variants']['mutation']:
                if 'build' in genomicEntry['variants']['mutation']['chromosome']:
                    return genomicEntry['variants']['mutation']['chromosome']['build']
    return ''


def getMutations_Result(genomicEntry):
    if 'result' in genomicEntry:
        return genomicEntry['result']
    return ''


def getMutations_InheritanceMode(genomicEntry):
    if 'variants' in genomicEntry:
        if 'notes' in genomicEntry['variants']:
            return findInheritanceModeInString(genomicEntry['variants']['notes'])
    return ''


def getHGVSCodeFromField(genomicEntry):
    if 'hgvs_variant_description' in genomicEntry:
        if genomicEntry['hgvs_variant_description']:
            return genomicEntry['hgvs_variant_description']
    if 'variants' in genomicEntry:
        if 'hgvs_variant_description' in genomicEntry['variants']:
            if genomicEntry['variants']['hgvs_variant_description']:
                return genomicEntry['variants']['hgvs_variant_description']

    return ''


def getMutations_HGVSCode(genomicEntry):
    """ Get two hgvs codes from a given genomic entry and order them by quality.
    They are obtained
    1) directly as a string from the item 'hgvs_variant_description'
    2) by combining information from various items describing the mutation
    If both methods lead to a similiar quality (both correct or both not fully
    correct but still meaningful) the second one is preferred.
    Either both returned strings or the one of lesser quality can be empty.
    """

    # check various keys for detailed genetic information and create hgvs key from them
    hgvsCodeCreated = getHGVSCodeFromDetails(genomicEntry)
    # check hgvs variant description item for a hgvs code string
    hgvsString = getHGVSCodeFromField(genomicEntry)

    # assess quality of hgvs codes and deal with them appropriatly
    # in case of similiar quality the created hgvs code is prioritized
    if isCorrectHGVSString(hgvsCodeCreated):
        primaryHGVSCode = hgvsCodeCreated
        alternativeHGVSCode = hgvsString
    elif isCorrectHGVSString(hgvsString):
        primaryHGVSCode = hgvsString
        alternativeHGVSCode = hgvsCodeCreated
    elif hgvsCodeCreated:
        primaryHGVSCode = hgvsCodeCreated
        alternativeHGVSCode = hgvsString
    else:
        # both might be empty
        primaryHGVSCode = hgvsString
        alternativeHGVSCode = hgvsCodeCreated

    return (primaryHGVSCode, alternativeHGVSCode)


def getHGVSCodeFromDetails(genomicEntry):
    """ Create a hgvs code string from various fields in the given genomic entry.

    The input information is not completed or validated in any way. If parts of
    the hgvs code are missing, the incomplete code is returned. If, however, the
    reference sequence is missing and additionally the variant description only
    consists of the level description ('c','g','p') the result is not considered
    meaningful and an empty string is returned instead of e.g. '.c:'.
    """

    if 'test_type' in genomicEntry:
        testType = genomicEntry['test_type']

        # TODO: include all test types
        # TODO: deal with missing test type
        # TODO: deal with missing variant info

        if testType in ['EXOME_SEQUENCING', 'WHOLE_GENE_SEQUENCING',
                        'MULTIGENE_PANEL', 'TARGETED_TESTING',
                        'SINGLE_GENE_SEQUENCING']:

            if 'variants' in genomicEntry:
                if 'variant_information' in genomicEntry['variants']:
                    variantInfo = genomicEntry['variants']['variant_information']

                    # variant type: coding dna level
                    if variantInfo == 'CDNA_LEVEL':
                        referenceSequenceHGVS = findReferenceSequenceDNA(
                            genomicEntry)
                        variantHGVS = findVariantDNA(genomicEntry)
                        if (referenceSequenceHGVS + variantHGVS):
                            return referenceSequenceHGVS + ':c.' + variantHGVS

                    # variant type: genomic dna level
                    elif variantInfo == 'GENOMIC_DNA_LEVEL':
                        referenceSequenceHGVS = findReferenceSequenceGenomic(
                            genomicEntry)
                        variantHGVS = findVariantDNA(genomicEntry)
                        if (referenceSequenceHGVS + variantHGVS):
                            return referenceSequenceHGVS + ':g.' + variantHGVS

                    # variant type: protein level
                    elif variantInfo == 'PROTEIN_LEVEL':
                        referenceSequenceHGVS = findReferenceSequenceProtein(
                            genomicEntry)
                        variantHGVS = findVariantProtein(genomicEntry)
                        if (referenceSequenceHGVS + variantHGVS):
                            return referenceSequenceHGVS + ':p.' + variantHGVS

    return ''


def findReferenceSequenceDNA(genomicEntry):

    # check for correct annotation
    if 'mutation' in genomicEntry['variants']:
        mutation = genomicEntry['variants']['mutation']
        if 'transcript' in mutation:
            return mutation['transcript']

    # if no correct annotation is found, check notes for reference sequences
    if 'notes' in genomicEntry['variants']:
        notesRefSeq = findReferenceSequenceDNAInString(
            genomicEntry['variants']['notes'])
        if notesRefSeq:
            return notesRefSeq

    # if no reference sequence is found, return an empty string
    return ''


def findReferenceSequenceGenomic(genomicEntry):
    if 'mutation' in genomicEntry['variants']:
        if 'chromosome' in genomicEntry['variants']['mutation']:
            if 'number' in genomicEntry['variants']['mutation']['chromosome']:
                return genomicEntry['variants']['mutation']['chromosome']['number']
    return ''


def findReferenceSequenceProtein(genomicEntry):
    # TODO: How could that even be entered through the new f2g interface?
    return ''


def findVariantDNA(genomicEntry):
    # check for correct annotation
    if 'mutation' in genomicEntry['variants']:
        mutation = genomicEntry['variants']['mutation']
        if 'mutation_type' in mutation:
            mutationType = mutation['mutation_type']

            # substitution
            if mutationType == 'SUBSTITUTION':
                location = originalBase = substitutedBase = ''
                if 'location' in mutation:
                    location = mutation['location']
                if 'original_base' in mutation:
                    originalBase = mutation['original_base']
                if 'substituted_base' in mutation:
                    substitutedBase = mutation['substituted_base']
                return location + originalBase + '>' + substitutedBase

            # duplication
            if mutationType == 'DUPLICATION':
                location = ''
                if 'location' in mutation:
                    location = mutation['location']
                    # TODO: remove error removal due to overlap with quality check?
                    # deal with errors such as location=3586dupC by deleting 'dup' and everything after it
                    if location.find('dup') != -1:
                        location = location[0:location.find('dup')]
                return location + 'dup'

            # insertion
            if mutationType == 'INSERTION':
                location = insertedBases = ''
                if 'location' in mutation:
                    location = mutation['location']
                if 'inserted_bases' in mutation:
                    insertedBases = mutation['inserted_bases']
                return location + 'ins' + insertedBases

            # deletion
            if mutationType == 'DELETION':
                location = ''
                if 'location' in mutation:
                    location = mutation['location']
                return location + 'del'

            # deletion insertion
            if mutationType == 'DELETION_INSERTION':
                location = insertedBases = ''
                if 'location' in mutation:
                    location = mutation['location']
                if 'inserted_bases' in mutation:
                    insertedBases = mutation['inserted_bases']
                return location + 'delins' + insertedBases
                # TODO: no mutation type present
    # TODO: no mutation present
    # TODO: check for info in notes

    # return default value
    return ''


def findVariantProtein(genomicEntry):
    # check for correct annotation
    if 'mutation' in genomicEntry['variants']:
        mutation = genomicEntry['variants']['mutation']
        if 'mutation_type' in mutation:
            mutationType = mutation['mutation_type']

            # substitution
            if mutationType == 'SUBSTITUTION':
                location = originalAminoAcid = substitutedAminoAcid = ''
                if 'first_amino_position' in mutation:
                    location = mutation['first_amino_position']
                if 'first_amino_acid' in mutation:
                    originalAminoAcid = mutation['first_amino_acid']
                if 'last_amino_acid' in mutation:
                    substitutedAminoAcid = mutation['last_amino_acid']
                return originalAminoAcid + location + substitutedAminoAcid

            # duplication, insertion, deletion, deletion insertion
                # TODO! There are no sample files (new format) to obtain the dictionary keys from
    # TODO: no mutation present
    # TODO: check for info in notes

    # return default value
    return ''


def findVariantGenomic(genomicEntry):
    # check for correct annotation
    if 'mutation' in genomicEntry['variants']:
        mutation = genomicEntry['variants']['mutation']
        if 'mutation_type' in mutation:
            mutationType = mutation['mutation_type']

            # substitution
            if mutationType == 'SUBSTITUTION':
                location = originalAminoAcid = substitutedAminoAcid = ''
                if 'first_amino_position' in mutation:
                    location = mutation['first_amino_position']
                if 'first_amino_acid' in mutation:
                    originalAminoAcid = mutation['first_amino_acid']
                if 'last_amino_acid' in mutation:
                    substitutedAminoAcid = mutation['last_amino_acid']
                return originalAminoAcid + location + substitutedAminoAcid

            # duplication, insertion, deletion, deletion insertion
                # TODO! There are no sample files (new format) to obtain the dictionary keys from
    # TODO: no mutation present
    # TODO: check for info in notes

    # return default value
    return ''


def findReferenceSequenceDNAInString(string):

    refSeqIdentifiers = ['NC', 'NG', 'NM']

    # if several identifiers are present the first number is returned
    for ident in refSeqIdentifiers:
        index = string.lower().find(ident.lower())
        if index != -1:
            if (string[index + 2] == '_' or string[index + 2] == '-'):
                begNumber = index + 3
            else:
                begNumber = index + 2
            if begNumber + 5 < len(string):
                if string[begNumber:begNumber + 6].isdigit():
                    return string[index:begNumber + 6]

    # return empty string if nothing is found
    return ''


def findMutationTypeInString(string):
    # TODO: maybe there are synonymous descriptions of these three mutation types?
    # TODO: maybe the mutation type is encoded in an item for certain test types ...
    # ... (definitely not the common ones, though)

    # note that mutation types are prioritized in that order
    mutationTypes = ['Epigenetic imprinting defect',
                     'Chromosomal', 'Monogenic']

    for mutationType in mutationTypes:
        if string.lower().find(mutationType.lower()) != -1:
            return mutationType

    return ''


def findInheritanceModeInString(string):
    # TODO: do something more elaborate
    inheritanceModes = ['Unknown', 'Autosomal Recessive', 'Autosomal Dominant (inherited)',
                        'X-Linked Dominant (inherited)', 'X-Linked Recessive (inherited)',
                        'Unknown, Autosomal Dominant (de novo)', 'X-Linked Dominant (de novo)',
                        'Mosaic', 'Autosomal Dominant (de novo)', 'Other Aberration Type',
                        'X-Linked Recessive (de novo)']

    for inheritanceMode in inheritanceModes:
        if string.lower().find(inheritanceMode.lower()) != -1:
            return inheritanceMode

    return ''


def isCorrectHGVSString(hgvsString):
    # TODO replace by hgvs library
    # This is not a proper verification of the HGVS standard
    # Actual behavior: Returns true, if the string contains exactly
    # one colon which separates a valid reference sequence and a valid
    # variant description. Their compatibility is not checked.

    parts = hgvsString.split(':')

    if (len(parts) == 2):
        reference = parts[0]
        variant = parts[1]
        if (isReferenceSequence(reference) and isVariantDescription(variant)):
            return 1
        else:
            return 0
    else:
        return 0


def isReferenceSequence(referenceString):
    refSeqIdents = {'NC', 'NG', 'NM', 'NR', 'NP'}
    # reference sequence has to have the following format: NC_000023.10
    # only reference sequences with version number are regarded as correct

    # The shortest reference sequences are the ones with a one digit
    # version number. They have exactly 11 characters.
    if (len(referenceString) >= 11):
        # Try all identifiers. When an identifier is found, a decision (positive or negative)
        # can be made in the same iteration. If it is not found, others still have to be tried.
        for ident in refSeqIdents:
            if referenceString.startswith(ident + '_'):
                if referenceString[3:9].isdigit() and referenceString[9] == '.' and referenceString[10:].isdigit():
                    return 1
                else:
                    return 0
            else:
                # do nothing and try other identifyers
                pass
        # return 0 because no identifier mathced the beginning of the string
        return 0

    else:
        return 0


def isVariantDescription(variantString):

    # exit if string is too short
    if (len(variantString) <= 2):
        return 0

    start = variantString[0:2]

    # TODO This is horribly incorrect, but works for now.
    if start == 'c.':
        return 1
    elif start == 'p.':
        return 1
    elif start == 'g.':
        return 1
    else:
        return 0


def getGenomicDataCorrectFormat(fileWrongFormat, path):
    genomicData = []

    genomicEntries = fileWrongFormat['genomic_entries']
    if genomicEntries:   # check whether genomicEntries is empty, because otherwise it is of type list which will cause an error in the next line
        for entry in genomicEntries:
            if not os.path.isfile(path + str(entry) + '.json'):
                continue
            genomicEntry = json.load(open(path + str(entry) + '.json'))
            # initialize entry for genomicData list
            genomicDataElement = {}
            # create content
            testInformation = getTestInformation(genomicEntry)
            mutations = getMutations(genomicEntry)
            # put together
            genomicDataElement['Test Information'] = testInformation
            genomicDataElement['Mutations'] = mutations
            genomicData.append(genomicDataElement)

    return genomicData


def inputIsEmpty(file):
    genomicEntries = file['genomic_entries']
    if genomicEntries:   # check whether genomicEntries is empty, because otherwise it is of type list which will cause an error in the next line
        for key, genomicEntry in genomicEntries.items():
            # TODO
            return 1
    else:
        return 1


def hasTestType(file):
    if 'genomic_entries' in file:
        if 'Mutation_1' in file['genomic_entries']:
            if 'test_type' in file['genomic_entries']['Mutation_1']:
                return 1
    return 0


def hasTwoMutations(file):
    try:
        mutation1 = file['genomic_entries']['Mutation_1']['variants']['mutation1']
        return 1
    except:
        return 0


def isTestTypeRelevant(file):
    """ Note that this function assumes the tet type is actually present in file"""
    # TODO This should be called above in getHGVSCodeFromDetails...

    relevantTestTypes = ['EXOME_SEQUENCING', 'WHOLE_GENE_SEQUENCING',
                         'MULTIGENE_PANEL', 'TARGETED_TESTING',
                         'SINGLE_GENE_SEQUENCING']

    testType = file['genomic_entries']['Mutation_1']['test_type']
    if testType in relevantTestTypes:
        return 1
    else:
        return 0


def isVariantsEmpty(file):
    if 'genomic_entries' in file:
        if 'Mutation_1' in file['genomic_entries']:
            if 'variants' in file['genomic_entries']['Mutation_1']:
                if file['genomic_entries']['Mutation_1']['variants']:
                    return 0
    return 1


def evaluateOutput(results, path):
    noGenomicData = []
    correctHGVS = []
    invalidNonemptyHGVS = []
    emptyHGVS = []
    # the following lists are subsets of emptyHGVS
    noTestType = []
    testTypeNotRelevant = []
    variantsEmpty = []
    twoMutations = []
    unknownIssues = []

    for item in results:
        caseID = item[0]

#        if caseID == '143386':
#            print('Here we go...')

        genomicData = item[1]

        # genetic data item present in dictionary?
        if genomicData:
            hgvsCode = genomicData[0]['Mutations']['HGVS-Code']
            file = json.load(open(path + '/' + caseID))

            # HGVS correct?
            if isCorrectHGVSString(hgvsCode):
                correctHGVS.append(caseID)
            # HGVS at least nonempty?
            elif hgvsCode:
                invalidNonemptyHGVS.append(caseID)
            # HGVS empty...
            else:
                emptyHGVS.append(caseID)

                # does test type cause problems?
                if not hasTestType(file):
                    noTestType.append(caseID)
                elif not isTestTypeRelevant(file):
                    testTypeNotRelevant.append(caseID)
                else:
                    # is file basically empty?
                    if isVariantsEmpty(file):
                        variantsEmpty.append(caseID)
                    # does it have two mutations?
                    elif hasTwoMutations(file):
                        twoMutations.append(caseID)
                    # UNKNOWN
                    else:
                        unknownIssues.append(caseID)
        else:
            noGenomicData.append(caseID)

    # outputs
    print('total:' + str(len(results)))
    print('correct:' + str(len(correctHGVS)))
    print('invalid nonempty:' + str(len(invalidNonemptyHGVS)))
    print('no genomic data:' + str(len(noGenomicData)))
    print('empty HGVS:' + str(len(emptyHGVS)))

    print('-- no test type:' + str(len(noTestType)))
    print('-- test type not relevant:' + str(len(testTypeNotRelevant)))
    print('-- variants empty:' + str(len(variantsEmpty)))
    print('-- two mutations:' + str(len(twoMutations)))
    print('-- unknown issues:' + str(len(unknownIssues)))

    return (correctHGVS, invalidNonemptyHGVS, emptyHGVS,
            noGenomicData, noTestType, testTypeNotRelevant,
            variantsEmpty, twoMutations, unknownIssues)


def getgenes(req, maplist):  # gets genes and creates genelist entries for syndrome with one omim id
    try:
        geneomimid = maplist[0]["phenotypeMap"]["mimNumber"]
        omimidsweb = maplist[0]["phenotypeMap"]["phenotypeMimNumber"]
        geneList.append(createentry(features, geneomimid, genedict))
    except KeyError:
        geneomimid = maplist[0]["phenotypeMap"]["mimNumber"]
        omimidsweb = maplist[0]["phenotypeMap"]["phenotypeMimNumber"]
        geneList.append(createentry(features, geneomimid, genedict))
    if(omimidsweb == idsinjson):
        print("Ok")
    else:
        with open("Exceptions_by_incompatible_omims.txt", "a") as exceptions:
            exceptions.write(str(fileName) + "\t" + syndromename + "\n")


def get_syndrome_omim(omim_str_array):
    omim = '0'
    for omim_str in omim_str_array:
        if '(' in omim_str and ')' in omim_str:
            omim_array = omim_str.split(' ')
            if len(omim_array) > 1:
                for i in range(0, len(omim_array)):
                    if "(" in omim_array[i] and ")" in omim_array[i]:
                        mapping_key_str = omim_array[i].split(
                            "(")[1].split(")")[0]
                        if mapping_key_str.isdigit() and int(mapping_key_str) < 4 and int(mapping_key_str) > 0:
                            mapping_key = int(mapping_key_str)
                            if omim_array[i - 1].isdigit():
                                if int(omim_array[i - 1]) > 100:
                                    omim = (omim_array[i - 1])
    return omim


def get_phenotype_gene_dict():
    count = 0
    pg_dict = {}
    with open("morbidmap.txt") as csvfile:
        r = csv.reader(csvfile, delimiter='\t')
        for row in r:
            if count > 3 and len(row) > 1:
                omim_str_array = row[0].split(",")
                gene = row[2]
                omim_syndrome_id = get_syndrome_omim(omim_str_array)
                if omim_syndrome_id != '0':
                    pg_dict.update({omim_syndrome_id: gene})
            count = count + 1
    return pg_dict


def rename_document_vcf(documents, case_id, vcf_path, new_vcf_path):
    vcf = {}
    vcf_file_list = os.listdir(vcf_path)
    # go to vcf folder to check file and copy to new data/PEDIA/vcfs
    
    # 46073 is special case
    if case_id == 46073:
        return vcf
    for doc in documents:
        if isinstance(doc, list):
            continue
        if doc['is_vcf'] == 1:
            vcf_name = ""
            for name in vcf_file_list:
                if str(case_id) in name:
                    vcf_name = name
            print(vcf_name)
            new_file = str(case_id) + '.vcf.gz'
            if vcf_name.endswith('.gz'):
                if new_file not in os.listdir(new_vcf_path):
                    cmd = 'cp ' + vcf_path + vcf_name + ' ' + new_vcf_path + new_file
                    os.system(cmd)
            elif vcf_name.endswith('.vcf'):
                if new_file not in os.listdir(new_vcf_path):
                    cmd = 'bgzip -c ' + vcf_path + vcf_name + ' > ' + new_vcf_path + new_file
                    os.system(cmd)
            elif vcf_name.endswith('.zip'):
                if new_file not in os.listdir(new_vcf_path):
                    cmd = 'unzip -p ' + vcf_path + vcf_name + ' | bgzip > ' + new_vcf_path + new_file
                    os.system(cmd)
            print(new_file)
            if vcf_name != "":
                vcf['original_filename'] = new_file
    return vcf
# ===============================
# ===== main script =============
# ===============================


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Mapping disorder to gene')
    parser.add_argument('-j', '--jsonsoriginal', help='path of original json folder')
    parser.add_argument('-m', '--mappedjsons', help='path of mapped json folder')
    parser.add_argument('-f', '--vcf', help='path of vcf file we want to copy to')
    args = parser.parse_args()

    path = args.jsonsoriginal
    case_path = args.jsonsoriginal + '/cases/'
    vcf_path = args.jsonsoriginal + '/vcfs/'
    new_vcf_path = args.vcf
    genomic_path = args.jsonsoriginal + '/genomics_entries/'
    newpath = args.mappedjsons + '/'
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    #results = []
    genedict = makegenedict()
    defaultfeatures = ["combined_score",
                       "feature_score", "gestalt_score", "has_mask"]
    syn_gene_dict = get_phenotype_gene_dict()
    for fileName in os.listdir(case_path):
        if fileName == '.gitignore':
            continue

        test = fileName
        # debugging for jsons file not yet mapped
        if test not in os.listdir(newpath) and os.path.isfile(case_path + test):
            mainurl = "https://api.omim.org/api/entry/search?"
            file_content = json.load(open(case_path + fileName))
            result = []  # HGVS result
            #result.append(fileName)
            result = getGenomicDataCorrectFormat(file_content, genomic_path)
            #results.append(result)
            geneList = []
            file_content["genomicData"] = result
            file_content['vcf'] = rename_document_vcf(file_content['documents'], file_content['case_id'], vcf_path, new_vcf_path)
            for syndrome in file_content["detected_syndromes"]:
                syndromename = syndrome["syndrome_name"]
                try:  # debugging
                    idsinjson = syndrome["omim_id"]
                    #print(idsinjson)
                    if (type(idsinjson) == str):
                        idsinjson = int(idsinjson)

                    if (type(idsinjson) == list):  # if syndrome has multiple omim ids in the json file
                        genes = []
                        features = {"syndrome_name": syndromename}
                        for feature in defaultfeatures:  # get features of the syndrome
                            if feature in syndrome:
                                features[feature] = syndrome[feature]
                        for omim in idsinjson:
                            if str(omim) in syn_gene_dict:
                                genes.append(syn_gene_dict[str(omim)])
                        genes = set(genes)  # discard duplicates
                        for gene in genes:
                            # create an entry with the features and gene id
                            entry = createentry(features, gene, genedict)
                            geneList.append(entry.copy())

                    elif (type(idsinjson) == int):  # for syndrome with only one omim id
                        features = {"syndrome_name": syndromename}
                        for feature in defaultfeatures:
                            if feature in syndrome:
                                features[feature] = syndrome[feature]

                        if str(idsinjson) in syn_gene_dict:
                            gene = syn_gene_dict[str(idsinjson)]
                        # create an entry with the features and gene id
                        entry = createentry(features, gene, genedict)
                        geneList.append(entry.copy())

                    else:  # catches syndromes without omim ids
                        with open("Exceptions_no_omimid_in_file.txt", "a"):
                            write(str(fileName) + "\t" + syndromename + "\n")
                    # adds the genelist to the json file
                    file_content["geneList"] = geneList
                    # saves json file to new location
                    newjson = open(newpath + str(fileName), "w")
                    json.dump(file_content, newjson)
                except Exception:  # catches all other exceptions
                    with open("Other_Exceptions.txt", "a") as exceptions:
                        exceptions.write(
                            str(fileName) + "\t" + syndromename + "\n")
        else:
            print("already done")

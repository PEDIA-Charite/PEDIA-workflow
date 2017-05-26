import json, requests, os


# GetOpt to read cli inputs
import sys, getopt

# CLI-Options

argv = sys.argv[1:]

import re

morbidmap = ''
mimfile = ''
phenomizerconfig = ''
outputfolder = ''
inputfile = ''
outputfile = ''

try:
	opts, args = getopt.getopt(argv,"h::",["help","mimfile=","morbidmap=","phenomizerconfig=","inputfile=","outputfile="])
except getopt.GetoptError as e:
    print(e)
    print('jsonToTable.py --jsonsoriginal --log --login ')
    sys.exit(2)

for opt, arg in opts:
    if opt in ("-h", "--help"):
        print('jsonToTable.py -i <input-folder> -o <output-file.tsv>')
        sys.exit(1)
    elif opt in ("--morbidmap"):
        morbidmap = arg
    elif opt in ("--mimfile"):
        mimfile = arg
    elif opt in ("--phenomizerconfig"):
        phenomizerconfig = arg
    elif opt in ("--outputfile"):
    	outputfile = arg
    elif opt in ("--inputfile"):
    	inputfile = arg

print 'mimfile:',mimfile
print 'morbidmap file:',morbidmap
print 'phenomizerconfig file:',phenomizerconfig
print 'input file:',inputfile
print 'output file:',outputfile


config = json.loads(open(phenomizerconfig).read())


def load_omim():
	global mim
	global omimID2gene
	global geneToGeneID
	geneToGeneID={}
	mim={}
	omimID2gene={}
	for line in open(morbidmap):
		fields=line[:-1].split('\t')
		if len(fields)>1:
			ID=fields[0].split(', ')[-1][:-4]
			if '}' not in ID:
				genes=fields[1]
				genes=genes.split(', ')
				#print(ID, genes)
				if ID in omimID2gene:
				    for gene in genes:
				        omimID2gene[ID].append(gene)
				if ID not in omimID2gene:
				    omimID2gene[ID]=genes

	for line in open(mimfile):
		if line[0]!='#':
			fields=line[:-1].split('\t')
			entry = fields[1]
			if ((entry == "gene") & (fields[2] != "")):
				mimid=fields[0]
				gene=fields[3]
				gene_id=int(fields[2])
				mim[gene_id]={"mimid": mimid,"gene": gene}
				geneToGeneID[gene] = gene_id



def annotate_pheno_and_boqa():   #load omim first



	file = inputfile
	with open(file) as json_data:
	    data = json.load(json_data)
	data['processing']=['python phenomize_jsons.py']

	#### Get phenomizer
	gene_pattern = re.compile("(.+) \((\d+)\)")
	scores={}
	prefix=config["phenomizer"]
	suffix="&numres=100"
	hpo=''
	for i in data['features']:
	    #print(i)
	    hpo+=(i+',')
	link=prefix +hpo[:-1] + suffix
	print(link)
	content = requests.get(link).text
	lines=content.split('\n')
	for line in lines:
		fields=line.split('\t')
		if len(fields)>3:
			phenoscore= 1-float(fields[0])
			genes=fields[4]
			if len(genes)>1:
			    genes=genes.split(', ')
			    for geneterm in genes:
					gene_id = int(gene_pattern.findall(geneterm)[0][1])
					gene = str(gene_pattern.findall(geneterm)[0][0])
					if gene_id not in mim:
						mim[gene_id]={"mimid": "?","gene": gene}
					if gene not in geneToGeneID:
						geneToGeneID[gene]=gene_id
					if gene_id not in scores:
						scores[gene_id] = phenoscore
					else:
						scores[gene_id]=max(scores[gene_id], phenoscore)

	usedIDs=set()
	for entry in data['geneList']:
	#                if data['geneList'][j]['gene_symbol'] not in scores:
	#                    data['geneList'][j]['pheno_score']=float(0)
		gene_id = int(entry['gene_id'])
		if gene_id in scores:
			entry['pheno_score']=scores[gene_id]
			usedIDs.add(gene_id)
	for gene_id in scores:
		if gene_id in usedIDs:
			continue
		entry={}
		entry['gene_omim_id']=mim[gene_id]['mimid']
		entry['gene_symbol']=mim[gene_id]['gene']
		entry['pheno_score']=scores[gene_id]
		entry['gene_id']=gene_id
		data['geneList'].append(entry)

	##### Get BOQA
	scores={}
	usedIDs = set()
	prefix=config["boqa"]
	link=prefix +hpo[:-1] + suffix
	print(link)
	content = requests.get(link).text
	lines=content.split('\n')
	for line in lines:
		fields=line.split('\t')
		if len(fields)>2:
			boqascore= float(fields[0])
			Omim=fields[2][10:]
			if Omim in omimID2gene:
				genes=omimID2gene[Omim]
				for gene in genes:
					if gene not in geneToGeneID:
						print "cannot get gene_id for gene"+gene
						continue
					gene_id = geneToGeneID[gene]
					if gene_id not in scores:
						scores[gene_id] = boqascore
					else:
						scores[gene_id]=max(scores[gene_id], boqascore)
	usedIDs=set()
	for entry in data['geneList']:
	#                if data['geneList'][j]['gene_symbol'] not in scores:
	#                    data['geneList'][j]['pheno_score']=float(0)
		gene_id = int(entry['gene_id'])
		if gene_id in scores:
			entry['boqa_score']=scores[gene_id]
			usedIDs.add(gene_id)
	for gene_id in scores:
		if gene_id in usedIDs:
			continue
		entry={}
		entry['gene_omim_id']=mim[gene_id]['mimid']
		entry['gene_symbol']=mim[gene_id]['gene']
		entry['boqa_score']=scores[gene_id]
		entry['gene_id']=gene_id
		data['geneList'].append(entry)

	with open(outputfile, 'w') as f:
		json.dump(data, f)

load_omim()
annotate_pheno_and_boqa()

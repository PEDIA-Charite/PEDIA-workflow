
# JSON, dirs and http requests
import json, requests, os


# GetOpt to read cli inputs
import sys, getopt

# use regular expressions
import re









##### command-line input ####

# CLI-Options
morbidmap = ''
mimfile = ''
phenomizerconfig = ''
outputfolder = ''
inputfile = ''
outputfile = ''

argv = sys.argv[1:]

try:
    opts, args = getopt.getopt(argv,"h::",["help","mimfile=","morbidmap=","phenomizerconfig=","inputfile=","outputfile="])
except getopt.GetoptError as e:
    print(e)
    print("phenomize_jsons.py --mimfile <mim2gene.txt> --morbidmap <morbidmap.txt> --phenomizerconfig <protected/phenomizer_server.json> --inputfile <input.json> --outputfile <output.json>")
    sys.exit(2)

for opt, arg in opts:
    if opt in ("-h", "--help"):
        print("phenomize_jsons.py --mimfile <mim2gene.txt> --morbidmap <morbidmap.txt> --phenomizerconfig <protected/phenomizer_server.json> --inputfile <input.json> --outputfile <output.json>")
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

print('Mimfile:' + mimfile)
print('Morbidmap file:' + morbidmap)
print('Phenomizerconfig file:' + phenomizerconfig)
print('Input file:' + inputfile)
print('Output file:' + outputfile)

####
# function to load mim and mobidmap file
####



def load_omim():
	# gene_id -> {mimid: <id>, gene: <symbol>}
	global mim
	# omim_syndrom_id -> [<symbol1>, <symbol2>,...]
	global omimID2gene
	# symbol -> gene_id
	global geneToGeneID

	# init maps
	geneToGeneID={}
	mim={}
	omimID2gene={}

	#load morbidmap and set omimID2gene
	omim_sydrom_id_pattern = re.compile(".+, (\d+) \(\d+\)")
	for line in open(morbidmap):
		if line.startswith("#") or line.startswith("{"): # skip headers or sucessibility
			continue
		fields=line[:-1].split('\t')

		#if len(fields)>1:
		search_result = omim_sydrom_id_pattern.findall(fields[0])
		if len(search_result) > 0:
			omim_syndrom_id = int(search_result[0])
			genes=fields[1].split(', ')
			#print(ID, genes)
			if omim_syndrom_id in omimID2gene:
			    for gene in genes:
			        omimID2gene[omim_syndrom_id].append(gene)
			if omim_syndrom_id not in omimID2gene:
			    omimID2gene[omim_syndrom_id]=genes

	# load mimfile and set mim and geneToGeneID
	for line in open(mimfile):
		if line.startswith("#"): # skip headers
			continue
		fields=line[:-1].split('\t')
		entry = fields[1]
		if ((entry == "gene") & (fields[2] != "")):
			mimid=fields[0]
			gene=fields[3]
			gene_id=int(fields[2])
			mim[gene_id]={"mimid": mimid,"gene": gene}
			geneToGeneID[gene] = gene_id


def gethpoIDs(data):
	return data['features']

###################
# annotate phenomizer
###################

def annotate_phenomizer(data, url_prefix, url_suffix):   #load omim first


	data['processing']=['python phenomize_jsons.py']

	#### Get phenomizer
	gene_pattern = re.compile("(.+) \((\d+)\)")
	scores={}
	url_hpos =','.join(gethpoIDs(data))
	link=url_prefix + url_hpos + url_suffix
	print(link)

	# recieve test
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

	return data


###################
# annotate boqa
###################
def annotate_boqa(data, url_prefix, url_suffix):   #load omim first

	##### geterate link
	scores={}
	usedIDs = set()
	url_hpos = ','.join(gethpoIDs(data));
	link = url_prefix + url_hpos + url_suffix
	print(link)

	## get context from link
	content = requests.get(link).text
	lines=content.split('\n')
	for line in lines:
		fields=line.split('\t')
		if len(fields)>2: # there could be lines with only 2 entries
			boqascore= float(fields[0])

			syndromDB = fields[2].split(":")[0]
			if (syndromDB == "OMIM"): # skip orphanet
				omim_syndrom_id = int(fields[2].split(":")[2])
				if omim_syndrom_id in omimID2gene:
					genes=omimID2gene[omim_syndrom_id]
					for gene in genes:
						if gene not in geneToGeneID:
							print("cannot get gene_id for gene " + gene)
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

	return data

load_omim()

## load config fille to login on phenomizer/boqa
config = json.loads(open(phenomizerconfig).read())

# load data to input
with open(inputfile) as json_data:
	data = json.load(json_data)

data = annotate_phenomizer(data, config['phenomizer'], "&numres=100")
data = annotate_boqa(data, config['boqa'], "&numres=100")

#write data in pouput
with open(outputfile, 'w') as f:
	json.dump(data, f)

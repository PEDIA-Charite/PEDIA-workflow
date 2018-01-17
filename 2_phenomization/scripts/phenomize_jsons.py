import json
import requests
import os
import sys
import re
from argparse import ArgumentParser
import requests
import csv

import pandas

'''
jSON Phenomization tools
---

'''

def retrieve_file(path, url, cached):
    '''
    Get file from local filesystem. If it doesnt exists or cached is false, save from remote url to the local path
    and thereafter read from local filesystem.
    '''
    if not os.path.exists(path) or not cached:
        r = requests.get(url)
        with open(path, 'wb') as f:
            for data in r.iter_content():
                f.write(data)

    return pandas.read_csv(path, delimiter='\t', comment='#')

def get_omim_files(mimdir='', api_key='', use_cached=True):
    '''
    Download omim files from official website resources or get them locally if they already exist

    Args:
        mimdir: Directory where omim files are saved and retrieved, defaults to the current working directory
        api_key: OMIM API key retrieve it from the official omim website
        use_cached: Whether already downloaded files are used, if False omim files are always downloaded

    Returns:
        Data stucture with mim2gene and morbidmap
    '''
    mim2gene = retrieve_file('mim2gene.txt', 'https://omim.org/static/omim/data/mim2gene.txt', use_cached)

    morbidmap_url = 'https://data.omim.org/downloads/{}/morbidmap.txt'.format(api_key)
    morbidmap = retrieve_file('morbidmap.txt', morbidmap_url, use_cached)

    return {
            'mim2gene' : mim2gene
            ,'morbidmap' : morbidmap
            }

class PhenomizerService(requests.Session):
    def __init__(self, url, user, password):
        '''
        Create a new phenomizer service instance.

        Params:
            url: Url of phenomizer service
        '''
        super().__init__()
        self.url = url
        self.user = user
        self.password = password

    def request_phenomize(self, hpo_ids):
        params = {
                'mobilequery' : 'true'
                ,'username' : self.user
                ,'password' : self.password
                ,'terms' : hpo_ids
                ,'numres' : 100
                }
        r = self.get(self.url, params=params)
        print(r.url)
        r.raise_for_status()
        return r

    def request_boqa(self, hpo_ids):
        params = {
                'username' : self.user
                ,'password' : self.password
                ,'mobilequery' : 'true'
                ,'doboqa' : 'true'
                ,'terms' : hpo_ids
                ,'numres' : 100
                }
        r = self.get(self.url, params=params)
        print(r.url)
        return r

def phenomize_argv():
    parser = ArgumentParser()
    parser.add_argument('inputdir')
    parser.add_argument('outputdir')
    args = parser.parse_args()
    return args.inputdir, args.outputdir


if __name__ == '__main__':
    # omim_key = os.getenv('OMIM_KEY')
    # cached = os.getenv('OMIM_CACHED')
    # cached = cached == 'TRUE'

    # if omim_key is None and not cached:
    #     print(
#'''Please set OMIM_KEY in your environment variables to retrieve OMIM morbidmap.
#Alternatively provide genemap and morbidmap in the script directory and set OMIM_CACHED to TRUE'''
    #     )
    #     sys.exit(1)

    # inputdir, outputdir = phenomize_argv()

    # m = get_omim_files(api_key=omim_key)

    phenomizer_url = os.getenv('PHENOMIZER_URL')
    phenomizer_user = os.getenv('PHENOMIZER_USER')
    phenomizer_pw = os.getenv('PHENOMIZER_PW')

    pheno = PhenomizerService(phenomizer_url, phenomizer_user, phenomizer_pw)
    r = pheno.request_phenomize(['HP:0000118'])
    print(r.text)


## ##### command-line input ####
## 
## # CLI-Options
## morbidmap = ''
## mimfile = ''
## phenomizerconfig = ''
## outputfolder = ''
## inputfile = ''
## outputfile = ''
## 
## argv = sys.argv[1:]
## 
## try:
## 	opts, args = getopt.getopt(argv,"h::",["help","mimfile=","morbidmap=","phenomizerconfig=","inputfile=","outputfile="])
## except getopt.GetoptError as e:
##     print(e)
##     print("phenomize_jsons.py --mimfile <mim2gene.txt> --morbidmap <morbidmap.txt> --phenomizerconfig <protected/phenomizer_server.json> --inputfile <input.json> --outputfile <output.json>")
##     sys.exit(2)
## 
## for opt, arg in opts:
##     if opt in ("-h", "--help"):
##         print("phenomize_jsons.py --mimfile <mim2gene.txt> --morbidmap <morbidmap.txt> --phenomizerconfig <protected/phenomizer_server.json> --inputfile <input.json> --outputfile <output.json>")
##         sys.exit(1)
##     elif opt in ("--morbidmap"):
##         morbidmap = arg
##     elif opt in ("--mimfile"):
##         mimfile = arg
##     elif opt in ("--phenomizerconfig"):
##         phenomizerconfig = arg
##     elif opt in ("--outputfile"):
##     	outputfile = arg
##     elif opt in ("--inputfile"):
##     	inputfile = arg
## 
## print 'Mimfile:',mimfile
## print 'Morbidmap file:',morbidmap
## print 'Phenomizerconfig file:',phenomizerconfig
## print 'Input file:',inputfile
## print 'Output file:',outputfile
## 
## ####
## # function to load mim and mobidmap file
## ####
## 
## 
## 
## def load_omim():
## 	# gene_id -> {mimid: <id>, gene: <symbol>}
## 	global mim
## 	# omim_syndrom_id -> [<symbol1>, <symbol2>,...]
## 	global omimID2gene
## 	# symbol -> gene_id
## 	global geneToGeneID
## 
## 	# init maps
## 	geneToGeneID={}
## 	mim={}
## 	omimID2gene={}
## 
## 	#load morbidmap and set omimID2gene
## 	omim_sydrom_id_pattern = re.compile(".+, (\d+) \(\d+\)")
## 	for line in open(morbidmap):
## 		if line.startswith("#") or line.startswith("{"): # skip headers or sucessibility
## 			continue
## 		fields=line[:-1].split('\t')
## 
## 		#if len(fields)>1:
## 		search_result = omim_sydrom_id_pattern.findall(fields[0])
## 		if len(search_result) > 0:
## 			omim_syndrom_id = int(search_result[0])
## 			genes=fields[1].split(', ')
## 			#print(ID, genes)
## 			if omim_syndrom_id in omimID2gene:
## 			    for gene in genes:
## 			        omimID2gene[omim_syndrom_id].append(gene)
## 			if omim_syndrom_id not in omimID2gene:
## 			    omimID2gene[omim_syndrom_id]=genes
## 
## 	# load mimfile and set mim and geneToGeneID
## 	for line in open(mimfile):
## 		if line.startswith("#"): # skip headers
## 			continue
## 		fields=line[:-1].split('\t')
## 		entry = fields[1]
## 		if ((entry == "gene") & (fields[2] != "")):
## 			mimid=fields[0]
## 			gene=fields[3]
## 			gene_id=int(fields[2])
## 			mim[gene_id]={"mimid": mimid,"gene": gene}
## 			geneToGeneID[gene] = gene_id
## 
## 
## def gethpoIDs(data):
## 	return data['features']
## 
## ###################
## # annotate phenomizer
## ###################
## 
## def annotate_phenomizer(data, url_prefix, url_suffix):   #load omim first
## 
## 
## 	data['processing']=['python phenomize_jsons.py']
## 
## 	#### Get phenomizer
## 	gene_pattern = re.compile("(.+) \((\d+)\)")
## 	scores={}
## 	url_hpos =','.join(gethpoIDs(data))
## 	link=url_prefix + url_hpos + url_suffix
## 	print(link)
## 
## 	# recieve test
## 	content = requests.get(link).text
## 	lines=content.split('\n')
## 	for line in lines:
## 		fields=line.split('\t')
## 		if len(fields)>3:
## 			phenoscore= 1-float(fields[0])
## 			genes=fields[4]
## 			if len(genes)>1:
## 			    genes=genes.split(', ')
## 			    for geneterm in genes:
## 					gene_id = int(gene_pattern.findall(geneterm)[0][1])
## 					gene = str(gene_pattern.findall(geneterm)[0][0])
## 					if gene_id not in mim:
## 						mim[gene_id]={"mimid": "?","gene": gene}
## 					if gene not in geneToGeneID:
## 						geneToGeneID[gene]=gene_id
## 					if gene_id not in scores:
## 						scores[gene_id] = phenoscore
## 					else:
## 						scores[gene_id]=max(scores[gene_id], phenoscore)
## 
## 	usedIDs=set()
## 	for entry in data['geneList']:
## 	#                if data['geneList'][j]['gene_symbol'] not in scores:
## 	#                    data['geneList'][j]['pheno_score']=float(0)
## 		gene_id = int(entry['gene_id'])
## 		if gene_id in scores:
## 			entry['pheno_score']=scores[gene_id]
## 			usedIDs.add(gene_id)
## 	for gene_id in scores:
## 		if gene_id in usedIDs:
## 			continue
## 		entry={}
## 		entry['gene_omim_id']=mim[gene_id]['mimid']
## 		entry['gene_symbol']=mim[gene_id]['gene']
## 		entry['pheno_score']=scores[gene_id]
## 		entry['gene_id']=gene_id
## 		data['geneList'].append(entry)
## 
## 	return data
## 
## 
## ###################
## # annotate boqa
## ###################
## def annotate_boqa(data, url_prefix, url_suffix):   #load omim first
## 
## 	##### geterate link
## 	scores={}
## 	usedIDs = set()
## 	url_hpos = ','.join(gethpoIDs(data));
## 	link = url_prefix + url_hpos + url_suffix
## 	print(link)
## 
## 	## get context from link
## 	content = requests.get(link).text
## 	lines=content.split('\n')
## 	for line in lines:
## 		fields=line.split('\t')
## 		if len(fields)>2: # there could be lines with only 2 entries
## 			boqascore= float(fields[0])
## 
## 			syndromDB = fields[2].split(":")[0]
## 			if (syndromDB == "OMIM"): # skip orphanet
## 				omim_syndrom_id = int(fields[2].split(":")[2])
## 				if omim_syndrom_id in omimID2gene:
## 					genes=omimID2gene[omim_syndrom_id]
## 					for gene in genes:
## 						if gene not in geneToGeneID:
## 							print "cannot get gene_id for gene " + gene
## 							continue
## 						gene_id = geneToGeneID[gene]
## 						if gene_id not in scores:
## 							scores[gene_id] = boqascore
## 						else:
## 							scores[gene_id]=max(scores[gene_id], boqascore)
## 	usedIDs=set()
## 	for entry in data['geneList']:
## 	#                if data['geneList'][j]['gene_symbol'] not in scores:
## 	#                    data['geneList'][j]['pheno_score']=float(0)
## 		gene_id = int(entry['gene_id'])
## 		if gene_id in scores:
## 			entry['boqa_score']=scores[gene_id]
## 			usedIDs.add(gene_id)
## 	for gene_id in scores:
## 		if gene_id in usedIDs:
## 			continue
## 		entry={}
## 		entry['gene_omim_id']=mim[gene_id]['mimid']
## 		entry['gene_symbol']=mim[gene_id]['gene']
## 		entry['boqa_score']=scores[gene_id]
## 		entry['gene_id']=gene_id
## 		data['geneList'].append(entry)
## 
## 	return data
## 
## load_omim()
## 
## ## load config fille to login on phenomizer/boqa
## config = json.loads(open(phenomizerconfig).read())
## 
## # load data to input
## with open(inputfile) as json_data:
## 	data = json.load(json_data)
## 
## data = annotate_phenomizer(data, config['phenomizer'], "&numres=100")
## data = annotate_boqa(data, config['boqa'], "&numres=100")
## 
## #write data in pouput
## with open(outputfile, 'w') as f:
## 	json.dump(data, f)

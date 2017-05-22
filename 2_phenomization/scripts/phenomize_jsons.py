import json, requests, os


# GetOpt to read cli inputs
import sys, getopt

# CLI-Options

argv = sys.argv[1:]

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
            mimid=fields[0]
            gene=fields[3]
            mim[gene]=mimid



def annotate_pheno_and_boqa():   #load omim first
    file = inputfile
    with open(file) as json_data:
        data = json.load(json_data)
    data['processing']=['python phenomize_jsons.py']
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
                for i in genes:
                    if i not in scores:
                        scores[i]=phenoscore

    for entry in data['geneList']:
#                if data['geneList'][j]['gene_symbol'] not in scores:
#                    data['geneList'][j]['pheno_score']=float(0)
        if entry['gene_symbol'] in scores:
            entry['pheno_score']=scores[entry['gene_symbol']]
            scores.pop(entry['gene_symbol'])
    for i in scores:
        entry={}
#                entry['combined_score']=float(0)
#                entry['feature_score']=float(0)
        entry['gene_omim_id']='?'
        if i in mim:
            entry['gene_omim_id']=mim[i]
            #print(mim[i])
        entry['gene_symbol']=i
#                entry['gestalt_score']=float(0)
        entry['pheno_score']=scores[i]
        data['geneList'].append(entry)
    scores={}
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
                for i in genes:
                    if i not in scores:
                        scores[i]=boqascore
    for j in range(len(data['geneList'])):
#                if data['geneList'][j]['gene_symbol'] not in scores:
#                    data['geneList'][j]['boqa_score']=float(0)
        if data['geneList'][j]['gene_symbol'] in scores:
            data['geneList'][j]['boqa_score']=scores[data['geneList'][j]['gene_symbol']]
            scores.pop(data['geneList'][j]['gene_symbol'])
    for i in scores:
        entry={}
#                entry['combined_score']=float(0)
#                entry['feature_score']=float(0)
        entry['gene_omim_id']='?'
        if i in mim:
            entry['gene_omim_id']=mim[i]
        entry['gene_symbol']=i
#                entry['gestalt_score']=float(0)
#                entry['pheno_score']=float(0)
        entry['boqa_score']=scores[i]
        data['geneList'].append(entry)

    with open(outputfile, 'w') as f:
         json.dump(data, f)

load_omim()
annotate_pheno_and_boqa()

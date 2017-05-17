import json, requests, os

location='current_serverstatus/'

def load_omim():
    global mim
    global omimID2gene
    mim={}
    omimID2gene={}
    for line in open('morbidmap.txt'):
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
    
    for line in open('mim2gene.txt'):
        if line[0]!='#':
            fields=line[:-1].split('\t')
            mimid=fields[0]
            gene=fields[3]
            mim[gene]=mimid
            


def annotate_pheno_and_boqa():   #load omim first
    for file in os.listdir(location):
        if file[-5:]=='.json':
            with open(location+file) as json_data:
                data = json.load(json_data)
            data['processing']=['python phenomize_jsons2.py']
            scores={}
            prefix="http://compbio.charite.de/phenomizer/phenomizer/PhenomizerServiceURI?mobilequery=true&username=mensah&password=martin123&terms="
            suffix="&numres=100"
            hpo=''
            for i in data['features']:
                #print(i)
                hpo+=(i+',')
            link=prefix +hpo[:-1] + suffix
            #print(link)
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
                                
            for j in range(len(data['geneList'])):
#                if data['geneList'][j]['gene_symbol'] not in scores:
#                    data['geneList'][j]['pheno_score']=float(0)
                if data['geneList'][j]['gene_symbol'] in scores:
                    data['geneList'][j]['pheno_score']=scores[data['geneList'][j]['gene_symbol']]
                    scores.pop(data['geneList'][j]['gene_symbol'])
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
            prefix="http://compbio.charite.de/phenomizer/phenomizer/PhenomizerServiceURI?mobilequery=true&username=mensah&password=martin123&doboqa=true&terms="
            link=prefix +hpo[:-1] + suffix
#            print(link)
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
                            
            with open(location+file, 'w') as f:
                 json.dump(data, f)

load_omim()
annotate_pheno_and_boqa()       
            
for file in os.listdir(location):
    if file[-5:]=='.json':
        with open(location+file) as json_data:
            data = json.load(json_data) 
            break
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 12:00:22 2017

@author: Tori
"""
import datetime as dt # Abspeicherung des Files nach Datum ( um Historie verfolgen zu können, eventuell noch Änderung)

import json # JSON öffnen, bearbeiten und speichern 7 open, change and save JSONs
import os, shutil, re # Directory-Informationen bekommen / get information of directory

# HGVS-Code überprüfen und übersetzen in Position / proofread HGVS-code and translate to position
import pyhgvs
import pyhgvs.utils as hgvs_utils
from pygr.seqdb import SequenceFileDB

'''
import hgvs.dataproviders.uta
import hgvs.parser
import hgvs.assemblymapper
import hgvs.validator
import hgvs.exceptions

hdp = hgvs.dataproviders.uta.connect()
hp = hgvs.parser.Parser()
'''

# change HGVS with mutalyzer
import urllib2
import re

# VCF-Tabelle einrichten / create table for multi-VCF
import csv
import pandas as pd

# GetOpt to read cli inputs
import sys, getopt

# CLI-Options

argv = sys.argv[1:]

jsonsoriginal = ''
logfile = ''
genomefile = ''
debugfolder = ''
jsoncurratedfolder = ''#
errordict= ''
mVCF= ''
configfile= ''
try:
	opts, args = getopt.getopt(argv,"h::",["help","mappedjsons=","log=","genomefile=","debugfolder=","genefile=","jsoncurrated=","errordict=","config=", "vcf="])
except getopt.GetoptError as e:
    print(e)
    print('JsonQC.py --jsonsoriginal --log ')
    sys.exit(2)

for opt, arg in opts:
    if opt in ("-h", "--help"):
        print('jsonToTable.py -i <input-folder> -o <output-file.tsv>')
        sys.exit(1)
    elif opt in ("--mappedjsons"):
        jsonsoriginal = arg
    elif opt in ("--log"):
        logfile = arg
    elif opt in ("--genomefile"):
    	genomefile = arg
    elif opt in ("--debugfolder"):
    	debugfolder = arg
    elif opt in ("--genefile"):
    	genefile = arg
    elif opt in ("--jsoncurrated"):
    	jsoncurratedfolder = arg
    elif opt in ("--errordict"):
    	errordict = arg
    elif opt in ("--vcf"):
    	mVCF = arg
    elif opt in ("--config"):
    	configfile = arg

print 'Downloadfolder jsons:',jsonsoriginal
print 'Gene file:',genefile
print 'Genome file:',genomefile
print 'Debug folder:',debugfolder
print 'JSON currated folder:',jsoncurratedfolder
print 'Errordict:',errordict
print 'ConfigFile:',configfile
print 'Multi-VCF:',mVCF

if not os.path.exists(debugfolder):
    os.makedirs(debugfolder)
# # Read genome sequence using pygr.
genome = SequenceFileDB(genomefile)
#
# # Read RefSeq transcripts into a python dict.
with open(genefile) as infile:
    transcripts = hgvs_utils.read_transcripts(infile)

# # Provide a callback for fetching a transcript by its name.
def get_transcript(name):
    return transcripts.get(name)
#
# location='C:/Users/Tori/Documents/Python Scripts/QualityCheck/'

now=dt.datetime.now()
date=now.strftime('%Y-%m-%d')
time=now.strftime('%H:%M:%S')

## Resultfile aktualisieren
def checkjsons(step,folder):
    # HGVS-Code überprüfen / check HGVS-code
    def check_hgvs(hgvs, submitterteam, submitter, file):
		try:
			chrom, offset, ref, alt = pyhgvs.parse_hgvs_name(str(hgvs), genome, get_transcript=get_transcript)
		except ValueError, e: #'falsche' HGVS-Codes überspringen und anzeigen
			print e,hgvs
			if file not in overview[submitterteam]['incorrect JSONs'].keys():
				overview[submitterteam]['incorrect JSONs'][file]={}
			if 'falscher HGVS-Code' not in overview[submitterteam]['incorrect JSONs'][file]:
				overview[submitterteam]['incorrect JSONs'][file]['falscher HGVS-Code']={}
				overview[submitterteam]['incorrect JSONs'][file]['submitter']=submitter
				overview[submitterteam]['incorrect JSONs'][file]['falscher HGVS-Code'][hgvs]=str(e)
			else:
			    overview[submitterteam]['incorrect JSONs'][file]['falscher HGVS-Code'][hgvs]=str(e)
		except NotImplementedError, e:
			print e,hgvs
			if file not in overview[submitterteam]['incorrect JSONs'].keys():
			    overview[submitterteam]['incorrect JSONs'][file]={}
			if 'falscher HGVS-Code' not in overview[submitterteam]['incorrect JSONs'][file]:
			    overview[submitterteam]['incorrect JSONs'][file]['falscher HGVS-Code']={}
			    overview[submitterteam]['incorrect JSONs'][file]['submitter']=submitter
			    overview[submitterteam]['incorrect JSONs'][file]['falscher HGVS-Code'][hgvs]=str(e)
			else:
			    overview[submitterteam]['incorrect JSONs'][file]['falscher HGVS-Code'][hgvs]=str(e)


     # JSON-File zu inkorrekten Files hinzufügen / appen JSON to list of incorrect files
    def append_incorrect(file, submitterteam, submitter, str):
        if file not in overview[submitterteam]['incorrect JSONs'].keys():
            overview[submitterteam]['incorrect JSONs'][file]={}
            overview[submitterteam]['incorrect JSONs'][file]['submitter']=submitter
        overview[submitterteam]['incorrect JSONs'][file][str]=True

    overview={}
    withvcf=[]

    for file in os.listdir(folder):
        if file[-5:]=='.json':
            with open(folder+'/'+file) as json_data:
                print file
                d=json.load(json_data)

                # Infos aus JSON extrahieren / get relevant information from JSON
                try:
                    submitterteam=d['submitter']['user_team']
                except KeyError:
                    print (file)
                    with open('nosubmitter.txt', 'a') as nosub:
                        nosub.write(file)
                        continue
                if submitterteam==None:
                    submitterteam=d['submitter']['user_name']
                submitter=d['submitter']['user_name']
                mail=d['submitter']['user_email']
                vcf=d['vcf']
                # Grundaufbau des Dictionaries mit Zaehlung der Cases und der VCFs / define dictionary
                if submitterteam not in overview.keys():
                    overview[submitterteam]={}
                    overview[submitterteam]['team members']={}
                    overview[submitterteam]['number of cases']=1
                    overview[submitterteam]['VCFs']=0
                    overview[submitterteam]['correct JSONs']={}
                    overview[submitterteam]['correct JSONs']['number of correct jsons']=0#[]
                    overview[submitterteam]['correct JSONs']['list of correct jsons']=[]
                    overview[submitterteam]['incorrect JSONs']={}
                    overview[submitterteam]['not monogenic']=[]
                else:
                    overview[submitterteam]['number of cases']=overview[submitterteam]['number of cases']+1
                if submitter not in overview[submitterteam]['team members'].keys():
                    overview[submitterteam]['team members'][submitter]=mail
                if vcf!='noVCF':
                    overview[submitterteam]['VCFs']=overview[submitterteam]['VCFs']+1
                    withvcf.append(file)
                # hoechster Gestaltscore im JSON / get highest gestalt score
                gscore_list=[]
                if isinstance(d['geneList'],list):
                    for gene in d['geneList']:
                        gscore=gene['gestalt_score']
                        gscore_list.append(float(gscore))
                if len(gscore_list) != 0:
                    max_gscore=max(gscore_list)
                else:
                    max_gscore = 0
                    append_incorrect(file, submitterteam, submitter, 'No gene list!')
                 # Aufzaehlung der JSONs mit jeweiligen Maengeln / list of incorrect JSONs with annotated error
                 ## keine eingetragenen Features / no annotated features
                if len(d['features'])==0:
                    append_incorrect(file, submitterteam, submitter, 'keine Features')
                 ## kein Bild hochgeladen / no image uploaded
                if max_gscore==0:
                    append_incorrect(file, submitterteam, submitter, 'kein Bild')
                 ## keine molekulare Diagnose / no molecular diagnosis
                if d['selected_syndromes']=='Non selected':
                    append_incorrect(file, submitterteam, submitter, 'keine Diagnose angegeben')
                 ## Mutationen falsch eingetragen / something wrong with hgvs
                 ### mehrere Mutationen in mehreren Genen durchgehen / more than one mutatuion
                if len(d['genomicData'])==0:
                    append_incorrect(file, submitterteam, submitter, 'keine Mutation eingetragen')
                for mutation in d['genomicData']:
                    if len(mutation['Mutations'])==0:
                        append_incorrect(file, submitterteam, submitter, 'einmal keine Mutation eingetragen')
                        continue
                     # eine Mutation pro Gen / one mutation per gene
                    if 'HGVS-code' in mutation['Mutations'].keys():#nochmal checken
                        hgvs=mutation['Mutations']['HGVS-code']
                        check_hgvs(hgvs, submitterteam, submitter, file)
                     # compound heterozygous
                    elif 'Mutation 1' in mutation['Mutations'].keys():
                        hgvslist=[]
                        for mutationnr, mutationdict in mutation['Mutations'].items():
                            if 'HGVS-code' in mutationdict.keys():
                                hgvs=mutationdict['HGVS-code']
                                hgvslist.append(hgvs)
                        for hgvs in hgvslist:
                            check_hgvs(hgvs, submitterteam, submitter, file)
                     # kein HGVS-Code enthalten, weil nicht eingetragen / no hgvs code given
                    elif mutation['Test Information']['Mutation Type']=='Monogenic':
                        if 'protein level' in mutation['Mutations'].keys():
                            append_incorrect(file, submitterteam, submitter, 'Mutation auf Proteinebene')
                        else:
                            if vcf=='noVCF':
                                append_incorrect(file, submitterteam, submitter, 'kein HGVS-Code angegeben')
                            else:
                                append_incorrect(file, submitterteam, submitter, 'VCF vorhanden, aber kein HGVS-Code angegeben')
                    else:
                        append_incorrect(file, submitterteam, submitter, 'kein HGVS-Code angegeben (2)')
                    if mutation['Test Information']['Mutation Type']!='Monogenic':
                        overview[submitterteam]['not monogenic'].append(file)
                if file not in overview[submitterteam]['incorrect JSONs'].keys():
                    overview[submitterteam]['correct JSONs']['number of correct jsons']=overview[submitterteam]['correct JSONs']['number of correct jsons']+1
                    overview[submitterteam]['correct JSONs']['list of correct jsons'].append(file)

    with open(debugfolder+'/result_'+date+'.json', 'w') as dicttojson:
        json.dump(overview, dicttojson)

    return overview

## Dateien aus lokalem Directory löschen / remove data from local directory
def purge(dir, pattern):
    for f in os.listdir(dir):
        if re.search(pattern, f):
            os.remove(os.path.join(dir, f))
def copy(dir, pattern, outdir):
    for f in os.listdir(dir):
        if re.search(pattern, f):
            if not (os.path.exists(os.path.join(outdir,f ))):
				print "COPY " + f
				with open(os.path.join(dir, f),'r') as json_data:
					data = json.load(json_data)
				for entry in data["geneList"]:
					if "gene_id" in entry:
						entry["gene_id"] = int(entry["gene_id"])
					if "gene_omim_id" in entry:
						entry["gene_omim_id"] = int(entry["gene_omim_id"])
				#write data in pouput
				with open(os.path.join(outdir,f ), 'w') as f:
					json.dump(data, f)

## Quality Check der JSONs
wantto=raw_input('Copy missing original to json corrected folder? Change gene_omim_id and gene_id to "int"! y or n ')

if wantto=='y':
    copy(jsonsoriginal,'^.*.json',jsoncurratedfolder)

## Bereits bekannte Fehler aus dem Fehlerwörterbuch anwenden / change HGVS which are already in errordict

wantto='y'#raw_input('Make a result file? (y OR n) ')

if wantto=='y':
     step=''
     overview=checkjsons(step,jsoncurratedfolder)

jsons=0
vcfs=0
submitter_count=0
cor_jsons=0

for submitter in overview.keys():
    submitter_count=submitter_count+1
    jsons=jsons+int(overview[submitter]['number of cases'])
    vcfs=vcfs+int(overview[submitter]['VCFs'])
    cor_jsons=cor_jsons+int(overview[submitter]['correct JSONs']['number of correct jsons'])


wantto=raw_input('Correct JSONs from errordict? (y OR n) ')

dicterrors=0

if wantto=='y':
    with open(errordict) as json_data:
        transcript_errors=json.load(json_data)

    with open (debugfolder+'/result_'+step+date+'.json') as json_data:
        result=json.load(json_data)
        for submitter in result.keys():
            for jsonname in result[submitter]['incorrect JSONs'].keys():
                if unicode('falscher HGVS-Code') in result[submitter]['incorrect JSONs'][jsonname].keys():
                    with open(jsoncurratedfolder+'/'+jsonname, 'r') as jsonfile:
                        jsondata=json.load(jsonfile)
                        for hgvs, error in result[submitter]['incorrect JSONs'][jsonname]['falscher HGVS-Code'].items():
                            #hgvs=str(hgvs.replace(" ", "")) #Unicode-Fehlermeldung entgehen
                            if hgvs in transcript_errors.keys():
                                for mutation in jsondata['genomicData']:
                                    if 'HGVS-code' in mutation['Mutations'].keys():
                                        if mutation['Mutations']['HGVS-code']==hgvs:
                                             mutation['Mutations']['HGVS-code']=transcript_errors[hgvs]
                                             dicterrors=dicterrors+1
                                    elif 'Mutation 1' in mutation['Mutations'].keys():
                                        for multimut, description in mutation['Mutations'].items():
                                            if mutation['Mutations'][multimut]['HGVS-code']==hgvs:
                                                mutation['Mutations'][multimut]['HGVS-code']=transcript_errors[hgvs]
                                                dicterrors=dicterrors+1
                        with open(jsoncurratedfolder+'/'+jsonname, 'w') as dicttojson:
                            json.dump(jsondata, dicttojson)

    step='errordict'
    overview=checkjsons(step,jsoncurratedfolder)


## Fehlende oder falsche Transkripte mit Mutalyzer berichtigen
wantto=raw_input('Correct JSONs with mutalyzer? (y OR n) ')

muterrors=0

if wantto=='y':

    with open(errordict, 'r') as dicttojson:
        transcript_errors=json.load(dicttojson)

    new_hgvs=[]
    corrected_transcripts={}

    with open (debugfolder+'/result_'+date+'.json') as json_data:
        result=json.load(json_data)
        for submitter in result.keys():
            for jsonfile in result[submitter]['incorrect JSONs'].keys():
                if unicode('falscher HGVS-Code') in result[submitter]['incorrect JSONs'][jsonfile].keys():
                   for key in result[submitter]['incorrect JSONs'][jsonfile]['falscher HGVS-Code'].keys():
                       if 'required' in result[submitter]['incorrect JSONs'][jsonfile]['falscher HGVS-Code'][key]:
                           url="https://mutalyzer.nl/position-converter?assembly_name_or_alias=GRCh37&description="+key
                           try:
                               page =urllib2.urlopen(url)
                               data=page.read()
                           except:
                               print 'Could not connect: ', jsonfile, key
                               data='empty'
                           if 'Found transcripts' in data:
                               data=data.split('variant region</h4>')[1].split('<br></pre>')[0]
                               transcript=key.split('.')[0]
                               transcriptversion=key.split(':')[0].split('.')[1]
                               positions=[m.start() for m in re.finditer(transcript, data)]
                               alt_transcripts=[]
                               for position in positions:
                                   trans=data[position:].split(':')[0]
                                   trans_vers=trans.split(transcript)[1].split(':')[0].split('.')[1]
                                   if trans_vers>transcriptversion:
                                       transcriptversion=trans_vers
                               transcriptversion=transcript+'.'+transcriptversion+':'+key.split(':')[1]
                               with open (jsoncurratedfolder+'/'+jsonfile) as json_data:
                                   d=json.load(json_data)
                                   for mutation in d['genomicData']:
                                       if 'HGVS-code' in mutation['Mutations'].keys():
                                           if mutation['Mutations']['HGVS-code']==key:
                                               mutation['Mutations']['HGVS-code']=transcriptversion
                                               #append_errordict()
                                           elif 'Mutation 1' in mutation['Mutations'].keys():
                                               for multimut, description in mutation['Mutations'].items():
                                                   if mutation['Mutations'][multimut]['HGVS-code']==key:
                                                        mutation['Mutations']['HGVS-code']=transcriptversion
                                                        #append_errordict
                               transcript_errors[key]=transcriptversion
                               muterrors=muterrors+1
                               new_hgvs.append(transcriptversion)
                               print jsonfile, 'vorher: ', key, 'nachher: ', transcriptversion
                               corrected_transcripts[key]=transcriptversion
                           elif 'We found these versions' in data:
                               #print jsonfile, data
                               newtranscript=data.split('We found these versions: ')[1].split('<p></p>')[0].split('</p>')[0]
                               newtranscript=newtranscript+':'+key.split(':')[1]
                               with open (jsoncurratedfolder+'/'+jsonfile) as json_data:
                                   d=json.load(json_data)
                                   for mutation in d['genomicData']:
                                       if 'HGVS-code' in mutation['Mutations'].keys():
                                           if mutation['Mutations']['HGVS-code']==key:
                                               mutation['Mutations']['HGVS-code']=newtranscript
                                               #append_errordict()
                                           elif 'Mutation 1' in mutation['Mutations'].keys():
                                               for multimut, description in mutation['Mutations'].items():
                                                   if mutation['Mutations'][multimut]['HGVS-code']==key:
                                                        mutation['Mutations']['HGVS-code']=newtranscript
                                                        #append_errordict()
                               transcript_errors[key]=newtranscript
                               muterrors=muterrors+1
                               print jsonfile,  'vorher: ', key, 'nachher: ', newtranscript
                               new_hgvs.append(newtranscript)
                               corrected_transcripts[key]=newtranscript
                           elif 'could not be found in our database (or is not a transcript).' in data:
                               print 'no transcript found: ', submitter, jsonfile, key
                           else:
                               print 'please check: ', jsonfile, key

    with open(errordict, 'w') as dicttojson:
        json.dump(transcript_errors, dicttojson)

    step='mutalyzer'
    overview=checkjsons(step,jsoncurratedfolder)

# Mutationen per Hand berichtigen / correct HGVS-codes manually
wantto=raw_input('Change HGVS-codes manually? (y OR n) ')

manerror=0

if wantto!='n':
    not_solved={}

    with open (debugfolder+'/result_'+date+'.json') as json_data:
        result=json.load(json_data)
        for submitter in result.keys():
            for jsonname in result[submitter]['incorrect JSONs'].keys():
                if unicode('falscher HGVS-Code') in result[submitter]['incorrect JSONs'][jsonname].keys():
                    with open(jsoncurratedfolder+'/'+jsonname, 'r') as jsonfile:
                        jsondata=json.load(jsonfile)
                        for hgvs, error in result[submitter]['incorrect JSONs'][jsonname]['falscher HGVS-Code'].items():
                            print hgvs, error
                            for mutation in jsondata['genomicData']:
                                if 'HGVS-code' in mutation['Mutations'].keys():
                                    if mutation['Mutations']['HGVS-code']==hgvs:
                                        print '\n\n', jsonname, '\n', hgvs, error, '\n', mutation.items()
                                        right_hgvs=raw_input('right HGVS:')
                                        if right_hgvs=='n':
                                                not_solved[jsonname]=hgvs
                                                continue
                                        else:
                                            mutation['Mutations']['HGVS-code']=right_hgvs
                                            transcript_errors[hgvs]=right_hgvs
                                            manerror=manerror+1
                                        with open(errordict, 'w') as dicttojson:
                                            json.dump(transcript_errors, dicttojson)
                                elif 'Mutation 1' in mutation['Mutations'].keys():
                                    for multimut, description in mutation['Mutations'].items():
                                        if mutation['Mutations'][multimut]['HGVS-code']==hgvs:
                                            print '\n\n', jsonname, '\n', hgvs, error, mutation.items()
                                            right_hgvs=raw_input('right HGVS:')
                                            if right_hgvs=='n':
                                                not_solved[jsonname]=hgvs
                                                continue
                                            else:
                                                mutation['Mutations'][multimut]['HGVS-code']=right_hgvs
                                                transcript_errors[hgvs]=right_hgvs
                                                manerror=manerror+1
                                            with open(errordict, 'w') as dicttojson:
                                                json.dump(transcript_errors, dicttojson)
                            with open(jsoncurratedfolder+'/'+jsonname, 'w') as dicttojson:
                                json.dump(jsondata, dicttojson)

    step='manually'
    overview=checkjsons(step,jsoncurratedfolder)

# korrigierte JSONs in MultiVCF / dump corrected JSONs to multi-VCF
wantto=raw_input('Dump to multiVCF? (y OR n) ')

if wantto =='y':

    # Liste mit korrekten JSONs
    with open(debugfolder+'/result_'+date+'.json') as json_data:
        results=json.load(json_data)

    jsonlist=[]
    jsonlist2=[]
    novcf=[]
    withvcf=[]

    for file in os.listdir(jsoncurratedfolder+'/'):
        if file[-5:]=='.json':
            with open (jsoncurratedfolder+'/'+file) as json_data:
                jsondata=json.load(json_data)
                for submitter in results.keys():
                    if file in results[submitter]['correct JSONs']['list of correct jsons']: #and jsondata['vcf']=='noVCF':
                        jsonlist.append(file)
                        if jsondata['vcf']=='noVCF':
                            novcf.append(file)
                        else:
                            withvcf.append(file)


    # multiVCF erstellen
    multivcf=pd.DataFrame(columns=['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT', 'NM'])
    vcfcounter=0
    x=0

#    for file in os.listdir(location+'current_serverstatus/'):
#        if 'json' in file and file in jsonlist:
    for file in jsonlist:
        with open(jsoncurratedfolder+'/'+file) as json_data:
            d=json.load(json_data)
            #if d['vcf']!='noVCF':
            #    vcfcounter=vcfcounter+1
            #else:
            caseID=d['case_id']
            hgvslist=[]
            #if len(d['genomicData'])!=0:
            for mutation in d['genomicData']:
                if 'HGVS-code' in mutation['Mutations'].keys():
                    hgvs=mutation['Mutations']['HGVS-code']
                    hgvslist.append(hgvs)
                elif 'Mutation 1' in mutation['Mutations'].keys():
                    for mutationnr, mutationdict in mutation['Mutations'].items():
                        if 'HGVS-code' in mutationdict.keys():
                            hgvs=mutationdict['HGVS-code']
                            hgvslist.append(hgvs)
                else:
                    print 'keine Mutationen: ', file
                    continue
                genotype=mutation['Test Information']['Genotype']
                if genotype=='Hemizygous':
                    genotype='1'
                elif genotype=='Homozygous':
                    genotype='1/1'
                elif genotype=='Heterozygous' or genotype=='Compound Heterozygous':
                    genotype='0/1'
                else:
                    genotype='./1'
                for hgvscode in hgvslist:
                    try:
                        chrom, offset, ref, alt = pyhgvs.parse_hgvs_name(
                            str(hgvscode), genome, get_transcript=get_transcript)
                        if hgvscode in multivcf['NM'].tolist():
                            index=multivcf['NM'].tolist().index(hgvscode)
                            multivcf.set_value(index, caseID, genotype)
                            if caseID not in jsonlist2:
                                jsonlist2.append(caseID)
                        else:
                            chromo=chrom.split('chr')[1]
                            multivcf.set_value(x, '#CHROM', str(chromo))
                            try:
                                multivcf.set_value(x,'sort',int(chromo))
                            except ValueError,e:
                                multivcf.set_value(x,'sort',30)
                            multivcf.set_value(x, 'NM', hgvscode)
                            multivcf.set_value(x, 'POS', offset)
                            multivcf.set_value(x, 'ID', '.')
                            multivcf.set_value(x, 'REF', str(ref))
                            multivcf.set_value(x, 'ALT', str(alt))
                            multivcf.set_value(x, 'QUAL', '.')
                            multivcf.set_value(x, 'FILTER', '.')
                            multivcf.set_value(x, 'INFO', 'HGVS="'+hgvscode+'"')
                            multivcf.set_value(x, 'FORMAT', 'GT')
                            multivcf.set_value(x, caseID, genotype)
                            x=x+1
                            if caseID not in jsonlist2:
                                jsonlist2.append(caseID)
                    except ValueError, e: #'falsche' HGVS-Codes überspringen und anzeigen
                        print 'Error:',file, hgvs, e
                        continue



    ##data_vcf sortieren
    print 'Sort DataFrame ...'
    multivcf=multivcf.sort_values(by=['sort', 'POS'])
    multivcf=multivcf.drop('sort',axis=1)
    multivcf=multivcf.drop('NM',axis=1)
    multivcf=multivcf.reset_index(drop=True)

    #leere Felder füllen
    multivcf=multivcf.fillna(value='0/0')

    multivcf.to_csv(mVCF+".tmp", sep='\t', index=False, header=True, quoting=csv.QUOTE_NONE)

    with open(mVCF, 'w') as outfile:
        outfile.write('##fileformat=VCFv4.1\n##INFO=<ID=HGVS,Number=1,Type=String,Description="HGVS-Code">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        with open(mVCF+".tmp",'r') as infile:
            for line in infile:
                outfile.write(line)

    os.remove(mVCF+".tmp")

    #Liste der Cases im MultiVCF
    with open (configfile, 'w') as jsoninmultivcf:
        jsoninmultivcf.write('SINGLE_SAMPLES:\n')
        for item in novcf:
            item=item.split('.')[0]
            jsoninmultivcf.write(" - %s\n" % item)
        jsoninmultivcf.write('VCF_SAMPLES:\n')
        for item in withvcf:
            item=item.split('.')[0]
            jsoninmultivcf.write(" - %s\n" % item)

## Ordnen der falschen JSONs nach Submitter / arrange incorrect JSONs per submitter
wantto=raw_input('Rearrange result file? (y OR n) ')

if wantto!='n':
    submitter_cases={}
    with open (debugfolder+'/result_'+date+'.json') as json_data:
        result=json.load(json_data)
        for submitterteam in result.keys():
            submitter_cases[submitterteam]={}
            for submitter in result[submitterteam]['team members'].keys():
                submitter_cases[submitterteam][submitter]={}
            for jsonfile in result[submitterteam]['incorrect JSONs'].keys():
                submitter=result[submitterteam]['incorrect JSONs'][jsonfile]['submitter']
                indeldups=['ins', 'del', 'dup']
                submitter_cases[submitterteam][submitter][jsonfile]={}
                submitter_cases[submitterteam][submitter][jsonfile].update(result[submitterteam]['incorrect JSONs'][jsonfile])
                del submitter_cases[submitterteam][submitter][jsonfile]['submitter']
                for error in result[submitterteam]['incorrect JSONs'][jsonfile].keys():
                    if 'falscher HGVS-Code' in error:
                        for key in result[submitterteam]['incorrect JSONs'][jsonfile]['falscher HGVS-Code'].keys():
                            if any(indeldup in key for indeldup in indeldups):
                                del submitter_cases[submitterteam][submitter][jsonfile]['falscher HGVS-Code']
                if len(submitter_cases[submitterteam][submitter][jsonfile])==0:
                    del submitter_cases[submitterteam][submitter][jsonfile]

    with open(debugfolder+'/resultpersubmitter'+date+'.json', 'w') as dicttojson:
        json.dump(submitter_cases, dicttojson)

## Dokumentation in Tabellenform

qc_jsons=0
qc_vcfs=0

with open(debugfolder+'/result_'+date+'.json') as jsondata:
    result=json.load(jsondata)
    for submitter in result.keys():
        qc_jsons=qc_jsons+int(result[submitter]['correct JSONs']['number of correct jsons'])
        for jsonfile in result[submitter]['correct JSONs']['list of correct jsons']:
            with open(jsoncurratedfolder+'/'+jsonfile) as json_data:
                vcfin=json.load(json_data)
                if vcfin['vcf']!='noVCF':
                    qc_vcfs=qc_vcfs+1


progress = pd.read_excel(open(logfile,'rb'), sheetname=0)
newest = [date, time, submitter_count, jsons, vcfs, cor_jsons, dicterrors, muterrors, manerror, qc_jsons, qc_vcfs]
progress.loc[len(progress)] = newest
progress.to_excel(logfile, index=False)

'''
#Phenomizer
wantto=raw_input('Phenomize? (y OR n) ')

if wantto!='n':
    os.system(location+'Phenomize/phenomize_jsons.py')
'''

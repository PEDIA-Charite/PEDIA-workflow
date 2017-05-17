# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 14:26:55 2017

@author: Tori
"""

# get data for Organigramm

import json # JSON öffnen, bearbeiten und speichern 7 open, change and save JSONs

location='C:/Users/Tori/Documents/Python Scripts/QualityCheck/current_serverstatus/'
resultfile='C:/Users/Tori/Documents/Python Scripts/QualityCheck/Results/result_manually100217.json'

jsons=0
vcfs=0
qc_jsons=0
qc_vcfs=0

with open(resultfile) as jsondata:
    result=json.load(jsondata)
    for submitter in result.keys():
        jsons=jsons+int(result[submitter]['number of cases'])
        vcfs=vcfs+int(result[submitter]['VCFs'])
        qc_jsons=qc_jsons+int(result[submitter]['correct JSONs']['number of correct jsons'])
        for jsonfile in result[submitter]['correct JSONs']['list of correct jsons']:
            with open(location+jsonfile) as json_data:
                vcfin=json.load(json_data)
                if vcfin['vcf']!='noVCF':
                    qc_vcfs=qc_vcfs+1
    #jsons
    
print 'JSON-Dateien: ', jsons, 'davon mit VCF: ', vcfs
print 'JSON-Dateien, die QC bestehen: ', qc_jsons, 'davon mit VCF: ', qc_vcfs
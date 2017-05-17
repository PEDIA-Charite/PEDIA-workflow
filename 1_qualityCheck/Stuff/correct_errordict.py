# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 14:59:42 2017

@author: Tori
"""

import json

empty=[]

with open ('C:/Users/Tori/Documents/Python Scripts/QualityCheck/hgvs_errordict.json') as json_data:
    d=json.load(json_data)
    for key, value in d.iteritems():
        if key==value:
            print key, value
            empty.append(key)
            
for item in empty:
    d.pop(item)
            
with open('C:/Users/Tori/Documents/Python Scripts/QualityCheck/hgvs_errordict.json', 'w') as dicttojson:
    json.dump(d, dicttojson)
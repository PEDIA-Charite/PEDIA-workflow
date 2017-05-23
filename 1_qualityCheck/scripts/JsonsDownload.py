# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 12:00:22 2017

@author: Tori
"""
import datetime as dt # Abspeicherung des Files nach Datum ( um Historie verfolgen zu können, eventuell noch Änderung)

import ftplib as ftp # Dateien vom Server holen / import JSONs from server
import json # JSON öffnen, bearbeiten und speichern 7 open, change and save JSONs
import os, shutil, re # Directory-Informationen bekommen / get information of directory

# GetOpt to read cli inputs
import sys, getopt

# CLI-Options

argv = sys.argv[1:]

output = ''
loginfile = ''
sample=''
try:
	opts, args = getopt.getopt(argv,"h::",["help","output=","login=","sample="])
except getopt.GetoptError as e:
    print(e)
    print('JsonsDownload.py --sample <sample> --output <output_file> --login <login_config_for_server>')
    sys.exit(2)

for opt, arg in opts:
    if opt in ("-h", "--help"):
        print('JsonsDownload.py --sample <sample> --output <file> --login <login_config_for_server>')
        sys.exit(1)
    elif opt in ("--output"):
        output = arg
    elif opt in ("--sample"):
        sample = arg
    elif opt in ("--login"):
        loginfile = arg

print 'Output file:',output
print 'Login file:',loginfile
print 'Sample:',sample

now=dt.datetime.now()
date=now.strftime('%Y-%m-%d')
time=now.strftime('%H:%M:%S')

## Dateien vom Server runterladen / get data from server
config = json.loads(open(loginfile).read())

json_server=ftp.FTP(config['url'])
json_server.login(config['login'],config['password'])

directory='/'
json_server.cwd(directory)

# download
filename = sample+ ".json"
file = open(output, 'wb')
print 'Downloading Sample ', sample, ' ...'
json_server.retrbinary('RETR '+filename, file.write)
file.close()
json_server.quit()

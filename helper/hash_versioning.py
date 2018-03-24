#!/usr/bin/env python3
import sys
import hashlib

CHUNK = 128

genhash = hashlib.md5()
inpath = sys.argv[1]
with open(inpath, "rb") as infile:
    while True:
        bytestr = infile.read(CHUNK)
        if not bytestr:
            break
        genhash.update(bytestr)
print(genhash.hexdigest())

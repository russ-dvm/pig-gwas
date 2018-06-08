#!/usr/bin/python
from __future__ import print_function
import requests, sys

server = "http://rest.ensembl.org"
ext = "/variation/sus_scrofa/"


rsIDs = open(sys.argv[1],'r')

for rsID in rsIDs:
    rsID = rsID.rstrip()
    extPlus = server+ext+rsID+"?"
    r = requests.get(extPlus, headers={ "Content-Type" : "application/json"})

    decoded = r.json()
    try:
        coords = decoded["mappings"][0]["location"]
        chrom = coords.split(":")[0]
        pos = coords.split(":")[1].split("-")[0]
        print(chrom, rsID, 0, pos, sep = "\t")
    except:
        print(chrom, rsID, 0, "FAIL", sep = "\t")


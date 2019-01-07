#!/bin/python

from __future__ import print_function
import sys

# Compare snps with rsIDs from Russ and from Jer

russ = open(sys.argv[1], 'r')
jer = open(sys.argv[2], 'r')


jer_list = []
russ_list = []

for snp in jer:
	snp = snp.rstrip()
	jer_list.append(snp)



for line in russ:
	line = line.rstrip()
	russ_snp = line.split("\t")[1]
	russ_list.append(russ_snp)


for item in russ_list:
	if item in jer_list:
		continue
	else:
		print(item, "in russ not in jer")
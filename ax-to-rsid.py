#!/bin/python

from __future__ import print_function
import sys

## Converts the AX ids to RSids. Requires two files: one) axids to convert and 2) axids with corresponding rsids. 
# Generate file1 - cut -f 1 xxx.map
# Generate file2 - cut -f 1,812 ref_files/all_snp_info.txt


ax = open(sys.argv[1])
ref = open(sys.argv[2])

def listify(x):
    x_list = []
    for line in x:
        line = line.rstrip()
        x_list.append(line)


    return x_list


def dictify(x):
    x_dict = {}
    for line in x:
        line = line.rstrip()
        lineFields = line.split("\t")
        try:
            x_dict[lineFields[0]] = lineFields[1]
        except:
            x_dict[lineFields[0]] = "no_match"
    return x_dict

ax_list = listify(ax)
ref_dict = dictify(ref)

for x in ax_list:
	if x in ref_dict:
		print(x, ref_dict[x], sep = "\t")
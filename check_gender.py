#!/usr/bin/python

from __future__ import print_function
import argparse


parser = argparse.ArgumentParser(description = "Checks that the gender computed by the SNP chip is the same as the known/submitted gender. True will print for those that are the same, false for those that are different, and NA if there is no corresponding ID in the computed gender file", epilog = "RS Fraser, 2018-03-16")

parser.add_argument("known", help = "A tab delimited file with pig ID and gender in separate columns.")
parser.add_argument("computed", help = "A tab delimited file with pig ID and the computed gender.")

args = parser.parse_args()

def dictify(x):
    y = {}
    try:
        for line in x:
            line = line.rstrip()
            lineField = line.split("\t")
            y[lineField[0]] = lineField[1]
    except:
        print("error. this is the data in the line that caused a problemo:\n", line)
        exit()
    return y


known = open(args.known,'r')
computed = open(args.computed, 'r')

known_dict = dictify(known)
comp_dict = dictify(computed)

for k,v in known_dict.iteritems():
    try:
        if comp_dict[k] == v:
            print(k, "TRUE", sep = "\t")
        elif comp_dict[k] == "unknown":
            print(k, "unknown_comp", sep = "\t")
        else:
            print(k, "FALSE", sep = "\t")
    except:
        print(k, "NA", sep = "\t")
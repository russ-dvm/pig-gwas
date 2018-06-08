#!/usr/bin/python

from __future__ import print_function
import sys

## written and run faster than it took grep to do the same damn thing
## basically this just spits out AX_ids that are missing from the final output but that are present on the snpchip.

full = open(sys.argv[1], 'r')
rest = open(sys.argv[2], 'r')

def listify(x):
    x_list = []
    for line in x:
        line = line.rstrip()
        x_list.append(line)


    return x_list


full_list = listify(full)
rest_list = listify(rest)

for item in full_list:
    if item in rest_list:
        continue
    else:
        print(item)



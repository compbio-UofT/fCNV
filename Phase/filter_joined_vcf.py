#!/usr/bin/python2

import argparse
from sys import exit
import random
import copy
import samParse as sp
from datetime import datetime

def is_within_intervals(num, intervals):
    for interval in intervals:
        if interval[0] <= num <= interval[1]:
            return True
    return False

def main():
    #parse ARGs
    parser = argparse.ArgumentParser(description='Filter joined VCF (M, P, F) and intersect with M Kitzman phasing. Read filenames: joined .vcf file and maternal phasing.')
    parser.add_argument('filenames', type=str, nargs='+', help='paths to 1) joined .vcf file; 2) Kitzman maternal phasing.')
    args = parser.parse_args()
    
    if len(args.filenames) != 2: exit("Unexpected number of arguments passed! Expecting 2 filenames.")
    
    #treat these as CONSTANTS!
    J = 0; PH = 1;
    ALL = [J, PH]
    
    #list of input files
    in_files = [open(args.filenames[i], "r" ) for i in ALL]
    
    #joint output file
    out_file = open(args.filenames[J]+".ftr.vcf", "w")

    kitz = {}
    for line in in_files[PH].readlines():
        line = line.rstrip('\n').split('\t')
        pos = line[2].split(':')[1]
        kitz[pos] = (line[5], line[6], line[7])
    
    #for x, y in kitz.items():
    #    print x, y
        
    overlap = missing = 0
    for line in in_files[J].readlines():
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 3:
            print >>out_file, line.rstrip('\n')
            continue
        pos = fields[1]
        m = fields[9].split(':')[0]
        if pos in kitz:
            #print kitz[pos], fields[3], fields[4], m
            overlap += 1
        else:
            if m == '0/0': print '0/0'
            if m == '.':
                overlap += 1 
            else:
                missing +=1
    
    print missing, overlap
    
if __name__ == '__main__':
    #import doctest
    #doctest.testmod()
    
    main()


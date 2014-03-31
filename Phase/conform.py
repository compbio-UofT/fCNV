#!/usr/bin/python2

import argparse
from sys import exit
import random
import copy
from datetime import datetime

def main():
    #parse ARGs
    parser = argparse.ArgumentParser(description='Filter joined VCF (M, P, F) and intersect with reference VCF. Read filenames: joined .vcf file and maternal phasing.')
    parser.add_argument('filenames', type=str, nargs='+', help='paths to 1) joined .vcf file; 2) reference phased VCF.')
    args = parser.parse_args()
    
    if len(args.filenames) != 2: exit("Unexpected number of arguments passed! Expecting 2 filenames.")
    
    #treat these as CONSTANTS!
    J = 0; R = 1;
    ALL = [J, R]
    
    #list of input files
    in_files = [open(args.filenames[i], "r" ) for i in ALL]
    
    #joint output file
    out_file = open(args.filenames[J][:-3]+"ftr.vcf", "w")

    ref = {}
    for line in in_files[R].readlines():
        line = line.rstrip('\n').split('\t')
        if line[0][0] == '#': continue
        pos = line[1]
        ref[pos] = (line[3], line[4])
    
    #for x, y in ref.items():
    #    print x, y
        
    okpos = wrpos = 0
    excluded = 0
    for line in in_files[J].readlines():
        fields = line.rstrip('\n').split('\t')
        if fields[0][0] == '#':
            print >>out_file, line.rstrip('\n')
            continue
        pos = fields[1]
        alleles = (fields[3], fields[4])
        if pos in ref:
            print >>out_file, line.rstrip('\n')
            if ref[pos] == alleles: #TODO: ?? exclude tri-allelic SNPs ?? 
                okpos += 1
            else:
                wrpos += 1 
        else:
            excluded += 1
    
    print 'alleles ok:', okpos, ' different:', wrpos
    print 'dropped positions:', excluded
    
if __name__ == '__main__':
    #import doctest
    #doctest.testmod()
    
    main()


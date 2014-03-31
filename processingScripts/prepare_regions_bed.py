#!/usr/bin/python2

import argparse
from sys import exit
import random
import copy
from datetime import datetime
import subprocess
import os

def main():
    #parse ARGs
    parser = argparse.ArgumentParser(description='Prepare BED files of valid regions for each chromosome.')
    parser.add_argument('filenames', type=str, nargs='+', help='paths to 1) dir with chr$NUM dirs; 2) centromeres list file.')
    args = parser.parse_args()
    
    if len(args.filenames) != 2: exit("Unexpected number of arguments passed! Expecting 1 parameters.")
    
    #treat these as CONSTANTS!
    BED = 0; CT = 1; #in_files
    ALL = [BED, CT]
    
    #list of input files
    in_files = [None for i in ALL]
    in_files[CT] = open(args.filenames[CT], "r" )
    res_path = args.filenames[BED]
    
    
    #read centromeres positions
    centromeres = dict()
    for line in in_files[CT].readlines():
        line = line.rstrip('\n').split('\t')
        if line[0] not in centromeres.keys(): centromeres[line[0]] = []
        centromeres[line[0]] += [(int(line[1]), int(line[2]))]
    
    for chrom, cent in centromeres.iteritems():
        if chrom == 'chrY': continue
        if len(centromeres[chrom]) not in [2,  3]: print "ERROR"
        
        out_file = open(res_path + "/" + chrom + "/regions.bed", "w")
        print res_path + "/" + chrom + "/regions.bed"
        
        print >>out_file, "{0}\t{1}\t{2}".format(chrom, cent[0][1], cent[1][0])
        if len(centromeres[chrom]) == 3: print >>out_file, "{0}\t{1}\t{2}".format(chrom, cent[1][1], cent[2][0])
        
        out_file.close()

    
if __name__ == '__main__':
    main()


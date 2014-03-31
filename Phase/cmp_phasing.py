#!/usr/bin/python2

import argparse
from sys import exit
import random
import copy
from datetime import datetime

def is_within_intervals(num, intervals):
    for interval in intervals:
        if interval[0] <= num <= interval[1]:
            return True
    return False

def main():
    #parse ARGs
    parser = argparse.ArgumentParser(description='Compare phased joined VCF (M, P, F) on intersect with M Kitzman phasing. Read filenames: joined .vcf file and maternal phasing.')
    parser.add_argument('filenames', type=str, nargs='+', help='paths to 1) phased joined .vcf file; 2) Kitzman maternal phasing.')
    args = parser.parse_args()
    
    if len(args.filenames) != 2: exit("Unexpected number of arguments passed! Expecting 2 filenames.")
    
    #treat these as CONSTANTS!
    J = 0; PH = 1;
    ALL = [J, PH]
    
    #list of input files
    in_files = [open(args.filenames[i], "r" ) for i in ALL]

    kitz = {}
    for line in in_files[PH].readlines():
        line = line.rstrip('\n').split('\t')
        pos = line[2].split(':')[1]
        kitz[pos] = (line[5], line[6], line[7])
    
    #for x, y in kitz.items():
    #    print x, y
        
    overlap = missing = 0
    nagr = agr = magr = 0
    aagr = []
    for line in in_files[J].readlines():
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 3:
            continue
        pos = fields[1]
        m = fields[9].split(':')[0]
        if pos in kitz:
            overlap += 1
            print kitz[pos], fields[3], fields[4], m
            hs = str(1-(kitz[pos][0]==fields[3]))+'|'+str(int(kitz[pos][1]==fields[4]))
            print "      ", hs, m, '                ', fields[9].split(':')[2]
            if hs==m: 
                agr+=1
                if nagr==0: continue
                magr = max(nagr, magr)
                aagr.append(nagr)
                nagr = 0
            else:
                nagr+=1
                if agr==0: continue
                magr = max(agr, magr)
                aagr.append(agr)
                agr = 0
            
            #print fields[9:12]
        else:
            #if m == '0|0': print '!!!!!!!! 0|0'
            if m == '.': print '!!!!!!!! .'
            missing +=1
                
    magr = max(agr, magr)
    aagr.append(agr)
    aagr.append(nagr)
    print "new_pos:", missing, " overlap:", overlap, "/", len(kitz), " max_agree:", magr
    print aagr
    
    count_in_small = num_of_small = 0
    for x in aagr:
        if x < 10: 
            count_in_small += x
            num_of_small += 1
    print num_of_small, count_in_small, sum(aagr)
    
    
if __name__ == '__main__':
    #import doctest
    #doctest.testmod()
    
    main()


#!/usr/bin/python2

import argparse
import sys
import random
import copy
from datetime import datetime
import matplotlib.pyplot as plt

def isWithinIntervals(num, intervals):
    for interval in intervals:
        if interval[0] <= num <= interval[1]:
            return True
    return False

def main():
    #parse ARGs
    parser = argparse.ArgumentParser(description='Compare phasings done by Kitzman@UW and Beagle.')
    parser.add_argument('filenames', type=str, nargs='+', help='paths to Kitzman@UW file and to phased .vcf file.')
    parser.add_argument('--centromeres', type=str, nargs=1, help='path centromeres list file.')
    args = parser.parse_args()
    
    if len(args.filenames) != 2: sys.exit("Unexpected number of arguments passed! Expecting 2 filenames.")
    
    #treat these as CONSTANTS!
    KT = 0; BG = 1; CT = 2;
    ALL = [KT, BG]
    
    #check if centromeres list is given
    filter_centromeres = 0
    if args.centromeres != None:
        filter_centromeres = 1
        ALL += [CT]
        args.filenames += args.centromeres
    
    #list of input files
    in_files = [open(args.filenames[i], "r" ) for i in ALL]
    
    #read centromeres positions
    if filter_centromeres:
        centromeres = dict()
        for line in in_files[CT].readlines():
            line = line.rstrip('\n').split('\t')
            if line[0] not in centromeres.keys(): centromeres[line[0]] = []
            centromeres[line[0]] += [(int(line[1]), int(line[2]))]
    
    #alleles and haplotypes for particular SNP positions
    lociKT = dict()
    lociBG = dict()
    processed_chr = ''
    
    #parse Beagle VCF file
    for line in in_files[BG].readlines():
        line = line.rstrip('\n')
        if line == '': continue
        #skip the header
        if len(line) > 0 and line[0] == '#': continue
        #split
        line = line.split('\t')
        
        snp = [0] * 4
        #get chromosome
        if processed_chr == '': processed_chr = line[0]
        if processed_chr != line[0]: print "WARNING: multiple chromosomes in the VCF", processed_chr, "|", line[0]
        #position
        snp[0] = int(line[1])
        #dbSNP rsid
        snp[1] = line[2]
        #ref and alt alleles
        snp[2] = (line[3], line[4])

        #parse out haplotype config info
        ht = map(int, line[9].split('/'))
        #hapA and hapB
        snp[3] = (snp[2][ht[0]], snp[2][ht[1]])
        
        lociBG[snp[0]] = tuple(snp)
    
    #parse Kitzman's haplotype file chr7:50354-231527	chr7_blk0000000	chr7:50354	C	A	A	C	rs6583339
    for line in in_files[KT].readlines():
        line = line.rstrip('\n')
        if line == '': continue
        #skip the header
        if len(line) > 0 and line[0] == 'b': continue
        #split
        line = line.split('\t')
        
        snp = [0] * 4
        #get chromosome
        coord = line[2].split(':')
        if coord[0] != processed_chr: continue
        #position
        snp[0] = int(coord[1])
        #dbSNP rsid
        snp[1] = line[7]
        #ref and alt alleles
        snp[2] = (line[3], line[4])

        #hapA and hapB
        snp[3] = (line[5], line[6])
        
        lociKT[snp[0]] = tuple(snp)
    
    #print lociBG.items()[4], '\n', lociKT.items()[4]
    
    segLengthsSet = []
    if filter_centromeres: centromere_regions = centromeres[processed_chr]
    #--------
    intersectSize1 = 0
    intersectMatchGT1 = 0
    intersectMatchHT1 = 0
    currentSegLen = 0; segCount = 0; segMax = -1; segMin = 1e10;
    revSegLen = 0; revSegCount = 0; revSegMax = -1; revSegMin = 1e10;
    prev_pos = 0
    skipped_in_centromere = 0
    missingStreak = 0
    #allPos = [(x, 0) for x in lociBG.keys()] + [(x, 1) for x in lociKT.keys()]
    for pos in sorted(lociBG.keys()):
        snp = lociBG[pos]
        #if snp[0] in lociKT.keys():
        try:
            snp2 = lociKT[snp[0]]
            
            if filter_centromeres and isWithinIntervals(pos, centromere_regions):
                skipped_in_centromere += 1
                continue
            
            intersectSize1 += 1
            
            
            if snp[2] == snp2[2]:
                intersectMatchGT1 += 1
                if snp[1] != snp2[1]: print snp[0], "rsid doesn't match", snp[1], snp2[1]
                
            if snp[3] == snp2[3]:
                intersectMatchHT1 += 1
                currentSegLen += 1
                if revSegLen > 0:
                    print "len+ %4d\tmiss: %3d\tdist: %4d" % (revSegLen, missingStreak, pos-prev_pos)
                    revSegCount += 1
                    segLengthsSet += [revSegLen]
                    revSegMax = max(revSegMax, revSegLen)
                    revSegMin = min(revSegMin, revSegLen)
                    revSegLen = 0
            else:
                revSegLen += 1
                if currentSegLen > 0:
                    print "len- %4d\tmiss: %3d\tdist: %4d" % (currentSegLen, missingStreak, pos-prev_pos)
                    segCount += 1
                    segLengthsSet += [currentSegLen]
                    segMax = max(segMax, currentSegLen)
                    segMin = min(segMin, currentSegLen)
                    currentSegLen = 0
                #print snp[0], "HT mismatch", snp[3], snp2[3]
            missingStreak = 0
        except KeyError, e:
            missingStreak += 1
        prev_pos = pos
      
    if currentSegLen != 0:
        segCount += 1
        segLengthsSet += [currentSegLen]
        segMax = max(segMax, currentSegLen)
        segMin = min(segMin, currentSegLen)
    if revSegLen != 0:
        revSegCount += 1
        segLengthsSet += [revSegLen]
        revSegMax = max(revSegMax, revSegLen)
        revSegMin = min(revSegMin, revSegLen)
        
    #--------
    '''
    intersectSize2 = 0
    intersectMatchGT2 = 0
    intersectMatchHT2 = 0
    for pos in sorted(lociKT.keys()):
        snp2 = lociKT[pos]
        if snp2[0] in lociBG.keys():
            intersectSize2 += 1
            snp = lociBG[snp2[0]]
            if snp2[2] == snp[2]:
                intersectMatchGT2 += 1
                if snp2[1] != snp[1]: print snp2[0], "rsid doesn't match", snp[1], snp2[1]
            if snp2[3] == snp[3]:
                intersectMatchHT2 += 1
            else:
                pass
                #print snp[0], "HT mismatch", snp[3], snp2[3]
    '''
                
    print "uniqLoci\tintersectSize\tintersectMatchGT\tintersectMatchHT"
    print "%d\t%d\t%d\t%d" % (len(lociBG)-intersectSize1, intersectSize1, intersectMatchGT1, intersectMatchHT1)
    #print "%d\t%d\t%d\t%d" % (len(lociKT)-intersectSize2, intersectSize2, intersectMatchGT2, intersectMatchHT2)
    print len(lociKT)-intersectSize1
    print ""
    print "skippedInCentromere:", skipped_in_centromere
    print "segCount:", segCount
    print "avgSegLength:", intersectMatchHT1/float(max(segCount,1))
    print "segMin:", segMin
    print "segMax:", segMax
    print ""
    print "revSegCount:", revSegCount
    print "avgRevSegLength:", (intersectMatchGT1-intersectMatchHT1)/float(max(revSegCount,1))
    print "revSegMin:", revSegMin
    print "revSegMax:", revSegMax
    print ""
    segLengthsSet.sort()
    print "N50:", segLengthsSet[len(segLengthsSet)/2]
    
    plot = plt.hist(segLengthsSet, bins=50, color='blue')
    plt.show()
    plot.savefig('plot1.png')
    
if __name__ == '__main__':
    #import doctest
    #doctest.testmod()
    
    main()


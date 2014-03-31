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
    parser = argparse.ArgumentParser(description='Prepare SNP data for FCNV. Read filenames: for M, P, and F .vcf files; plasma, M, and P .sam files; and for centromeres list.')
    parser.add_argument('filenames', type=str, nargs='+', help='paths to 1) .vcf files with M, P, F SNPs; 2) reads in SAM format for plasma, M, and P samples; 3) centromeres list file.')
    args = parser.parse_args()
    
    if len(args.filenames) != 7: exit("Unexpected number of arguments passed! Expecting 7 filenames.")
    
    #treat these as CONSTANTS!
    M = 0; P = 1; F = 2; PLASMA = 3; MR = 4; PR = 5; CT = 6;
    ALL = [M, P, F, PLASMA, MR, PR, CT]
    
    #list of input files
    in_files = [open(args.filenames[i], "r" ) for i in ALL]
    
    #list of output files
    out_files = [None for i in [M, P, F, PLASMA]]
    out_files[M] = open("M_alleles.txt", "w")
    out_files[P] = open("P_alleles.txt", "w")
    out_files[F] = open("F_alleles.txt", "w")
    out_files[PLASMA] = open("plasma_samples.txt", "w")
    date = datetime.now().strftime('%m-%d-%H-%M')
    out_pos_file = open("positions" + date + ".txt", "w")
    
    #read centromeres positions
    centromeres = dict()
    for line in in_files[CT].readlines():
        line = line.rstrip('\n').split('\t')
        if line[0] not in centromeres.keys(): centromeres[line[0]] = []
        centromeres[line[0]] += [(int(line[1]), int(line[2]))]
        
    
    #allele counts in plasma samples for particular positions
    data = dict()
    loci = dict()
    processed_chr = ''
    
    print "  Getting union of SNP positions and corresponding list of alleles"
    #read SNPs from M, P, F vcf files
    snps = [[] for i in [M, P, F]]
    for i in [M, P, F]:
        #skip the header
        line = in_files[i].readline()
        while len(line) > 0 and line[0] == '#': line = in_files[i].readline()
        #split
        snps[i] = line.split('\t')
    
    #get genotypes for all positions in UNION of M and P SNP positions
    while len(snps[M])>2 or len(snps[P])>2: #while there is a SNP positions in M or P
        #get the position of SNP that occure first
        for i in [M, P, F]:
            if snps[i][0] == '': #if an input files is already at EOF
                snps[i][0] = 'chrZZ'
                snps[i].append(1e15)
            else: #convert to int
                snps[i][1] = int(snps[i][1])
        #chromosome
        min_chr = min(snps[M][0], snps[P][0])
        if processed_chr == '': processed_chr = min_chr
        if processed_chr != min_chr: print "WARNING: multiple chromosomes in the input", processed_chr, "|", min_chr
        #position
        min_pos = 1e15
        for i in [M, P]:
            if min_chr == snps[i][0] and snps[i][1] < min_pos:
                min_pos = snps[i][1]
        
        #get alleles
        alleles = [['.', '.'] for x in [M, P, F]]
        for i in [M, P, F]:
            #if there is a SNP in the data at this position, use it
            if min_chr == snps[i][0] and min_pos == snps[i][1]:
                alleles[i] = [snps[i][3][0], snps[i][4][0]]

            
        #if there is no info on some allele, impute the reference allele
        ref_allele = alleles[M][0]
        if ref_allele == '.': ref_allele = alleles[P][0]
        for i in [M, P, F]:
            for a in [0, 1]:
                if alleles[i][a] == '.': alleles[i][a] = ref_allele

        #check for homozygous alternative sites in Fetal VCF
        for i in [F]:
            #if there is a SNP in the data at this position
            if min_chr == snps[i][0] and min_pos == snps[i][1]:
                info = snps[i]
                if len(info) <= 2: continue
                #parse out genotype config info
                gt = info[9].split(':')[0]
                #if homozygous alternative
                if gt[0] == '1':
                    alleles[i][0] = snps[i][4][0]
                if gt[2] == '1':
                    alleles[i][1] = snps[i][4][0]
                if gt[0] == '0':
                    alleles[i][0] = snps[i][3][0]
                if gt[2] == '0':
                    alleles[i][1] = snps[i][3][0] 
                    
        #organize the haplotypes in M, P (phased VCF files)
        for i in [M, P]:
            #if there is a SNP in the data at this position
            if min_chr == snps[i][0] and min_pos == snps[i][1]:
                info = snps[i]
                if len(info) <= 2: continue
                #parse out haplotype config info
                ht = map(int, info[9].split('/'))
                #get the configuration
                phased_alleles = [alleles[i][ht[0]], alleles[i][ht[1]]]
                alleles[i] = phased_alleles
        
        #take note that for this position we need to get allele counts in plasma samaples
        loci[min_pos] = alleles
        sp.add_pos(min_pos, data)
            
        #read input: next SNP
        for i in [M, P, F]:
            if min_chr >= snps[i][0] and min_pos >= snps[i][1]:
                snps[i] = in_files[i].readline().split('\t')

        #END WHILE
    
    print "  Aligning the reads"
    #set up datastructures for counting allele support in diffrenct SAM files
    posInfo = [dict() for i in ALL]
    for R in [PLASMA, MR, PR]:
        posInfo[R] = copy.deepcopy(data)
        
    #fetch the reads in plasma SAM file and get counts for the positions originally specified in 'data'
    for R in [PLASMA, MR, PR]:
        while True:
            line = in_files[R].readline()
            if not line: break
            if len(line) > 0 and line[0] == '@': continue
            sp.pile_up(sp.mapping_parser(line), posInfo[R])    
    
    print "  Writing output"
    skipped_in_centromere = 0
    skipped_low = 0
    centromere_regions = centromeres[processed_chr]
    #print info / compute stats for each SNP position
    for pos in sorted(data.keys()):
        alleles = loci[pos]
        #if alleles[M][0] != alleles[M][1]:
        
        #print the plasma allele counts 
        nuc_counts = posInfo[PLASMA][pos]
        tmp = []
        for nuc in 'ACGT': #to make sure they are in the right order
            try:
                tmp.append(str(nuc_counts[nuc]))
            except KeyError:
                tmp.append('0')
                
        #if the plasma coverage is too low, skip this position        
        if sum(map(int, tmp)) < 20: 
            print pos, "- low overall coverage", sum(map(int, tmp))
            skipped_low += 1
            continue
        #if in centromere region, skip
        if is_within_intervals(pos, centromere_regions):
            skipped_in_centromere += 1
            continue
        
        print >>out_files[PLASMA], ' '.join(tmp)
        
        #output M, P, F alleles at this SNP locus
        for i, r in [(M, MR), (P, PR)]:
            a1 = alleles[i][0]
            a2 = alleles[i][1]
            count_a1 = 0
            count_a2 = 0
            try: count_a1 = posInfo[r][pos][a1]
            except: print i, pos, a1, posInfo[r][pos], alleles[i]
            try: count_a2 = posInfo[r][pos][a2]
            except: print i, pos, a2, posInfo[r][pos], alleles[i]
            
            if a1 == a2:
                count_a1 /= 2.
                count_a2 /= 2.
            
            print >>out_files[i], a1, a2, count_a1, count_a2
            
        print >>out_files[F], alleles[F][0], alleles[F][1], 3
        print >>out_pos_file, pos, "- M:", alleles[M], " P:", alleles[P], " F:", alleles[F]
     
    print "Low overall coverage positions ignored:", skipped_low
    print "Ignored positions in centromere regions:", skipped_in_centromere   
    
    
if __name__ == '__main__':
    #import doctest
    #doctest.testmod()
    
    main()


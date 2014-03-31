#!/usr/bin/python2

import argparse
import samParse as sp
import copy

def main():
    #parse ARGs
    parser = argparse.ArgumentParser(description='Analyze mixture and allele ratios for SNP positions in union(M, P). Read filenames for M, P .vcf files. Further M, F .sam files with reads that together form plasma reads.')
    parser.add_argument('filenames', type=str, nargs='+', help='paths to .vcf files with M, P SNPs and *sorted* M, F reads in SAM format')
    args = parser.parse_args()
    
    if len(args.filenames) != 4: die("Unexpected number of arguments passed! Expecting 4 filenames.")
    
    #treat these as CONSTANTS!
    M = 0; P = 1; MR = 2; FR = 3;
    ALL = [M, P, MR, FR]
    
    #list of input files
    in_files = [open(args.filenames[i], "r" ) for i in ALL]
    
    #list of output files
    #out_files = [None for i in ALL]
    #out_files[M] = open("M_alleles.txt", "w")

    
    #union of maternal and paternal SNP positions
    data = dict()
    loci = dict()
    
    #read SNPs from M, P, F vcf files
    snps = [[] for i in [M, P]]
    for i in [M, P]:
        #skip the header
        line = in_files[i].readline()
        while len(line) > 0 and line[0] == '#': line = in_files[i].readline()
        #split
        snps[i] = line.split('\t')
    
    #union the maternal and paternals SNPs positions
    while len(snps[M])>2 or len(snps[P])>2: #while there is a SNP positions in M or P
        #get the position of SNP that occure first
        for i in [M, P]:
            if snps[i][0] == '': #if an input files is already at EOF
                snps[i][0] = 'chrZZ'
                snps[i].append(1e15)
            else: #convert to int
                snps[i][1] = int(snps[i][1])
        #chromosome
        min_chr = min(snps[M][0], snps[P][0])
        #position
        min_pos = 1e15
        for i in [M, P]:
            if min_chr == snps[i][0] and snps[i][1] < min_pos:
                min_pos = snps[i][1]
        
        #get alleles
        alleles = [['.', '.'] for x in [M, P]]
        for i in [M, P]:
            #if there is a SNP in the data at this position, use it
            if min_chr == snps[i][0] and min_pos == snps[i][1]:
                alleles[i] = [snps[i][3], snps[i][4]]

            
        #if there is no info on some allele, impute the reference allele
        ref_allele = alleles[M][0]
        if ref_allele == '.': ref_allele = alleles[P][0]
        for i in [M, P]:
            for a in [0, 1]:
                if alleles[i][a] == '.': alleles[i][a] = ref_allele

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
        sp.add_pos(min_pos, data)
        #loci[min_pos] = alleles
        #print min_pos, ": M:", alleles[M], " P:", alleles[P], " F:", alleles[F]
            
        #read input: next SNP
        for i in [M, P]:
            if min_chr >= snps[i][0] and min_pos >= snps[i][1]:
                snps[i] = in_files[i].readline().split('\t')
    
    #fetch the maternal and fetal portion of plasma reads from SAM files,
    # and get counts for the positions specified in 'data'
    posInfo = [dict() for i in ALL]
    posInfo[MR] = copy.deepcopy(data)
    posInfo[FR] = copy.deepcopy(data)
    
    for R in [MR, FR]:
        while True:
            line = in_files[R].readline()
            if not line: break
            if len(line) > 0 and line[0] == '@': continue
            sp.pile_up(sp.mapping_parser(line), posInfo[R])

    #compute and print the stats
    for pos in sorted(data.keys()):
        MR_nuc_counts = posInfo[MR][pos]
        FR_nuc_counts = posInfo[FR][pos]
        try:
            local_mix_ratio = float(sum(FR_nuc_counts.values())) / (sum(FR_nuc_counts.values()) + 9*sum(MR_nuc_counts.values()))
        except ZeroDivisionError:
            local_mix_ratio = 0
        
        print pos, local_mix_ratio,
        '''
        for i, NC in enumerate([MR_nuc_counts, FR_nuc_counts]):
            tmp = []
            for nuc in 'ACGT': #to make sure they are in the right order
                try:
                    tmp.append(NC[nuc])
                except KeyError:
                    tmp.append(0)
            try:    
                summ = float(sum(tmp))
                tmp = [tmp[i]/summ for i in range(len(tmp))]
            except ZeroDivisionError:
                pass
                #tmp = [0]
            #print sorted(tmp),
            try:
                ind1 = 'ACGT'.index(loci[pos][i][0].upper())
                ind2 = 'ACGT'.index(loci[pos][i][1].upper())
                print abs(tmp[ind1] - tmp[ind2]), 
            except IndexError:
                print "X", loci[pos][i][0].upper(), loci[pos][i][1].upper(), tmp,
                print 0, 
        '''    
        print " "
        
    
if __name__ == '__main__':
    #import doctest
    #doctest.testmod()
    
    main()


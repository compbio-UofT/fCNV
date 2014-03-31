#!/usr/bin/python2

import argparse
from sys import exit
import subprocess
import copy
import samParse as sp


prev_lines = ['', '']
def getlineFromFile(in_files, out_files, idn):
    global prev_lines
    
    if prev_lines[idn] != '':
        print >>out_files[idn], prev_lines[idn],
    
    line = in_files[idn].readline()
    while len(line) > 0 and line[0] == '#':
        #ignore the header
        print >>out_files[idn], line,
        line = in_files[idn].readline()
    
    
    split_line = line.split('\t')
    if len(split_line)>=5 and len(split_line[4]) > 1:
        split_line[4] = split_line[4][0]
    line = '\t'.join(split_line)
    
    prev_lines[idn] = line
    return split_line
    

def main():
    global prev_lines
    #parse ARGs
    parser = argparse.ArgumentParser(description='Filter SNP positions by call quality and min. coverage. Awaits filenames for M, P .vcf files, and M, P .sam files.')
    parser.add_argument('filenames', type=str, nargs='+', help='paths to .vcf files with M, P SNPs and to corresponding .sam files')
    args = parser.parse_args()
    
    if len(args.filenames) != 4: exit("Unexpected number of arguments passed! Expecting 4 filenames.")
    
    #treat these as CONSTANTS!
    M = 0; P = 1; MR = 2; PR = 3;
    ALL_VCF = [M, P]
    ALL = [M, P, MR, PR]
    
    '''
    Get union of the M and P SNP positions
    '''
    print "  Stage 1"
    #list of input files
    in_files = [open(args.filenames[i], "r" ) for i in ALL]
    
    #dictionary to store the union of SNP positions
    data = dict()
    loci = dict()
    
    #read the first SNP from M, P vcf files
    snps = [[] for i in ALL_VCF]
    for i in ALL_VCF:
        #skip the header
        line = in_files[i].readline()
        while len(line) > 0 and line[0] == '#': line = in_files[i].readline()
        #split
        snps[i] = line.split('\t')
    
    #get list of all positions in UNION of M and P SNP positions
    while len(snps[M])>2 or len(snps[P])>2: #while there is a SNP positions in M or P
        #get the position of SNP that occure first
        for i in ALL_VCF:
            if snps[i][0] == '': #if an input files is already at EOF
                snps[i][0] = 'chrZZ'
                snps[i].append(1e15)
            else: #convert to int
                snps[i][1] = int(snps[i][1])
        #chromosome
        min_chr = min(snps[M][0], snps[P][0])
        #position
        min_pos = 1e15
        for i in ALL_VCF:
            if min_chr == snps[i][0] and snps[i][1] < min_pos:
                min_pos = snps[i][1]
        
        #get alleles
        alleles = [['.', '.'] for x in [M, P]]
        for i in [M, P]:
            #if there is a SNP in the data at this position, use it
            if min_chr == snps[i][0] and min_pos == snps[i][1]:
                alleles[i] = [snps[i][3][0], snps[i][4][0]]

        #if there is no info on some allele, impute the reference allele
        ref_allele = alleles[M][0]
        if ref_allele == '.': ref_allele = alleles[P][0]
        for i in [M, P]:
            for a in [0, 1]:
                if alleles[i][a] == '.': alleles[i][a] = ref_allele

        #check for homozygous alternative sites
        for i in [M, P]:
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
        
        #take note that for this position we need to get allele counts in plasma samaples
        loci[min_pos] = alleles
        sp.add_pos(min_pos, data)
        
        #read input: next SNP
        for i in ALL_VCF:
            if min_chr >= snps[i][0] and min_pos >= snps[i][1]:
                snps[i] = in_files[i].readline().split('\t')
        #END WHILE
        
        
    '''
    Get coverage information of the SNP positions from corresponding .sam files
    '''
    print "  Stage 2"
    #fetch allele support for the UNION positions in maternal and paternal reads
    #set up datastructures for counting allele support in diffrenct SAM files
    posInfo = [dict() for i in ALL]
    for R in [MR, PR]:
        posInfo[R] = copy.deepcopy(data)
    
    #fetch the reads in plasma SAM file and get counts for the positions originally specified in 'data'   
    for R in [MR, PR]:
        while True:
            line = in_files[R].readline()
            if not line: break
            if len(line) > 0 and line[0] == '@': continue
            sp.pile_up(sp.mapping_parser(line), posInfo[R])
                
    
    '''
    Filter the SNP positions according to call quality and coverage
    '''
    print "  Stage 3"
    #reopen VCF files
    for f in in_files: f.close()
    in_files = [open(args.filenames[i], "r" ) for i in ALL_VCF]
    #list of output files
    out_files = [open(args.filenames[i][:-3]+"ftr.vcf", "w") for i in ALL_VCF]
    
    #read the first SNP from M, P vcf files
    snps = [[] for i in ALL_VCF]
    for i in ALL_VCF:
        snps[i] = getlineFromFile(in_files, out_files, i)
        
    #positions ignored from the union of M and P
    ignored_pos = 0
    
    while len(snps[M])>2 or len(snps[P])>2: #while there is a SNP positions in M or P
        #get the position of SNP that occure first
        for i in ALL_VCF:
            if snps[i][0] == '': #if an input files is already at EOF
                snps[i][0] = 'chrZZ'
                snps[i].append(1e15)
            else: #convert to int
                snps[i][1] = int(snps[i][1])
        #chromosome
        min_chr = min(snps[M][0], snps[P][0])
        #position
        min_pos = 1e15
        for i in ALL_VCF:
            if min_chr == snps[i][0] and snps[i][1] < min_pos:
                min_pos = snps[i][1]
        
        #get genotype call quality
        callQ = [0. for i in ALL_VCF]
        for i in ALL_VCF:
            #if there is a SNP in the data at this position
            if min_chr == snps[i][0] and min_pos == snps[i][1]:
                info = snps[i]
                if len(info) <= 2: continue
                #parse out quality info
                callQ[i] = float(info[5])
        qualityOK = bool(callQ[M] >= 75 or callQ[P] >= 75)
        
        #get coverage info
        alleles = loci[min_pos]
        coverage = [0 for i in [M, P]]
        count_sum = [0 for i in [M, P]]
        for i in [M, P]:
            a1 = alleles[i][0]
            a2 = alleles[i][1]
            count_a1 = 0
            count_a2 = 0
            try: count_a1 = posInfo[i+2][min_pos][a1]
            except: print i, min_pos, a1, posInfo[i+2][min_pos], alleles[i]
            try: count_a2 = posInfo[i+2][min_pos][a2]
            except: print i, min_pos, a2, posInfo[i+2][min_pos], alleles[i]
            
            count_sum[i] = sum(posInfo[i+2][min_pos].values())
            coverage[i] = count_a1 + count_a2 #posInfo[i+2][min_pos][a1] + posInfo[i+2][min_pos][a2]
            """ #using mpileup to get the coverage info
            cmd = 'samtools mpileup -r %(chr)s:%(pos)d-%(pos)d __%(gnm)s.part.bam' % {'chr':min_chr, 'pos':min_pos, 'gnm':'MP'[i]}
            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = process.communicate()
            if process.returncode == 0:
                fields = out.split('\t')
                if len(fields) < 4:
                    #there is no coverage in the BAM file for this position
                    coverage[i] = 0
                    #print "!!! :", out, "|", err, "|", cmd
                else:
                    coverage[i] = int(fields[3])
                #print out, '>>coverage>>', coverage[i]
            else:
                print err
            """
        coverageOK = bool(coverage[M] >= 15 and coverage[P] >= 15)
        contaminationOK = False
        try: contaminationOK = bool(float(coverage[M])/count_sum[M] >= 0.9 and float(coverage[P])/count_sum[P] >= 0.9)
        except: contaminationOK = False
            
        if coverageOK and not contaminationOK:
            print min_pos, "- contamination M:", coverage[M], posInfo[M+2][min_pos], alleles[M], "P:", coverage[P], posInfo[P+2][min_pos], alleles[P]
        
        #MOK = bool(callQ[M] >= 75 or coverage[M] >= 15)
        #POK = bool(callQ[P] >= 75 or coverage[P] >= 15)
        #if not (MOK and POK): #ignore positions that are not good enough
        if not (qualityOK and coverageOK and contaminationOK): #ignore positions that are not good enough
            ignored_pos += 1
            for i in ALL_VCF:
                #if there is a SNP in the data at this position, skip it
                if min_chr == snps[i][0] and min_pos == snps[i][1]:
                    #print min_pos, callQ[M], callQ[P], coverage[M], coverage[P],  prev_lines[i],
                    prev_lines[i] = ''      
            
        
        #read input: next SNP
        for i in ALL_VCF:
            if min_chr >= snps[i][0] and min_pos >= snps[i][1]:
                snps[i] = getlineFromFile(in_files, out_files, i)
        
        #END WHILE
    
    print "Low quality positions ignored in the region:", ignored_pos  
    
if __name__ == '__main__':
    #import doctest
    #doctest.testmod()
    
    main()


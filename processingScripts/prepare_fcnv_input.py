#!/usr/bin/python2
#example usage: time prepare_fcnv_input.py mp.phase.vcf __plasma.dup100kb.sam __M.part.sam __P.part.sam /dupa-filer/laci/centromeres >log_prepare_fcnv_input.txt 2>>log_prepare_fcnv_input.txt &

import argparse
from sys import exit
import random
import copy
from datetime import datetime
import subprocess
import os

tmp_dir = '/tmp/'
if os.environ['TMPDIR']: tmp_dir = os.environ['TMPDIR']

def is_within_intervals(num, intervals):
    for interval in intervals:
        if interval[0] <= num <= interval[1]:
            return True
    return False

def main():
    #parse ARGs
    parser = argparse.ArgumentParser(description='Prepare SNP support data for FCNV. Read filenames: for joined M&P phased .vcf file; plasma, M, and P .bam files; and for centromeres list.')
    parser.add_argument('filenames', type=str, nargs='+', help='paths to 1) .vcf file with phased M & P SNPs; 2) reads in SAM format for plasma, M, and P samples; 3) valid regions BED file 4) centromeres list file; 5) result path.')
    args = parser.parse_args()
    
    if len(args.filenames) != 7: exit("Unexpected number of arguments passed! Expecting 7 parameters.")
    
    #treat these as CONSTANTS!
    MP = 0; PLR = 1; MR = 2; PR = 3; BED = 4; CT = 5; RES_PATH=6; #in_files
    ALL = [MP, PLR, MR, BED, PR, CT]
    M = 0; P = 1; #maternal, paternal
    ALDOC = 0; GT = 1; #out_files: allele DOC and ground truth
    
    #list of input files
    in_files = [None for i in ALL]
    in_files[MP] = open(args.filenames[MP], "r" )
    in_files[CT] = open(args.filenames[CT], "r" )
    res_path = args.filenames[RES_PATH]
    
    plasma_id = args.filenames[PLR][:-4].replace(':', '-').replace('.sort', '')
    plasma_path = '/'.join(plasma_id.split('/')[:-1])
    if len(plasma_path) != 0: plasma_path += '/'
    plasma_id = plasma_id.split('/')[-1]
    tmp_pos_file_name = tmp_dir + "/__tmp" + plasma_id + "_snp_pos.txt"
    tmp_pos_file = open(tmp_pos_file_name, "w")
    
    #parse CNV type and position from the plasma .bam file name
    try:
        plasmaFN = args.filenames[PLR].split('/')[-1]
        reg_begin = int(plasmaFN.split(':')[1].split('-')[0])
        reg_end = int(plasmaFN.split(':')[1].split('-')[1])
        
        maternal_origin = False
        paternal_origin = False
        if plasmaFN.split('-')[0].lower().find('m') != -1:
            maternal_origin = True
        elif plasmaFN.split('-')[0].lower().find('p') != -1:
            paternal_origin = True
        
        #[(0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1)]
        cnv_state_id = 9    
        if plasmaFN.find('duplicate') != -1:
            if maternal_origin: cnv_state_id = 6
            elif paternal_origin: cnv_state_id = 4
        if plasmaFN.find('delete') != -1:
            cnv_state_id = 0
    except:
        reg_begin = 0
        reg_end = 0
        cnv_state_id = 9
    
    #read centromeres positions
    centromeres = dict()
    for line in in_files[CT].readlines():
        line = line.rstrip('\n').split('\t')
        if line[0] not in centromeres.keys(): centromeres[line[0]] = []
        centromeres[line[0]] += [(int(line[1]), int(line[2]))]
    
    #allele counts in plasma samples for particular positions
    loci = dict()
    processed_chr = ''
    skipped_in_centromere = 0
    
    print "  Getting union of SNP positions and corresponding list of alleles " + datetime.now().strftime('%m-%d-%H-%M')
    #read SNPs from M, P, F vcf files
    snps = [[] for i in [M, P]]
    
    #get genotypes for all positions in UNION of M and P SNP positions
    while True:
        line = in_files[MP].readline()
        #skip if part of the header
        while len(line) > 0 and line[0] == '#': line = in_files[MP].readline()
        if not line: break
        
        fields = line.rstrip('\n').split('\t')
        if processed_chr == '': processed_chr = 'chr'+fields[0]
        if processed_chr != 'chr'+fields[0]:
            print "WARNING: multiple chromosomes in the input", processed_chr, "|", 'chr'+fields[0]
        pos = int(fields[1])
        ref = fields[3]
        alt = fields[4]
        
        #get M and P haplotypes
        snps[M] = fields[9].split(':')[0].split('|')
        for x in [0, 1]:
            if snps[M][x] == '0': snps[M][x] = ref
            else: snps[M][x] = alt
            
        snps[P] = fields[10].split(':')[0].split('|')
        for x in [0, 1]:
            if snps[P][x] == '0': snps[P][x] = ref
            else: snps[P][x] = alt
        
        #if in centromere region, skip
        centromere_regions = centromeres[processed_chr]
        if is_within_intervals(pos, centromere_regions):
            skipped_in_centromere += 1
            continue
        
        #get +/-20MB window begin and end
        wbegin = reg_begin - 20000000
        wend = reg_end + 20000000
        tel1_end = centromeres[processed_chr][0][1]
        cent_beg = centromeres[processed_chr][1][0]
        cent_end = centromeres[processed_chr][1][1]
        tel2_beg = centromeres[processed_chr][2][0]
        #check centromiers for wbegin 
        if wbegin <= tel1_end: 
            wend = wend + tel1_end - wbegin
            wbegin = tel1_end
        if cent_beg <= wbegin <= cent_end:
            wend = wend + cent_end - wbegin
            wbegin = cent_beg
            
        #check centromiers for wend
        if wend >= tel2_beg:
            wbegin = wbegin - (wend-tel2_beg)
            wend = tel2_beg
        if cent_beg <= wend <= cent_end:
            wbegin = wbegin - (wend-cent_beg)    
            wend = cent_beg
        
        if not (wbegin <= pos <= wend): continue
        
        #take note that for this position we need to get allele counts in plasma samaples
        alleles = (snps[M], snps[P], (ref, alt))
        loci[pos] = alleles
        print >>tmp_pos_file, processed_chr, pos
    #END WHILE
    tmp_pos_file.close()
    
    print "  Piling up the reads " + datetime.now().strftime('%m-%d-%H-%M')
    #call samtools mpileup to get allele counts for positions in 'loci'
    #cdir = os.getcwd() + '/'
    tmp_vcf_prefix = tmp_dir + '/__tmp' + plasma_id
    pile_prefix = res_path + '/' + plasma_id
    cmd = "span_samtools.sh /filer/hg19/hg19.fa {0} {1} {2} {3} {4} {5} {6} {7} {8}".format(tmp_pos_file_name, \
        args.filenames[PLR], args.filenames[MR], args.filenames[PR], tmp_vcf_prefix, wbegin, wend, processed_chr, pile_prefix)
    os.system(cmd)
    
    posInfo = [dict() for i in range(4)]
    for R in [PLR, MR, PR]:
        tmp_vcf_name = tmp_vcf_prefix + "." + str(R) + ".vcf"
        vcf_file = open(tmp_vcf_name, 'r')
        while True:
            line = vcf_file.readline()
            if not line: break
            if len(line) > 0 and line[0] == '#': continue
            line = line.rstrip('\n').split('\t')
            
            ac = {'A':0, 'C':0, 'G':0, 'T':0}
            pos = int(line[1])
            ref = loci[pos][2][0]
            alt = loci[pos][2][1]
            if ref != line[3]: print "SOMETHING IS WRONG WITH TEMP VCF POSITIONS" 
            for x in line[7].split(';'):
                if len(x) > 3 and x[0:3]=='DP4':
                    counts = map(int, x[4:].split(','))
                    ac[ref] = counts[0] + counts[1]
                    ac[alt] = counts[2] + counts[3]
            posInfo[R][pos] = ac
            
        if len(loci) != len(posInfo[R]): print "DIFFERENT NUMBER OF POSITIONS IN TEMP VCF " + str(R) 
        vcf_file.close()    
        os.remove(tmp_vcf_name)
    os.remove(tmp_pos_file_name)
        
    print "  Writing output " + datetime.now().strftime('%m-%d-%H-%M')
    #list of output files
    out_files = [None for i in [ALDOC, GT]]
    out_files[ALDOC] = open(res_path + '/' + plasma_id + ".alleles_doc.txt", "w")
    out_files[GT] = open(res_path + '/' + plasma_id + ".target.txt", "w")
    print >>out_files[ALDOC], '#POS\tA\tC\tG\tT\tM_hapA\tM_hapB\tDP_hapA\tDP_hapB\tP_hapA\tP_hapB\tDP_hapA\tDP_hapB'
    
    skipped_low_doc = 0
    forced_to_ignore = 0
    #print info / compute stats for each SNP position
    for pos in sorted(loci.keys()):
        alleles = loci[pos]
        
        #print the plasma allele counts 
        if pos not in posInfo[PLR]:
            print "FORCED TO IGNORE POS:", pos
            forced_to_ignore += 1
            continue
        nuc_counts = posInfo[PLR][pos]
        tmp = []
        for nuc in 'ACGT': #to make sure they are in the right order
            try:
                tmp.append(str(nuc_counts[nuc]))
            except KeyError:
                tmp.append('0')
                
        #if the plasma coverage is too low, skip this position        
        if sum(map(int, tmp)) < 30: 
            print pos, "- low overall coverage", sum(map(int, tmp))
            skipped_low_doc += 1
            continue
        
        out_str = str(pos) + '\t' + '\t'.join(tmp)
        
        #output M, P alleles at this SNP locus
        error = False
        for i, r in [(M, MR), (P, PR)]:
            a1 = alleles[i][0]
            a2 = alleles[i][1]
            count_a1 = 0
            count_a2 = 0
            try:
                count_a1 = posInfo[r][pos][a1]
            except KeyError:
                error = True
                break
                #print i, pos, a1, posInfo[r][pos], alleles[i]
            try:
                count_a2 = posInfo[r][pos][a2]
            except KeyError:
                error = True
                break
                #print i, pos, a2, posInfo[r][pos], alleles[i]
            
            if a1 == a2:
                count_a1 /= 2.
                count_a2 /= 2.
            
            out_str += '\t{0}\t{1}\t{2}\t{3}'.format(a1, a2, count_a1, count_a2)
        if error: continue
        
        print >>out_files[ALDOC], out_str
        
        if reg_begin <= pos and pos <= reg_end:
            print >>out_files[GT], '{0}\t{1}\t{2}\t{3}'.format(pos, 'N', 'N', cnv_state_id)
        else:
            print >>out_files[GT], '{0}\t{1}\t{2}\t{3}'.format(pos, 'N', 'N', 3)
    
    out_files[ALDOC].close()
    out_files[GT].close() 
    print "Low overall coverage positions ignored:", skipped_low_doc
    print "Forced to ignore due to missing in plasma:", forced_to_ignore
    print "Ignored positions in centromere regions:", skipped_in_centromere   
    print "DONE " + datetime.now().strftime('%m-%d-%H-%M')
    
if __name__ == '__main__':
    main()


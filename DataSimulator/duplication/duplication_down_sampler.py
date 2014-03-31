#!/usr/bin/pypy
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import argparse
import random

UNMAPPED = 0x4

def mapping_parser(m):
    '''
    Parse a read in SAM format, return a dictionary with filled in fields of interest.
    '''
    if isinstance(m, str):
        m = m.strip().split('\t')
        d = {}
        d['flag'] = int(m[1])   # flags
        d['chr'] = m[2]         # chr
        d['pos'] = int(m[3])    # pos
        d['mapq'] = int(m[4])   # mapping quality
        d['cigar'] = m[5]       # cigar string
        d['seq'] = m[9]         # sequence
        d['qual'] = m[10]       # sequencing quality

    return d

def main():
    #duplication_down_sampler.py  $filtered_res_file $numReads $readLength $tmp_pileup_file $plasmaFetusRate $region
    parser = argparse.ArgumentParser(description='Down sample reads from the target haplotype to DOC for mixing with plasma reads.')
    parser.add_argument('hapReadsFile', type=str, nargs=1, help='path to SAM file with filtered target haplotype reads')
    parser.add_argument('numReads', type=int, nargs=1, help='number of reads in the hapReadsFile')
    parser.add_argument('readLength', type=int, nargs=1, help='read length')
    parser.add_argument('plasmaPileup', type=str, nargs=1, help='path to file with piled up plasma reads for the target region')
    parser.add_argument('fetalRate', type=float, nargs=1, help='fetal fraction in plasma')
    parser.add_argument('region', type=str, nargs=1, help='coordinates of target region for duplication')
    
    
    args = parser.parse_args()

    hap_reads_file = open(args.hapReadsFile[0], "r")
    plasma_doc_file = open(args.plasmaPileup[0], "r")
    
    numReads = int(args.numReads[0])
    readLength = int(args.readLength[0])
    ffmix = float(args.fetalRate[0])
    region = map(int, args.region[0].split(':')[1].split('-'))
    
    #read the plasma piled up info and compute prefix sums
    prefix_sum = [0] * (region[1] - region[0] + 2342)
    prefix_count = [0] * (region[1] - region[0] + 2342)
    offset = -1
    last = 0
    for line in plasma_doc_file:
        row = map(int, line.split(' '))
        if offset == -1:
            offset = row[0]
        row[0] -= offset
        
        prefix_sum[row[0]] = prefix_sum[last] + row[1]
        prefix_count[row[0]] = prefix_count[last] + 1
        last = row[0]
    plasma_doc_file.close()
    
    #average DOC of target haplotype in filtered reads
    hapDOC = numReads * readLength / float(region[1] - region[0])
    
    #read the SAM file and down sample the reads
    for line in hap_reads_file:
        read = mapping_parser(line)

        # If the read is not aligned, ignore it
        if read['flag'] & UNMAPPED != 0: continue
        
        begin = max(0, read['pos'] - 300 - offset)
        end = begin + readLength + 300 
        
        try:
            while prefix_sum[begin]==0: begin += 1
            while end > begin and prefix_sum[end]==0: end -= 1
        except:
            if not plasmaDOC: plasmaDOC = 70
        
        try:
            plasmaDOC = (prefix_sum[end]-prefix_sum[begin]) / (prefix_count[end]-prefix_count[begin])
        except:
            if not plasmaDOC: plasmaDOC = 70
        
        targetDOC = plasmaDOC * ffmix / 2.
        outputRate = targetDOC / hapDOC
        if random.random() < outputRate:             
            print line,
    

if __name__ == '__main__':
    main()


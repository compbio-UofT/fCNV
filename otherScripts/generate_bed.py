#!/usr/bin/python2

import argparse
from sys import exit
import copy
   

def main():
    #parse ARGs
    parser = argparse.ArgumentParser(description='Generate BED track files for visualization of data and predictions.')
    parser.add_argument('--chr', type=str, nargs=1, help='chromosome number')
    parser.add_argument('--filenames', type=str, nargs='+', help='paths to .txt files with positions summary, M, P, F alleles counts and predictions.')
    args = parser.parse_args()
    
    #if len(args.filenames) != 5: exit("Unexpected number of arguments passed! Expecting 5 filenames.")
        
    chrNo = args.chr[0]
    
    #treat these as CONSTANTS!
    S = 0; M = 1; P = 2; F = 3; A = 4;
    ALL = [S, M, P, F, A] #, M, P, F, A]
    
    #list of input files
    in_files = [open(args.filenames[i], "r" ) for i in ALL]
    
    #list of output files
    ST = 0; DT = 1; AT = 2;
    out_files = [None for i in [ST, DT, AT]]
    out_files[ST] = open("track_summary.txt", "w")
    out_files[DT] = open("track_details.txt", "w")
    out_files[AT] = open("track_annotation.txt", "w")
    
    
    prev_ip = -1
    prev_breakpoint = 1
    colors = ['255,0,0', '255,128,0', '255,255,0', '0,255,0', '0,255,255', '0,0,255', '128,128,128']
    
    summary_lines = in_files[S].readlines()
    for sline in summary_lines:
        sline = sline.split('-')
        pos = int(sline[0])
        summary = sline[1].strip(" ").rstrip('\n')
        summary = summary.replace("'", "")
        summary = summary.replace(": ", ":")
        
        
        #predictions and anotations
        #1 (51, 0, 0, 0) ('A', 'A') [5.75, 5.75] ('A', 'A') [7.25, 7.25] 4 4 [(6, -5.30695), (5, -5.30695), (4, -5.30695), (3, -5.30695), (2, -5.30695), (1, -5.30695), (0, -5.30695)]
        line = in_files[A].readline().rstrip('\n')
        plasma_counts = 'plasma: ' + line.split(')')[0].split('(')[1]
        annot = line.split('[')[-2].split(']')[1].strip(' ').rstrip(' ')
        
        details = line.split('[')[-1].rstrip(']')
        out_files[DT].write("chr%s\t%d\t%d\t%s\t0\t+\t%d\t%d\t255,0,0\n" % (chrNo, pos-1, pos, details, pos-1, pos))  
        
        ip = int(annot[0])
        if prev_ip == -1: prev_ip = ip
        if ip != prev_ip:
            out_files[AT].write("chr%s\t%d\t%d\t%s\t0\t+\t%d\t%d\t%s\n" % (chrNo, prev_breakpoint-1, pos-1, "IP: "+str(prev_ip), prev_breakpoint-1, pos-1, colors[prev_ip]))
            prev_ip = ip
            prev_breakpoint = pos
        
        out_files[ST].write("chr%s\t%d\t%d\t%s\t0\t+\t%d\t%d\t255,0,0\n" % (chrNo, pos-1, pos, 'IP: '+annot, pos-1, pos))
        out_files[ST].write("chr%s\t%d\t%d\t%s\t0\t+\t%d\t%d\t0,0,0\n" % (chrNo, pos-1, pos, plasma_counts, pos-1, pos))  
        
        
        #maternal and paternal alleles&counts
        for (g, gLetter) in [(P, 'P'), (M, 'M')]:
            line = in_files[g].readline()
            line = line.rstrip('\n').split(' ')
            if line[0] != line[1]:
                txt = gLetter + ': ' + line[0] + '=' + line[2] + ' ' + line[1] + '=' + line[3]
            else:
                txt = gLetter + ': ' + line[0] + '=' + str(int(float(line[2]) + float(line[3])))
            out_files[ST].write("chr%s\t%d\t%d\t%s\t0\t+\t%d\t%d\t0,0,0\n" % (chrNo, pos-1, pos, txt, pos-1, pos))

        out_files[ST].write("chr%s\t%d\t%d\t%s\t0\t+\t%d\t%d\t255,0,0\n" % (chrNo, pos-1, pos, summary, pos-1, pos))
    
    out_files[AT].write("chr%s\t%d\t%d\t%s\t0\t+\t%d\t%d\t%s\n" % (chrNo, prev_breakpoint-1, prev_breakpoint+50, "IP: "+str(prev_ip), prev_breakpoint-1, prev_breakpoint+50, colors[prev_ip]))
    
if __name__ == '__main__':
    #import doctest
    #doctest.testmod()
    
    main()


#!/usr/bin/pypy
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import argparse
from random import randint 
from tempfile import mkstemp
from shutil import move
from os import remove, close

def main():
    parser = argparse.ArgumentParser(description='Finding the origin of deletions')
    parser.add_argument('vcfFile', type=str, nargs=1, help='path to the .vcf  phasing file')
    parser.add_argument('targetFile', type=str, nargs=1, help='path to the .target.txt file')
    parser.add_argument('--targetHap', type=str, help='target haplotype')
    parser.add_argument('--aBegin', type=int, help='begin of the region')
    parser.add_argument('--aEnd', type=int, help='end of the region')
    args = parser.parse_args()
    
    vcfFile = open(args.vcfFile[0], "r")
    targetFile = open(args.targetFile[0], "r")
    
    tmpFileName = '/'.join(args.targetFile[0].split('/')[:-1])
    if len(tmpFileName) != 0: tmpFileName += '/'
    tmpFileName += '__' + args.targetFile[0].split('/')[-1] + '.tmp'
    tmpFile = open(tmpFileName, "w")
    
    #if the target haplotype and deletion position are not specified explicitly,
    #try to parse it from the targetFile name
    if not (args.targetHap and args.aBegin and args.aEnd):
        try:
            name = args.targetFile[0].split('/')[-1]
            begin = int(name.split('-')[2])
            end = int(name.split('-')[3])
            hap = name[0].upper()
            tHaplotype=ord(hap)-ord('A')
        except:
            print "Error: target haplotype and/or position can't be parsed from .target.txt file name"
            return 1
        if tHaplotype not in [0, 1]:
            print "Error: target haplotype and/or position can't be parsed from .target.txt file name"
            return 1
    else:
        begin = args.aBegin[0]
        end = args.aEnd[0]
        tHaplotype = ord(args.targetHap[0])-ord('A')
        
    counts = {'PA':0, 'PB':0, 'MA':0, 'MB':0}
    ans = {'PA':0, 'PB':0, 'MA':0, 'MB':0}
    ct = 0
    haplos = ['PA', 'PB', 'MA', 'MB']

    for row in vcfFile:
        parts = row.strip().split('\t')
        
        if parts[0][0]=='#' :
            continue
        
        if ct==100 or (ct > 0 and int(parts[1])>end):
            mx=-1
            ct=0
            for haplo in haplos:
                if counts[haplo] == mx:
                    match.append(haplo)
                if counts[haplo] > mx:
                    mx = counts[haplo]
                    match = [haplo]
                counts[haplo] = 0
            lucky = randint(0,len(match)-1)
            ans[match[len(match)-1]] += 1
            print '['+match[len(match)-1]+']',
        
        if int(parts[1])<begin or int(parts[1])>end:
            continue
        
        ct += 1
        maternal = parts[9].strip().split(':')[0].strip().split('|')
        paternal = parts[10].strip().split(':')[0].strip().split('|')
        fetal = parts[11].strip().split(':')[0].strip().split('|')
        if fetal[tHaplotype]==paternal[0]:
            counts['PA'] += 1
        if fetal[tHaplotype]==paternal[1]:
            counts['PB'] += 1
        if fetal[tHaplotype]==maternal[0]:
            counts['MA'] += 1
        if fetal[tHaplotype]==maternal[1]:
            counts['MB'] += 1
            
    print ' '

#    print 'MA: ', ans['MA']
#    print 'MB: ', ans['MB']
#    print 'PA: ', ans['PA']
#    print 'PB: ', ans['PB']

    if ans['PA']+ans['PB'] > ans['MA']+ans['MB']:
        subs='2'
    else:
        subs='0'

    for row in targetFile:
        parts=row.strip().split('\t')
        if parts[3]!='3':
            parts[3]=subs
        for i in range(3):
            tmpFile.write(parts[i]+'\t')
        tmpFile.write(parts[3]+'\n')

    remove(args.targetFile[0])
    move(tmpFileName,args.targetFile[0])


if __name__=='__main__':
    main()
        

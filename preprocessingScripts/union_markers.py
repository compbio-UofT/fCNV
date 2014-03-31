#!/usr/bin/python2

import argparse
from sys import exit

def main():
    parser = argparse.ArgumentParser(description='Takes two .markers files, returns markers from the first file agumented by alleles from the second file.')
    parser.add_argument('filenames', type=str, nargs='+', help='paths to .markers files')
    args = parser.parse_args()
    
    if len(args.filenames) != 2: exit("No enough arguments: missing file names!")

    filename1 = args.filenames[0]
    filename2 = args.filenames[1]
    
    f1 = open(filename1, 'r')
    f2 = open(filename2, 'r')

    mem = {}
    #read the markers from the second file and store the alleles to dictionary
    while True:
        line = f2.readline()
        if not line: break
        fields = line.rstrip("\n").split("\t")
        mem[fields[0]] = tuple(fields[2:4])

    while True:
        line = f1.readline().rstrip("\n")
        if not line: break
        fields = line.split("\t")
        try:
            malleles = mem[fields[0]]
            toAppend = ""
            if malleles[0] != fields[2]: 
                toAppend += malleles[0]
            if malleles[1] != fields[3]:
                toAppend += malleles[1]
            if toAppend != "":
                for x in toAppend:
                    if x!=fields[2] and x!=fields[3]:
                        line += "\t" + x
            print line
            
        except KeyError:
            print line

if __name__ == "__main__":
    import doctest
    #doctest.testmod()
    
    main()



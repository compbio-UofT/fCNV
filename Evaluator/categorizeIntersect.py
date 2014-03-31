#!/usr/bin/python2

import argparse

k=1000
binSizes = [0, 100, 300, 500, 1000, 5000, 10000] 
numBins = 7

binSizes= map(lambda x: x*k, binSizes)

def findIndex(sz):
    difs=map(lambda x: abs(sz-x), binSizes)
    return difs.index(min(difs))

def main():
    parser = argparse.ArgumentParser(description='Categorizes intersect of CNVnator calls with prediction by size of overlap.')
    parser.add_argument('intersectFile', type=str, nargs=1, help='path to bed intersect file')
    args = parser.parse_args()
    
    intersect_file=open(args.intersectFile[0], "r")
    
    predictions = {}
    
    for line in intersect_file:
        line = line.rstrip().split('\t')
        if not line[0].startswith('chr'): continue
        
        call_id = (line[0], line[1], line[2])
        if call_id not in predictions: predictions[call_id] = 0
        predictions[call_id] += int(line[-1])
        
    
    categories = [0]*numBins
    for call_id, overlap in predictions.iteritems():
        #print call_id, overlap, findIndex(overlap), binSizes[findIndex(overlap)]
        if findIndex(overlap) == findIndex(int(call_id[2])-int(call_id[1])):
            categories[findIndex(overlap)] += 1
    
    for sz in binSizes:
        print "{0}k\t".format(sz/k),
    print ""
    for num in categories:
        print "{0}\t".format(num),
    print ""
    

if __name__ == '__main__':
    main()
    
    

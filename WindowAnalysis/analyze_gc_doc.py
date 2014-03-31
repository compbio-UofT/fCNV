#!/usr/bin/pypy

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import argparse
import math

def main():
    parser = argparse.ArgumentParser(description='DOC statistics per windows of different size.')
    parser.add_argument('filenames', type=str, nargs='+', help='paths to: 1) text file with "pos coverage" per row info; 2) reference seq in fasta format')
    args = parser.parse_args()
    if len(args.filenames) != 2: exit("Unexpected number of arguments passed! Expecting 2 filenames.")
    
    pile_file = open(args.filenames[0], "r")
    ref_file = open(args.filenames[1], "r")
    
    gc_sum = [0]*70000000
    prefix_sum = [0]*70000000
    prefix_count = [0]*70000000
    
    #get GC content prefix sums from the reference 
    gen_pos = 0
    while True:
        line = ref_file.readline().strip().upper()
        if len(line) == 0: break
        if line[0] == '>': continue
        for i in range(len(line)):
            gc_sum[gen_pos] = gc_sum[max(gen_pos - 1, 0)]
            if line[i] in 'GC': gc_sum[gen_pos] += 1
            gen_pos += 1 
    
    #read the coverage data
    last = 0
    for position in pile_file:
        row=map(int, position.strip().split(' '))
        prefix_sum[row[0]] = prefix_sum[last] + row[1]
        prefix_count[row[0]] = prefix_count[last] + 1
        last = row[0]
        #print prefix_sum[last],
        #print prefix_count[last]


    for l in [10, 100, 1000, 10000, 100000, 1000000]:
        #create bins by GC content ratio
        bin_count = 40
        mean = [0.]*bin_count
        variance = [0.]*bin_count
        win_count = [0]*bin_count
        
        #get the DOC means for GC-ratio bins
        for i in range(len(prefix_sum)):
            if prefix_sum[i]==0 or prefix_sum[i+l]==0:
                continue
            tmp_mean = (prefix_sum[i+l]-prefix_sum[i]) / (prefix_count[i+l]-prefix_count[i]) * l / 200.
            
            gc_ratio = (gc_sum[i+l] - gc_sum[i]) / float(l)
            index = min(int(math.floor(gc_ratio * bin_count)), bin_count-1)
            
            mean[index] += tmp_mean
            win_count[index] += 1
        
        for i in range(bin_count):
            try: mean[i] /= win_count[i]
            except: mean[i] = 0.

        #calculate the variance of windows' DOC for each GC ratio
        for i in range(len(prefix_sum)):
            if prefix_sum[i]==0 or prefix_sum[i+l]==0:
                continue
            doc = (prefix_sum[i+l]-prefix_sum[i]) / (prefix_count[i+l]-prefix_count[i]) * l / 200.

            gc_ratio = (gc_sum[i+l] - gc_sum[i]) / float(l)
            index = min(int(math.floor(gc_ratio * bin_count)), bin_count-1)
            
            variance[index] += (doc - mean[index])**2
            
        for i in range(bin_count):
            try: variance[i] /= win_count[i]
            except: variance[i] = 0.
            
        #print the stats
        print '-------- length=', l
        print 'GC:',
        for i in range(15, bin_count-15): print "\t(%.2f,%.2f)" % (i*1./bin_count, (i+1)*1./bin_count),
        print ''
        print 'mean:',
        for i in range(15, bin_count-15): print "\t%.5f" % mean[i],
        print ''
        print 'dev:',
        for i in range(15, bin_count-15): print "\t%.5f" % math.sqrt(variance[i]),
        print ''
        print 'var:',
        for i in range(15, bin_count-15): print "\t%.5f" % variance[i],
        print ''
        print 'wc:',
        for i in range(15, bin_count-15): print "\t%.5f" % win_count[i],
        print ''
            


if __name__ == '__main__':
    main()


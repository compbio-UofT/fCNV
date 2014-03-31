#!/usr/bin/pypy
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import argparse

k=1000
binSizes = [0, 100, 300, 500, 1000, 5000, 10000] 
numBins = 7

truePos= {'MDup': [0]*numBins, 'MDel': [0]*numBins,'PDup': [0]*numBins,'PDel': [0]*numBins, 'sum':[0]*numBins}
falsePos={'MDup': [0]*numBins, 'MDel': [0]*numBins,'PDup': [0]*numBins,'PDel': [0]*numBins, 'sum':[0]*numBins} 
total= {'MDup': [0]*numBins, 'MDel': [0]*numBins,'PDup': [0]*numBins,'PDel': [0]*numBins, 'sum':[0]*numBins}

binSizes= map(lambda x: x*k, binSizes)

def findIndex(sz):
    difs=map(lambda x: abs(sz-x), binSizes)
    return difs.index(min(difs))

def main():
    parser = argparse.ArgumentParser(description='This script filters the given reads of a certain region of plasma so that it simulates a deletion in this region. It requires the read file, snips file, the rate of fetus DNA in the plasma and the target haplotype of fetus for deletion')
    parser.add_argument('summaryFile', type=str, nargs=1, help='the address to the summary file')
    parser.add_argument('dgtFile', type=str, nargs=1, help='path to dataset ground truth file')
    args = parser.parse_args()

    summary_file=open(args.summaryFile[0], "r")
    
    dgt_file=open(args.dgtFile[0], "r")
    DGT={}
    for line in dgt_file:
        line = line.split(" ")
        DGT[line[0]] = int(line[1])
    
    nameMap={}
    
    stateToType=['MDel', '', 'PDel', '', 'PDup', '', 'MDup']
    
    rawHeader = ''
    for line in summary_file:
        if line.find(".txt")!=-1:
            #if len(rawHeader) > 0 and rawHeader not in nameMap.keys():
            #    print '>>UNCATEGORIZED: ', rawHeader,
            #    #nameMap[rawHeader]='PDel'
            rawHeader=line
            fname = line.split(".")[0]
            if fname not in DGT:
                print '>>UNCATEGORIZED: ', rawHeader,
                continue
            if (DGT[fname]==2):
                nameMap[rawHeader]='PDel'
            if (DGT[fname]==0):
                nameMap[rawHeader]='MDel'
            if (DGT[fname]==6):
                nameMap[rawHeader]='MDup'
            if (DGT[fname]==4):
                nameMap[rawHeader]='PDup'
                
#        if line.find('Real State:')!=-1:
#            parts=line.split(' ')
#            if (rawHeader.find('cvrg')==-1 and int(parts[2])==2) or (rawHeader.find('cvrg')!=-1 and int(parts[2])==0):
#                nameMap[rawHeader]='PDel'
#            if (rawHeader.find('cvrg')==-1 and int(parts[2])==0):
#                nameMap[rawHeader]='MDel'
#            if (rawHeader.find('cvrg')==-1 and int(parts[2])==6) or (rawHeader.find('cvrg')!=-1 and int(parts[2])==2 and rawHeader.find('IM1')!=-1):
#                nameMap[rawHeader]='MDup'
#            if (rawHeader.find('cvrg')==-1 and int(parts[2])==4) or (rawHeader.find('cvrg')!=-1 and int(parts[2])==2 and rawHeader.find('IP1')!=-1):
#                nameMap[rawHeader]='PDup'
    
    if len(rawHeader) > 0 and rawHeader not in nameMap.keys():
        print '>>UNCATEGORIZED: ', rawHeader,
        #nameMap[rawHeader]='PDel'

    summary_file.seek(0)
    for line in summary_file:
        #if a new result has started
        if line.find(".txt")!=-1:
            rawHeader=line
            header= line.split('-')
            if (header[4].startswith('delete')):
                simSZ=int(header[3])-int(header[2])
            else:
                simSZ=int(header[4])-int(header[3])    
                     
        if rawHeader not in nameMap:
            if rawHeader[-4:]!='!@#$':
                print '>>ERROR: ', rawHeader,
                rawHeader += '!@#$'
            continue
        
        if line.find('predicted')!=-1:
            parts=line.split(' ')
            reg= map(int, parts[0].split('-'))
            sz=reg[1]-reg[0]
            if int(parts[2])!=int(parts[5]) and ((rawHeader.find('cvrg')==-1 and (int(parts[5])==3)) or (rawHeader.find('cvrg')!=-1 and int(parts[5])==1)):
                falsePos[stateToType[int(parts[2])]][findIndex(sz)]+=1
                falsePos['sum'][findIndex(sz)]+=1
        
        if (line.startswith('Recall')):
            tmp_rc = int(line.split('\t')[1])
            if line.split('\t')[2].startswith("(NAN)"):
                tmp_rc = 0
            truePos[nameMap[rawHeader]][findIndex(simSZ)]+= tmp_rc
            truePos['sum'][findIndex(simSZ)]+= tmp_rc
        
        if (line.startswith('Recall')):
            total[nameMap[rawHeader]][findIndex(simSZ)]+=1
            total['sum'][findIndex(simSZ)]+=1
#        if (line.startswith('Recall')):
#            if int(line.split('\t')[1])==0 and simSZ==1000*k:
#                print rawHeader


    print '*',
    print "\t",
    print 'Size:',
    print "\t",
    print "\t",
    for i in range(1, len(binSizes)-1, 2):
            print "{0}k+{1}k\t".format(binSizes[i]/k, binSizes[i+1]/k),
    print 'total',

    #precision and recall:
    for key in sorted(total.keys()):
        print ''
        print key,
        print '\tRecall:',
        print "\t",
        for i in range(1, len(binSizes)-1, 2):
            if total[key][i]==0:
                print 1,
                print "\t",
            else: 
                print '%0.0f%%\t' % (float(truePos[key][i]+truePos[key][i+1])/(total[key][i]+total[key][i+1]) *100),
                print str(truePos[key][i]+truePos[key][i+1])+'/'+str(total[key][i]+total[key][i+1]),
                print "\t",
        #        print "Recall for ",
        #        print num,
        #        print ":: ",
        if sum(total[key])==0:
            print 1,
            print "\t",
        else:
            #print '%0.4f' % (float(sum(truePos[key]))/sum(total[key])),
            print str(sum(truePos[key]))+'/'+str(sum(total[key])),
            print "\t",

        print ''
        print key,
        print '\tPrecision:',
        print "\t",
        for i in range(1, len(binSizes)-1, 2):
            #print float(truePos[key][i]+truePos[key][i+1])/(truePos[key][i]+falsePos[key][i]+truePos[key][i+1]+falsePos[key][i+1]),
            if float(truePos[key][i]+falsePos[key][i]+truePos[key][i+1]+falsePos[key][i+1]) == 0:
                print "NA\t",
            else:
                print "%0.0f%%\t" % (float(truePos[key][i]+truePos[key][i+1])/(truePos[key][i]+falsePos[key][i]+truePos[key][i+1]+falsePos[key][i+1]) *100),
            
            print str(truePos[key][i]+truePos[key][i+1])+'/'+str(truePos[key][i]+falsePos[key][i]+truePos[key][i+1]+falsePos[key][i+1]),
            print "\t",
        print str(sum(truePos[key]))+'/'+str(sum(truePos[key])+sum(falsePos[key])),
        print "\t",
        print ''

    print "---------------------------"
    print '*',
    print "\t",
    print 'Size:',
    print "\t",
    print "\t",
    for i in range(len(binSizes)):
            print "{0}k\t".format(binSizes[i]/k),
    print 'total',
    #precision and recall:
    for key in sorted(total.keys()):
        print ''
        print key,
        print '\tRecall:',
        print "\t",
        for (i,num) in enumerate(binSizes):
            if total[key][i]==0:
                print 1,
                print "\t",
            else: 
                #print '%0.4f' % (float(truePos[key][i])/total[key][i]),
                print str(truePos[key][i])+'/'+str(total[key][i]),
                print "\t",
        #        print "Recall for ",
        #        print num,
        #        print ":: ",
        if sum(total[key])==0:
            print 1,
            print "\t",
        else:
            #print '%0.4f' % (float(sum(truePos[key]))/sum(total[key])),
            print str(sum(truePos[key]))+'/'+str(sum(total[key])),
            print "\t",

        print ''
        print key,
        print '\tPrecision:',
        print "\t",
        for (i,num) in enumerate(binSizes):
            #print float(truePos[i])/(truePos[i]+falsePos[i]),
            #if float(truePos[key][i]+falsePos[key][i]) == 0: print "1\t",
            #else: print "%.4f\t" % (truePos[key][i] / float(truePos[key][i]+falsePos[key][i])),
            print str(truePos[key][i])+'/'+str(truePos[key][i]+falsePos[key][i]),
            print "\t",
        print str(sum(truePos[key]))+'/'+str(sum(truePos[key])+sum(falsePos[key])),
        print "\t",
        print ''

    #fmeasure
#    for key in sorted(total.keys()):
#        print ''
#        print key,
#        print '\tF-score:',
#        print "\t",
#        for (i,num) in enumerate(binSizes):
#            if total[key][i] == 0: recall = 1
#            else: recall = truePos[key][i]/float(total[key][i])
#            if truePos[key][i]+falsePos[key][i] == 0: precis = 1
#            else: precis = truePos[key][i]/float(truePos[key][i]+falsePos[key][i])
#            
#            if recall + precis == 0: fscore = 0
#            else: fscore = 2 * (recall * precis) / float( recall + precis)
#            print "%.3f\t" % (fscore),

#        if sum(total[key])==0:
#            print 1,
#            print "\t",
#        else:
#            if sum(total[key]) == 0: recall = 1
#            else: recall = sum(truePos[key])/float(sum(total[key]))
#            if sum(truePos[key])+sum(falsePos[key]) == 0: precis = 1
#            else: precis =sum(truePos[key])/float(sum(truePos[key])+sum(falsePos[key]))
#            
#            if recall + precis == 0: fscore = 0
#            else: fscore = 2 * (recall * precis) / float( recall + precis)
#            print "%.3f\t" % (fscore),
#            
#    print ''   

 

if __name__ == '__main__':
    main()


#!/usr/bin/python2

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import argparse

normal=3

sameCopyCount={0:[0, 2], 2:[0, 2], 3:[3], 4:[4, 6], 6:[4, 6]}

def main():
    parser = argparse.ArgumentParser(description='Evaluates FCNV prediction recall and precision. Takes an annotation file produced by FCNV.')
    parser.add_argument('filenames', type=str, nargs=1, help='path to the annotation file')
    parser.add_argument('thereshold', type=int, nargs=1, help='threshold for min length of false positives to be considered')
    args = parser.parse_args()
    thresh=args.thereshold[0]

    results_file=open(args.filenames[0], "r")
    results=[]

    for line in results_file:
            line = line.strip().split(' ')
            results.append({"pos": int(line[0]),
                            "real": int(line[1]),
                            "pred": int(line[2]),})
    normal=int(results[0]['real'])

    foundSnips=0
    regionSnips=0
    for res in results:
        if res["real"]!=normal:
            regionSnips=regionSnips+1
            if res["real"] in sameCopyCount[res["pred"]]:
                foundSnips=foundSnips+1
    #print foundSnips
    #print regionSnips

    if regionSnips>0:
        recall=float(foundSnips)/regionSnips
    else:
        recall='NAN'

    
    

    foundSnips=0
    regionSnips=0
    for res in results:
       if res["pred"]!=normal:
           regionSnips=regionSnips+1
           if res["real"]==res["pred"]:
               foundSnips=foundSnips+1
               
    if regionSnips>0:
        precision=float(foundSnips)/regionSnips
    else:
        precision='NAN'
        
    #print foundSnips
    #print regionSnips
    

    #print "F-Measure: "+str(fmeasure)

    length=0
    lastCH='!'
    print_buffer=''
    for res in results:
        if res["pred"]!=res["real"]:
            if length==0:
                print_buffer += str(res["pred"])+": "
            length=length+1
        if length!=0 and res["pred"]==res["real"]:
            print_buffer += str(length)+'\n'
            length=0
    if length!=0:
        print_buffer += str(length)+'\n'
            
    length=0
    badCalls=0
    lastCH='!'
    count=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    for index, res in enumerate(results):
        if res["pred"]==lastCH:
            length=length+1
            count[res["real"]]=count[res["real"]]+1
        if res["pred"]!=lastCH:
            majority=0
            ctMax=count[0]
            i=0
            for ct in count:
                if ct>ctMax:
                    ctMax=ct
                    majority=i
                i=i+1

            if (length!=0):
                if length!=ctMax or majority!=lastCH:
                    print_buffer += str(results[index-length]["pos"])+'-'+str(results[index-1]["pos"])+' predicted: '+str(lastCH)+' ['+str(length)+"] real: "+str(majority)+" ["+str(ctMax)+"]\n"
            if length>thresh and lastCH!=normal and majority==normal:
                badCalls+=1
            
            count=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            length=1
            lastCH=res["pred"]
            count[res["real"]]=count[res["real"]]+1

    mRecall=0
    if (recall!='NAN' and recall>0.5):
        mRecall=1
    if (recall=='NAN'):
        mRecall='NAN'
    
    if recall == 'NAN':
        mRecall=1
    
    if badCalls == 0 and (mRecall==0 or mRecall=='NAN'):
        mPrecision=1
    else:
        mPrecision=float(mRecall)/(badCalls+mRecall)

    print "Recall:   \t" + str(mRecall) + "\t(" + str(recall) + ')'
    print "Precision:\t" + str(mPrecision) + "\t(" +str(precision) + ')'

    #print "recall: "+str(recall)
    #print "precision: "+str(precision)

    if len(print_buffer)!=0:
        print print_buffer,
    
        
    inTheRegion=False
    for res in results:
        if inTheRegion==False and res["real"]!=normal:
            count=0
            lastch=res['pred']
            inTheRegion=True
            print 'Real State:',res["real"]
            print 'Inside the region'
        if inTheRegion and res["real"]==normal:
            print '[('+str(lastch)+') - '+str(count)+']',
            break
        if inTheRegion:
            if res['pred']==lastch:
                count=count+1
            else:
                print '[('+str(lastch)+') - '+str(count)+']',
                count=1
                lastch=res['pred']
            
            
            


if __name__ == '__main__':
    main()

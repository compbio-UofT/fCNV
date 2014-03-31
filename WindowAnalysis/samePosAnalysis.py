#!/usr/bin/pypy
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import argparse

def main():
    parser = argparse.ArgumentParser(description='This script calculates the mean and variance for difference of coverage between maternal and plasma sequencing data. It needs path to maternal and plasma .piled.txt files as arguments')
    parser.add_argument('maternalFile', type=str, nargs=1, help='Path to .piled.txt file for maternal reads')
    parser.add_argument('plasmaFile', type=str, nargs=1, help='Path to .piled.txt file for plasma reads')
    args = parser.parse_args()

    maternal_file=open(args.maternalFile[0], "r")
    plasma_file=open(args.plasmaFile[0], "r")
    prefix_value_maternal=[0]*70000000
    prefix_value_plasma=[0]*70000000
    prefix_count_maternal=[0]*70000000
    prefix_count_plasma=[0]*70000000
    mToPlasma=67.256440892/29.5912297825

    

    last=0
    chert=0
    for position in maternal_file:
        #chert+=1
        #if chert==10000:
        #    break
        row=map(int,position.strip().split(' '))
        prefix_value_maternal[row[0]]=prefix_value_maternal[last]+row[1]
        prefix_count_maternal[row[0]]=prefix_count_maternal[last]+1
        last=int(row[0])
        #print prefix_value[last] ,

    last=0
    chert=0
    for position in plasma_file:
        #chert+=1
        #if chert==10000:
        #    break
        row=position.strip().split(' ')
        prefix_value_plasma[int(row[0])]=prefix_value_plasma[last]+int(row[1])
        prefix_count_plasma[int(row[0])]=prefix_count_plasma[last]+1
        last=int(row[0])
       #print prefix_count[last]
    
    for l in [10, 100, 1000, 10000, 100000, 1000000]:
        mean=0.
        ct=0
        for i in range(0,int(len(prefix_value_maternal))):
            if prefix_value_maternal[i]==0 or prefix_value_maternal[i+l]==0 or prefix_value_plasma[i]==0 or prefix_value_plasma[i+l]==0:
                continue
            maternal=(prefix_value_maternal[i+l]-prefix_value_maternal[i])/(prefix_count_maternal[i+l]-prefix_count_maternal[i])
            plasma=(prefix_value_plasma[i+l]-prefix_value_plasma[i])/(prefix_count_plasma[i+l]-prefix_count_plasma[i])
            maternal*=mToPlasma
            mean+=plasma-maternal
            ct+=1
        mean/=ct
       # if (prefix_value[i+l]-prefix_value[i])/(prefix_count[i+l]-prefix_count[i]) < 30.:
       #     print (prefix_value[i+l]-prefix_value[i])/(prefix_count[i+l]-prefix_count[i])
        #print float(prefix_value[i+l]-prefix_value[i])/l

        variance=0
        for i in range(0,int(len(prefix_value_maternal))):
            if prefix_value_maternal[i]==0 or prefix_value_maternal[i+l]==0 or prefix_value_plasma[i]==0 or prefix_value_plasma[i+l]==0:
                continue
            maternal=(prefix_value_maternal[i+l]-prefix_value_maternal[i])/(prefix_count_maternal[i+l]-prefix_count_maternal[i])
            plasma=(prefix_value_plasma[i+l]-prefix_value_plasma[i])/(prefix_count_plasma[i+l]-prefix_count_plasma[i])
            maternal*=mToPlasma
            variance+=(plasma-maternal-mean)**2

        variance/=ct
        print '--------'
        print 'length= '+str(l)
        print 'mean= '+str(mean)
        print 'variance= '+str(variance)

if __name__ == '__main__':
    main()


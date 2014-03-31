# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import argparse

def main():
    parser = argparse.ArgumentParser(description='This script calculates the mean and variance for coverage of different regions with the same size. It needs path to the .piled.txt file of the target reads as argument')
    parser.add_argument('PileUpFile', type=str, nargs=1, help='Path to .piled.txt file for the target reads')
    args = parser.parse_args()

    pile_file=open(args.filenames[0], "r")
    prefix_value=[0]*70000000
    prefix_count=[0]*70000000

    last=0
    chert=0
    for position in pile_file:
   #     chert+=1
   #     if chert==10000:
   #         break
        row=position.strip().split(' ')
        prefix_value[int(row[0])]=prefix_value[last]+int(row[1])
        prefix_count[int(row[0])]=prefix_count[last]+1
        last=int(row[0])
        #print prefix_value[last] ,
        #print prefix_count[last]


    for l in [10, 100, 1000, 10000, 100000, 1000000]:
        mean=0.
        ct=0
        for i in range(0,int(len(prefix_value))):
            if prefix_value[i]==0 or prefix_value[i+l]==0:
                continue
            mean+=(prefix_value[i+l]-prefix_value[i])/(prefix_count[i+l]-prefix_count[i])
            ct+=1
        mean/=ct
       # if (prefix_value[i+l]-prefix_value[i])/(prefix_count[i+l]-prefix_count[i]) < 30.:
       #     print (prefix_value[i+l]-prefix_value[i])/(prefix_count[i+l]-prefix_count[i])
        #print float(prefix_value[i+l]-prefix_value[i])/l


        variance=0
        for i in range(0,int(len(prefix_value))):
            if prefix_value[i]==0 or prefix_value[i+l]==0:
                continue
            doc= (prefix_value[i+l]-prefix_value[i])/(prefix_count[i+l]-prefix_count[i])
            variance+=(doc-mean)**2

        variance/=ct
        print '--------'
        print 'length= '+str(l)
        print 'mean= '+str(mean)
        print 'variance= '+str(variance)







if __name__ == '__main__':
    main()


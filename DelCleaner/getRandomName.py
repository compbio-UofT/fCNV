#!/usr/bin/pypy
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import argparse
from random import randint 

def main():
    parser = argparse.ArgumentParser(description='generating name for deletion file')
    parser.add_argument('size', type=int, nargs=1, help='cnv size')
    args = parser.parse_args()
    if randint(0,1)==1:
        name="A-chr1-"
    else:
        name="B-chr1-"
    if randint(0,1)==1:
        num=randint(2400000,110600000)
    else:
        num=randint(148000000,233700000)
    name+=str(num)+"-"+str(num+args.size[0])
    name+="-delete.target.txt"
    print name
        
 

if __name__=='__main__':
    main()

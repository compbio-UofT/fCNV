#!/bin/bash

/usr/bin/time -p pypy ~/Ubuntu\ One/FetalCNV/fcnv/fcnv.py $2 > $1 2>> $1

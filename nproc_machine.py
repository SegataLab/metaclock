#!/usr/bin/env python

from multiprocessing import Pool

def f(x):
    return {x:[x*x]}

if __name__ == '__main__':
    with Pool(processes=4) as pool:
    	print(pool.map(f,range(100)))

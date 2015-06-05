"""A solution to the Rosalind problem "Mendel's First Law."
http://rosalind.info/problems/iprb/
"""
import sys

def choose2(i):
    return i * (i - 1.) / 2.

def propability(k, m, n):
    return 1. - (choose2(n) + choose2(m) / 4. + m * n / 2.) / choose2(k + m + n)     
    
with open(sys.argv[1]) as f:
        k, m, n = map(int, f.read().split())    
	      print propability(k, m, n)

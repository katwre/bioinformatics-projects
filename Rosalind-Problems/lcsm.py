#!/usr/bin/env python
# coding=utf-8

#python 2.x
from Bio import SeqIO #biopython

"""
Finding a Shared Motif
======================
***Problem

A common substring of a collection of strings is a substring of every member of the collection.
We say that a common substring is a longest common substring if there does not exist a longer
common substring. For example, "CG" is a common substring of "ACGTACGT" and "AACCGGTATA", but 
it is not as long as possible; in this case, "GTA" is a longest common substring of "ACGTACGT" 
and "AACCGTATA".Note that the longest common substring is not necessarily unique; for a simple
example, "AA" and "CC" are both longest common substrings of "AACC" and "CCAA".

Given: A collection of k (k≤100) DNA strings of length at most 1 kbp each in FASTA format.
Return: A longest common substring of the collection. (If multiple solutions exist, you may 
return any single solution.)

***Sample Dataset

>Rosalind_1
GATTACA
>Rosalind_2
TAGACCA
>Rosalind_3
ATACA
Sample Output

AC
"""

def search_substrings(short_string):
  
  o = []
  k = len(short_string)

  while k>0:
    x = set()
    for i in xrange(0,len(short_string)-k+1):
      x.add(short_string[i:i+k])
    o += x
    k-=1

  return o
      
      
def if_substr_in_all(substr,l):
  
    for i in l:
        if substr not in i:  return False
    return True



def lcsm(records):
    """the longest common consecutive substring"""

    records = sorted(records, key=lambda i: len(i))
    short_string = records[0]
    other_strings = records[1:]
    substrings = search_substrings(short_string)

    for substr in substrings:
        if if_substr_in_all(substr,other_strings):
            return substr
    
    return ""
    

if __name__ == "__main__":
     
      handle = open("./datasets/rosalind_lcsm.txt", "rU")
      records = [str(i.seq) for i in SeqIO.parse(handle, "fasta")]
      handle.close()

      print(lcsm(records))

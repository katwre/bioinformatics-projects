#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Problem - https://rosalind.info/problems/revc/
Get reverse complement of a DNA string.
"""

def reverse_complement(l):
    complements = {"A":"T","T":"A","G":"C","C":"G"}
    return ''.join([complements[i] for i in l][::-1])


if __name__ == "__main__":

    # read an input file
    with open("./datasets/rosalind_revc.txt") as f:
        dna_seq = f.read().splitlines()[0]
    
    # print reverse complement sequence
    revc_dna_seq = reverse_complement(dna_seq)
    print(revc_dna_seq)





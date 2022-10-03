#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Problem
https://rosalind.info/problems/gc/

Given: At most 10 DNA strings in FASTA format (of length at most 1 kbp each).

Return: The ID of the string having the highest GC-content, followed by the GC-content of that string. Rosalind allows for a default error of 0.001 in all decimal answers unless otherwise stated; please see the note on absolute error below.
"""
from Bio import SeqIO

def gc_content(seq):
    g = seq.count("G")
    c = seq.count("C")
    seq_len = len(seq)
    return( (g+c) / seq_len * 100 )


if __name__ == "__main__":
    
    gc_dict = {}
    for record in SeqIO.parse("./datasets/rosalind_gc.txt", "fasta"):
        gc_dict[record.id] = gc_content(record.seq)
    id_max = max(gc_dict, key = lambda k: gc_dict[k])
    
    print( id_max )
    print( gc_dict[id_max] )
    

        

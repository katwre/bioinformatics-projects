#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Problem - https://rosalind.info/problems/rna/

An RNA string is a string formed from the alphabet containing 'A', 'C', 'G', and 'U'.
Given a DNA string tt corresponding to a coding strand, its transcribed RNA string uu is formed by replacing all occurrences of 'T' in tt with 'U' in uu.
Given: A DNA string tt having length at most 1000 nt.
Return: The transcribed RNA string of tt.

Sample Dataset
GATGGAACTTGACTACGTAAATT
Sample Output
GAUGGAACUUGACUACGUAAAUU
"""

# read an input file
with open("./datasets/rosalind_rna.txt") as f:
    dna_seq = f.read().splitlines()

rna_seq = dna_seq[0].replace("T","U") 
print(rna_seq)


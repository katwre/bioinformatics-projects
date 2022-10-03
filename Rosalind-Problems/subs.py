#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Problem - https://rosalind.info/problems/subs/

Given: Two DNA strings s and t (each of length at most 1 kbp).
Return: All locations of t as a substring of s.
"""


def find_motif(s,t):

    positions = []
    for i in range(len(s)):
        if s[i]==t[0]:
            if s[ i:(i+len(t)) ] == t:
                positions.append(i+1)
    return(positions)



if __name__ == "__main__":
    
    with open('./datasets/rosalind_subs.txt','r') as file:
        content = file.read()

    DNA, subDNA = content.splitlines()
    positions = find_motif(DNA, subDNA)
    print( " ".join([str(x) for x in positions])  )



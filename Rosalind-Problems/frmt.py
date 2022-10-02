#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Problem

GenBank can be accessed here. A detailed description of the GenBank format can be found here. A tool, from the SMS 2 package, for converting GenBank to FASTA can be found here.

Given: A collection of n (n≤10) GenBank entry IDs.
Return: The shortest of the strings associated with the IDs in FASTA format.


Sample Dataset

FJ817486 JX069768 JX469983


Sample Output

>JX469983.1 Zea mays subsp. mays clone UT3343 G2-like transcription factor mRNA, partial cds
ATGATGTATCATGCGAAGAATTTTTCTGTGCCCTTTGCTCCGCAGAGGGCACAGGATAATGAGCATGCAA
GTAATATTGGAGGTATTGGTGGACCCAACATAAGCAACCCTGCTAATCCTGTAGGAAGTGGGAAACAACG
GCTACGGTGGACATCGGATCTTCATAATCGCTTTGTGGATGCCATCGCCCAGCTTGGTGGACCAGACAGA
GCTACACCTAAAGGGGTTCTCACTGTGATGGGTGTACCAGGGATCACAATTTATCATGTGAAGAGCCATC
TGCAGAAGTATCGCCTTGCAAAGTATATACCCGACTCTCCTGCTGAAGGTTCCAAGGACGAAAAGAAAGA
TTCGAGTGATTCCCTCTCGAACACGGATTCGGCACCAGGATTGCAAATCAATGAGGCACTAAAGATGCAA
ATGGAGGTTCAGAAGCGACTACATGAGCAACTCGAGGTTCAAAGACAACTGCAACTAAGAATTGAAGCAC
AAGGAAGATACTTGCAGATGATCATTGAGGAGCAACAAAAGCTTGGTGGATCAATTAAGGCTTCTGAGGA
TCAGAAGCTTTCTGATTCACCTCCAAGCTTAGATGACTACCCAGAGAGCATGCAACCTTCTCCCAAGAAA
CCAAGGATAGACGCATTATCACCAGATTCAGAGCGCGATACAACACAACCTGAATTCGAATCCCATTTGA
TCGGTCCGTGGGATCACGGCATTGCATTCCCAGTGGAGGAGTTCAAAGCAGGCCCTGCTATGAGCAAGTC
A

"""

from Bio import Entrez, SeqIO

def shortest_entry(entry_ids):
    Entrez.email = "katwre@gmail.com"
    print([", ".join(entry_ids)])
    handle = Entrez.efetch(db="nucleotide", id=[", ".join(entry_ids)], rettype="fasta")
    records = list(SeqIO.parse(handle, "fasta"))
    print(min(records, key=lambda s: len(s.seq)).format("fasta"))

if __name__ == "__main__":
    with open("./datasets/rosalind_frmt.txt", "r") as f:
        entry_ids = f.readline().strip().split()
    shortest_entry(entry_ids)





#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Problem

The novel-motif finding tool MEME can be found here (http://meme-suite.org/tools/meme).

Given: A set of protein strings in FASTA format that share some motif with minimum length 20.

Return: Regular expression for the best-scoring motif.
"""

from Bio import SeqIO
import subprocess
from Bio.motifs import meme
import numpy as np


def run_motif(fasta_file, outputdir):
    cmd = ["/root/meme/bin/meme", fasta_file, "-minw", "20", "-oc", outputdir]
    print(cmd)
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    out, err = process.communicate()
    return out


def extract_regex(motif):
    charar = np.array( [motif[0].instances])
    output_list = []
    for i in range(charar.shape[2]):
        nuc = list(set(charar[0][:,i]))
        if len(nuc) == 1:
            output_list.append(nuc[0])
        else:
            output_list.append("[" + "".join(nuc) + "]")
    return( "".join(output_list) )



if __name__ == "__main__":
    
    input_file = "./datasets/rosalind_meme.txt"

    # run meme
    meme_outdir = "/root/Rosalind-problems/bioinformatics-armory/meme_out/"
    meme_out = run_motif(input_file, meme_outdir)
    meme_out = meme_outdir + "meme.xml"
    record = meme.read(meme_out)
    
    # extract regex from the meme output
    print(extract_regex(record))





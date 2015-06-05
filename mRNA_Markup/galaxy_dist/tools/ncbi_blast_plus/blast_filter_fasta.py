#!/usr/bin/env python

import sys
from Bio import SeqIO

def blast_filter_fasta(blast_file,in_file, positive_w, negative_w):
  """Filter a FASTA file using tabular output, e.g. from BLAST."""

  #Read tabular BLAST file and record all queries with hit(s)
  ids = set()
  blast_handle = open(blast_file, "rU")  
  for line in blast_handle:
      ids.add(line.split("\t")[0])
  blast_handle.close()

  #Write filtered FASTA file based on IDs from BLAST file
  reader = open(in_file, "rU")
  positive_writer = open(positive_w, "w")
  negative_writer = open(negative_w, "w")

  for record in SeqIO.parse(reader, "fasta"):
    if record.name in ids:
      SeqIO.write(record, positive_writer, "fasta")
    else:
      SeqIO.write(record, negative_writer, "fasta")
  
  reader.close()    
  positive_writer.close()
  negative_writer.close()







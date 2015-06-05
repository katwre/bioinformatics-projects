#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
from Bio import SeqIO



inputfile = sys.argv[1]
input_data_path = sys.argv[2]
outputfile1 = sys.argv[3]
outputfile2 = sys.argv[4]



def seqnames_from_museqboxfile(filename):
  """parse sequences id's from MuSeqBox output"""
  
  plik_msb = open(filename,"r")

  j=0
  nazwy_seq = []
  for i in plik_msb.readlines():

    if i[0]=="-" and j==0:
	j+=1
	continue
    if i[0]=="-" and j>0: 
      break
    
    if j>0 and i[0] != "\n" and i[0] != " ":
      nazwy_seq.append( i.split()[0] )

  return nazwy_seq

  
def MuSeqBox_Partition(museqboxfile_path, inputfile_path):

  museqboxfile = open(museqboxfile_path,"r")
  
  l = seqnames_from_museqboxfile(museqboxfile_path)
  
  sth_mRNA = []
  not_sth_mRNA = []
  
  for record in SeqIO.parse(open(inputfile_path, "rU"), "fasta"):

    if record.name.split("|")[-1] in l:
      sth_mRNA.append(record)
    else:
      not_sth_mRNA.append(record)
  
  return (sth_mRNA, not_sth_mRNA) 
  
  
  
(VC_mRNA, not_VC_mRNA) =  MuSeqBox_Partition(inputfile, input_data_path)

output_handle1 = open(outputfile1, "w")
SeqIO.write(VC_mRNA, output_handle1, "fasta")
output_handle1.close()

output_handle2 = open(outputfile2, "w")
SeqIO.write(not_VC_mRNA, output_handle2, "fasta")
output_handle2.close()  
  
  

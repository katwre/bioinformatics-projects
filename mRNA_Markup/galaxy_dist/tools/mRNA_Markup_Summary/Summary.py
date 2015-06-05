#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
from Bio import SeqIO


outputfile = sys.argv[1]
input_data,VC,BC,RA,FL,not_FL,PC,not_PC,AA,CD = sys.argv[2:]


def l(arg):
  """count sequences in FASTA file"""
  return len(list(SeqIO.parse(open(arg, "rU"), "fasta")))
  

summary_text = """

mRNAmarkup Report

  Number of input sequences:                            %d
  Number of potential vector-contaminated sequences:    %d (file: VC_mRNA.fas)
  Number of potential bacterial-contaminated sequences: %d (file: BC_mRNA.fas)
  Number of sequences matching the ReferenceDB:         %d (file: RA_mRNA.fas)
    Number of potential full-length coding sequences:     %d (file: FL_mRNA.fas)
    Non-qualifying sequences:                             %d (file: not_FL_mRNA.fas)
      Number of potential chimeric sequences:               %d (file: PC_mRNA.fas)
      Non-qualifying sequences:                             %d (file: not_PC_mRNA.fas)
  Number of sequences matching the AllProteinDB:        %d (file: AA_mRNA.fas)
  Number of sequences matching the ProteinDomainDB:     %d (file: CD_mRNA.fas)
  Number of remaining sequences:                        %d (file: remaining-mRNA)
  
""" % (l(input_data),l(VC),l(BC),l(RA),l(FL),l(not_FL),l(PC),l(not_PC),l(AA),l(CD),l(input_data)-l(VC)-l(BC)-l(RA)-l(PC)-l(AA)-l(CD))


summary = open(outputfile,"w")
summary.write(summary_text)
summary.close() 

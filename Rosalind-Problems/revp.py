"""A solution to the Rosalind problem "Locating Restriction Sites."
http://rosalind.info/problems/revp/

restriction enzyme cuts between i and i+1
"""
import sys
from Bio import SeqIO #biopython

        
def _reverse_complement(l):
    complements = {"A":"T","T":"A","G":"C","C":"G"}
    return ''.join([complements[i] for i in l][::-1])

    
def _palindroms(seq):
    
    results=[]

    for i in xrange(1,len(seq)):
      
      for j in xrange(1, 6):
	if i>=j and i+j+2<len(seq):
	    
	    if seq[i-j:i+1]==_reverse_complement(seq[i+1:i+j+2]):
	      results.append((i-j+1,len(seq[i-j:i+j+2])))
	      
    return results
     


  
if __name__ == "__main__":
     
      seq = list(SeqIO.parse(open("./datasets/rosalind_revp.txt", "rU"), "fasta"))[0].seq

      palindroms = sorted(_palindroms(str(seq)), key=lambda i: i[0])
      for i in palindroms:
	      print i[0],i[1]

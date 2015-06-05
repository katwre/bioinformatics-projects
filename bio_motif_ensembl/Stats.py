
from Bio import Motif
from scipy import stats 


class Stats(object):

    def __init__(self,MSAs):
      """MSAs : a list of MSA's, in particular list of class AlingIO objects"""

      self.MSAs = MSAs


    def number_of_occurrences(self,motif,fpr=0.001,precision=10**4):
      """returns number of occurrences of a given motif. Motif is represented as class Motif object"""
      
      sd = Motif.ScoreDistribution(motif,precision=precision)
      threshold = sd.threshold_fpr(fpr)
      n=0
      for msa in self.MSAs:
        for s in msa:
          for pos,score in motif.search_pwm(s.seq,threshold=threshold):
            if pos >= 0:
              n +=1
      return n

    def pomocnicza(self,seq,motif):
      i=0
      wynik = 0
      while i < len(seq):
  while i< len(seq) and seq[i]=="-":
	  i+=1
	j=i
	while j< len(seq) and seq[j]!="-":
	  j=j+1
	
	
	wynik+=max(0,j-i-len(motif)+1)
	i=j
      
      return wynik
      
      
      
    def _Binomial_distribution(self,motif,number_of_occurrences,p=0.001,nmbr_of_gaps=0):
      """returns probability of p or more success in Bernoulli process"""

      n = sum([i.get_alignment_length() -len(motif)+1 for i in self.MSAs]) 
      k = number_of_occurrences                                            
      return stats.binom.cdf(n-k,n,1-p)
      
      
    def overrepresented_motif(self,motif,fpr=0.001,alfa=0.05,nmbr_of_gaps=0):
      """motif is represented as class Motif object"""
         
      k = self.number_of_occurrences(motif=motif,fpr=fpr)
      p=self._Binomial_distribution(motif=motif,number_of_occurrences=k,p=fpr)
      if p < alfa:
            return p,motif
            
            
            
            
            
            
            
            
            

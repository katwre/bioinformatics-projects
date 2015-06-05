#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
The Felsenstein's algorithm implementation
"""
import math

NUCLEOTIDES = "ACTG"


def children(w, alpha):
  """Returns children of node w"""
  assert(w in alpha)
    
  for i in xrange(len(alpha)):
    if alpha[i]==w:
      return i, i+1


def distance(w1, w2, alpha, t):
  """Returns distance between node w1 and w2"""
  for i in xrange(len(alpha)):
    if w1 in [i,alpha[i]] and w2 in [i,alpha[i]]:
      return t[i]


def jukes_cantor_model(alfa, n1, n2, t):
  """JC69 model of DNA evolution"""
  
  if n1==n2:
    return 1/4. + 3/4. * math.exp(-4.*alfa*t)
  else:
    return 1/4. - 1/4. * math.exp(-4.*alfa*t)


def Pr(x,y,T):
  """Return propability of substitution nucelotide x to y. T=length of edge."""
  return jukes_cantor_model(1, x, y, T)



def PrL(k,a,u,x,t,alpha):
  """Returns propability of subtree observed from node k (k is the root) at site u in the sequence,
  where in k is nucleotide a.
  ( P(L_k | a) )"""
  
  if k not in alpha: # if k in a leaf
    if x[k][u]==a: 	return 1
    else:		return 0
  else: # if k is the interior node
    
    s=0.
    b, c = children(k, alpha)
    for i in NUCLEOTIDES:
      for j in NUCLEOTIDES:
	s = s + PrL(b,i,u,x,t,alpha) * PrL(c,j,u,x,t,alpha) * Pr(i, a, distance(k,b,alpha,t)) * Pr(j, a, distance(k,c,alpha,t))
	
  return s


def PrLBis(u, x, t, alpha):
  """Returns propability of tree in u position"""
  s = 0.
  q = 1/4.
  root = 2*len(x) - 2
  
  for a in NUCLEOTIDES:
    s  = s + q * PrL(root, a, u, x, t, alpha) # o is the root
  return s


def TreeLikelihood(x, t, alpha):
  """Returns tree likelihood, where:
     x = table of sequences (alignment),
     t = lenght of edges,
     alfa = relationship between son and father"""
      
  likelihood = 1.
  for u in xrange(len(x[0])):
    likelihood = likelihood * PrLBis(u, x, t, alpha)
  
  return likelihood


def main():
  
      x = ["ACC","CCT","AAC","GGG"] #sequences of leaves
      alpha = [4, 4, 5, 5, 6, 6, -1.] #tree with root; ( (0, 1)4, (2, 3)5 )6
      #alpha[1..2n-1], where alpha[i] is father of node i; 1-n leaves, n+1,...,2m-1 interior nodes
      t = [1., 2.5, 1., 1., 3., 4.5, -1.] #length of edges

      print TreeLikelihood(x, t, alpha)



if __name__ == "__main__":
      main()




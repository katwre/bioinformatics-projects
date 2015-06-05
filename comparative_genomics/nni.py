#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
The NNI algorithm implementation
"""

class Node:

  def __init__(self,n):
    self.n = n
    self.father = None

    if type(n)==str:
      self.if_leaf = True
    else:
      self.if_leaf = False
      
      self.l = Node(n[0])
      self.l.father = self
      
      self.r = Node(n[1])
      self.r.father = self


  def sibling(self):
    if self.father != None:
      f=self.father
      if self==f.r:
        return f.l
      else:
        return f.r

  def nodes(self):
    	if self.if_leaf:
      		return []
   	else:
		if self.father:
      			return [self] + self.l.nodes() + self.r.nodes()
		else:
			return self.l.nodes() + self.r.nodes()
   
  def NNI(self):
	a=self.sibling()
	f=self.father
	b=self.l
	c=self.r
	if f.r==self:
		self.l=a
		a.father=self
		self.r=b
		b.father=self
		self.father.l=c
		c.father=self.father
	else:
		self.l=a
		a.father=self
		self.r=b
		b.father=self
		self.father.r=c
		c.father=self.father

  def __str__(self):
    if self.if_leaf:
      return self.n
    else:
      return "(" + self.l.__str__() + "," + self.r.__str__() + ")"
    
  def __repr__(self):
    return self.__str__()
  
  
  
  
if __name__ == "__main__":

  g = Node( ('a',('c',('d','e'))) )
  for x in g.nodes():
	  print g
	  x.NNI()
	  print g
	  x.NNI()
	  print g
	  x.NNI()
	  print g
	  print "end of changes"

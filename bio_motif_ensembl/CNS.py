

class CNS(object):
    """Class representing conserved noncoding sequence"""
    
    
    def __init__(self,MSA):
        self.MSA = MSA
        self.length_MSA = MSA.get_alignment_length()
        
      
    def _msa_to_binary_table():
      
      table=[]
      
      length_of_column = len(self.MSA[:])
      A,T,C,G = "A"*length_of_column,"T"*length_of_column,"C"*length_of_column,"G"*length_of_column

      for i in xrange(self.length_MSA):
  if "-" not in self.MSA[:,i] and ( self.MSA[:,i]==A or self.MSA[:,i]==T or self.MSA[:,i]==C or self.MSA[:,i]==G):
	  table.append(1)
	else:
	  table.append(0)
	  
      return table
    
    
    def _fragment_identity(self,i,j,cumsum):
      return float(cumsum[j] - cumsum[i]) / float(j-i)
      
    def _cumulative_sum(self,vector):
      total = 0
      yield total
      for i in vector:
	total += i
	yield total
      
      
    def _modif_naive_method(self,binary_table,block_length,ident_threshold):
	
	stack=[(0,len(binary_table),len(binary_table))] #(start,end,len(block))

	blocks_ok=[]
	cumsum = list(self._cumulative_sum(binary_table))
	
	
	while stack:
	  
	  i,j,ln = stack.pop()
	  
	  while ln >= block_length: #looping length of blocks
	    
	    found=False
	    for k in xrange(i,j-ln+1): #looping position of block
	      
	      ident = self._fragment_identity(k,k+ln,cumsum)
	      if ident >= ident_threshold:
		
		blocks_ok.append((k,k+ln,ident))
		stack.append((i,k,ln-1))
		stack.append((k+ln,j,ln))
		found=True
		break
		
	    if found:
		#when block is ok
		break
	    else:
		#if block with length ln is not found, search for block with length ln-1
		ln -= 1
 
	return blocks_ok
    
      
      
    def get_conserved_blocks_of_MSA(self,block_length,ident):
      
      assert(type(block_length)==int and block_length > 0)
      assert(type(ident)==float)
      
      binary_table = self._msa_to_binary_table()
      blocks =  self._modif_naive_method(binary_table,block_length,ident)

      MSAs = []
      for block in blocks:
	MSAs.append( (self.MSA[:,block[0]:block[1]],block[0],block[1],block[2]) ) # (seq of MSA, start position, end position)
      
      return MSAs
            
  
  
      
  

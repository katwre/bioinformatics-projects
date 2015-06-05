"""
@description: De Bruijn graph implementation in genome assembly
"""

class Node:
  
  def __init__(self, s):
    self.n_in = 0
    self.n_out = 0
    self.s = s
  
  def __str__(self):
    return self.s
  
  def __repr__(self):
    return self.__str__()
  
  def __eq__(self, node):
    return self.s == node.s
      
  def __hash__(self):
    return hash(self.s) #Two objects with the same value have the same hash value.
  
  def isBalanced(self):
    return self.n_in == self.n_out
  
  def isSemiBalanced(self):
    return self.n_in + 1 == self.n_out or self.n_in - 1 == self.n_out




class deBruijnGraph: 
   
    
  def __init__(self, s, k):
    """de Bruijn multigraph with given string iterator and k-mer length k 
    """
    self.k = k
    self.string = s

    self.edges, self.nodes = [], {} #[(node_from, node_to)], {name:Node}
    self.build_graph() 

  def generate_k_mers(self):
    for i in xrange(len(self.string)-self.k+1):
      yield self.string[i:i+self.k]  
    
   
  def add_node(self, km): 
    if km not in self.nodes.keys():
	self.nodes[km] = Node(km)
   

  def build_graph(self):
    """Build MultiGraph"""
        
    for i in self.generate_k_mers():
      
      km1, km2 = i[:-1], i[1:]
      if not km1==km2: #without edges that lead from and to the same node, e.g. node1 -> node1
	  self.add_node(km1)
	  self.add_node(km2)
	  n1, n2 = self.nodes[km1], self.nodes[km2]
	  
	  self.edges.append( (n1, n2) )
	  n1.n_out += 1
	  n2.n_in += 1  
	  
	  self.nodes[km1], self.nodes[km2] = n1, n2
  
  def isEulerian(self):
      """Eulerian graph is a directed, connected graph that has at most 2 semi-balanced
	 nodes and all other nodes are balanced"""
      semi_balanced, balanced = 0, 0
      for k,v in self.nodes.iteritems():
	if v.isBalanced():	balanced +=1
	if v.isSemiBalanced():	semi_balanced +=1
      return semi_balanced <=2 and balanced==len(self.nodes) - semi_balanced
    
  def eulerian_walk(self):
      """For Eulerian graph, Eulerian walk can be found in O(|E|) time. |E| is # edges."""
      
      if self.isEulerian:
	
	tour = []
	graph = list(self.edges) #copy
	
	current_vertex = graph[0][0]
	tour.append(current_vertex)

	while len(graph) > 0:
	    #print(graph, current_vertex)
	    for edge in graph:
		if current_vertex in edge:
		    if edge[0] == current_vertex:
			current_vertex = edge[1]
		    else:
			current_vertex = edge[0]

		    graph.remove(edge)
		    tour.append(current_vertex)
		    break
		else:
		    # Edit to account for case no tour is possible
		    return False
	return tour
  
  
  def to_dot(self):
      """Return string with graphviz representation
      """
      ss = "digraph deBruijnGraph {\n" #directed graph
      for e in self.edges:
	ss += '"%s" -> "%s";\n' % (e[0], e[1])
      return ss + "}"
   
  def to_png(self, png_path = "deBruijnGraph.png"):
      """Return file that contains graph visualization in the PNG format
      """
      import subprocess, os
      dot_path = "./deBruijnGraph.dot"
      f = open(dot_path, "w")
      f.write(self.to_dot())
      f.close()
      subprocess.call("dot -Tpng %s -o %s" % (dot_path, png_path), shell=True)
      os.remove(dot_path)
   
   
G = deBruijnGraph("AAABB", 3)  
#G = deBruijnGraph("a_long_long_long_time", 5) 
print G.eulerian_walk()
#print G.edges
#print G.nodes
#for i in G.nodes.values():
  #print i, i.n_in, i.n_out, i.isBalanced()
G.to_png()

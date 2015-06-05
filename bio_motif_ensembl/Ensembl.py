import MySQLdb,re

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC



ENSEMBL_HOST = "ensembldb.ensembl.org"
ENSEMBL_USER = "anonymous"
ENSEMBL_PASSWD = ""




class EnsemblException(Exception):
  def __init__(self, value):
      self.value = value
  def __str__(self):
      return repr(self.value) 
        
class ConnectionError(EnsemblException):
  pass    
class RelationshipError(EnsemblException):
  pass     
class SpieciesError(EnsemblException):
  pass
class ContigError(EnsemblException):
  pass  

  
  
  
class Ensembl_Connection(object):
  
  def __init__(self,host,user,passwd):
    
    self.host = host
    self.user = user
    self.passwd = passwd
    
    self.con = None
  
  def connect_with_ensembl(self,release=None,port=None):
    """type release OR port"""
  
    if release!=None and port==None:
      assert (type(release)==int)
      if release < 48:
  port1  = 3306
      else:
	port1 = 5306
    elif release==None and port!=None:
      assert(type(port)==int)
      port1 = port
    
    #connecting with db
    try:
	self.con = MySQLdb.connect(host=self.host, user=self.user, passwd=self.passwd, port = port1)
    except Exception as e:
	raise ConnectionError(e)

    return self.con.cursor(), self.con
  
  
  def close(self):
    self.con.close() 

  
  

class Genome(object):
  """Class representing information contained in the Core database."""
  
  def __init__(self, release, genome_name):
    
    
    assert(type(release)==int and release > 0)
    assert(genome_name!=None)
    
    self.genome_name = genome_name

    connection = Ensembl_Connection(ENSEMBL_HOST,ENSEMBL_USER,ENSEMBL_PASSWD)
    cursor,con = connection.connect_with_ensembl(release=release)
    self.cursor,self.con = cursor,con
    
    dgenome_name = genome_name.replace(" ","_").lower()
    query = "show databases like '"+dgenome_name+"_core_"+str(release)+"%';"
    
    
    cursor.execute(query)
    assert cursor.rowcount==1
    database = cursor.fetchone()[0]

    if len(database)==0:
       raise SpieciesError("Given genome's name doesn't exist in the Ensembl database.")
      
    con.select_db(database)

    
    
    
  def close(self):
      """method closes connection with the Ensembl database."""
      
      self.con.close()  
   
  def __repr__(self):
      return "Genome(Name='%s')" % (self.genome_name)
      
      
  def _to_SeqRecord(self,seq, name ="<unknown name>", status="<unknown status>", biotype="<unknown type>",gene_id='<unknown id>',description='<unknown description>',genome_name="<unknown genome's name>",start = None,end=None,strand = None,source='<unknown source>',chromosome=None,ref="",ref_db=""):
      """
      id : stable_id of gene from the Ensembl database.
      ref: primary accession number.
      ref_db: external database's name.
      """
      assert(type(seq)==Seq)
      assert(type(start)==int)
      assert(type(end)==int)
      assert(type(chromosome)==int or type(chromosome)==str)     
     
      length = end-start #sequence starts at position 1
      a = {"genome": genome_name,
	    "start": start,
	    "end":end,
	    "length":length,
	    "chromosome":chromosome,
	    "strand": strand,
	    "source":source,
	    "biotype":biotype,
	    "status":status,
	    "ref":ref,
	    "ref_db":ref_db,
	    
	    }    
      f = SeqFeature(FeatureLocation(start, end), type=biotype, strand=strand,ref=ref,ref_db=ref_db)
      record = SeqRecord(seq,
                   id=gene_id,
                   name=name,
                   description=description,
                   annotations=a,
                   features = [f],
                   )
      return record
      
      
      
  def _get_seq_from_contigs(self,coordinates_and_seqs,strand):
      """coordinates_and_seqs : a list of tuples like:
	 (contig_start, contig_end, region_start, region_end, strand, sequence, chromosome,position_start,position_end).
         strand : sequence region strand: 1 - forward; -1 - reverse."""
      
      for c in coordinates_and_seqs:
	c[5] = c[5][(int(c[0])-1):int(c[1])]
      
      
      assert(len(set([(c[2],c[3]) for c in coordinates_and_seqs]))==1)
      assert(all([len(c[5])==int(c[1])-int(c[0])+1 for c in coordinates_and_seqs])) # all -> czy wszystkie elementy listy to true
      
      
      coordinates_and_seqs = sorted(coordinates_and_seqs, key=lambda contig: contig[7]) #order by contigs start position
      
      #Check if contigs overlap
      if len(coordinates_and_seqs) > 1:
	for i in xrange(len(coordinates_and_seqs)-1):
	  assert((int(coordinates_and_seqs[i+1][7]) - int(coordinates_and_seqs[i][8])) >= 1)
	  

      if len(coordinates_and_seqs) > 1:
	for i in xrange(len(coordinates_and_seqs)-1):
	  #if there is a gap between contigs
	  if (int(coordinates_and_seqs[i+1][7]) - int(coordinates_and_seqs[i][8])) > 1:
	    n = int(coordinates_and_seqs[i+1][7]) - int(coordinates_and_seqs[i][8]) - 1
	    coordinates_and_seqs[i][8] = coordinates_and_seqs[i][8] + n
	    coordinates_and_seqs[i][5] += n * "N"  

      for contig in coordinates_and_seqs:
	
	#all contigs to the same strand
	if contig[4]!=str(strand):
	  contig[5]=Seq(contig[5],IUPAC.IUPACAmbiguousDNA()).complement()
	  contig[4]=str(strand)
	else:
	  contig[5]=Seq(contig[5],IUPAC.IUPACAmbiguousDNA())
      
	
	#trimming contigs to the sequence region
	if contig[7] < contig[2]:
	 cut = int(contig[2]) - int(contig[7])
	 contig[5] = contig[5][cut:]
	 contig[7] = contig[2]
	if contig[8] > contig[3]:
	 cut = int(contig[8]) - int(contig[3])
	 contig[5] = contig[5][:-cut]
	 contig[8] = contig[3] 


      if int(coordinates_and_seqs[0][7]) - int(coordinates_and_seqs[0][2]) > 0:
	
	n = int(coordinates_and_seqs[0][7]) - int(coordinates_and_seqs[0][2])
	coordinates_and_seqs[0][7] = coordinates_and_seqs[0][7] - n
	coordinates_and_seqs[0][5] = n * "N" + coordinates_and_seqs[0][5]

      if int(coordinates_and_seqs[-1][3]) - int(coordinates_and_seqs[-1][8]) > 0:

	n = int(coordinates_and_seqs[-1][3]) - int(coordinates_and_seqs[-1][8])
	coordinates_and_seqs[-1][8] = coordinates_and_seqs[-1][8] + n
	coordinates_and_seqs[-1][5] = coordinates_and_seqs[-1][5] + n * "N"
	
	
      #I assume that contigs covered our region are next to each other and they don't overlap
      seq=Seq("",IUPAC.IUPACAmbiguousDNA())
      for c in coordinates_and_seqs:
	seq += c[5]

      return seq	
      
      
      
      
  def _get_seq_by_stable_id(self,stable_id, strand):
      """stable_id : gene's id in the Ensembl database.
         strand : sequence region strand: 1 - forward; -1 - reverse."""
      
      assert(strand==1 or strand==-1)
      assert(type(stable_id)==str)

      query = " select s.cmp_start, s.cmp_end,g.seq_region_start, g.seq_region_end, s.ori, d.sequence, r.name,s.asm_start, s.asm_end  \
	       from assembly s, seq_region r, coord_system o,seq_region r2, gene g, dna d   \
	       where r.seq_region_id=g.seq_region_id and r.seq_region_id=s.asm_seq_region_id and  \
	       s.cmp_seq_region_id=r2.seq_region_id and r2.coord_system_id=o.coord_system_id and   \
	       o.name='contig' and g.stable_id='%s' and   \
	       ((g.seq_region_end >= s.asm_start and g.seq_region_end <= s.asm_end) or (g.seq_region_start >= s.asm_start and g.seq_region_start <= s.asm_end) or (g.seq_region_start<= s.asm_start and g.seq_region_end>= s.asm_end))  \
	       and d.seq_region_id=s.cmp_seq_region_id; " % (str(stable_id))

      self.cursor.execute(query)
      results = self.cursor.fetchall()
      coordinates_and_seqs = [list(i) for i in results]
      
      return self._get_seq_from_contigs(coordinates_and_seqs,strand=strand)
     
     
     
     
     
  def get_gene_by_stable_id(self,stable_id):
      """using stable_id of gene obtained from the Ensembl database. Returns SeqRecord object of gene coordinates."""
    
      assert(type(stable_id)==str)
      
      query = "select g.biotype, g.status, g.seq_region_start, g.seq_region_end, g.seq_region_strand,g.description,g.source,r.name,x.display_label,x.dbprimary_acc, e.db_name \
	      from gene g,xref x,external_db e,seq_region r \
	      where g.stable_id='%s' and r.seq_region_id=g.seq_region_id \
	      and g.display_xref_id=x.xref_id and x.external_db_id=e.external_db_id;" % stable_id

      self.cursor.execute(query)
      assert self.cursor.rowcount==1 
      result = self.cursor.fetchone()
      
      biotype,status,start,end,strand,description,source,chromosome = str(result[0]),str(result[1]),int(result[2]),int(result[3]),int(result[4]),str(result[5]),str(result[6]),str(result[7])
      name,ref,ref_db = str(result[8]),str(result[9]),str(result[10])
      seq = self._get_seq_by_stable_id(stable_id=stable_id,strand=strand)
      
      s = self._to_SeqRecord(seq=seq, name = name,biotype=biotype,status=status,gene_id=stable_id,description=description,genome_name=self.genome_name,start = start,end=end,strand = strand,source=source,chromosome=chromosome,ref=ref,ref_db=ref_db)
      return s
         
      
  def get_genes_by_description(self,description):
      """Using description of gene obtained from the Ensembl database. Searching for a given pattern in descriptions of genes. 
      Name pattern passed to the SQL like statement (i.e. '%','_',... are used as special characters). See http://dev.mysql.com/doc/refman/5.6/en/string-comparison-functions.html#operator_like.
      example:
      'myosin'   no results
      'myosin%'  e.g.: myosin IE [Source:HGNC Symbol;Acc:7599]
      '%myosin%' e.g.: tropomyosin 1 (alpha) [Source:HGNC Symbol;Acc:12010]
      """
      assert(type(description)==str)
      query = "select stable_id from gene where description like '"+description+"';"
      self.cursor.execute(query)
      results = self.cursor.fetchall()
      if self.cursor.rowcount==1:
	genes = [results[0]]
      else:
	genes = [list(i)[0] for i in results]
      
      g=[]
      for stable_id in genes:
	pattern = "LRG"
	r = re.search(pattern,stable_id)
	if not r: #without LRG sequence (see http://www.lrg-sequence.org/LRG/)
	  t = self.get_gene_by_stable_id(stable_id)
	  g.append(t)
	  
      return g

      
  def get_gene_by_name(self,name):
      """using name of gene obtained from the Ensembl database. Returns SeqRecord object of gene coordinates."""
      assert(type(name)==str)
      query = "select gene.stable_id from gene, object_xref, xref, external_db \
		where gene.gene_id = object_xref.ensembl_id  and   object_xref.ensembl_object_type = 'Gene' \
		and object_xref.xref_id =   xref.xref_id and xref.external_db_id = external_db.external_db_id \
		and xref.display_label = '%s' group by gene.stable_id ;" % name
      self.cursor.execute(query)
      results = self.cursor.fetchall()
      if self.cursor.rowcount==1:
	results = [results[0]]
      else:
	results = [list(i)[0] for i in results]

      g=[]
      for stable_id in results:
	g.append(self.get_gene_by_stable_id(stable_id[0]))

      return g
      

  def get_region(self,coord_name,coord_start,coord_end,coord_strand,region_id=None):
      """coord_name : a string identifying a chromosome.
	 coord_start : the sequence's start position.
	 coord_end : the sequence's end position.
	 coord_strand : the sequence's region strand: 1 - forward; -1 - reverse."""
      assert(type(coord_name)==str)
      assert(type(coord_start)==int)
      assert(type(coord_end)==int)
      assert(coord_strand==1 or coord_strand==-1)
      assert(coord_start < coord_end)
      
      coord_start,coord_end = str(coord_start),str(coord_end)

      query = " select s.cmp_start, s.cmp_end,"+coord_start+","+coord_end+", s.ori, d.sequence,r.name,s.asm_start, s.asm_end from assembly s, dna d, seq_region r where   \
	  (("+coord_end+" >= s.asm_start and "+coord_end+" <= s.asm_end) or ("+coord_start+" >= s.asm_start and "+coord_start+" <= s.asm_end) or   \
	  ("+coord_start+"<= s.asm_start and "+coord_end+">= s.asm_end)) and r.seq_region_id =s.asm_seq_region_id and r.name='"+coord_name+"'   \
	  and d.seq_region_id = s.cmp_seq_region_id; "

	  
      coord_start,coord_end = int(coord_start),int(coord_end)
      
      
      self.cursor.execute(query)
      results = self.cursor.fetchall()
      coordinates_and_seqs = [list(i) for i in results]

      seq = self._get_seq_from_contigs(coordinates_and_seqs,strand=coord_strand)
      
      assert(len(seq)>0)

      if region_id!=None:
	return self._to_SeqRecord(seq=seq,genome_name=self.genome_name,start=coord_start,end=coord_end,strand=coord_strand,chromosome=coord_name, gene_id=region_id)
      return self._to_SeqRecord(seq=seq,genome_name=self.genome_name,start=coord_start,end=coord_end,strand=coord_strand,chromosome=coord_name)

    
 
 
 
    
    
class Compara(object):
    """Class representing information contained in the Compara database."""
    
    def __init__(self, release):
      
      assert(type(release)==int and release > 0)

      connection = Ensembl_Connection(ENSEMBL_HOST,ENSEMBL_USER,ENSEMBL_PASSWD)
      cursor,con = connection.connect_with_ensembl(release=release)
      self.cursor,self.con = cursor,con

      con.select_db("ensembl_compara_"+str(release))
      
     
    def close(self):
      """closes connection with the Ensembl database."""
      
      self.con.close()
      
      
      
    def available_relationships(self):
      """returns list of available relationships in the Ensembl database."""
      
      query = "select distinct description from homology;"
      self.cursor.execute(query)
      results = self.cursor.fetchall()
      return [i[0] for i in results]
      
      
      
    def available_spieces(self):
      """returns list of available spiecies in the Ensembl database."""
      
      query = "select name from genome_db;"
      self.cursor.execute(query)
      results = self.cursor.fetchall()
      return [i[0].replace("_"," ").title() for i in results]
    
    
	
    def get_related_genes(self,species = [], ref_gene_id = None ,relationship = None):
      """
      species : a list of at least 1 species name.
      ref_gene_id : stable_id of gene on the reference genome.
      relationship : relationship between genes (reference and those given in spiecies). A list of possible relationships can be obtained using the available_relationships method.
      """

      assert(len(species)>=1)
      assert(type(ref_gene_id)==str)
      assert(type(relationship)==str)
       
      for sp in species:
	if sp.title() not in self.available_spieces():
	  raise SpieciesError("Invalid name of spiecies")
      if relationship not in self.available_relationships():
	raise RelationshipError("Invalid name of relationship")
      

      species = [i.replace(" ","_") for i in species]
      spieces_names = ",".join(map(repr,species))
      
      query = " select m2.stable_id,g.name from member m,homology_member hm, \
	      homology h,homology_member hm2, member m2, genome_db g \
	      where m.stable_id='%s' and \
	      m.member_id=hm.member_id and h.homology_id=hm.homology_id and \
	      h.description='%s' and hm2.homology_id=hm.homology_id \
	      and m2.member_id=hm2.member_id and m.stable_id!=m2.stable_id and \
	      m2.genome_db_id=g.genome_db_id and g.name in (%s); " % (ref_gene_id,relationship,spieces_names)

      self.cursor.execute(query)
      results = self.cursor.fetchall()

      rowcount = self.cursor.rowcount
      if rowcount==0: #e.g. when list of given species contains the same name's of species
	return []
      elif rowcount==1:
	return [(results[0][0],results[0][1].replace("_"," ").title())]
      else:
	d=[]
	for i in results:
	  d.append((i[0],i[1].replace("_"," ").title()))
	return d
    
    

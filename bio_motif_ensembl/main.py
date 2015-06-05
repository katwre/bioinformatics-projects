
import subprocess as sub
import StringIO

from Ensembl import *
from CNS import *
from Stats import *

from Bio import SeqIO,AlignIO
from Bio.Align.Applications import MuscleCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC



command_MEME = "./meme_4.8.1/scripts/meme"
command_AlignACE="./alignace2004/AlignACE"
command_muscle="./muscle3.8.31_i86linux64"

seq_region = 10000
ident = 0.70
block_length = 100

strand = 1


human = Genome(genome_name="Homo Sapiens",release=68)
genes = human.get_genes_by_description("myosin%")


compara = Compara(release=68)
rat = Genome(genome_name='Rattus Norvegicus',release=68)
mouse = Genome(genome_name='Mus Musculus',release=68)





# 1. Znalezienie ortologow, sekwencji o 10 kb oddalonych od TSS'ow i ich multiuliniowienie

for gene in genes:
  
  orthologs = compara.get_related_genes(['Rattus Norvegicus','Mus musculus'],gene.id,"ortholog_one2one")
  
  if len(orthologs)==2:
    
    gene = human.get_gene_by_stable_id(gene.id) #gen referencyjny
    h = human.get_region(coord_name=gene.annotations['chromosome'],coord_start=gene.annotations['start'] - seq_region, coord_end=gene.annotations['start'] + seq_region - 1,coord_strand=strand,region_id="1")
    seqs=[h]
    
    for o in orthologs:
      
      if o[1]=='Rattus Norvegicus':

  r = rat.get_gene_by_stable_id(o[0])
	r = rat.get_region(coord_name=r.annotations['chromosome'],coord_start=r.annotations['start'] - seq_region, coord_end=r.annotations['start'] + seq_region - 1,coord_strand=strand,region_id="2")
	seqs.append(r)

      if o[1]=='Mus Musculus':
      
	m = mouse.get_gene_by_stable_id(o[0])
	m = mouse.get_region(coord_name=m.annotations['chromosome'],coord_start=m.annotations['start'] - seq_region,coord_end=m.annotations['start'] + seq_region - 1,coord_strand=strand,region_id="3")
	seqs.append(m)
    
    if len(seqs)==3: 
    
      output_handle = open("tmp_files/"+gene.name+".fasta", "wb")
      SeqIO.write(seqs, output_handle, "fasta")
      output_handle.close()
      
      muscle_cline = MuscleCommandline(command_muscle, input="tmp_files/"+gene.name+".fasta", out="tmp_files/"+gene.name+".aln")
      muscle_cline()
     
    
    
    
    
# 2. Znalezienie konserwowanych blokow sekwencji


j=0
for gene in genes:
  
    try:
      
      MSA = AlignIO.read("tmp_files/" + gene.name + ".aln", "fasta",alphabet=IUPAC.ambiguous_dna)
      cns = CNS(MSA)
      bloki = cns.get_conserved_blocks_of_MSA(block_length,ident)

      
      for b in bloki:
	output_handle = open("tmp_segments/"+gene.name+".fasta", "w")
	SeqIO.write(b[0]._records, output_handle, "fasta")
	output_handle.close()
	j+=1
      
    except:
      pass




# 3. Znaleznie nadreprezentowanych motywow w blokach sekwencji


conn = sub.Popen(["ls", "./temp_segments"],stdout=sub.PIPE, stdin=sub.PIPE,stderr=sub.PIPE)
temp = conn.stdout.read()
x=StringIO.StringIO(temp)
msa_list = [i.rstrip() for i in list(x)]


plik = open("MSA.fasta","w")
for i in msa_list:
  align = AlignIO.read("./tmp_segments/"+str(i), "fasta")
  plik.write(align.format("fasta"))
plik.close()  


motifs = []


#MEME
def meme_motifs(CNS_name,command_MEME,q):
  
  motifs = []

  conn = sub.Popen([command_MEME,"-dna",CNS_name, "-text","-maxw","10"],stdout=sub.PIPE, stdin=sub.PIPE,stderr=sub.PIPE)
  temp = conn.stdout.read()
  x=StringIO.StringIO(temp)
  records = list(Motif.parse(x,"MEME"))
  w=0
  for i in records:
    i.name = str(q)+"_meme_"+str(w)
    w+=1
  
  return records

meme = alignace_motifs("tmp_segments/MSA.fasta",command_MEME,0)
motifs.extend(meme)  
  
  
  
  
#ALIGNACE  
def alignace_motifs(CNS_name,command_AlignACE,q):
 
  motifs = []

  conn = sub.Popen([command_AlignACE,"-i",CNS_name,"-numcols","10","-gcback","0.525"],stdout=sub.PIPE, stdin=sub.PIPE,stderr=sub.PIPE)
  temp = conn.stdout.read()
  x=StringIO.StringIO(temp)
  records=list(Motif.parse(x,"AlignAce"))
  w=0
  for i in records:
    i.name = str(q)+"_alignace_"+str(w)
    w+=1

  return records

alignace = alignace_motifs("tmp_segments/MSA.fasta",command_AlignACE,0)
motifs.extend(alignace)  
  
  
  
#JASPAR  
def from_jaspar_site(motifs):
  
  import urllib
  import StringIO
  
  
  link = "http://jaspar.genereg.net/html/DOWNLOAD/jaspar_CORE/non_redundant/by_tax_group/vertebrates/FlatFileDir/"
  
  m=[]
  for motif in motifs:
    try:
      data = urllib.urlopen(link+motif).read()
    except:
      print "Motif's id %s does not exist in jaspar database" % motif
      motifs.remove(motif)
      continue
    x=StringIO.StringIO(data)
    k = Motif.read(x,"jaspar-pfm")
    k.name = motif
    m.append(k)
    
  if m==[]:
    print "Connection error or invalid names of motifs"
    exit()
  return m 

jaspar_motifs_names = ['MA0002.2.pfm', 'MA0003.1.pfm', 'MA0004.1.pfm', 'MA0006.1.pfm', 'MA0007.1.pfm', 'MA0009.1.pfm', 'MA0014.1.pfm', 'MA0017.1.pfm', 'MA0018.2.pfm', 'MA0019.1.pfm', 'MA0024.1.pfm', 'MA0025.1.pfm', 'MA0027.1.pfm', 'MA0028.1.pfm', 'MA0029.1.pfm', 'MA0030.1.pfm', 'MA0031.1.pfm', 'MA0032.1.pfm', 'MA0033.1.pfm', 'MA0035.2.pfm', 'MA0036.1.pfm', 'MA0037.1.pfm', 'MA0038.1.pfm', 'MA0039.2.pfm', 'MA0040.1.pfm', 'MA0041.1.pfm', 'MA0042.1.pfm', 'MA0043.1.pfm', 'MA0046.1.pfm', 'MA0047.2.pfm', 'MA0048.1.pfm', 'MA0050.1.pfm', 'MA0051.1.pfm', 'MA0052.1.pfm', 'MA0055.1.pfm', 'MA0056.1.pfm', 'MA0057.1.pfm', 'MA0058.1.pfm', 'MA0059.1.pfm', 'MA0060.1.pfm', 'MA0061.1.pfm', 'MA0062.2.pfm', 'MA0063.1.pfm', 'MA0065.2.pfm', 'MA0066.1.pfm', 'MA0067.1.pfm', 'MA0068.1.pfm', 'MA0069.1.pfm', 'MA0070.1.pfm', 'MA0071.1.pfm', 'MA0072.1.pfm', 'MA0073.1.pfm', 'MA0074.1.pfm', 'MA0075.1.pfm', 'MA0076.1.pfm', 'MA0077.1.pfm', 'MA0078.1.pfm', 'MA0079.2.pfm', 'MA0080.2.pfm', 'MA0081.1.pfm', 'MA0083.1.pfm', 'MA0084.1.pfm', 'MA0087.1.pfm', 'MA0088.1.pfm', 'MA0089.1.pfm', 'MA0090.1.pfm', 'MA0091.1.pfm', 'MA0092.1.pfm', 'MA0093.1.pfm', 'MA0095.1.pfm', 'MA0098.1.pfm', 'MA0099.2.pfm', 'MA0100.1.pfm', 'MA0101.1.pfm', 'MA0102.2.pfm', 'MA0103.1.pfm', 'MA0104.2.pfm', 'MA0105.1.pfm', 'MA0106.1.pfm', 'MA0107.1.pfm', 'MA0108.2.pfm', 'MA0109.1.pfm', 'MA0111.1.pfm', 'MA0112.2.pfm', 'MA0113.1.pfm', 'MA0114.1.pfm', 'MA0115.1.pfm', 'MA0116.1.pfm', 'MA0117.1.pfm', 'MA0119.1.pfm', 'MA0122.1.pfm', 'MA0124.1.pfm', 'MA0125.1.pfm', 'MA0130.1.pfm', 'MA0131.1.pfm', 'MA0132.1.pfm', 'MA0133.1.pfm', 'MA0135.1.pfm', 'MA0136.1.pfm', 'MA0137.2.pfm', 'MA0138.2.pfm', 'MA0139.1.pfm', 'MA0140.1.pfm', 'MA0141.1.pfm', 'MA0142.1.pfm', 'MA0143.1.pfm', 'MA0144.1.pfm', 'MA0145.1.pfm', 'MA0146.1.pfm', 'MA0147.1.pfm', 'MA0148.1.pfm', 'MA0149.1.pfm', 'MA0150.1.pfm', 'MA0151.1.pfm', 'MA0152.1.pfm', 'MA0153.1.pfm', 'MA0154.1.pfm', 'MA0155.1.pfm', 'MA0156.1.pfm', 'MA0157.1.pfm', 'MA0158.1.pfm', 'MA0159.1.pfm', 'MA0160.1.pfm', 'MA0161.1.pfm', 'MA0162.1.pfm', 'MA0163.1.pfm', 'MA0164.1.pfm', 'MA0258.1.pfm', 'MA0259.1.pfm', 'MA0442.1.pfm']
jaspar_motifs = from_jaspar_site([jaspar_motifs_names[0]]) 
motifs.extend(jaspar_motifs)






stats = Stats(MSAs)
for motif in motifs:

  d = stats.overrepresented_motif(motif=motif,fpr=0.0001)
  if d:
     plik = open("./nadreprezentowane_motywy/"+d[1].name+".fasta","w")
     plik.write(d[1]._to_fasta())
     plik.close()
     d[1].weblogo("./nadreprezentowane_motywy/"+d[1].name+".png")
  


#!/usr/bin/python
# -*- coding: utf-8 -*-

from optparse import OptionParser
import sys


version = "0.01"
usage = "usage: %prog [options]\n \
Example:\n \
python merge_cif.py -r example_cif_files/14ni_4f_ref.cif  -e ./example_cif_files/14ni_4f_exp.cif -o ./example_cif_files/14ni_4f_out.cif"

usage_description = """Merges two files in CIF format. Tested with python 2.7.
Files in CIF format may consist of one or more data blocks which begin with the declaration data_, but
this program binds them to only one block named data_merged.
For more information obout CIF format see:
http://www.iucr.org/__data/iucr/cif/standard/cifstd1.html"""

epilog = u"merge_cif.py Copyright (C) 2013 Katarzyna Wręczycka"

help_interactive = "If not provided will be asked for interactively."


def _is_valid(options):
    return options.__dict__.values() == [None]*3

def parse_options():
    """
    Handles options parsing and getting filenames of
    reference and experimental file in CIF format. 
    Returns a pair (options, args) in the same format as
    OptionParser.parse_args from optparse module.
    """
    parser = OptionParser(usage, description=usage_description, version=version,
                        epilog=epilog)
    parser.add_option("-r", "--ref", help=help_interactive)
    parser.add_option("-e", "--exp", help=help_interactive)
    parser.add_option("-o", "--out", help="output filename.")
    (options, args) = parser.parse_args()
    
    if _is_valid(options):
        parser.print_help()
        sys.exit(1)
    if options.ref is None:
        sys.stdout.write("Enter reference file in CIF format: ")
        options.ref = sys.stdin.readline()[:-1]
    if options.exp is None:
        sys.stdout.write("Enter experimental file in CIF format: ")
        options.ref = sys.stdin.readline()[:-1]
    return (options, args)







def parse_variable(current_line,fl,result):
  """Parses all variables with their values (those having many lines or only one line at the same time
     and not belonging to 'loop_')
     
  So variables can look like:
  _audit_contact_author.name  'G.G.Dodson'
  or
  _audit_creation_date
  ;
  'Wed Dec 12 18:00:41 2012'
  ;
  """

  variable_name = fl[current_line].split()[0]
  start = current_line

 
  current_line += 1
  while current_line<len(fl) and (fl[current_line].split()==[] or ((not _is_variable(fl[current_line])) and (not _is_loop(fl[current_line])))):
    current_line += 1
  result.append((variable_name,fl[start:current_line]))
  return current_line
  

def parse_loop(current_line,fl,result):
  """Parses all loops with their values, for example:
  loop_
  _struct_sheet.id
  _struct_sheet.number_strands
  A     5
  B     2
  C     2
  """

  variable_name = "loop_"
  start = current_line

  current_line += 1
  while current_line < len(fl) and (fl[current_line].split()==[] or _is_variable(fl[current_line])):
    current_line += 1
  while current_line < len(fl) and (fl[current_line].split()==[] or (not _is_variable(fl[current_line]) and not _is_loop(fl[current_line]))):
    current_line += 1
    
  result.append((variable_name,fl[start:current_line]))
  return current_line


  

def _is_loop(line):
  return line.split()[0]=="loop_"
  
def _is_variable(line):
  return line.split()[0][0]=="_"
 
 
 
 
 
def _data(f):
  """calculates the first occurrence of pattern "data_" in ref file (for example "data_xcalibur")
  """
  for i, v in enumerate(f):
    if v[:5]=="data_":
      return i+1
  
def parse_file(filename):
  
  ref_file = open(filename,"r").readlines()
  current_line = _data(ref_file)
  result = []

  while current_line != len(ref_file):

    if ref_file[current_line][:5]=="data_" or ref_file[current_line]=="\n":
      current_line+=1

    elif _is_variable(ref_file[current_line]):
      current_line = parse_variable(current_line,ref_file,result)
    
    elif _is_loop(ref_file[current_line]):
      current_line = parse_loop(current_line,ref_file,result)
    

  return result



 
def merge_and_save(ref, exp, out):

    r = dict(i for i in ref if i[0]!="loop_")
    e = dict(i for i in exp if i[0]!="loop_")
    r.update(e)
    
    outfile = open(out,"w")
    outfile.write("data_merged\n")
    
    w = []
    for i in ref:
      if i[0]!="loop_":
	outfile.write("".join(r[i[0]]))
      else:
	outfile.write("".join(i[1]))
      w.append(i[0])
      
    for i in exp:
      if i[0]=="loop_":
	outfile.write("\n"+"".join(i[1]))
      elif i[0]!="loop_" and i[0] not in w:
	outfile.write("".join(i[1]))
	
    outfile.close()

    
     
      
def main():
  
    (options,args) = parse_options()    
    ref = parse_file(options.ref)
    exp = parse_file(options.exp)

    merge_and_save(ref, exp, options.out)
    
    
    
if __name__ == '__main__':
    main()
    
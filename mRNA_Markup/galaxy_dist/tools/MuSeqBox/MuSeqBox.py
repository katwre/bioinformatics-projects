#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import subprocess as sub

museqbox_path = "/home/user/MUSEQBOX/bin/MuSeqBox" 



##Parse Command Line
infile = sys.argv[1]
outfile = sys.argv[2]
arguments = sys.argv[3:]


def _args4museqbox(args):
  
  tlist = []
  for i in xrange(len(args)):
    
    if args[i]=="-n":
      if args[i+1]!="-":
	tlist.extend(["-n",args[i+1]])
	
    if args[i]=="-s":
      if args[i+1]!="-":
	tlist.extend(["-s",args[i+1]])
	
    if args[i]=="-F":
      if args[i+1]!="-":
	tlist.extend(["-F",args[i+1],args[i+2],args[i+3],args[i+4],args[i+5],args[i+6]])
	
    if args[i]=="-M":
      if args[i+1]!="-":
	tlist.extend(["-M",args[i+1],args[i+2]])
	
    if args[i]=="-q":
      tlist.append("-q")
      
  return " ".join(tlist) 
  
arguments1 = _args4museqbox(arguments) 
  
  
  
##Call MuSeqBox

k = museqbox_path + " -i " + infile + " " + arguments1 + " -o " + outfile
sub.Popen(k,stdout=sub.PIPE, shell=True) 



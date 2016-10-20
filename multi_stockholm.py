#!/usr/bin/python

#This script takes many fasta alignment files (MSA's) and writes them
#all to a single stockholm file suitable for analysis in X-rate.  

#Usage similar to 'cp' command:
#
#~$multi_stockholm.py MSA1 MSA2 MSA3... outfile.stock
#
#OR
#
#~$multi_stockholm.py *.fas outfile.stock (all files ending in .fas)

import sys
import os

from Bio.Alphabet import generic_protein
from Bio import AlignIO

#THIS CODE ACTUALLY DOES THE WORK
########################################################################
msa_files = []

for msa_file in sys.argv[1:-1]:
	msa_files.append(AlignIO.read(msa_file, "fasta"))
#AlignIO.write(msa_files, sys.argv[-1], "stockholm") #user-proofed below
########################################################################



#ALL OF THIS OTHER CODE JUST KEEPS YOU FROM RUINING YOUR LIFE
########################################################################
if os.path.isfile(sys.argv[-1]):
	user_response = raw_input("Warning! destination file %s already exists. Overwrite (Y/n)?" % sys.argv[-1])
	
	if user_response  in ["y","Y","yes","Yes","YES"]:
		AlignIO.write(msa_files, sys.argv[-1], "stockholm")
	
	elif user_response in ["n","N","no","No","NO"]:
		print "Aborting...."
		raise SystemExit
		
	else:
		print "Invalid option. Aborting...."
		raise SystemExit
		
else:
	AlignIO.write(msa_files, sys.argv[-1], "stockholm")
########################################################################

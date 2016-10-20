#!/usr/bin/python



import sys
import os

from Bio.Alphabet import generic_protein
from Bio import AlignIO

msa_files = []

for msa_file in sys.argv[1:-1]:
	msa_files.append(AlignIO.read(msa_file, "fasta"))

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

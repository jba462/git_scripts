#!/usr/bin/python

import sys

fasta_file = sys.argv[1]

with open(fasta_file, 'r') as inp:
	for line in inp:
		if line.startswith(">"):
			print ">" + line.split("|")[1]
		else:
			print line.rstrip("\n")

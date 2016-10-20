#!/usr/bin/python

import sys

map_file = sys.argv[1]

tree_file = sys.argv[2]

with open(map_file, 'r') as inp:
	maps = [line.rstrip('\r\n') for line in inp]
	#print maps
with open(tree_file, 'r') as inp:
	as_list = [line.rstrip('\r\n') for line in inp]
	#print len(as_list)
for record in maps:
	id_only = record.split(" ")[0]
	#print id_only
	name_only = " ".join(record.split(" ")[1:])
	#print name_only
	for index, line in enumerate(as_list):
		as_list[index] = as_list[index].replace(name_only,(name_only + " " + id_only))
		print as_list[index]
with open(tree_file + ".mod", 'w') as outp:
	for line in as_list:	
		outp.write(line + "\n")

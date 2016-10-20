#!/usr/bin/python
#
#to run:
#
#~$ss_dataframe_maker.py aa_alignment 'path_to_horiz_files' 'path_to_ss_files'
#
#NOTES: Encapsulate your paths in single quotes; alignment must be in Fasta format

import sys
sys.path.append("/home/joseph/bin") #uses the jbamsa library

aa_alignment = sys.argv[1]
path_to_horiz_file = sys.argv[2]
path_to_ss_file = sys.argv[3]
path_to_ss2_file = sys.argv[4]

import glob

from jbamsa import FasMSA


alignment_data = FasMSA(aa_alignment)
headers = alignment_data.headers
aligned_sequences = alignment_data.alignedSequences

###
#NOTE: THIS IS THE FUNCTION THAT EXCRACTS THE SEQUENCE IDENTIFIER FROM
#      THE EACH SEQUENCE HEADER IN YOU ORIGINAL ALIGNED FASTA FILE. IF
#      YOU WANT TO USE THIS SCRIPT, YOU MAY HAVE TO MODIFY THIS FUNCTION
#      IN ORDER FOR IT TO WORK PROPERLY.
def get_identifier(sequence_header):
	return header.split("_")[1]
###

with open(aa_alignment + "_ss_df_full.txt", 'w') as outp:
	outp.write("#Names" + "\t" + "\t".join([str(i) for i in range(1,(len(aligned_sequences[0]) + 1))]) + "\n")
	
for index,header in enumerate(headers):
	#print header
	residues = []
	predictions = []
	confidence_values = []
	propensity_C = []
	propensity_H = []
	propensity_E = []
	propensity_C2 = []
	propensity_H2 = []
	propensity_E2 = []
	
	#print path_to_horiz_file.rstrip("/") + "/" + get_identifier(header) + ".horiz"
	#print glob.glob(path_to_horiz_file.rstrip("/") + "/")# + get_identifier(header) + ".horiz")
	if len(glob.glob(path_to_horiz_file.rstrip("/") + "/" + get_identifier(header) + ".horiz")) < 1:
		print "Error: horiz file not found for %s" % header
		raise SystemExit
	with open(glob.glob(path_to_horiz_file.rstrip("/") + "/" + get_identifier(header) + ".horiz")[0], 'r') as inp:
		horiz_file_lines = [line.rstrip("\r\n") for line in inp]
	
	for line in horiz_file_lines:
		if line.startswith("  AA: "):
			residues += list(line.replace("  AA: ",""))
		if line.startswith("Pred: "):
			predictions += list(line.replace("Pred: ",""))
		if line.startswith("Conf: "):
			confidence_values += list(line.replace("Conf: ",""))

	if len(glob.glob(path_to_ss_file.rstrip("/") + "/" + get_identifier(header) + ".ss")) < 1:
		print "Error: ss file not found for %s" % header
		raise SystemExit

	with open(glob.glob(path_to_ss_file.rstrip("/") + "/" + get_identifier(header) + ".ss")[0], 'r') as inp:
		ss_file_lines = [line.rstrip("\r\n") for line in inp]

	for line in ss_file_lines:
		propensity_C.append(line.split()[3])
		propensity_H.append(line.split()[4])
		propensity_E.append(line.split()[5])

	if len(glob.glob(path_to_ss2_file.rstrip("/") + "/" + get_identifier(header) + ".ss2")) < 1:
		print "Error: ss2 file not found for %s" % header
		raise SystemExit

	
	with open(glob.glob(path_to_ss2_file.rstrip("/") + "/" + get_identifier(header) + ".ss2")[0], 'r') as inp:
		ss2_file_lines = [line.rstrip("\r\n") for line in inp if not line.startswith("#") and line.split()]
	
	for line in ss2_file_lines:
		propensity_C2.append(line.split()[3])
		propensity_H2.append(line.split()[4])
		propensity_E2.append(line.split()[5])	
	
	if residues == list(aligned_sequences[index].replace("-","")):
		with open(aa_alignment + "_ss_df_full.txt", 'a') as outp:
			p_index = 0
			outp.write(header + "\t")
			for s_index,site in enumerate(list(aligned_sequences[index])):
				if site == "-":
					outp.write("NaN\t")
				else:
					outp.write(','.join([residues[p_index],
										predictions[p_index],
										confidence_values[p_index],
										propensity_C[p_index],
										propensity_H[p_index],
										propensity_E[p_index],
										propensity_C2[p_index],
										propensity_H2[p_index],
										propensity_E2[p_index]
										]) + "\t")
					p_index += 1
			outp.write("\n")
	
	else:
		print "ERROR: %s does not match the original sequence" % header
		raise SystemExit
		

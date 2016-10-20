#!/usr/bin/python

import sys

input_file = sys.argv[1]

coverage_per_site = sys.argv[2]


#--------------------------------------------------------------------------------------------
#############################################################################################
#############################################################################################
#
#The following setion (#----- to #------) can be pasted into other scripts (make sure to import the "sys" module
#and set sys.argv[1] to "input_file"). It takes a FASTA MSA as input and returns several lists (see below) that are useful
#in general statistical analyses.
#
#PLEASE COPY THE ABOVE INSTRUCTIONS AS WELL
#
#  Joseph
#
#############################################################################################
#This function linearizes sequences in an input file (removes text-wrapping) and
#returns a list of sequences WITHOUT their headers, BUT in the same
#order as in the original alignment.
#
#Input: a complete list of lines from an alignment file (including headers) with newlines (\r, \n) removed
#
def linearize(alignment_as_list):
	
	#Index changes to 0 (first item in list) after the loop below reads the ">" character
	#on the first line (see below for exlpanation).
	index = -1

	sequences = []
	
	#Iterate over the list of lines (input): 
	for line in alignment_as_list:
		
		#If the line is a header:
		if line.startswith(">"):
			
			#Add 1 to the "index" variable. On the first iteration, this changes "index"
			#from -1 to 0, so it points to the first item in "sequences" which is 
			#"created" just below.
			#
			#On subsequent iterations, this line moves "index" to the next position in 
			#the list and the line below will, again, create an empty string for
			#"index" to point to.
			index += 1
			
			#Append an empty string to "sequences".  This gives "index" something to 
			#point to before any strings are added (see above).
			sequences.append("")
		
		#If the line is NOT a header:
		else:
			
			#If "index" is NOT pointing to an empty string (meaning some portion of this sequence
			#has already been added):
			if sequences[index] != "":
				
				#Add the current line in "alignment_as_list" to "index" in "sequences"
				sequences[index] = "".join([sequences[index], line])
			
			#If the "index" in "sequences" is currently empty:
			else:
				
				#Set it equal to the current line from "alignment_as_list"
				sequences[index] = line
	
	#OUTPUT is a list of "unwrapped" sequences, in original order, from the alignment file.
	return sequences
############################################################################################


#Open and read in the FASTA MSA
with open(input_file, 'r') as inp:
	
	#Convert the file to a list of lines (with newline characters remmoved)
	file_as_list = [line.rstrip("\r\n") for line in inp]
	
	#Create a list of all headers in the file
	header_list = [line for line in file_as_list if line.startswith(">")]
	
	#Create a list of linearized, aligned  sequences from the original MSA, 
	#using the function above and "file_as_list" as input.
	aligned_sequence_list = [line for line in linearize(file_as_list)]
	
	#Create a list of unaligned (totally de-gapped) sequences.
	unaligned_sequence_list = [line.replace("-","") for line in linearize(file_as_list)]
	
	#Create a list of end-trimmed (internal gaps left alone) sequences.
	end_trimmed_sequence_list = [line.strip("-") for line in linearize(file_as_list)]

############################################################################################
############################################################################################
#-------------------------------------------------------------------------------------------




#ALL FUNCTIONS FOR THIS SCRIPT ARE DEFINED BELOW:
################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





#This computes proportion of non-gap characters in each column and returns a list of proportions
#for each column.
############################################################################################
def coverage_per_site_plot(alignment):
	
	coverage_per_site = []
	
	for index in range(len(alignment[0])):
		
		residues = 0
		
		gaps = 0
		
		for sequence in alignment:
			
			if sequence[index] == "-":
				
				gaps += 1
				
			else:
				
				residues += 1
				
		coverage_per_site.append(float(float(residues) / (float(residues) + float(gaps))))
		
	return coverage_per_site
#################################################################################################




#This returns a single number: the proportion of vertical columns in an alignment that contain 
#no gaps.
#################################################################################################
def gap_free_columns_per_alignment(alignment):

	gap_free_columns = 0

	for index in range(len(alignment[0])):
		
		column = [sequence[index] for sequence in alignment]
		
		if "-" not in column:
			
			gap_free_columns += 1
			
	return float(float(gap_free_columns) / float(len(alignment[0])))
################################################################################################




#This returns the proportion of non-gap (DNA/AA) characters in the whole alignment.	
################################################################################################	
def non_gap_chracters_per_alignment(alignment):			
		
		all_sites = [site for sequence in aligned_sequence_list for site in sequence]
		
		non_gap_sites = [site for site in all_sites if site != "-"]
		
		return float(float(len(non_gap_sites) / float(len(all_sites))))
################################################################################################




#This returns the shortest OR longest sequence in a list.
################################################################################################
def min_max_length(sequence_list, min_max):
	
	if min_max == "min":
			
		return len(min(sequence_list, key=len))

	if min_max == "max":
			
		return len(max(sequence_list, key=len))
###############################################################################################



#This returns the smallest OR largest alignment coverage value (sequence length / alignment length)
##############################################################################################
def min_max_alignment_coverage(sequence_list, min_max):
	
	if min_max == "min":
	
		return float(len(min(sequence_list, key=len))) / float(len(aligned_sequence_list[0]))
	
	elif min_max == "max":
		
		return float(len(max(sequence_list, key=len))) / float(len(aligned_sequence_list[0]))
##############################################################################################





#Compute the pairwise distance between two aligned sequences (relative to total alignment length)
#############################################################################################	
def pairwise_distance_alignment_length(sequence1, sequence2):
	
	mismatches = 0
	
	for index, state in enumerate(sequence1):
		
		if sequence1[index] != sequence2[index]:
			
			mismatches += 1

	return float(float(mismatches) / float(len(sequence1)))
#############################################################################################






#Compute pairwise distance in the same way as the program "distmat" embedded in EMBOSS
#############################################################################################
def pairwise_distance_like_emboss(sequence1, sequence2):
	
	shared_sites = 0
	
	mismatches = 0
	
	for index, state in enumerate(sequence1):
		
		if sequence1[index] != "-" and sequence2[index] != "-":
			
			shared_sites += 1
			
			if sequence1[index] != sequence2[index]:
				
				mismatches += 1
	
	try:
				
		return float(float(mismatches) / float(shared_sites))
	
	except:
		
		return 0.0
############################################################################################





#Compute a full distance matrix using the EMBOSS-like method from above. 
############################################################################################
def distance_matrix(alignment):
	matrix = []
	
	for sequence in alignment:
		
		row = []
		
		#"Other" really means "all sequences in the alignment INCLUDING 
		#the one the loop is currently working on."  So the function 
		#returns a square matrix including all self-to-self comparisons.
		for other_sequence in alignment:
			
			row.append(pairwise_distance_like_emboss(sequence,other_sequence))

		matrix.append(row)

	return matrix

#DEBUG - print distance matrix to screen
########################################

#for row in distance_matrix(aligned_sequence_list):

#	print " ".join(row)
############################################################################################






#Find the minimum OR maximum sequence identity in the matrix calculated above.
############################################################################################
def min_max_identity(alignment, min_max):
	
	matrix = distance_matrix(alignment)
	
	values = [value for row in matrix for value in row]
	
	#DEBUG
	#print values
	#print [(1.0-value) for value in values]
	
	if min_max == "max":
	
		return (1.0 - min(values))
		
	elif min_max == "min":
		
		return (1.0 - max(values))
############################################################################################



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
############################################################################################
#END OF FUNCTIONS




if coverage_per_site == "coverage_per_site":
	
	print "|".join([str(stat) for stat in coverage_per_site_plot(aligned_sequence_list)])
	
	#for site_stat in coverage_per_site_plot(aligned_sequence_list):
		
	#	print site_stat

else:
		
	print input_file + " " + str(len(aligned_sequence_list)) + " " + str(len(aligned_sequence_list[0])) + " " + str(min_max_length(unaligned_sequence_list, "min")) + " " + str(min_max_length(unaligned_sequence_list, "max")) + " " + str(min_max_length(end_trimmed_sequence_list, "min")) + " " + str(min_max_length(end_trimmed_sequence_list, "max")) + " " + str(min_max_alignment_coverage(unaligned_sequence_list, "min")) + " " + str(min_max_alignment_coverage(unaligned_sequence_list, "max")) + " " + str(min_max_identity(aligned_sequence_list, "min")) + " " + str(min_max_identity(aligned_sequence_list, "max")) + " " + str(gap_free_columns_per_alignment(aligned_sequence_list)) + " " + str(non_gap_chracters_per_alignment(aligned_sequence_list))

#!/usr/bin/python

class FasMSA(object): #Accepts fasta-formatted mutliple sequence alignment (MSA) files
	
	def __init__(self,name):
		self.name = name
		with open(self.name, 'r') as inp:
			file_as_generator = (line.rstrip("\r\n") for line in inp) #generator object containing pointers to file lines (low memory usage)
			################################################################
			index = -1							#Index changes to 0 (first item in list) after the loop below reads the ">" character										
			sequences = []						#on the first line (see below for exlpanation).	
			headers = []
			for line in file_as_generator:  	#Iterate over the list of lines (input): 
				if line.startswith(">"): 		#If the line is a header:	
					index += 1					#Add 1 to the "index" variable. On the first iteration, this changes "index"
					sequences.append("")		#from -1 to 0, so it points to the first item in "sequences" which is 
					headers.append(line[1:])	#"created" just below.
												#On subsequent iterations, this line moves "index" to the next position in 
												#the list and the line below will, again, create an empty string for
												#"index" to point to.
												#Append an empty string to "sequences".  This gives "index" something to 
												#point to before any strings are added (see above).		
				else: 							#If the line is NOT a header:

					sequences[index] += line #Add the current line in "alignment_as_list" to "index" in "sequences"				
			sequences = tuple(sequences)		#Make "sequences" immutable
			headers = tuple(headers)
			################################################################
		
		self.headers = headers
		
		self.alignedSequences = sequences #see above
		
		self.numberOfSequences = len(self.alignedSequences)
		
		self.numberOfSites = len(self.alignedSequences[0])
		
		self.msaDict = {header: self.alignedSequences[index] for index, header in enumerate(self.headers)}
		
		self.rawSequences = tuple([sequence.replace("-","") for sequence in self.alignedSequences])
		

		
		
	def __pairwise_distance_like_emboss(self, sequence1, sequence2):
		shared_sites = 0
		mismatches = 0	
		for index, state in enumerate(sequence1):
			if sequence1[index] != "-" and sequence2[index] != "-":
				shared_sites += 1	
				if sequence1[index] != sequence2[index]:
					mismatches += 1
		if shared_sites != 0:
			return float(float(mismatches) / float(shared_sites))
		else:
			return 1.0

	def distance_matrix(self):
		matrix = []
		for sequence in self.alignedSequences:
			row = []
			#"Other" really means "all sequences in the alignment INCLUDING 
			#the one the loop is currently working on."  So the function 
			#returns a square matrix including all self-to-self comparisons.
			for other_sequence in self.alignedSequences:
				row.append(self.__pairwise_distance_like_emboss(sequence,other_sequence))
			matrix.append(row)
		return matrix
		
	def species_frequency_dict(self):
		
		sp_freq_dict = {}
		
		#populate sp_freq_dict
		###################################################################
		for header in self.headers:
			if header.split("_")[0] not in sp_freq_dict:
				sp_freq_dict[header.split("_")[0]] = 1
			else:
				sp_freq_dict[header.split("_")[0]] += 1		
		return sp_freq_dict
		###################################################################
	
	def species_count(self):
		
		return len(self.species_frequency_dict())
		
	def cluster_attributes(self):
		
		sp_dict = self.species_frequency_dict()
		
		attributes = [
						"species_specific",
						"non_redundant",
						"mammal_specific",
						"arthropod_specific",
						"all_species_included"
					]
		
		if len(sp_dict) > 1:
			attributes[0] = "not_" + attributes[0]
		
		if len(sp_dict) != len(self.alignedSequences):
			attributes[1] = "not_" + attributes[1]
		
		if not set([key for key in sp_dict]).issubset(["HUMAN","MOUSE","HORSE","AILME","MONDO","ORNAN"]):
			attributes[2] = "not_" + attributes[2]
			
		if not set([key for key in sp_dict]).issubset(["ANOGA","AEDAE","DROME","APIME","DAPPU"]):
			attributes[3] = "not_" + attributes[3]
			
		if len(sp_dict) < 25:
			attributes[4] = "not_" + attributes[4]
		
		return attributes








class BioMatrix(object):
	
	def __init__(self,name):
		
		self.name = name
		
		headers = []
		data = []
		
		with open(self.name, 'r') as inp:
			
			for line in inp:
				line = line.strip("\r\n").split("\t")
				headers.append(line[0])
				data.append(tuple(line[1:]))
		self.headline = headers[0] + "\t" + "\t".join(data[0])		
		self.headers = tuple(headers[1:])
		self.data = tuple(data[1:])
	
	def __rmsd(self,profile1, profile2):
        
		squared_errors = []
			
		for index, value in enumerate(profile1):
				
			squared_error = (float(value) - float(profile2[index]))**2
				
			squared_errors.append(squared_error)

		return float(sum(squared_errors)) / float(len(squared_errors))
            
            		
	def rmsd_matrix(self):
		
		matrix = []
		
		for row in self.data:
			
			mat_row = []
				
			for other_row in self.data:
					
				mat_row.append(__rmsd(float(row), float(other_row)))
					
				matrix.append(mat_row)
					
		return matrix
            
            

			
			
		
	def returnColumns(self):
		
		columns = []
		
		for index, site in enumerate(self.data[0]):
			
			column = tuple([entry[index] for entry in self.data])
			
			columns.append(column)
			
		return tuple(columns)
		
		
		

				


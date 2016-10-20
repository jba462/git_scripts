#!/home/joseph/anaconda2/bin/python

import argparse

parser=argparse.ArgumentParser()

parser.add_argument("input_tree")

parser.add_argument("data_file")

parser.add_argument("--attach_heatmap", 
					"-hm", 
					action="store_true", 
					help="creates a figure containing a tree with a \
					heatmap attached")

parser.add_argument("--attach_alignment", 
					"-am", 
					action="store_true", 
					help="creates a figure containing a tree with an alignment attached")

parser.add_argument("--pfam_full", 
					"-pf", 
					action="store_true", 
					help="creates a figure containing a tree with a \
					pfam domain architecture map attached")
					
parser.add_argument("--pfam_colors", 
					"-pc", 
					type=str, 
					help="Read domain name colors from an input file \
					and add them to the color dictionary")
					


parser.add_argument("--pfam_data_width",  
					type=int,
					default = 4300,
					help="adjusts the width of each domain motif in the\
					 final figure")

parser.add_argument("--hmap_mean",
					"-mn", 
					type=float,
					default=0.5,
					help="set the gradient inflection point for the \
					heatmap")

parser.add_argument("-i",
					"--image", 
					type=str, 
					help="generate an image of the tree/heatmap figure \
					(specify filetype using extension)")

args = parser.parse_args()

input_tree = args.input_tree

data_file = args.data_file

from ete3 import PhyloTree, ClusterTree, ProfileFace, TreeStyle, Tree, TextFace, BarChartFace, SeqMotifFace, NodeStyle

from ete3.treeview.faces import add_face_to_node

def convert_to_dataframe(data):
	
	with open(data, 'r') as inp:
		converted_data_file = [] #holds dataframe as list of lines
		lineCount = 1
		
		for line in inp:
			
			if lineCount == 1:
				line = line.replace("LABELS","#Names") #only first line! 
			
			lineCount += 1
			
			for missing_value_char in ["\tX","\t-","\t?"]:
				line = line.replace(missing_value_char, "\tNaN")
			
			converted_data_file.append(line)
	#print "".join(converted_data_file[:-1])
	return "".join(converted_data_file)
					
def make_binary_values_visible(data):
	return data.replace("\t0", "\t0.1").replace("\t1","\t0.9")

def remove_quotations(tree):
	
	with open(tree, 'r') as inp:
		
		converted_tree_file = [((line).rstrip('\r\n')).replace("'","") for line in inp]
		
	return converted_tree_file[0]
	
if args.attach_heatmap:
	
	heat_figure = ClusterTree(remove_quotations(input_tree), text_array=make_binary_values_visible(convert_to_dataframe(data_file)))
	matrix_max = 1.0
	matrix_min = -0.01
	matrix_mean = args.hmap_mean #default is 0.5
	h_map = ProfileFace(matrix_max, matrix_min, matrix_mean, 39580, 50, 'heatmap', colorscheme=2)
	h_map.margin_left = 30
	h_map.margin_bottom = -0.5
	
	def my_layout(node):
		
		if node:	
			ns = NodeStyle()
			ns["hz_line_width"]=5
			ns["vt_line_width"]=5
			node.set_style(ns)
			
		if node.is_leaf():
			nameface = TextFace(node.name)
			add_face_to_node(nameface, node, 0, aligned=True)
			add_face_to_node(h_map, node, 1, aligned=True)
			node.img_style["size"]=0
		
		elif node.is_root():
			node.dist = 0.01
			node.img_style["size"]=10

	ts = TreeStyle()	
	ts.show_leaf_name=False
	ts.scale=200#len(heat_figure.get_leaf_names())*10
	ts.layout_fn = my_layout
	
	if args.image:
		heat_figure.render(args.image, w=800, units="mm", tree_style = ts)
	
	else:
		heat_figure.show(tree_style = ts)

if args.attach_alignment:
	alignment_figure = PhyloTree(input_tree, alignment=args.data_file, alg_format="fasta")
	
	if args.image:
		alignment_figure.render(args.image, w=100, units="mm")
	
	else:
		alignment_figure.show()
		
if args.pfam_full:
	
	tree = Tree(input_tree, format=1)
	number_of_sequences = len(tree.get_leaf_names())
	print number_of_sequences
	pfam_data = convert_to_dataframe(data_file).split("\n")[1:-1]
	print len(pfam_data)
	domain_dict = {} #header:data dictionary constructed below	
	for pfam_prediction in pfam_data:
		#print pfam_prediction
		domain_dict[pfam_prediction.split("\t")[0]] = pfam_prediction.split("\t")[1:]
	#print domain_dict
	default_colors = ["red", "orange", "yellow", "green", "blue", "indigo", "violet"]
	all_predictions = []
	
	for key in domain_dict:
		#print key + " " + str(len(domain_dict[key]))
		all_predictions += domain_dict[key]
	
	#print all_predictions
	
	domain_counter = 0
	domain_color_dict = {}
	if args.pfam_colors:
		
		with open(args.pfam_colors, 'r') as inp:
			#inp.readlines()
			for line in inp:
				if len(line.split("\t")) == 2:
					domain_color_dict[line.split("\t")[0]] = line.rstrip("\r\n").split("\t")[1]
	
	for item in all_predictions:
		if item not in ["NaN", "0"] and item not in domain_color_dict:
			domain_color_dict[item] = default_colors[domain_counter % len(default_colors)]
			domain_counter += 1
	domain_color_dict["NaN"] = "white"
	domain_color_dict["0"] = "grey"
	#domain_color_dict["1"] = "purple"
	#domain_color_dict["2"] = "orange"
	print domain_color_dict
	
	def build_pfam_motif_face(data):
		output_lists = []
		data_index = 0
		pattern_index = 0
		motif_index = 0
		
		while data_index < len(data):
			
			if pattern_index == 0:
				output_lists.append([data_index, 
									data_index + 1, 
									"[]", 
									0, 
									args.pfam_data_width / number_of_sequences, 
									"transparent",#domain_color_dict[data[data_index]], 
									domain_color_dict[data[data_index]], 
									None])
				pattern_index += 1
				data_index += 1
				
			elif data[data_index] == data[data_index - 1]:
				output_lists[motif_index][1] += 1
				pattern_index += 1
				data_index += 1
			#Write another 'elif' statement here to deal with ending values that do not match their previous values	
			elif data_index == len(data) - 1:
				output_lists.append([data_index, 
									data_index + 1, 
									"[]", 
									0, 
									args.pfam_data_width / number_of_sequences, 
									"transparent",#domain_color_dict[data[data_index]], 
									domain_color_dict[data[data_index]], 
									None])
				pattern_index += 1
				data_index += 1
			else:
				data_index += 1
				pattern_index = 0
				motif_index += 1
		
		#print output_lists		
		return output_lists
				
				

	def pfam_layout(node):
		
		if node:	
			ns = NodeStyle()
			ns["hz_line_width"]=5
			ns["vt_line_width"]=5
			node.set_style(ns)		
	
		if node.is_leaf():
			
			seq = "A" * len(domain_dict[node.name])
			#print len(seq)
			#print seq
			nameface = TextFace(node.name)
			add_face_to_node(nameface, node, 0, aligned=True) #display name
			if node.name.startswith("NEMVE_A7S1E2") or node.name.startswith("NEMVE_A7RG38"):
				print build_pfam_motif_face(domain_dict[node.name]) 
			#print build_pfam_motif_face(domain_dict[node.name])
			motif_face = SeqMotifFace(seq, build_pfam_motif_face(domain_dict[node.name]), scale_factor=1, seq_format = "compactseq")
			#motif_face.border.width = 10
			#motif_face.border.color = "pink"
			#motif_face.inner_border.width = 2
			#motif_face.inner_border.color = "orange"
			node.add_face(motif_face, column=1, position = "aligned") #display pfam domains
			
			node.img_style["size"]=0	
		elif node.is_root():
			node.dist = 0.01
			node.img_style["size"]=10

	ts = TreeStyle()	
	ts.show_leaf_name=False
	ts.scale=200#len(heat_figure.get_leaf_names())*10
	ts.layout_fn = pfam_layout
	
	if args.image:
		tree.render(args.image, w=800, units="mm", tree_style = ts)
	
	else:
		tree.show(tree_style = ts)


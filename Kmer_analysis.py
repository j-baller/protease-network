import os;
import sys;
import argparse
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import random
import math
import compare_peptides as cp #Local Module
 

parser = argparse.ArgumentParser(description=\
  'This program accepts fasta or fastq files and analyzes kmer distribution of a fixed size')

parser.add_argument(dest="in_filename",\
  help="specifies the input file (fasta or fastq)")
parser.add_argument('-b', '--background', action="store", dest="in_background", type=str, default="", required=False,\
  help="specifies the input file (fasta or fastq) for the background distribution") 
parser.add_argument('--fastq',action="store_true",dest="is_fastq",default=False ,required=False,\
  help="flag for fastq input, currently does not handle new lines within a fastq sequence (standard)")
parser.add_argument('--random',action="store_true",dest="randomize_flag",default=False ,required=False,\
  help="flag for randomization control")
parser.add_argument('-k',action="store",type=int,dest="kmer_size",required=False, default=4,\
  help="kmer size to use")
parser.add_argument('-c',action="store",type=float,dest="min_edge",required=False, default=.32,\
  help="minimal edge weight to keep")
parser.add_argument('-n',action="store",type=int,dest="min_node",required=False, default=1,\
  help="minimal node weight to keep")
parser.add_argument('-m',action="store",type=int,dest="max_node_count",required=False, default=5000,\
  help="maximum number of nodes afer enrichment step. Nodes are pruned based on the least kmer count first")
parser.add_argument('-d',action="store",type=int,dest="min_degree",required=False, default=3,\
  help="minimal node degree to keep")
parser.add_argument('-e',action="store",type=float,dest="min_enrichment",required=False, default=2,\
  help="Proportional fold enrichment of node weight (after pruning by minimum node weight. Only has an effect if a background distribution is provided.")
parser.add_argument("-l", "--left_flank", help="Amino Acid sequence of Left Flank", default='PQP')
parser.add_argument("-r", "--right_flank", help="Amino Acid sequence of Right Flank", default='PQP')
parser.add_argument("--replace_X", help="replace X's with flanking sequences",action="store_true",dest="replace_flag",default=False)
parser.add_argument('-x', "--enrichment_only", help="only perform enrichment checking, skip network filtering",action="store_true",dest="enrich_flag",default=False)
parser.add_argument('-z', "--control_zero_option", help="How to handle cases where the control pool has zeros when determining enrichment",required=False,dest="control_zero",default="average",choices=['average','avg','pseudo'])
parser.add_argument('-s', '--normalize_off',action="store_false",dest="norm_flag",default=True, \
  help="Shuts per comparision normalization off, returns a score in cale of actual amino acid matrix")
parser.add_argument('-a', '--align',action="store_true",dest="align_flag",default=False, \
  help="Return maximum score obtainable by sliding kmers relative to each other rather than assuming they are already aligned. Edges are padded with Xs")
cmd_args = parser.parse_args()

#### Development Code Snippets ############
#nx.connected_component_subgraphs(G) Can be used to get a list of subgraphs in graph.

#Get number of subgraphs
#[len(c) for c in net.connected_component_subgraphs(e) if len(c) > 10]

#Copy Network with Trimming:
# def trim_edges(g, weight=1):
        # g2=net.Graph()
        # for f, to, edata in g.edges(data=True):
                # if edata['weight'] > weight:
                        # g2.add_edge(f,to,edata)
        # return g2



##################################
#### Hardcoded Values ############ #Min/Max Score Should no longer be used, derived directly from loaded matrix
#max_score=cmd_args.kmer_size*17 #based on maximum positive in PAM250 matrix. This is the best possible score,
#min_score=cmd_args.kmer_size*-8
# rescaling is then, for score S, edge_weight =(S-min_score)/(max_score-min_score), making a score of 1 == perfect and a score of 0 equal to worst possible.
AA_poss = 20
##################################
input_name = cmd_args.in_filename.split(".")[0]

def replace_Xs(peptide,replace_flag):
	front_list = [i for i, ltr in enumerate(peptide) if ltr == 'x' and i < len(peptide)/2]
	back_list = [i for i, ltr in enumerate(peptide) if ltr == 'x' and i > len(peptide)/2]
	curr_str = [ltr for ltr in peptide if ltr != 'x']
	if replace_flag:
		return [cmd_args.left_flank[(len(cmd_args.left_flank)-len(front_list)):], "".join(curr_str), cmd_args.right_flank[:len(back_list)]]
	else:
		return [peptide[:len(front_list)], "".join(curr_str), peptide[(len(peptide)-len(back_list)):]]

def load_sequences(curr_filename):
	total_kmers = 0
	in_handle = open(curr_filename, 'r')
	master_dict = {}
	if cmd_args.is_fastq: #if is_fastq
		raise ValueError("Currently missing fastq module")	
		for line in in_handle:
			if line_count%4 == 1:
				curr_seq = line.rstrip().lower()
				f_str, mid_str, l_str = replace_Xs(curr_seq, cmd_args.replace_flag)
				if cmd_args.randomize_flag:
					mid_str = ''.join(random.sample(mid_str, len(mid_str)))
				curr_seq = f_str+mid_str+l_str
				curr_seq = curr_seq.lower()
				pos_zip = zip(range(0,len(curr_seq)-cmd_args.kmer_size+1), range(cmd_args.kmer_size, len(curr_seq)+1))
				for coord in pos_zip:
					reduc_denom = max(len(f_str)-coord[0],0) - max(len(f_str)-coord[1],0) - max(coord[0]- (len(curr_seq)-len(l_str)),0) + max(coord[1]- (len(curr_seq)-len(l_str)),0)
					master_dict[curr_seq[coord[0]:coord[1]]] = master_dict.get(curr_seq[coord[0]:coord[1]], 0) + 1/(AA_poss**reduc_denom)
					total_kmers += 1/(AA_poss**reduc_denom)
		in_handle.close()
		
	else: #is_fasta
		curr_seq = ""
		for line in in_handle:
			if line[0] == ">":
				if curr_seq != "":
					f_str, mid_str, l_str = replace_Xs(curr_seq, cmd_args.replace_flag)
					if cmd_args.randomize_flag:
						mid_str = ''.join(random.sample(mid_str, len(mid_str)))
					curr_seq = f_str+mid_str+l_str
					curr_seq = curr_seq.lower()
					pos_zip = zip(range(0,len(curr_seq)-cmd_args.kmer_size+1), range(cmd_args.kmer_size, len(curr_seq)+1))
					for coord in pos_zip:
						reduc_denom = max(len(f_str)-coord[0],0) - max(len(f_str)-coord[1],0) - max(coord[0]- (len(curr_seq)-len(l_str)),0) + max(coord[1]- (len(curr_seq)-len(l_str)),0)
						sav_kmer = curr_seq[coord[0]:coord[1]]
						if "_" in sav_kmer:
							print(line, file=sys.stderr)
							print(curr_seq, file=sys.stderr)
							raise ValueError()
						master_dict[sav_kmer] = master_dict.get(sav_kmer, 0) + 1/(AA_poss**reduc_denom)
						total_kmers += 1/(AA_poss**reduc_denom)
			else:
				curr_seq = curr_seq + line.rstrip().lower()
		curr_seq = line.rstrip().lower()
		f_str, mid_str, l_str = replace_Xs(curr_seq, cmd_args.replace_flag)
		if cmd_args.randomize_flag:
			mid_str = ''.join(random.sample(mid_str, len(mid_str)))
		curr_seq = f_str+mid_str+l_str
		curr_seq = curr_seq.lower()
		pos_zip = zip(range(0,len(curr_seq)-cmd_args.kmer_size+1), range(cmd_args.kmer_size, len(curr_seq)+1))
		for coord in pos_zip:
			reduc_denom = max(len(f_str)-coord[0],0) - max(len(f_str)-coord[1],0) - max(coord[0]- (len(curr_seq)-len(l_str)),0) + max(coord[1]- (len(curr_seq)-len(l_str)),0)
			sav_kmer = curr_seq[coord[0]:coord[1]]
			if "_" in sav_kmer:
				print(line, file=sys.stderr)
				print(curr_seq, file=sys.stderr)
				raise ValueError()
			master_dict[sav_kmer] = master_dict.get(sav_kmer, 0) + 1/(AA_poss**reduc_denom)
			total_kmers += 1/(AA_poss**reduc_denom)
		in_handle.close()
	return((master_dict, total_kmers))

def plot_network(t_graph, filename_detail, fig_x=16,fig_y=22):
	#Prepare plot of kmer network
	plt.clf()
	plt.figure(figsize=(fig_x,fig_y))
	pos = nx.spring_layout(m_graph) 

	#get max node weight
	node_weights = max(list(nx.get_node_attributes(t_graph,'weight').values()))
	
	#Draw Graph Elements
	nx.draw_networkx_edges(t_graph, pos,width=[math.log(edata['weight']-cmd_args.min_edge+.001)*2.5 for u,v,edata in t_graph.edges(data=True)],alpha=.3,edge_color='m')
	nx.draw_networkx_nodes(t_graph, pos,with_labels=True,node_size=[ndata['weight']/node_weights*100 for n,ndata in t_graph.nodes(data=True)])
	nx.draw_networkx_labels(t_graph,pos,fontsize=14)
	
	plt.savefig("".join(["network",filename_detail])+".png",format="png")


if cmd_args.randomize_flag:
	output_id=input_name+"_"+str(cmd_args.kmer_size)+"kmer_size"+str(cmd_args.min_edge)+"min_edge"+str(cmd_args.min_node)+"min_node"+str(cmd_args.min_degree)+"min_degree"+str(cmd_args.min_enrichment)+"fold_enrichment"+str(random.randint(0, 10000000))+"_rand"
else:
	output_id=input_name+"_"+str(cmd_args.kmer_size)+"kmer_size"+str(cmd_args.min_edge)+"min_edge"+str(cmd_args.min_node)+"min_node"+str(cmd_args.min_degree)+"min_degree"+str(cmd_args.min_enrichment)+"fold_enrichment"+str(random.randint(0, 10000000))+"_real"
print("Output Name: "+output_id)

master_dict = {}
m_graph = nx.Graph()

#Define Distance Matrix for comparing strings.
#score_mat = PAM250_special()
pd = cp.scoring_distance(cmd_args.norm_flag)

print("Input path resolved to:", os.path.abspath(cmd_args.in_filename))
main_dict, total_main = load_sequences(os.path.abspath(cmd_args.in_filename))

#Convert Primary Dictionary to Node_set, filter out kmers with fewer than a minimum occurance cutoff
node_set = [(i,{'weight':j}) for i,j in main_dict.items() if j >= cmd_args.min_node]

#Create Nodes in Network with a weight characteristic
m_graph.add_nodes_from(node_set)
node_weights = nx.get_node_attributes(m_graph,'weight').values()

#Print plot showing kmer_distribution
fig=plt.hist([math.log10(x) for x in list(node_weights)])
#fig = plt.gcf()
plt.savefig("kmer_frequency_histo"+output_id+".png",format="png")

#If background is defined load the data, check enrichment levels
#Remove all nodes that have an enrichment lower than allowable enrichment.
if cmd_args.in_background != "":
	print("Input path resolved to:", os.path.abspath(cmd_args.in_background))
	background_dict, total_back = load_sequences(os.path.abspath(cmd_args.in_background))

	back_set = [(i,{'weight':j}) for i,j in background_dict.items()]
	b_graph = nx.Graph()
	b_graph.add_nodes_from(back_set)
	
	remove_list = []
	enrich_attribute_dict = {}
	for n in m_graph:
		if cmd_args.control_zero in ['average','avg']:
			if n in b_graph:		
				#Filter on enrichment of each kmer node. Enrichment is the proportional occurance in the main set, divided by the proportional occurance in the background set.
				temp_enrich = (m_graph.node[n]["weight"]/total_main)/(b_graph.node[n]["weight"]/total_back)
				if temp_enrich < cmd_args.min_enrichment:
					remove_list.append(n)
				else:
					enrich_attribute_dict[n] = temp_enrich
			else:
				#Set background to average Kmer proportion
				temp_enrich = (m_graph.node[n]["weight"]/total_main)/(1.0/len(back_set))
				if temp_enrich < cmd_args.min_enrichment:
					remove_list.append(n)
				else:
					enrich_attribute_dict[n] = temp_enrich
		elif cmd_args.control_zero in ['pseudo']:
			#Filter on enrichment of each kmer node. Enrichment is the proportional occurance in the main set, divided by the proportional occurance in the background set. Adding Pseudocount of 1 to the background and foreground
			temp_enrich = ((m_graph.node.get(n,{'weight':0})["weight"]+1)/total_main)/((b_graph.node.get(n,{'weight':0})["weight"]+1)/total_back)
			if temp_enrich < cmd_args.min_enrichment:
				remove_list.append(n)
			else:
				enrich_attribute_dict[n] = temp_enrich
			
	m_graph.remove_nodes_from(remove_list)
	#Add enrich attribute to nodes
	nx.set_node_attributes(m_graph, 'enrich', enrich_attribute_dict)			
del back_set			
	

past_set = []



print("Total Number of filtered Nodes", str(len(m_graph.node)))

#Compute distances between nodes
for i in m_graph.nodes_iter(): 
	for j in past_set:
		#print(i,j,[score_mat[i[c].upper(),j[c].upper()] for c in range(cmd_args.kmer_size)])
		#Distance calculation between nodes
		try:
			if cmd_args.align_flag:
				edge_weight = pd.compare_peptides_mp(i,j)
			else:
				edge_weight = pd.compare_peptides_sp(i,j)
		except KeyError:
			print(str(i)+" "+str(j), file=sys.stderr)
			print(list(m_graph.nodes_iter()), file=sys.stderr)
			raise
		#edge_weight = (sum((score_mat[i[c].upper(),j[c].upper()] for c in range(cmd_args.kmer_size)))-min_score)/(max_score-min_score) #Old edge weight code, now in external module
		if edge_weight >= cmd_args.min_edge:
			m_graph.add_edge(i,j,weight=edge_weight) #Insert Edges if greater than minimum
		
	past_set.append(i)
	if len(past_set) % 1000 == 0:
		print(len(past_set))

#If a background is defined, print out the nodes that meet the enrichment criteria but are being cut due to degree.
del_out = open("kmer_deleted"+input_name+"_"+output_id+".txt", 'w')
node_weight = nx.get_node_attributes(m_graph, 'weight')
if cmd_args.in_background != "":
	for n,d in m_graph.degree_iter():
		if d < cmd_args.min_degree:
			for i in range(node_weight[n]):
				print(n,file=del_out)

#Delete nodes failing degree criteria
m_graph.remove_nodes_from([n for n,d in m_graph.degree_iter() if d < cmd_args.min_degree])

#Plot Graph ####
plot_network(m_graph, output_id)

#Output Kmers that pass all criteria, replicated for the number of times they occur (WebLogo input)
kmer_out = open("kmer_output"+input_name+"_"+output_id+".txt", 'w')
for n,ndata in m_graph.nodes(data=True):
	for i in range(ndata['weight']):
		print(n,file=kmer_out)
kmer_out.close()

#Output information about all Kmer nodes that pass criteria
kmer_out = open("node_info"+input_name+"_"+output_id+".txt", 'w')
for n,ndata in m_graph.nodes(data=True):
	print(n+"\t"+str(ndata['weight'])+"\t"+str(ndata['enrich'])+"\t"+str(m_graph.degree(n))+"\t"+str(m_graph.degree(n,weight="weight")),file=kmer_out)
kmer_out.close()

#print(m_graph.edges(data=True))

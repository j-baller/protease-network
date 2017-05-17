import os;
import argparse
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import random


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
parser.add_argument('-d',action="store",type=int,dest="min_degree",required=False, default=1,\
  help="minimal node degree to keep")
parser.add_argument('-e',action="store",type=float,dest="min_enrichment",required=False, default=2,\
  help="Proportional fold enrichment of node weight (after pruning by minimum node weight. Only has an effect if a background distribution is provided.")
parser.add_argument("-l", "--left_flank", help="Amino Acid sequence of Left Flank", default='PQP')
parser.add_argument("-r", "--right_flank", help="Amino Acid sequence of Right Flank", default='PQP')
parser.add_argument("--replace_X", help="replace X's with flanking sequences",action="store_true",dest="replace_flag",default=False)
cmd_args = parser.parse_args()


#### Hardcoded Values ############
max_score=cmd_args.kmer_size*17 #based on maximum positive in PAM250 matrix. This is the best possible score,
min_score=cmd_args.kmer_size*-8
# rescaling is then, for score S, edge_weight =(S-min_score)/(max_score-min_score), making a score of 1 == perfect and a score of 0 equal to worst possible.
##################################
def replace_Xs(peptide,replace_flag):
	front_list = [i for i, ltr in enumerate(peptide) if ltr == 'x' and i < len(peptide)/2]
	back_list = [i for i, ltr in enumerate(peptide) if ltr == 'x' and i > len(peptide)/2]
	curr_str = [ltr for ltr in peptide if ltr != 'x']
	if replace_flag:
		return [cmd_args.left_flank[(len(cmd_args.left_flank)-len(front_list)):], "".join(curr_str), cmd_args.right_flank[:len(back_list)]]
	else:
		return [peptide[:len(front_list)], "".join(curr_str), peptide[(len(peptide)-len(back_list)):]]

input_name = cmd_args.in_filename.split(".")[0]
	

class PAM250(object):
    """The PAM250 scoring matrix class."""

    def __init__(self):
        """Initialize the scoring matrix."""
        with open(os.path.join(os.path.dirname(__file__), 'PAM250wX.txt')) as input_data:
            items = [line.strip().split() for line in input_data.readlines()]
            self.scoring_matrix = {(item[0], item[1]):int(item[2]) for item in items}

    def __getitem__(self, pair):
        """Returns the score of the given pair of protein."""
        return self.scoring_matrix[pair[0], pair[1]]
        
class PAM250_special(object):
	"""The PAM250 scoring matrix class."""
	def __init__(self):
		"""Initialize the scoring matrix."""
		with open(os.path.join(os.path.dirname(__file__), 'PAM250wX.txt')) as input_data:
			items = [line.strip().split() for line in input_data.readlines()]
			self.scoring_matrix = {}
			for item in items:
				if item[0] in self.scoring_matrix:
					self.scoring_matrix[item[0]][item[1]] = int(item[2])
				else:
					self.scoring_matrix[item[0]] = {item[1]:int(item[2])}
	def __getitem__(self, pair):
		"""Returns the score of the given pair of protein."""
		char_max = max(list(self.scoring_matrix[pair[0]].values())+list(self.scoring_matrix[pair[1]].values())) 
		char_min = min(list(self.scoring_matrix[pair[0]].values())+list(self.scoring_matrix[pair[1]].values()))
		return ((self.scoring_matrix[pair[0]][pair[1]]-char_min)/(char_max-char_min))*25-8


if cmd_args.randomize_flag:
	output_id=input_name+"_"+str(cmd_args.kmer_size)+"kmer_size"+str(cmd_args.min_edge)+"min_edge"+str(cmd_args.min_node)+"min_node"+str(cmd_args.min_degree)+"min_degree"+str(cmd_args.min_enrichment)+"fold_enrichment"+str(random.randint(0, 10000000))+"_rand"
else:
	output_id=input_name+"_"+str(cmd_args.kmer_size)+"kmer_size"+str(cmd_args.min_edge)+"min_edge"+str(cmd_args.min_node)+"min_node"+str(cmd_args.min_degree)+"min_degree"+str(cmd_args.min_enrichment)+"fold_enrichment"+str(random.randint(0, 10000000))+"_real"
print("Output Name: "+output_id)

master_dict = {}
m_graph = nx.Graph()

def load_sequences(curr_filename):
	total_kmers = 0
	in_handle = open(curr_filename, 'r')
	master_dict = {}
	if cmd_args.is_fastq:
		for line in in_handle:
			if line_count%4 == 1:
				curr_seq = line.rstrip().lower()
				f_str, mid_str, l_str = replace_Xs(curr_seq, cmd_args.replace_flag)
				if cmd_args.randomize_flag:
					mid_str = ''.join(random.sample(mid_str, len(mid_str)))
				curr_seq = f_str+mid_str+l_str
				while len(curr_seq) >= cmd_args.kmer_size:
					master_dict[curr_seq[0:cmd_args.kmer_size]] = master_dict.get(curr_seq[0:cmd_args.kmer_size], 0) + 1
					curr_seq = curr_seq[1:]
					total_kmers +=1
		in_handle.close()
		
	else: #is_fasta
		curr_seq = ""
		for line in in_handle:
			if line[0] == ">":
				f_str, mid_str, l_str = replace_Xs(curr_seq, cmd_args.replace_flag)
				if cmd_args.randomize_flag:
					mid_str = ''.join(random.sample(mid_str, len(mid_str)))
				curr_seq = f_str+mid_str+l_str
				while len(curr_seq) >= cmd_args.kmer_size:
					master_dict[curr_seq[0:cmd_args.kmer_size]] = master_dict.get(curr_seq[0:cmd_args.kmer_size], 0) + 1
					curr_seq = curr_seq[1:]
					total_kmers +=1
				curr_seq = ""
			else:
				curr_seq = curr_seq + line.rstrip().lower()
		f_str, mid_str, l_str = replace_Xs(curr_seq, cmd_args.replace_flag)
		if cmd_args.randomize_flag:
			mid_str = ''.join(random.sample(mid_str, len(mid_str)))
		curr_seq = f_str+mid_str+l_str
		while len(curr_seq) >= cmd_args.kmer_size:
					master_dict[curr_seq[0:cmd_args.kmer_size]] = master_dict.get(curr_seq[0:cmd_args.kmer_size], 0) + 1
					curr_seq = curr_seq[1:]
					total_kmers +=1
		in_handle.close()
	return((master_dict, total_kmers))
	
print("Input path resolved to:", os.path.abspath(cmd_args.in_filename))
main_dict, total_main = load_sequences(os.path.abspath(cmd_args.in_filename))

if cmd_args.in_background != "":
	print("Input path resolved to:", os.path.abspath(cmd_args.in_background))
	background_dict, total_back = load_sequences(os.path.abspath(cmd_args.in_background))
	
score_mat = PAM250_special() ## Selection of distance metric
node_set = [(i,{'weight':j}) for i,j in main_dict.items() if j >= cmd_args.min_node]
m_graph.add_nodes_from(node_set)
node_weights = nx.get_node_attributes(m_graph,'weight').values()

#If background is defined, check enrichment levels
if cmd_args.in_background != "":
	back_set = [(i,{'weight':j}) for i,j in background_dict.items()]
	b_graph = nx.Graph()
	b_graph.add_nodes_from(back_set)
	
	remove_list = []
	for n in m_graph:
		if n in b_graph:
			#print((m_graph.node[n]["weight"]),(b_graph.node[n]["weight"]))
			if (m_graph.node[n]["weight"]/total_main)/(b_graph.node[n]["weight"]/total_back) < cmd_args.min_enrichment:
				remove_list.append(n)
		else: # if not in background, assume present at average rate
			#print(n,"not in background, substituting average occurance")
			if (m_graph.node[n]["weight"]/total_main)/(1.0/len(back_set)) < cmd_args.min_enrichment:
				remove_list.append(n)
	m_graph.remove_nodes_from(remove_list)
				
del back_set			
	

past_set = []
fig=plt.hist(list(nx.get_node_attributes(m_graph,'weight').values()))
#fig = plt.gcf()
plt.savefig("kmer_frequency_histo"+output_id+".png",format="png")

print("Total Number of Nodes", str(len(m_graph.node)))

for i in m_graph.nodes_iter(): #Insert Edges
	for j in past_set:
		#print(i,j,[score_mat[i[c].upper(),j[c].upper()] for c in range(cmd_args.kmer_size)])
		edge_weight = (sum((score_mat[i[c].upper(),j[c].upper()] for c in range(cmd_args.kmer_size)))-min_score)/(max_score-min_score)
		if edge_weight >= cmd_args.min_edge:
			m_graph.add_edge(i,j,weight=edge_weight)
		
	past_set.append(i)
	if len(past_set) % 1000 == 0:
		print(len(past_set))

m_graph.remove_nodes_from([n for n,d in m_graph.degree_iter() if d < cmd_args.min_degree])
node_weights = list(nx.get_node_attributes(b_graph,'weight').values())

plt.clf()
plt.figure(figsize=(16,22))
pos = nx.spring_layout(m_graph) 

nx.draw_networkx_edges(m_graph, pos,width=[edata['weight']*edata['weight']*10 for u,v,edata in m_graph.edges(data=True)],alpha=.3,edge_color='m')
nx.draw_networkx_nodes(m_graph, pos,with_labels=True,node_size=[ndata['weight']/max(node_weights)*100 for n,ndata in m_graph.nodes(data=True)])
nx.draw_networkx_labels(m_graph,pos,fontsize=14)
plt.savefig("".join(["network",output_id])+".png",format="png")

kmer_out = open("kmer_output"+input_name+"_"+output_id+".txt", 'w')
for n,ndata in m_graph.nodes(data=True):
	for i in range(ndata['weight']):
		print(n,file=kmer_out)

#print(m_graph.edges(data=True))

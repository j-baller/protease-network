import os
import argparse
from collections import Counter
import random
import statistics
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

##CMD parser options:  #####################################################
parser = argparse.ArgumentParser()
parser.add_argument("input_file", help="Input fastq file representing inframe peptides")
parser.add_argument("-p", "--cleavage_peptides", help="Fasta file, contains centered cleavage sites", default="")
parser.add_argument("-b", "--peptide_background", help="Fasta file, background sites for cleavage. Ignored if cleavage-peptides is blank", default="")
parser.add_argument("-o","--output_name", help="Output filename root",default="AA_output")
parser.add_argument("-v", "--verbose", help="increase output verbosity", action="count", default=0)
parser.add_argument("-d", "--drop_stop", help="Drop peptides containing a stop codon", action="store_true")
parser.add_argument("-m", "--min_correct", help="Minimum number of correct codons", type=int, default=0)
parser.add_argument("-l", "--left_flank", help="Amino Acid sequence of Left Flank", default='PQP')
parser.add_argument("-r", "--right_flank", help="Amino Acid sequence of Right Flank", default='PQP')
parser.add_argument("-c", "--PWM_score_cutoff", help="Determines minimal score for a site to be included as a cleavage site", default=-99999,type=float)
parser.add_argument("--PWM_prior", help="pseudocount to add to all aminoacids in all positions", default=.05,type=float)
cmd_args = parser.parse_args()
#############################################################################################

##Function to generate an AminoAcid Dictionary:  ###########################
def AA_lookup_gen(trans_type="default"):
	nts = ['t', 'c', 'a', 'g']
	codons = [a+b+c for a in nts for b in nts for c in nts]
	if trans_type=="default": 
		amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
	if trans_type=="seelig_artifical": 
		amino_acids = 'fFllSsssyY**cC*WlllLpppPHhqQRrrriIiMTtttnNKkssrrVvvvAaaadDEeGggg'
	codon_table = dict(zip(codons, amino_acids))
	return(codon_table)
amino_acid_list = ['a','r','n','d','c','e','q','g','h','i','l','k','m','f','p','s','t','w','y','v']
#############################################################################################


def Process_Handle_into_PWM(curr_handle): #Read in a file of aligned, equal length, proteins and generate a position weight matrix. This was for the planned protein data.
	init_flag=False
	count = cmd_args.PWM_prior*len(amino_acid_list) #To account for Pseudocount
	for line in curr_handle:				
		curr_seq = line.rstrip().lower()
		if curr_seq[0] != ">":
			count += 1
			if not init_flag:
				position_matrix = [dict((i, cmd_args.PWM_prior) for i in amino_acid_list) for j in range(len(curr_seq))]
				init_flag = True
			for i,curr_chr in enumerate(curr_seq):
				if curr_chr == 'x':
					for j in amino_acid_list:
						position_matrix[i][j] += 1/len(amino_acid_list)
				else:
					position_matrix[i][curr_seq[i]]  += 1
	return((position_matrix,count))
	
def Process_Handle_into_rand_PWM(background_list, num_seq): #Read in a file of aligned, equal length, proteins and generate a position weight matrix
	total_reps=10000
	position_matrix = [dict((i, [cmd_args.PWM_prior]*total_reps) for i in amino_acid_list) for j in range(len(background_list[0]))]
	for k in range(total_reps):
		curr_sample = random.sample(background_list, int(num_seq))
		for line in curr_sample:				
			for i,curr_chr in enumerate(line):
				if curr_chr == 'x':
					for j in amino_acid_list:
						position_matrix[i][j][k] =  position_matrix[i][j][k] + 1/len(amino_acid_list)
				else:
					position_matrix[i][curr_chr][k]  = position_matrix[i][curr_chr][k] + 1
	return(position_matrix)
		

#Load AA lookup table that uses lowercase to denote Amino Acids that should not exist in the pool.	
AA_lookup = AA_lookup_gen(trans_type="seelig_artifical")

#Define input data
print("Input path resolved to:", os.path.abspath(cmd_args.input_file))
in_handle = open(os.path.abspath(cmd_args.input_file), 'r')

#Read in fastq file with pre trimmed dna sequences. Sequences are pushed into "peptide_list"
line_count=0
peptide_list = []
nt_list = []
for line in in_handle:
	if line_count%4 == 1:
		curr_seq = line.rstrip().lower()
		if len(curr_seq)%3 != 0:
			raise RuntimeError("Current DNA sequence is not divisable my 3, File:",os.path.abspath(cmd_args.input_file),"Line:",line_count)
		out_peptide = ""
		out_nt = []
		while len(curr_seq) > 0:
			out_peptide += AA_lookup[curr_seq[0:3]]
			out_nt.append(curr_seq[0:3])
			curr_seq = curr_seq[3:]
		peptide_list.append(out_peptide)
		nt_list.append(out_nt)
	line_count+=1	
in_handle.close()


list_of_correct = []
#If option is set to drop stop codons, count the occurence of '*' in each string and purge from NT list and Pepetide lists those containing any.
if cmd_args.drop_stop:
	nt_list = [nt for peptide,nt in zip(peptide_list,nt_list) if peptide.count('*') == 0]
	peptide_list = [peptide for peptide in peptide_list if peptide.count('*') == 0]
#Count the number of 'correct' peptides. 
list_of_correct = (sum(1 for c in peptide if c.isupper()) for peptide in peptide_list)
error_counts = Counter(list_of_correct)

out_handle = open(cmd_args.output_name+".csv", "w")
for vals in error_counts.items():
	print(str(vals[0])+","+str(vals[1]),file=out_handle)
out_handle.close()


out_str = ""
# if cmd_args.cleavage_peptides != "": #If peptides are provided, subselect for cleavage sites
	# in_pep = open(os.path.abspath(cmd_args.cleavage_peptides), 'r')
	# Obs_matrix, num_seq = Process_Handle_into_PWM(in_pep)
	# out_str += "w_cleave"
	
# test_out = open("test.pwm.csv", "w")
# if cmd_args.cleavage_peptides != "" and cmd_args.peptide_background != "":
	# in_back = open(os.path.abspath(cmd_args.peptide_background), 'r')
	# background_list = []
	# for line in in_back:
		# if line[0] != ">":
			# curr_seq = line.rstrip().lower()
			# background_list.append(curr_seq)
	# Obs_matrix_rand = Process_Handle_into_rand_PWM(background_list, num_seq-cmd_args.PWM_prior*20)
	
	# print("\t".join(amino_acid_list), file=test_out) # test
	# for i in range(len(Obs_matrix)):
		# out_list = [] # test
		# out_list2 = [] # test
		# for j in amino_acid_list:
			# rand_mean = statistics.mean(Obs_matrix_rand[i][j])
			# out_list.append(str(rand_mean)) # test
			# out_list2.append(str(Obs_matrix[i][j])) # test
			# cert_adj = abs((sum(k >= Obs_matrix[i][j] for k in Obs_matrix_rand[i][j])/len(Obs_matrix_rand[i][j]) -.5)*2)
			# Obs_matrix[i][j] = math.log((Obs_matrix[i][j])/rand_mean,2)*cert_adj
		# print("\t".join(out_list), file=test_out) # test
		# print("\t".join(out_list2), file=test_out) # test
	# out_str += "_w_back"
# test_out.close()

# if cmd_args.cleavage_peptides != "":
	# PWM_out = open(cmd_args.output_name+out_str+".pwm.csv", "w")
	# for j in amino_acid_list:
		# out_list = [str(j)]
		# for i in range(len(Obs_matrix)):
			# out_list.append(str(Obs_matrix[i][j]))
		# print("\t".join(out_list), file=PWM_out)
	# PWM_out.close()

out_handle = open(cmd_args.output_name+out_str+".pt.fasta", "w")
list_of_correct = (sum(1 for c in peptide if c.isupper()) for peptide in peptide_list)

score_list = []
count=0
for peptide in peptide_list:
	corr = next(list_of_correct)
	if corr >= cmd_args.min_correct: #Only keep sequences with fewer errors than threshold			
		if cmd_args.cleavage_peptides != "":
			peptide = cmd_args.left_flank+peptide+cmd_args.right_flank
			max_vals = (0,-999999999)
			for i in range(len(peptide)-len(Obs_matrix)):
				curr_score = 0
				for j in range(len(Obs_matrix)):
					curr_score += Obs_matrix[j][peptide[i+j].lower()]
				if curr_score > max_vals[1]:
					max_vals = (i,curr_score)
			max_vals[0]
			peptide = 'X'*len(cmd_args.left_flank)+peptide[len(cmd_args.left_flank):(len(cmd_args.left_flank)+len(Obs_matrix))]+'X'*len(cmd_args.right_flank)
			peptide = peptide[max_vals[0]:(max_vals[0]+len(Obs_matrix))]
		if cmd_args.cleavage_peptides != "":
			score_list.append(max_vals[1]/len(Obs_matrix))
		if cmd_args.cleavage_peptides == "" or cmd_args.PWM_score_cutoff <= (max_vals[1])/len(Obs_matrix):
			peptide = 'X'*len(cmd_args.left_flank)+peptide+'X'*len(cmd_args.right_flank)
			print(">Peptide_"+str(count),file=out_handle)
			print(peptide,file=out_handle)
	count+=1
out_handle.close()

past_set = []

if cmd_args.cleavage_peptides != "":
	fig=plt.hist(score_list)
	#fig = plt.gcf()
	plt.savefig(cmd_args.output_name+"_PWM_max_score_hist.png",format="png")

#Not currently being filtered on PWM
out_handle = open(cmd_args.output_name+".filtered.nt.fasta", "w")
list_of_correct = (sum(1 for c in peptide if c.isupper()) for peptide in peptide_list)
count=0
for nt in nt_list:
	corr = next(list_of_correct)
	if corr >= cmd_args.min_correct:
		
		print(">Fragment_"+str(count),file=out_handle)
		print("".join(nt),file=out_handle)
	count+=1
out_handle.close()


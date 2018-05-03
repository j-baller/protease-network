import os
import sys
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import Counter

##CMD parser options:  #####################################################
parser = argparse.ArgumentParser()
parser.add_argument(dest="in_filename",\
  help="specifies the input file (fasta)")
parser.add_argument("-k","--query_kmers", nargs="+", help="Query sequences to be scanned for", required=True)
parser.add_argument("-o","--output_name", help="Output root, if blank, write to terminal",default="")
parser.add_argument("-l", "--left_flank", help="Amino Acid sequence of Left Flank", default='PQP')
parser.add_argument("-r", "--right_flank", help="Amino Acid sequence of Right Flank", default='PQP')
parser.add_argument("--dont_replace_X", help="replace X's with flanking sequences",action="store_false",dest="replace_flag",default=True)
parser.add_argument("-d",help="Allowed distance from center, 0 will return no result if the parity of the peptide and kmer are different", dest="dist_center", default=.5,type=float)

cmd_args = parser.parse_args()
###########################################################################################

def replace_Xs(peptide,replace_flag):
	front_list = [i for i, ltr in enumerate(peptide) if ltr == 'x' and i < len(peptide)/2]
	back_list = [i for i, ltr in enumerate(peptide) if ltr == 'x' and i > len(peptide)/2]
	curr_str = [ltr for ltr in peptide if ltr != 'x']
	if replace_flag:
		return [cmd_args.left_flank[(len(cmd_args.left_flank)-len(front_list)):], "".join(curr_str), cmd_args.right_flank[:len(back_list)]]
	else:
		return [peptide[:len(front_list)], "".join(curr_str), peptide[(len(peptide)-len(back_list)):]]

if cmd_args.output_name == "":
	out_handle = sys.stdout
else:
	out_handle = open(os.path.abspath(cmd_args.output_name+"full_seq.tab"), 'w')

def scan_sequences(curr_filename):

	in_handle = open(curr_filename, 'r')
	curr_seq = ""
	out_arr = [[] for i in cmd_args.query_kmers] 
	for line in in_handle:
		if line[0] == ">":
			#most lines
			f_str, mid_str, l_str = replace_Xs(curr_seq, cmd_args.replace_flag)
			curr_seq = f_str+mid_str+l_str
			seq_cent = len(curr_seq)/2
			for n,query_kmer in enumerate(cmd_args.query_kmers):
				q_cent = len(query_kmer)/2  
				if abs(curr_seq.lower().find(query_kmer.lower())+q_cent - seq_cent) <= cmd_args.dist_center:
					if curr_seq.strip() == "":
						print("Attempted pass of empty line passed to next stage")
						print("q_cent="+str(q_cent))
						print("seq_cent="+str(seq_cent))
						print("line="+line.strip())
					else:
						out_arr[n].append((line.strip(), curr_seq.lower()))
			curr_seq = ""
			hold_title = line
		else:
			curr_seq = curr_seq + line.rstrip().lower()
	#final lines
	f_str, mid_str, l_str = replace_Xs(curr_seq, cmd_args.replace_flag)
	curr_seq = f_str+mid_str+l_str
	if abs(curr_seq.lower().find(query_kmer.lower())+q_cent - seq_cent) <= cmd_args.dist_center:
		out_arr[n].append((line.strip(), curr_seq.lower()))
	in_handle.close()
	return(out_arr)
	
out_arrs = scan_sequences(os.path.abspath(cmd_args.in_filename))
print("Query\tSequence\tFrequency\tPass_Freq\tHydro\tPass_Hydro",file=out_handle)
for n,out_arr in enumerate(out_arrs):
	if len(out_arr) == 0:
		print(cmd_args.query_kmers[n]+"\t-\t0\t0\t-\t"+str(0),file=out_handle)
		continue
	seq_dict = {}
	for it in out_arr:
		seq_dict[it[1]] = seq_dict.get(it[1],0) + 1
	max_cnt = max(seq_dict.values())
	filtered_arr = [i for i in seq_dict.keys() if seq_dict[i] == max_cnt ]	

	#load external table
	table_handle = open(os.path.abspath("./Hydro_tab.txt"),'r')
	tab_dict = {}
	for line in table_handle:
		spl_arr = line.strip().split("\t")
		tab_dict[spl_arr[2].lower()]= float(spl_arr[3])

	hydro_arr = [(pep, sum(tab_dict[aa] for aa in pep)) for pep in filtered_arr]
	hydro_min = min([j for i,j in hydro_arr])
	hydro_arr = [i for i,j in hydro_arr if j==hydro_min] 

	print(cmd_args.query_kmers[n]+"\t"+hydro_arr[0]+"\t"+str(max_cnt)+"\t"+str(len(filtered_arr))+"\t"+str(hydro_min)+"\t"+str(len(hydro_arr)),file=out_handle)


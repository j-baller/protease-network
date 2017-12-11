import os
import sys
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

##CMD parser options:  #####################################################
parser = argparse.ArgumentParser()
parser.add_argument(dest="in_filename",\
  help="specifies the input file (fasta)")
parser.add_argument("-k","--query_kmer", help="Query sequence to be scanned for", required=True)
parser.add_argument("-o","--output_name", help="Output root, if blank, write to terminal",default="")
parser.add_argument("-l", "--left_flank", help="Amino Acid sequence of Left Flank", default='PQP')
parser.add_argument("-r", "--right_flank", help="Amino Acid sequence of Right Flank", default='PQP')
parser.add_argument("--replace_X", help="replace X's with flanking sequences",action="store_true",dest="replace_flag",default=False)
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

def scan_sequences(curr_filename):
	if cmd_args.output_name == "":
		out_handle = sys.stdout
	else:
		out_handle = open(os.path.abspath(cmd_args.output_name+"_"+cmd_args.query_kmer+".fa"), 'w')
	in_handle = open(curr_filename, 'r')
	curr_seq = ""
	out_arr = []
	for line in in_handle:
		if line[0] == ">":
			#most lines
			f_str, mid_str, l_str = replace_Xs(curr_seq, cmd_args.replace_flag)
			curr_seq = f_str+mid_str+l_str
			if cmd_args.query_kmer.lower() in curr_seq.lower():
				print(hold_title.strip(), file=out_handle)
				print(curr_seq.strip(), file=out_handle)
			curr_seq = ""
			hold_title = line
		else:
			curr_seq = curr_seq + line.rstrip().lower()
	#final lines
	f_str, mid_str, l_str = replace_Xs(curr_seq, cmd_args.replace_flag)
	curr_seq = f_str+mid_str+l_str
	if cmd_args.query_kmer.lower() in curr_seq:
		print(hold_title.strip(), file=out_handle)
		print(curr_seq.strip(), file=out_handle)
	in_handle.close()
	return(out_arr)
	
out_arr = scan_sequences(os.path.abspath(cmd_args.in_filename))

	


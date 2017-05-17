import os
import argparse

##CMD parser options:  #####################################################
parser = argparse.ArgumentParser()
parser.add_argument("input_file", help="Input csv file representing peptides")
cmd_args = parser.parse_args()
#############################################################################################

in_handle = open(cmd_args.input_file, 'r')
count=0
for line in in_handle:
	curr_list = line.rstrip().split(",")
	if count == 0:
		out_handles = []
		for i in curr_list:
			print(i)
			out_handles.append(open(i+".fasta", 'w'))
	else:
		for i,curr_str in enumerate(curr_list):
			if curr_str != "":
				print(">peptide_"+str(i)+str(count),file=out_handles[i])
				print(curr_str,file=out_handles[i])
	count +=1




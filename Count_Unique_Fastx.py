import os;
import argparse
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from collections import Counter

parser = argparse.ArgumentParser(description=\
  'This program accepts fasta or fastq files and analyzes kmer distribution of a fixed size')
parser.add_argument(dest="in_filenames",\
  help="specifies the input file (fasta or fastq)", nargs="+")
parser.add_argument('--fastq',action="store_true",dest="is_fastq",default=False ,required=False,\
  help="flag for fastq input, currently does not handle new lines within a fastq sequence (standard)")
parser.add_argument('-k',action="store",type=int,dest="kmer_size",required=False, default=14,\
  help="kmer size to use")
parser.add_argument("--replace_X", help="replace X's with flanking sequences",action="store_true",dest="replace_flag",default=False)
parser.add_argument('-b', '--background', action="store", dest="in_background", type=str, default="", required=False,\
  help="specifies the input file (fasta or fastq) for the background distribution") 
parser.add_argument("-c", help="The minimum number of times a sequence needs to show up to be output", type=int, dest="min_count", required=False, default=5)
parser.add_argument("-l", "--left_flank", help="Amino Acid sequence of Left Flank", default='')
parser.add_argument("-r", "--right_flank", help="Amino Acid sequence of Right Flank", default='')

cmd_args = parser.parse_args()

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
		for line in in_handle:
			if line_count%4 == 1:
				curr_seq = line.rstrip().lower()
				f_str, mid_str, l_str = replace_Xs(curr_seq, cmd_args.replace_flag)
				if cmd_args.randomize_flag:
					mid_str = ''.join(random.sample(mid_str, len(mid_str)))
				curr_seq = f_str+mid_str+l_str
				curr_seq = curr_seq.lower()
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
				curr_seq = f_str+mid_str+l_str
				curr_seq = curr_seq.lower()
				while len(curr_seq) >= cmd_args.kmer_size:
					master_dict[curr_seq[0:cmd_args.kmer_size]] = master_dict.get(curr_seq[0:cmd_args.kmer_size], 0) + 1
					curr_seq = curr_seq[1:]
					total_kmers +=1
				curr_seq = ""
			else:
				curr_seq = curr_seq + line.rstrip().lower()
			f_str, mid_str, l_str = replace_Xs(curr_seq, cmd_args.replace_flag)
			
			curr_seq = f_str+mid_str+l_str
			curr_seq = curr_seq.lower()
			while len(curr_seq) >= cmd_args.kmer_size:
						master_dict[curr_seq[0:cmd_args.kmer_size]] = master_dict.get(curr_seq[0:cmd_args.kmer_size], 0) + 1
						curr_seq = curr_seq[1:]
						total_kmers +=1
		in_handle.close()
	return((master_dict, total_kmers))

print("Background path resolved to:", os.path.abspath(cmd_args.in_background))
back_dict, total_back = load_sequences(os.path.abspath(cmd_args.in_background))
out_handle=open(cmd_args.in_background+"."+str(cmd_args.kmer_size)+".comparisons.csv", 'w')
out_handle_2 = open(cmd_args.in_background+"."+str(cmd_args.kmer_size)+".duplication.csv", 'w')
out_handle_3 = open(cmd_args.in_background+"."+str(cmd_args.kmer_size)+".all_seq.csv","w")
for key, value in back_dict.items():
	if cmd_args.min_count <= value:
		print(str(key)+","+str(value), file=out_handle_3)
out_handle_3.close()


count_dict = Counter(list(back_dict.values()))
for key, value in count_dict.items():
	print(str(key)+","+str(value), file=out_handle_2)
out_handle_2.close()


print("sample_name,background_uniq_kmers,sample_uniq_kmer,overlap_uniq_kmer",file=out_handle)
for in_filename in cmd_args.in_filenames:
	print("Input path resolved to:", os.path.abspath(in_filename))
	main_dict, total_main = load_sequences(os.path.abspath(in_filename))
	print(in_filename+","+str(len(back_dict))+","+str(len(main_dict))+","+str(sum([1 for x in main_dict if x in back_dict])), file=out_handle)
	out_handle_2 = open(in_filename+"."+str(cmd_args.kmer_size)+".duplication.csv", 'w')
	count_dict = Counter(list(main_dict.values()))
	for key, value in count_dict.items():
		print(str(key)+","+str(value), file=out_handle_2)
	out_handle_2.close()
	out_handle_3 = open(in_filename+"."+str(cmd_args.kmer_size)+".all_seq.csv","w")
	for key, value in main_dict.items():
        	if cmd_args.min_count <= value:
        	        print(str(key)+","+str(value), file=out_handle_3)
	out_handle_3.close()
#python3 Count_Unique_Fastx.py Fxa1R1.pep.pt.fasta Fxa20R1.pep.pt.fasta Fxa400R1.pep.pt.fasta ADAM171R1.pep.pt.fasta ADAM1750R1.pep.pt.fasta Library.pep.pt.fasta -b ResinLibrary.pep.pt.fasta



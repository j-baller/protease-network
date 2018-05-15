import sys
sys.path.insert(0, "/Users/Joshua Baller/Documents/Seelig/protease-sites/")
import compare_peptides as cp
from subprocess import call
import argparse

parser = argparse.ArgumentParser(description=\
  'This program accepts fasta or fastq files and analyzes kmer distribution of a fixed size')

parser.add_argument(dest="in_kmers",\
  help="specifies the input file (kmer_output)")
parser.add_argument('-s', action="store", dest="score_mat", type=str, default="", required=True, choices=['even','group1','group2'],\
  help="Choose between 'even', 'group1', 'group2'. Group1 and Group2 are two seperate divisions of the aminoacids (Group2 groups aromatics)") 
parser.add_argument('-r','--root',action="store",dest="root_path",default=False ,required=False,\
  help="root name for output")
parser.add_argument('-c',action="store",dest="cleave_path",default='' ,required=False,\
  help="Location of cleavage file if one is provided")
cmd_args = parser.parse_args()

if cmd_args.score_mat == 'even':
	sd = cp.scoring_distance(normalize=False,matrix_dict=cp.even_scoring)
if cmd_args.score_mat == 'group1':
	sd = cp.scoring_distance(normalize=False,matrix_dict=cp.group_scoring_v1)
if cmd_args.score_mat == 'group2':
	sd = cp.scoring_distance(normalize=False,matrix_dict=cp.group_scoring_v2)


kmer_dict = cp.read_kmer_out_to_dict(cmd_args.in_kmers)
if cmd_args.cleave_path != '':
	clv_mat = cp.Cleavage_matrix(cmd_args.cleave_path)
	clust = cp.Clustering(kmer_dict,sd,known_cleavage=clv_mat)
	clust.write_verbose_flat(cmd_args.root_path+"_"+cmd_args.score_mat+"_cleave_tree/")
else:
	clust = cp.Clustering(kmer_dict,sd)
	clust.write_verbose_flat(cmd_args.root_path+"_"+cmd_args.score_mat+"_tree/")


#clust.write_history(cmd_args.root_path+"_"+cmd_args.score_mat+".pickle")
del clust



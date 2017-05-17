import os
import random
import argparse


def mutation_float(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return x
    
parser = argparse.ArgumentParser()
parser.add_argument("output_name", help="Output filename root")
parser.add_argument("-m", "--mutation", help="Mutation Rate, fraction of DNA sites", default=0, type=mutation_float)
parser.add_argument("-c", "--count", help="Number of Sequences", default=100,type=int)
parser.add_argument("-l", "--length", help="Codons in each sequence", default=8, type=int)
cmd_args = parser.parse_args()

building_blocks = ['GCT','TGC','GAC','GAA','TTC','GGT','CAT','ATC','AAA','CTG', 'ATG','AAC', 'CCG', 'CAG', 'CGT', 'TCT', 'ACT', 'GTT', 'TAC']
nt_set = {'a','t','c','g'}

out_handle = open(cmd_args.output_name+"_random_"+str(cmd_args.length)+"_AA_"+str(cmd_args.mutation)+"rate.fastq", 'w')
for i in range(cmd_args.count):
	seq_list = list("".join(random.sample(building_blocks, cmd_args.length)))
	for j in range(len(seq_list)):
		if random.random() < cmd_args.mutation:
			seq_list[j] = random.sample(nt_set - {seq_list[j].lower()}, 1)[0]
	print("@Pseudo_Read_"+str(i),file=out_handle)
	print("".join(seq_list),file=out_handle)
	print("+",file=out_handle)
	print("".join(['I']*cmd_args.length*3),file=out_handle)
	
out_handle.close()
	
	

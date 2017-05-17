import os
import argparse

##CMD parser options:  #####################################################
parser = argparse.ArgumentParser()
parser.add_argument("input_fastq", help="Input fastq file representing inframe peptides")
parser.add_argument("-o","--output_name", help="Output filename root",default="AA_output")
parser.add_argument("-v", "--verbose", help="increase output verbosity", action="count", default=0)
parser.add_argument("-d", "--drop_stop", help="Drop peptides containing a stop codon", action="store_true")
parser.add_argument("-m", "--min_correct", help="Minimum number of correct codons", type=int, default=0)
cmd_args = parser.parse_args()
#############################################################################################

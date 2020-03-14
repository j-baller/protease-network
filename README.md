# protease-site-detection
## Background
These scripts provide the different steps for identifying overrepresented motifs within the provided sequences. Script are written for Python3.

## Primary Scripts
1. **ProcessingScript.py** - Accepts pretrimmed fastq containing only cleavage target sequences (nt). Outputs a fasta file with cleavage target sequences converted to AA and sequences not meeting quality minimums filtered out. Options available through --help allow the user to set filter conditions. Examples include removal of all sequences with a stop codon, and removal of sequences containing other improper codons. 
1. **Kmer_Analysis.py** - Accepts fasta/fastq sequence files as input, such as those output from the ProcessingScript.py. Outputs are diverse and represent different stages of processing. Options can be listed with --help. Key options include -k to indicate the kmer size to use for decomposing the input sequences. -c, -n, -m, -d and -e to filter the resultant kmers. This will return the primary kmers of interest. Also responsible for generating randomized datasets with --random flag.

## Helper Scripts
* **Extract_Fullseqs.py** - Pull all the full sequences that contain a specific kmer.
* **ImageComposite.py** - Accepts a directory path to a folder of WebLogo images. Makes a composite image from these individual images.

## Data Tables
* **Entropy_calc.txt** - Pre computed table of entropy values.
* **Hydro_tab.txt** - Hydrophobicity table for filtering peptides for mass spec analysis.
* **PAM250wX.txt** - PAM250 matrix for use as a distance matrix in the kmer analysis. Xs were added to allow for special treatment fixed sequences flanking the peptides.

## Function Library
* **compare_peptides.py** - Library of functions called by other scripts.

## Required Environment/Software
* Weblogo2 - Weblogo is used to visualize motifs
* Python3 - The provided scripts were developed and tested in the following python3 environment. The scripts are likely to work in other environments but were not tested for such.
  * Python 3.7.1
  * GCC 7.3.0, Anaconda on linux (CentOS 7)
  * Packages include: Matplotlib and NetworkX
* Shell scipts for testing are compatible with sh/bash
* System requirements are highly dependent on kmer count and complexity. Analysis from the associated paper used 4 cores, 200gb RAM and 72 hours max walltime to generate kmer networks 4 at a time 

## Outputs
To console:
* Input path resolved to: 
  * Resolved path to listed input file
* Total Seq:
  * Number of sequences in input file
* Total background seq: 
  * Number of sequences in the background seq file
* Total Nodes:
  * Number of unique kmers identified. Covers both input and background
* Total kmers:
  * Sum of all input kmer weights. This includes frequency of repeat kmers and downweighting due to fixed flanking sequences.
* Total background kmers:
  * Sum of all background kmer weights. This includes frequency of repeat kmers and downweighting due to fixed flanking sequences.
* First Pass Kmers:
  * Sum of all input kmer weights after filtering on Enrichment vs background and minimum frequency. This includes frequency of repeat kmers and downweighting due to fixed flanking sequences.
* First Pass Nodes:
  * Number of unique kmers identified after filtering on Enrichment vs background and minimum frequency. Covers both input and background
* Second Pass Kmers:
  * Sum of all input kmer weights after filtering on network degree (based on minimum edge weight). This includes frequency of repeat kmers and downweighting due to fixed flanking sequences.
* Second Pass Nodes:
  * Number of unique kmers identified after filtering on network degree (based on minimum edge weight). Covers both input and background
* Output Name
  * root used for naming all output files.
  
File outputs
Square bracket denote variable fields based on parameters and script.
* Kmer_Analysis.py
  * kmer_delete_[File1]_vs_[File2]_[kSize]kmer_size[edge_filter]min_edge[min_nodes]min_node[min_deg]min_degree[min_enrich]fold_enrichment[random_num]_[real or rand]_[dist_mat].txt
    * Kmers filtered out based on parameters supplied
  * kmer_frequency_histo_[File1]_vs_[File2]_[kSize]kmer_size[edge_filter]min_edge[min_nodes]min_node[min_deg]min_degree[min_enrich]fold_enrichment[random_num]_[real or rand]_[dist_mat].png
    * Frequency histogram showing the distribution of frequencies for the kmers.
  * kmer_output_[File1]_vs_[File2]_[kSize]kmer_size[edge_filter]min_edge[min_nodes]min_node[min_deg]min_degree[min_enrich]fold_enrichment[random_num]_[real or rand]_[dist_mat].txt
    * Kmers that passed all filtering based on parameters supplied
  * node_info_[File1]_vs_[File2]_[kSize]kmer_size[edge_filter]min_edge[min_nodes]min_node[min_deg]min_degree[min_enrich]fold_enrichment[random_num]_[real or rand]_[dist_mat].txt
    * More in depth information about each kmer that passed filtering. Columns are: Kmer sequence, Weighted frequency of occurence, Enrichment fraction, Node Degree, Weighted Node Degree
  * network_[File1]_vs_[File2]_[kSize]kmer_size[edge_filter]min_edge[min_nodes]min_node[min_deg]min_degree[min_enrich]fold_enrichment[random_num]_[real or rand]_[dist_mat].png
    * A plot of kmers in the network after filtering on all paramters. Line width denotes edge weight.
  * [File1]_vs_[File2]_[kSize]kmer_size[edge_filter]min_edge[min_nodes]min_node[min_deg]min_degree[min_enrich]fold_enrichment[random_num]_[real or rand]_[dist_mat].log
    * A copy of the run log for each sample

## License

Copyright 2020 Regents of the University of Minnesota

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


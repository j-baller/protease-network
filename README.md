# protease-site-detection
## Background
These scripts provide the different steps for identifying overrepresented motifs within the provided sequences. Script are written for Python3.

## Primary Scripts
1. **ProcessingScript.py** - Accepts pretrimmed fastq containing only cleavage target sequences (nt). Outputs a fasta file with cleavage target sequences converted to AA and sequences not meeting quality minimums filtered out. Options available through --help allow the user to set filter conditions. Examples include removal of all sequences with a stop codon, and removal of sequences containing other improper codons. 
1. **Kmer_Analysis.py** - Accepts fasta/fastq sequence files as input, such as those output from the ProcessingScript.py. Outputs are diverse and represent different stages of processing. Options can be listed with --help. Key options include -k to indicate the kmer size to use for decomposing the input sequences. -c, -n, -m, -d and -e to filter the resultant kmers. This will return the primary kmers of interest. 

## Helper Scripts
* **Extract_Fullseqs.py** - Pull all the full sequences that contain a specific kmer.
* **ImageComposite.py** - Accepts a directory path to a folder of WebLogo images. Makes a composite image from these individual images.

## Function Library
* **compare_peptides.py** - Library of functions called by other scripts.



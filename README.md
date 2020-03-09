# protease-site-detection
## Background
These scripts provide the different steps for identifying overrepresented motifs within the provided sequences. Script are written for Python3.

## Primary Scripts
1. **ProcessingScript.py** - Accepts pretrimmed fastq containing only cleavage target sequences (nt). Outputs a fasta file with cleavage target sequences converted to AA and sequences not meeting quality minimums filtered out. Options available through --help allow the user to set filter conditions. Examples include removal of all sequences with a stop codon, and removal of sequences containing other improper codons. 
1. **Kmer_Analysis.py** - Accepts fasta/fastq sequence files as input, such as those output from the ProcessingScript.py. Outputs are diverse and represent different stages of processing. Options can be listed with --help. Key options include -k to indicate the kmer size to use for decomposing the input sequences. -c, -n, -m, -d and -e to filter the resultant kmers. This will return the primary kmers of interest. Also responsible for generating randomized datasets with --random flag.

## Helper Scripts
* **Extract_Fullseqs.py** - Pull all the full sequences that contain a specific kmer.
* **ImageComposite.py** - Accepts a directory path to a folder of WebLogo images. Makes a composite image from these individual images.

## Function Library
* **compare_peptides.py** - Library of functions called by other scripts.

## Required Environment/Software
* Weblogo2 - Weblogo is used to visualize motifs
* Python3 - The provided scripts were developed and tested in a python3 environment.
  * Python 3.7.1
  * GCC 7.3.0, Anaconda on linux
  * Packages include: Matplotlib and NetworkX
* Shell scipts for testing are compatible with sh/bash
* System requirements are highly dependent on kmer count and complexity. Analysis from the associated paper used 4 cores, 200gb RAM and 72 hours max walltime to generate kmer networks 4 at a time 

## License

Copyright 2020 Regents of the University of Minnesota

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


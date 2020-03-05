#!/bin/bash -el
#PBS -l nodes=1:ppn=10,mem=200gb,walltime=72:00:00 -q ram256g
#PBS -m abe

module load python3

cd ~/Seelig

kmer_set()
{
	casefile=$1
	controlfile=$2
	k=$3
	n=$4
	e=$5
	d=$6
	c=$7
	python3 Kmer_analysis.py $casefile -b $controlfile -k $k -n $n -e $e -d $d -c $c -a -z "pseudo" --replace_X --debug -p 'even' 
	#python3 Kmer_analysis.py $casefile -b $controlfile -k $k -n $n -e $e -d $d -c $c --random --replace_X &
	#python3 Kmer_analysis.py $casefile -b $controlfile -k $k -n $n -e $e -d $d -c $c --random --replace_X &
	#python3 Kmer_analysis.py $casefile -b $controlfile -k $k -n $n -e $e -d $d -c $c --random --replace_X &
	#python3 Kmer_analysis.py $casefile -b $controlfile -k $k -n $n -e $e -d $d -c $c --random --replace_X &
	#python3 Kmer_analysis.py $casefile -b $controlfile -k $k -n $n -e $e -d $d -c $c --random --replace_X &
	#python3 Kmer_analysis.py $casefile -b $controlfile -k $k -n $n -e $e -d $d -c $c --random --replace_X &
	#python3 Kmer_analysis.py $casefile -b $controlfile -k $k -n $n -e $e -d $d -c $c --random --replace_X &
	#python3 Kmer_analysis.py $casefile -b $controlfile -k $k -n $n -e $e -d $d -c $c --random --replace_X &
	#python3 Kmer_analysis.py $casefile -b $controlfile -k $k -n $n -e $e -d $d -c $c --random --replace_X &
	#python3 Kmer_analysis.py $casefile -b $controlfile -k $k -n $n -e $e -d $d -c $c --random --replace_X &

}
# Input, Background, kmer_size, minimum_count, fold_enrichment, min_degree, similarity_cutoff 

kmer_set "Fxa20R1.pep.pt.fasta" "ResinLibrary.pep.pt.fasta" 6 7 6 4 0.55 &
kmer_set "Fxa1R1.pep.pt.fasta"  "ResinLibrary.pep.pt.fasta" 6 10 10 4 0.55 &
kmer_set "ADAM171R1.pep.pt.fasta"  "ResinLibrary.pep.pt.fasta" 6 7 6 4 0.55 &
kmer_set "ADAM1750R1.pep.pt.fasta"  "ResinLibrary.pep.pt.fasta" 6 7 6 4 0.55 &
wait
kmer_set "Fxa400R1.pep.pt.fasta"  "ResinLibrary.pep.pt.fasta" 6 7 6 4 0.55 &
kmer_set "SpeB1R1.pep.pt.fasta"  "ResinLibrary.pep.pt.fasta" 6 7 6 4 0.55 &
kmer_set "SpeB20R1.pep.pt.fasta"  "ResinLibrary.pep.pt.fasta" 6 7 6 4 0.55 &
kmer_set "SpeB400R1.pep.pt.fasta"  "ResinLibrary.pep.pt.fasta" 6 7 6 4 0.55 &
wait

#kmer_set "Fxa20R1.pep.pt.fasta" "Library.pep.pt.fasta" 6 10 8 4 0.55 &
#kmer_set "Fxa1R1.pep.pt.fasta"  "Library.pep.pt.fasta" 6 10 10 4 0.55 &
#kmer_set "ADAM171R1.pep.pt.fasta"  "Library.pep.pt.fasta" 6 10 8 4 0.55 &
#kmer_set "ADAM1750R1.pep.pt.fasta"  "Library.pep.pt.fasta" 6 10 6 4 0.55 &
#wait

#kmer_set "Fxa400R1.pep.pt.fasta"  "Library.pep.pt.fasta" 6 10 5 4 0.55 &
#kmer_set "SpeB1R1.pep.pt.fasta"  "Library.pep.pt.fasta" 6 10 10 4 0.55 &
#kmer_set "SpeB20R1.pep.pt.fasta"  "Library.pep.pt.fasta" 6 10 12 4 0.55 &
#kmer_set "SpeB400R1.pep.pt.fasta"  "Library.pep.pt.fasta" 6 10 6 4 0.55 &
#wait
#
#kmer_set "Fxa20R2.pep.pt.fasta" "ResinLibrary.pep.pt.fasta" 6 10 8 4 0.55 &
#kmer_set "Fxa1R2.pep.pt.fasta"  "ResinLibrary.pep.pt.fasta" 6 10 10 4 0.55 &
#kmer_set "ADAM171R2.pep.pt.fasta"  "ResinLibrary.pep.pt.fasta" 6 10 8 4 0.55 &
#kmer_set "ADAM1750R2.pep.pt.fasta"  "ResinLibrary.pep.pt.fasta" 6 10 6 4 0.55 &
#wait
#kmer_set "Fxa400R2.pep.pt.fasta"  "ResinLibrary.pep.pt.fasta" 6 10 5 4 0.55 &
#kmer_set "SpeB1R2.pep.pt.fasta"  "ResinLibrary.pep.pt.fasta" 6 10 10 4 0.55 &
#kmer_set "SpeB20R2.pep.pt.fasta"  "ResinLibrary.pep.pt.fasta" 6 10 12 4 0.55 &
#kmer_set "SpeB400R2.pep.pt.fasta"  "ResinLibrary.pep.pt.fasta" 6 10 6 4 0.55 &
#wait
#


#171R1.pep.pt.fasta                Fxa20R2.pep.pt.fasta              SpeB20R1.pep.pt.fasta
#A1750R1.pep.pt.fasta               Fxa400R1.pep.pt.fasta             SpeB20R1.pt.fasta
#AA_outputw_cleave_w_back.pt.fasta  Fxa400R2.pep.pt.fasta             SpeB20R2.pep.pt.fasta
#ADAM171R1.pep.pt.fasta             Library.pep.pt.fasta              SpeB20R2.pt.fasta
#ADAM171R2.pep.pt.fasta             Library.pt.fasta                  SpeB400R1.pep.pt.fasta
#ADAM1750R1.pep.pt.fasta            ResinLibrary.pep.pt.fasta         SpeB400R2.pep.pt.fasta
#ADAM1750R2.pep.pt.fasta            ResinLibrary.pt.fasta             SpeB400R2.pt.fasta
#Fxa1R1.pep.pt.fasta                Selectedw_cleave_w_back.pt.fasta  Unselectedw_cleave_w_back.pt.fasta
#Fxa1R2.pep.pt.fasta                SpeB1R1.pep.pt.fasta
#Fxa20R1.pep.pt.fasta               SpeB1R2.pep.pt.fasta

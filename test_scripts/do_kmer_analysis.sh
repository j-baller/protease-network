
module load python3

kmer_set()
{
	casefile=$1
	controlfile=$2
	k=$3
	n=$4
	e=$5
	d=$6
	c=$7
	python3 ../Kmer_analysis.py $casefile -b $controlfile -k $k -n $n -e $e -d $d -c $c -a -z "pseudo" --replace_X --debug -p "even" 
	python3 ../Kmer_analysis.py $casefile -b $controlfile -k $k -n $n -e $e -d $d -c $c --random --replace_X &
}
# Input, Background, kmer_size, minimum_count, fold_enrichment, min_degree, similarity_cutoff 

kmer_set "../test_data/SpeB20_samp.fasta" "../test_data/immobilized_lib_samp.fasta" 4 2 2 3 0.45 

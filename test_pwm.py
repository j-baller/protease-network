import sys
sys.path.insert(0, "/Users/Joshua Baller/Documents/Seelig/protease-sites/")
import compare_peptides as cp
sd = cp.scoring_distance()

pwm1 = cp.PWM('ACDR',sd)
(pwm1+pwm1)._curr_PWM

pwm2 = cp.PWM('CDRA',sd)
(pwm1+pwm2)._curr_PWM
(pwm1+pwm2+pwm2)._curr_PWM

pwm3 = cp.PWM('TTAC',sd,rep=5)
(pwm1+pwm2+pwm2+pwm3)._curr_PWM

#clust = cp.Clustering({'ACDR':3,'CDRA':5,'TTAC':4, 'TTTC':2, 'PQDR':1},sd)
#clust.write_alignment('/home/support/jballer/Seelig/test_logo_dat.txt')
#print(clust.get_pwm()._curr_PWM)

kmer_dict = cp.read_kmer_out_to_dict('/home/support/jballer/Seelig/kmer_outputADAM171R1_ADAM171R1_6kmer_size0.65min_edge10min_node4min_degree8.0fold_enrichment1295751_real.txt')
clust = cp.Clustering(kmer_dict,sd)
clust.write_alignment('/home/support/jballer/Seelig/test_logo_dat.txt')
#print(clust.get_pwm()._curr_PWM)

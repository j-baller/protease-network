import sys
sys.path.insert(0, "/Users/Joshua Baller/Documents/Seelig/protease-sites/")
import compare_peptides as cp
from subprocess import call

sd = cp.scoring_distance(normalize=True)

pwm1 = cp.PWM('ACDR',sd)
(pwm1+pwm1)._curr_PWM

pwm2 = cp.PWM('CDRA',sd)
(pwm1+pwm2)._curr_PWM
(pwm1+pwm2+pwm2)._curr_PWM

pwm3 = cp.PWM('TTAC',sd,rep=5)
(pwm1+pwm2+pwm2+pwm3)._curr_PWM

clust = cp.Clustering({'ACDR':3,'CDRA':5,'TTAC':4, 'TTTC':2, 'PQDR':1},sd)
clust.write_alignment('/home/support/jballer/Seelig/test_logo_dat_simple.txt')
print(clust.get_pwm()._curr_PWM)
call(['./WebLogo/weblogo/weblogo', '-Fpdf', '-slarge', '-Aprotein'],stdin=open('test_logo_dat_simple.txt'),stdout=open('test_logo.pdf','w'))

kmer_dict = cp.read_kmer_out_to_dict('/home/support/jballer/Seelig/kmer_outputADAM171R1_ADAM171R1_6kmer_size0.65min_edge10min_node4min_degree8.0fold_enrichment1295751_real.txt')
clust = cp.Clustering(kmer_dict,sd)
clust.write_verbose_flat('/home/support/jballer/Seelig/ADAM171_tree/')
del clust

#kmer_dict = cp.read_kmer_out_to_dict('/home/support/jballer/Seelig/kmer_outputFxa400R1_Fxa400R1_6kmer_size0.65min_edge10min_node4min_degree5.0fold_enrichment7851219_real.txt')
#clust = cp.Clustering(kmer_dict,sd)
#clust.write_alignment('/home/support/jballer/Seelig/test_logo_Fxa400_dat.txt')
#call(['./WebLogo/weblogo/weblogo', '-Fpdf', '-slarge', '-Aprotein'],stdin=open('test_logo_Fxa400_dat.txt'),stdout=open('Fxa400_logo.pdf','w'))
#del clust

#kmer_dict = cp.read_kmer_out_to_dict('/home/support/jballer/Seelig/kmer_outputSpeB1R1_SpeB1R1_6kmer_size0.65min_edge10min_node4min_degree10.0fold_enrichment3265074_real.txt')
#clust = cp.Clustering(kmer_dict,sd)
#clust.write_alignment('/home/support/jballer/Seelig/test_logo_SpeB1_dat.txt')
#call(['./WebLogo/weblogo/weblogo', '-Fpdf', '-slarge', '-Aprotein'],stdin=open('test_logo_SpeB1_dat.txt'),stdout=open('SpeB1_logo.pdf','w'))
#del clust

#kmer_dict = cp.read_kmer_out_to_dict('/home/support/jballer/Seelig/kmer_outputSpeB400R1_SpeB400R1_6kmer_size0.65min_edge10min_node4min_degree6.0fold_enrichment8762309_real.txt')
#clust = cp.Clustering(kmer_dict,sd)
#clust.write_alignment('/home/support/jballer/Seelig/test_logo_SpeB400_dat.txt')
#call(['./WebLogo/weblogo/weblogo', '-Fpdf', '-slarge', '-Aprotein'],stdin=open('test_logo_SpeB400_dat.txt'),stdout=open('SpeB400_logo.pdf','w'))

#del clust


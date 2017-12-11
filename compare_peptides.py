import os;
import math;
import itertools
from collections import Counter 
from operator import itemgetter
import copy

#For testing
#import sys
#sys.path.insert(0, "/Users/Joshua Baller/Documents/Seelig/protease-sites/")
#import compare_peptides as cp
#sd = cp.scoring_distance()
#cp.PWM('ACDR',sd)

def read_kmer_out_to_dict(filename):
	curr_handle = open(filename, 'r')
	out_dict = Counter()
	for line in curr_handle:
		line = line.rstrip()
		out_dict[line] += 1
	return(out_dict)
		
		
	

class Clustering:
	def __init__(self, kmer_dict,scoring_class):
		self.kmers = kmer_dict
		self.score = scoring_class
		self._full_PWM = self._build_hiercluster()
	def _build_hiercluster(self):
		item_list = [PWM(i,self.score,j) for i,j in self.kmers.items()] # Convert kmer input into a list of individual PWMs
		pair_iter =itertools.combinations(item_list,2) 
		score_list = [(i,j,i*j) for i,j in pair_iter]
		while len(item_list) > 1:
			print(len(item_list))
			max_item = max(score_list,key=itemgetter(2))
			score_list = [(i,j,k) for i,j,k in score_list if i not in (max_item[0], max_item[1]) and j not in (max_item[0], max_item[1])] # strip nodes being removed
			item_list = [i for i in item_list if i not in (max_item[0], max_item[1])]
			new_node = max_item[0]+max_item[1]
			score_list = score_list + [(new_node, i, new_node*i) for i in item_list]
			item_list.append(new_node)
		return(item_list[0])
		
	def get_pwm(self):
		return(self._full_PWM)
		
	def write_alignment(self,out_file):
		out_handle = open(out_file, 'w')
		PWM_elem = (t.elements() for t in self._full_PWM._curr_PWM)
		PWM_elem = ("".join(z) for z in zip(*PWM_elem))
		for line in PWM_elem:
			print(line, file=out_handle)
		
		
			
		
		
class PWM:
	def __init__(self,in_str,score_obj,rep=1):
		self._curr_PWM = [Counter({char:rep}) for char in in_str.upper()]
		self.score_obj = score_obj
		self._depth = rep
	
	def __copy__(self):
		new_pwm = PWM('',self.score_obj)
		new_pwm._depth = self._depth
		new_pwm._curr_PWM = copy.deepcopy(self._curr_PWM)
		return(new_pwm)	
	
	def length(self):
		return len(self._curr_PWM)
	
	def _compare_PWM_sp(self, pwm1, pwm2, phase=0):
		idx1 = list(range(0,len(pwm1)))
		idx2 = list(range(0,len(pwm2)))
		if phase > 0:
			idx2 = [-1]*phase + idx2
		elif phase < 0:
			idx1 = [-1]*abs(phase) + idx1
		if len(idx1) > len(idx2):
			idx2 = idx2 + [-1]*(len(idx1) - len(idx2))
		elif len(idx2) > len(idx1):
			idx1 = idx1 + [-1]*(len(idx2) - len(idx1))
		n_iter = len(idx1) #equal to len(pep2) due to 'X' padding
		
		#Still needs work below here
		align_score =0
		for p in range(n_iter):
			if idx1[p] == -1 and idx2[p] == -1: #This case should probably never happen
				align_score += self.score_obj['X','X']
				print("Double end gap occurred")
			elif idx1[p] == -1:
				align_score += sum([self.score_obj[AA,'X']*cnt for AA, cnt in pwm2[idx2[p]].items()])/sum(pwm2[idx2[p]].values())
			elif idx2[p] == -1:
				align_score += sum([self.score_obj[AA,'X']*cnt for AA, cnt in pwm1[idx1[p]].items()])/sum(pwm1[idx1[p]].values())
			else:
				align_score += sum((self.score_obj[tup1[0],tup2[0]]*tup1[1]*tup2[1] for tup1, tup2 in itertools.product(pwm1[idx1[p]].items(), pwm2[idx2[p]].items())))/(sum(pwm1[idx1[p]].values())*sum(pwm2[idx2[p]].values()))
		return((align_score,(idx1, idx2)))
		#edge_weight = (sum((self[pep1[c].upper(),pep2[c].upper()]-self.min_score)/(self.max_score-self.min_score) for c in range(n_iter)))/n_iter

	def _compare_PWM_mp(self, pwm1,pwm2):
		phases = len(pwm1)+len(pwm2)-1
		min_phase = -1*len(pwm2)+1
		out = []
		#Total phases are (n+m)-1, 
		for p in range(min_phase, phases+min_phase):
			out.append(self._compare_PWM_sp(pwm1,pwm2,phase=p))
		res_max = max(out,key=itemgetter(0))
		return((res_max[0],res_max[1]))
		
	def __mul__(self, other):
		if type(other) == str:
			return(self*PWM(other, self.score_obj))
		elif type(other) == PWM:
			if self.score_obj != other.score_obj:
				raise TypeError('Scoring of PWMs is only defined for PWMs using the same scoring object')
			m_score, m_idxs = self._compare_PWM_mp(self._curr_PWM, other._curr_PWM)
			return(m_score)
		else:
			print((self,other))
			raise TypeError('Scoring of PWMs is only defined for Strings and PWMs')
			
	def __rmul__(self, other):
		return(self.__mul__(other))
		
	def __add__(self, other):
		if type(other) == str:
			return(self+PWM(other, self.score_obj))
		elif type(other) == PWM:
			if self.score_obj != other.score_obj:
				raise TypeError('Addition of PWMs is only defined for PWMs using the same scoring object')
			m_score, m_idxs = self._compare_PWM_mp(self._curr_PWM, other._curr_PWM)
			m_idx1, m_idx2 = m_idxs
			out_PWM = copy.copy(self)
			for pos, idxs in enumerate(zip(m_idx1,m_idx2)):
				if idxs[0] == -1 and idxs[1] == -1:
					print("Double end gap occurred")
					out_PWM._curr_PWM.insert(pos, Counter({'-':(self._depth+other._depth)}))
				elif idxs[0] == -1:
					out_PWM._curr_PWM.insert(pos, Counter({'-':self._depth})+other._curr_PWM[idxs[1]])
				elif idxs[1] == -1:
					out_PWM._curr_PWM[pos] = self._curr_PWM[idxs[0]] + Counter({'-':other._depth})
				else:
					out_PWM._curr_PWM[pos] = self._curr_PWM[idxs[0]] + other._curr_PWM[idxs[1]]
			out_PWM._depth = self._depth+ other._depth
			return(out_PWM)
		else:
			raise TypeError('Addition of PWMs is only defined for Strings and PWMs')
	def __radd__(self,other):
		return(self.__add__(other))
		
import sys;
   
class scoring_distance(object):
	"""The PAM250 scoring matrix class."""
	def __init__(self, normalize = False, debug=False):
		"""Initialize the scoring matrix."""
		with open(os.path.join(os.path.dirname(__file__), 'PAM250wX.txt')) as input_data:
			items = [line.strip().split() for line in input_data.readlines()]
			self.scoring_matrix = {}
			for item in items:
				if item[0] in self.scoring_matrix:
					self.scoring_matrix[item[0]][item[1]] = int(item[2])
				else:
					self.scoring_matrix[item[0]] = {item[1]:int(item[2])}
		self.max_score = max(max(self.scoring_matrix.values(), key=lambda v: max(v.values())).values())
		self.min_score = min(min(self.scoring_matrix.values(), key=lambda v: min(v.values())).values())
		self.per_pos_norm = normalize
		self._debug=debug
	def __getitem__(self, pair):
		"""Returns the score of the given pair of AA."""
		if self.per_pos_norm==True:
			if self.debug:
				print(pair[0])
				print(pair[1])
			try:
				char_max = max(list(self.scoring_matrix[pair[0]].values())+list(self.scoring_matrix[pair[1]].values())) #Capture the Maximum score either AA could get
				char_min = min(list(self.scoring_matrix[pair[0]].values())+list(self.scoring_matrix[pair[1]].values())) #Capture the Minimum score either AA could get
			except KeyError:
				print("Key Error occured while attempting to normalize AA matrix, likely causes are invalid comparison matrix file or invalid peptide input", file=sys.stderr)
				print("First item is "+str(pair[0])+", second item is "+str(pair[0])+".", file=sys.stderr)
				raise 
			return ((self.scoring_matrix[pair[0]][pair[1]]-char_min)/(char_max-char_min))*(self.max_score-self.min_score)+self.min_score #Per position normalization. Norm is 0 to 1, *25-8 resets range to -8 to 17
		else:
			try:
				return self.scoring_matrix[pair[0]][pair[1]]
			except KeyError:
				print("Key Error occured while attempting to extract a particular pair from the AA matrix, likely causes are invalid comparison matrix file or invalid peptide input", file=sys.stderr)
				print("First item is "+str(pair[0])+", second item is "+str(pair[1])+".", file=sys.stderr)
				print("First level of scoring matrix", file=sys.stderr)
				print(self.scoring_matrix[pair[0]], file=sys.stderr)
				raise 
			
	def compare_peptides_sp(self,pep1, pep2, phase=0):
		if phase > 0:
			pep2 = '-'*phase + pep2
		elif phase < 0:
			pep1 = '-'*abs(phase) + pep1
		if len(pep1) > len(pep2):
			pep2 = pep2 + '-'*(len(pep1) - len(pep2))
		elif len(pep2) > len(pep1):
			pep1 = pep1 + '-'*(len(pep2) - len(pep1))
		n_iter = len(pep1) #equal to len(pep2) due to 'X' padding

		edge_weight = (sum((self[pep1[c].upper(),pep2[c].upper()]-self.min_score)/(self.max_score-self.min_score) for c in range(n_iter)))/n_iter
		if self._debug: print(pep1.upper(), pep2.upper(), n_iter)
		if self._debug: print(edge_weight)

		try:
			edge_weight = (sum((self[pep1[c].upper(),pep2[c].upper()]-self.min_score)/(self.max_score-self.min_score) for c in range(n_iter)))/n_iter
		except KeyError:
			print("Key Error occured while attempting to compute Single Phase edge weight.", file=sys.stderr)
			print("Peptides being analyzed on error are: "+str(pep1)+" "+str(pep2), file=sys.stderr)
			raise
		if self.debug: print(pep1.upper(), pep2.upper(), n_iter)
		if self.debug: print(edge_weight)
		return(edge_weight)
	
	def compare_peptides_mp(self,pep1, pep2):
		return(self._compare_peptide_core_mp(pep1,pep2)[0])
	
	def _compare_peptide_core_mp(self,pep1,pep2):
		phases = len(pep1)+len(pep2)-1
		min_phase = -1*len(pep2)+1
		out = []
		#Total phases are (n+m)-1, 
		for p in range(min_phase, phases+min_phase):
			out.append(self.compare_peptides_sp(pep1,pep2,phase=p))
		if self._debug: print(list(zip(range(min_phase, phases+min_phase), out)))
		return(max(out), [min_phase+i for i, j in enumerate(out) if j == max(out)])
			 
		
if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description=\
  'This program accepts two peptides as arguments and computes the distance between them')
	parser.add_argument(dest="peptide_1",\
  help="The first peptide to compare")
	parser.add_argument(dest="peptide_2",\
  help="The second peptide to compare")
	parser.add_argument('-n', '--normalize',action="store_true",dest="norm_flag",default=False, \
  help="Specifies if a per position normalization will be performed")
	parser.add_argument('-a', '--align',action="store_true",dest="align_flag",default=False, \
  help="Specifies if multiple alignments should be evaluated")
	parser.add_argument('-d', '--debug',action="store_true",dest="debug_flag",default=False, \
  help="Prints debug information")   
	cmd_args = parser.parse_args()
	
	pd = scoring_distance(cmd_args.norm_flag, cmd_args._debug_flag)

	if cmd_args.align_flag:
		print(pd.compare_peptides_mp(cmd_args.peptide_1,cmd_args.peptide_2))
	else:
		print(pd.compare_peptides_sp(cmd_args.peptide_1,cmd_args.peptide_2))

import os;
import math;
import itertools
from collections import Counter 
from operator import itemgetter
import copy
import pickle
import errno
from subprocess import call

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
		item_list = [PWM(i,self.score,j,history=True) for i,j in self.kmers.items()] # Convert kmer input into a list of individual PWMs
		pair_iter =itertools.combinations(item_list,2) # Generate all of the combinations of Kmer comparisions
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
		self.get_pwm().write_alignment(out_file)
	
	def write_history(self,out_file):
		pickle.dump(self._full_PWM._hist_list, open(out_file, 'wb'))
	
	def write_verbose_dirs(self,out_dir):
		if not self._full_PWM._history:
			raise TypeError("Clustering.write_verbose requires that the underlying PWM retained history in order to produce meaningful results")
		hist_list = self._full_PWM._hist_list
		def try_makedir(mk_loc):
			try:
				os.makedirs(mk_loc)
			except OSError as e:
				if e.errno != errno.EEXIST:
					raise
		try_makedir(out_dir)
		try_makedir(out_dir+"/align_root")
		def outer_rec(br,base_dir):
			if type(br) == list and len(br) == 4:
				try_makedir(base_dir+"/0")
				try_makedir(base_dir+"/1")
				curr_join = outer_rec(br[0],base_dir+"/0") + outer_rec(br[1],base_dir+"/1")
				curr_join.write_logo(base_dir+"/partial_logo")
				return(curr_join)
			else:
				br.write_logo(base_dir+"/partial_logo")
				return(br)		
		outer_rec(hist_list,out_dir+"/align_root")
	
	def write_verbose_flat(self,out_dir,leaves=False):
		if not self._full_PWM._history:
			raise TypeError("Clustering.write_verbose requires that the underlying PWM retained history in order to produce meaningful results")
		hist_list = self._full_PWM._hist_list
		def try_makedir(mk_loc):
			try:
				os.makedirs(mk_loc)
			except OSError as e:
				if e.errno != errno.EEXIST:
					raise
		try_makedir(out_dir)
		def outer_rec(br,base_dir,idx=0):
			if type(br) == list and len(br) == 4:
				curr_join = outer_rec(br[0],base_dir+"_0",idx+br[3])+outer_rec(br[1],base_dir+"_1",idx+br[3])
				curr_join.write_logo(base_dir+"_")
				return(curr_join)
			else:
				if leaves:
					br.write_logo(base_dir+"_")
				return(br)		
		outer_rec(hist_list,out_dir+"/logo_out")
		
class Empirical_Error:
	def __init__(self):
		self._uniqueAA=20
		ent_handle = open("/home/support/jballer/Seelig/Entropy_calc.txt", 'r')
		self.entropy_dict = {}
		for line in ent_handle:
			samp_size,samp_ent = line.split(" ")
			self.entropy_dict[int(samp_size)] = float(samp_ent)
	
	def get_error(self, N):
		if N <= 100:
				E = self.entropy_dict[N]
				return(E)
		else:
			Hg = -1*math.log2(1/self._uniqueAA) #Assumes equal occurance of the 20 AAs, otherwise Sum p*log2(p) over each aminoacid (p is proportion)
			AE = Hg - (self._uniqueAA - 1)/(2*math.log(2)*N)
			return(AE)

	
class PWM:
	err_obj = Empirical_Error()
	def __init__(self,in_str,score_obj,rep=1,history=False):
		if len(in_str) > 0:
			self._curr_PWM = [Counter({char:rep}) for char in in_str.upper()] # Getting null strings, need to investigate
			self._calc_entropy()
		else:
			self._entropy=0
		self.score_obj = score_obj
		self._depth = rep
		self._history = history
		self._hist_list = self
	
		
	def _calc_entropy(self):
		accum = 0
		for c in self._curr_PWM:
			total_samples = sum((aa[1] for aa in c.items() if aa[0] not in ('-','X')))
			accum = accum + self.__class__.err_obj.get_error(total_samples) + sum(((aa[1]/total_samples)*math.log2(aa[1]/total_samples) for aa in c.items() if aa[0] not in ('-','X')))
		if len(self._curr_PWM) == 0:
			print('NULL')
			self._entropy = 0
		else:
			self._entropy = accum/len(self._curr_PWM)
				
			
	
	def __copy__(self):
		new_pwm = PWM('',self.score_obj)
		new_pwm._depth = self._depth
		new_pwm._curr_PWM = copy.deepcopy(self._curr_PWM)
		new_pwm._history = self._history
		new_pwm._hist_list = self._hist_list
		new_pwm._calc_entropy()
		return(new_pwm)	
	
	def length(self):
		return len(self._curr_PWM)
	
	def write_alignment(self,out_file):
		if self._history:
			out_handle = open(out_file, 'w')
			def outer_rec(br,offset=0):
				if isinstance(br, PWM):
					PWM_elem1 = (t.elements() for t in br._curr_PWM)
					PWM_elem1 = ("".join(z) for z in zip(*PWM_elem1))
					if offset > 0:
						prefix = '-'*offset
						PWM_elem1 = (prefix+s for s in PWM_elem1)
					if br.length() + offset < self.length():
						suffix = '-'*(self.length() - br.length() - offset)
						PWM_elem1 = (s+suffix for s in PWM_elem1)
					for line in PWM_elem1:
						print(line, file=out_handle)
				else:
					if br[3] >= 0:
						outer_rec(br[0],offset+br[3])
						outer_rec(br[1],offset)
					else:
						outer_rec(br[0],offset)
						outer_rec(br[1],offset+abs(br[3]))
			outer_rec(self._hist_list)
			out_handle.close()
		else:
			out_handle = open(out_file, 'w')
			PWM_elem = (t.elements() for t in self._curr_PWM)
			PWM_elem = ("".join(z) for z in zip(*PWM_elem))
			for line in PWM_elem:
				print(line, file=out_handle)
				
	def write_logo(self, out_file_root, weblogo_exec='/panfs/roc/groups/2/support/jballer/Seelig/WebLogo/weblogo/weblogo'):
		l_filepath = out_file_root+"ent_"+str(self._entropy)+"_cnt_"+str(sum(self._curr_PWM[0].values()))
		self.write_alignment(l_filepath+".txt")
		call([weblogo_exec, '-Feps', '-slarge', '-Aprotein'],stdin=open(l_filepath+".txt"),stdout=open(l_filepath+".eps",'w'))
		
	
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
		return(self._score_and_align(other)[0])

	def _score_and_align(self,other):
		if type(other) == str:
			return(self*PWM(other, self.score_obj))
		elif type(other) == PWM:
			if self.score_obj != other.score_obj:
				raise TypeError('Scoring of PWMs is only defined for PWMs using the same scoring object')
			m_score, m_idxs = self._compare_PWM_mp(self._curr_PWM, other._curr_PWM)
			m_idx = len(list(itertools.takewhile(lambda x: x==-1, m_idxs[0])))-len(list(itertools.takewhile(lambda x: x==-1, m_idxs[1])))
			return([m_score,m_idx])
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
			if out_PWM._history:
				out_PWM._hist_list = [self._hist_list, other._hist_list] + self._score_and_align(other)
			else:
				out_PWM._hist_list = out_PWM
			out_PWM._calc_entropy()
			return(out_PWM)
		else:
			raise TypeError('Addition of PWMs is only defined for Strings and PWMs')
	def __radd__(self,other):
		return(self.__add__(other))
		
import sys;
   
group_scoring_v1={'groups':[['R','H','K'],['D','E'],['S','T','C','N','Q'], ['G','A','V','I','L','M'],['F'],['Y'],['W'],['P']], 'AA_X':-1,'AA_hyphen':-2,'X_X':0,'hyphen_hyphen':-1,'AA_match':2,'AA_within_grp':1, 'AA_between_grp':-2, 'hyphen_X':-1}
group_scoring_v2={'groups':[['R','H','K'],['D','E'],['S','T','C','N','Q'], ['G','A','V','I','L','M'],['F','Y','W'],['P']], 'AA_X':-1,'AA_hyphen':-2,'X_X':0,'hyphen_hyphen':-1,'AA_match':2,'AA_within_grp':1, 'AA_between_grp':-2, 'hyphen_X':-1}
even_scoring={'groups':[['R'],['H'],['K'],['D'],['E'],['S'],['T'],['C'],['N'],['Q'], ['G'],['A'],['V'],['I'],['L'],['M'],['P'],['F'],['Y'],['W']], 'AA_X':-1,'AA_hyphen':-2,'X_X':0,'hyphen_hyphen':-1,'AA_match':2,'AA_within_grp':1, 'AA_between_grp':-2, 'hyphen_X':-1}


class scoring_distance(object):
	"""The PAM250 scoring matrix class."""
	def __init__(self, normalize = False, debug=False, matrix_name='PAM250wX.txt',  matrix_dict=None ):
		"""Initialize the scoring matrix."""
		"""Scoring matrices should include all amino acids as well as 'X' and '-'. 'X' is used for network alignment, '-' is used for motif alignment and hierarchical clustering """
		if matrix_dict is None:	
			with open(os.path.join(os.path.dirname(__file__), matrix_name)) as input_data:
				items = [line.strip().split() for line in input_data.readlines()]
				self.scoring_matrix = {}
				for item in items:
					if item[0] in self.scoring_matrix:
						self.scoring_matrix[item[0]][item[1]] = int(item[2])
					else:
						self.scoring_matrix[item[0]] = {item[1]:int(item[2])}
		else:
			print("Loading Alternative Matrix")
			self.scoring_matrix = {}
			matrix_dict['groups'].extend([['-'],['X']])
			for i in matrix_dict['groups']:
				for j in matrix_dict['groups']:
					for ij in itertools.product(i,j):
						if i == j: # if amino acid group matches
							if ij[0] == ij[1]: #If actual amino acid matches
								if ij[0] == '-':
									curr_score = matrix_dict['hyphen_hyphen']
								elif ij[0] == 'X':
									curr_score = matrix_dict['X_X']
								else:
									curr_score = matrix_dict['AA_match']
							else: #if amino acid group matches, but not the same amino acid
								curr_score = matrix_dict['AA_within_grp']
						else: #if in different amino acid groups
							if ij[0] == '-' or ij[1] == '-': #if compared to hyphen
								if ij[0] == 'X' or ij[1] == 'X': #Unless you happen to be an 'X'
									curr_score = matrix_dict['hyphen_X']
								else:
									curr_score = matrix_dict['AA_hyphen']
							elif ij[0] == 'X' or ij[1] == 'X': #If compared to X (compared to hyphen coverd in previous case)
								curr_score = matrix_dict['AA_X']
							else: #All between group AAs handled here
								curr_score = matrix_dict['AA_between_grp']
						if ij[0] in self.scoring_matrix:
							self.scoring_matrix[ij[0]][ij[1]] = curr_score
						else:
							self.scoring_matrix[ij[0]] = {ij[1]:curr_score}
			#print(self.scoring_matrix)					
							
		self.max_score = max(max(self.scoring_matrix.values(), key=lambda v: max(v.values())).values())
		self.min_score = min(min(self.scoring_matrix.values(), key=lambda v: min(v.values())).values())
		self.per_pos_norm = normalize
		self._debug=debug
		
	def __getitem__(self, pair):
		"""Returns the score of the given pair of AA."""
		if self.per_pos_norm==True:
			if self._debug:
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
		if self._debug: print(pep1.upper(), pep2.upper(), n_iter)
		if self._debug: print(edge_weight)
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

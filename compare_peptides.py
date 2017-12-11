import os;
import math;
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
		self.debug=debug
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
			return self.scoring_matrix[pair[0]][pair[1]]
			
	def compare_peptides_sp(self,pep1, pep2, phase=0):
		if phase > 0:
			pep2 = 'X'*phase + pep2
		elif phase < 0:
			pep1 = 'X'*abs(phase) + pep1
		if len(pep1) > len(pep2):
			pep2 = pep2 + 'X'*(len(pep1) - len(pep2))
		elif len(pep2) > len(pep1):
			pep1 = pep1 + 'X'*(len(pep2) - len(pep1))
		n_iter = len(pep1) #equal to len(pep2) due to 'X' padding
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
		phases = len(pep1)+len(pep2)-1
		min_phase = -1*len(pep2)+1
		out = []
		#Total phases are (n+m)-1, 
		for p in range(min_phase, phases+min_phase):
			out.append(self.compare_peptides_sp(pep1,pep2,phase=p))
		if self.debug: print(list(zip(range(min_phase, phases+min_phase), out)))
		return(max(out))
			 
		
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
	
	pd = scoring_distance(cmd_args.norm_flag, cmd_args.debug_flag)

	if cmd_args.align_flag:
		print(pd.compare_peptides_mp(cmd_args.peptide_1,cmd_args.peptide_2))
	else:
		print(pd.compare_peptides_sp(cmd_args.peptide_1,cmd_args.peptide_2))

# Combine multiple images into one.

#

import os
import glob
import math
import argparse
import warnings
import sys

from PIL import Image,ImageDraw,ImageFont

Image.MAX_IMAGE_PIXELS = 10000000000

parser = argparse.ArgumentParser(description=\
 'This folder generates a composite image from the collection of images in a folder')
parser.add_argument(dest="in_folder",\
  help="specifies the input folder location")
parser.add_argument('-d',action="store",type=int,dest="max_depth",required=False,\
  help="kmer size to use")
parser.add_argument('-x',action="store",type=int,dest="motif_x",required=False, default=600,\
  help="Sets the pixel width per image")
parser.add_argument('-y',action="store",type=int,dest="motif_y",required=False, default=300,\
  help="Sets the pixel height per image")
parser.add_argument('-e',action="store",type=float,dest="min_ent",required=False, default=0,\
  help="Minimum entropy per figure")
parser.add_argument('-c',action="store",type=int,dest="min_count",required=False, default=0,\
  help="Minimum count per figure")
parser.add_argument('-p', action="store", type=float,dest="min_prop",required=False, default=0,\
  help="Minimum proportion per figure")
parser.add_argument('-s', action="store", type=float,dest="min_sim",required=False,default=-999,\
  help="Minimum similarity per jump")
parser.add_argument('-o',action="store",dest="output_root",required=False, default="composite_image",\
  help="Root of output filename")
parser.add_argument('-l',action="store_true",dest="listmode_flag",required=False, default=False,\
  help="Outputs a list of images meeting criteria rather than a tree. Overrides -d flag, setting it to infinite")
parser.add_argument('-v',action="store_true",dest="verbose_flag",required=False, default=False,\
  help="Enables more extensive information at runtime")
parser.add_argument('--cleave-only', action="store_true", dest="cleave_only_flag", required=False, default=False,\
  help="Filter to only motifs with cleavage sites")
cmd_args = parser.parse_args()


size_mult=2.5
if cmd_args.max_depth is None:
	if cmd_args.listmode_flag:
		depth_max=9999
	else:
		depth_max=5
else:
	if cmd_args.listmode_flag:
		warnings.warn("Max_depth (-d) is overwritten by Listmode (-l)")
		depth_max=9999
	else:
		depth_max=cmd.max_depth

working_dir = os.path.expanduser(cmd_args.in_folder)

print(working_dir)

#Find all eps files in listed directory
file_list = glob.glob(working_dir+'/*.eps')
if len(file_list) == 0:
	sys.exit("No .eps Files in Listed Directory")
text_list = [os.path.splitext(thisFile)[0]+".txt" for thisFile in file_list]
param_list = [os.path.splitext(thisFile)[0]+"param.txt" for thisFile in file_list]
if len(file_list) != len(text_list) or len(file_list) != len(param_list):
	sys.exit("Missing a .txt or param.txt file for one of the provided .eps files")

tail_size = 1
#For each file identified, split parameters in filename, forms list of lists
mid_list = [os.path.basename(in_str).split('_')[2:] for in_str in file_list]

#for each file generate a list from 0 to the depth of the file in the tree
spl_list = [in_ls[0:(len(in_ls)-tail_size)] for in_ls in mid_list]
#For each file generate a list with the remainder of the information from the file, entropy, count and distance (pending) information.
inf_list = []
for thisFile in param_list:
	tempHandle = open(thisFile)
	temp_dict = {}
	for line in tempHandle:
		line = line.strip()
		line_spl = line.split(",")
		temp_dict[line_spl[0]] = line_spl[1]
	inf_list.append(temp_dict)
	tempHandle.close()

#inf_list = [in_ls[(len(in_ls)-10):len(in_ls)] for in_ls in mid_list]
m_dep = max([len(sub_list) for sub_list in spl_list])+1

#print(spl_list)
#print(inf_list)

m_dep = min(depth_max, m_dep)
print("Max Depth = "+str(m_dep))

unit_w = cmd_args.motif_x
unit_h = cmd_args.motif_y

w=int(unit_w*size_mult)
h=int(unit_h*size_mult)



cnt_lst= []
for ind, pos_list in enumerate(spl_list):
	cnt_lst.append(float(inf_list[ind]['count'].split('.')[0]))
max_cnt = max(cnt_lst)

full_x=math.pow(2,(m_dep-1))*unit_w*size_mult
max_y_obs=0
min_y_obs=1
min_x_obs=full_x
max_x_obs=0
good_count=0


#first pass
for ind, pos_list in enumerate(spl_list):
	ent = float(inf_list[ind]['entropy'])
	cnt = float(inf_list[ind]['count'])
	sim = float(inf_list[ind]['similariy']) #FIX TYPO
	wid = float(inf_list[ind]['width'])
	ks = float(inf_list[ind]['kmers'])
	clv = inf_list[ind]['cleavage']
	

	x_pos = .5
	dep = 1
	frac_diff = 1/math.pow(2,(dep))
	
	for i in pos_list:
		dep = dep+1
		frac_diff = 1/math.pow(2,(dep))
		if i == '0':
			x_pos=x_pos - frac_diff
		else:
			x_pos=x_pos + frac_diff
			
	if dep > m_dep:
		continue
		
	y_pos = (dep-1)/(m_dep)
	
	
	if ent < cmd_args.min_ent or cnt < cmd_args.min_count or sim < cmd_args.min_sim or cnt/max_cnt < cmd_args.min_prop:
		continue

	if clv is None and cmd_args.cleave_only_flag:
		continue
		
	if y_pos > max_y_obs:
		max_y_obs=y_pos
	if y_pos < min_y_obs:
		min_y_obs=y_pos
	if int((x_pos*full_x)+w/2) > max_x_obs:
		max_x_obs=int((x_pos*full_x)+w/2)
	if int((x_pos*full_x)-w/2) < min_x_obs:
		min_x_obs=int((x_pos*full_x)-w/2)
	
	good_count += 1

#determine palette size
if cmd_args.listmode_flag:
	mar_offset=200
	x_dim = int(w)
	y_dim = int(h*good_count)
else:
	mar_offset=0
	x_dim = int((max_x_obs- min_x_obs))
	y_dim = int(h*(m_dep)*(max_y_obs-min_y_obs))

print("xdim="+str(x_dim)+", ydim="+str(y_dim))
	
	
#Create Palette
im = Image.new("RGB", (x_dim+mar_offset, y_dim),color=(255,255,255))
draw = ImageDraw.Draw(im)
font = ImageFont.truetype("arial.ttf", 50)

	
#print pass
print("Number of images passing all filters is: "+str(good_count))
if(good_count == 0):
	sys.exit("No Good Motifs Left")

def round_to_n(x,n=4):
	if x == 0:
		return(0.0)
	try:
		return round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1))
	except ValueError:
		print(["function args:",x,n])
		raise ValueError
	except:
		raise

out_txt_handle =  open(working_dir+"/"+cmd_args.output_root+".composite.txt",'w')
image_count=0
for ind, pos_list in enumerate(spl_list):
	ent = float(inf_list[ind]['entropy'])
	cnt = float(inf_list[ind]['count'])
	sim = float(inf_list[ind]['similariy']) #FIX TYPO
	wid = float(inf_list[ind]['width'])
	ks = float(inf_list[ind]['kmers'])
	clv = inf_list[ind]['cleavage']
	
	
	x_pos = .5
	dep = 1
	frac_diff = 1/math.pow(2,(dep))
	
	for i in pos_list:
		dep = dep+1
		frac_diff = 1/math.pow(2,(dep))
		if i == '0':
			x_pos=x_pos - frac_diff
		else:
			x_pos=x_pos + frac_diff
			
	if dep > m_dep:
		continue
		
	y_pos = (dep-1)/(m_dep)
	
	if cmd_args.verbose_flag: print(x_pos,y_pos)
	if dep < m_dep and not cmd_args.listmode_flag:
		#print(((int(x_pos*x_dim-frac_diff*x_dim/2+unit_w/2), int(y_pos*y_dim+unit_h*3/2), int(x_pos*x_dim+frac_diff*x_dim/2-unit_w/2), int(y_pos*y_dim+unit_h*3/2)), (int(x_pos*x_dim), int(y_pos*y_dim), int(x_pos*x_dim), in.t(y_pos*y_dim+unit_h*3/2))))
		draw.line((int(x_pos*x_dim-frac_diff*x_dim/2+unit_w/2), int(y_pos*y_dim+unit_h*3/2), int(x_pos*x_dim+frac_diff*x_dim/2-unit_w/2), int(y_pos*y_dim+unit_h*3/2)),fill=0,width=5)
		draw.line((int(x_pos*x_dim), int(y_pos*y_dim), int(x_pos*x_dim), int(y_pos*y_dim+unit_h*3/2)),fill=0,width=5)
	
	if ent < cmd_args.min_ent or cnt < cmd_args.min_count or sim < cmd_args.min_sim or cnt/max_cnt < cmd_args.min_prop:
		continue

	if clv is None and cmd_args.cleave_only_flag:
		continue


	
	if image_count >= good_count: print("exceeded expected good motif count, some images may be missing")

	print(text_list[ind], file=out_txt_handle)
	uniq_kmer = set()
	for curr_line in open(text_list[ind],'r'):
		temp_line = curr_line.strip()
		uniq_kmer.add(temp_line)
		print(temp_line, file=out_txt_handle)

 
	if not cmd_args.verbose_flag: print(x_pos,y_pos)
	img = Image.open(file_list[ind])
	img.load(scale=size_mult)
	img = img.resize((w, h), Image.BILINEAR)

	if cmd_args.listmode_flag:
		im.paste(img, (int((.5*x_dim)-w/2+mar_offset), int((image_count/good_count)*y_dim)) )
		try:	
			draw.text((int((.5*x_dim)-w/2), int((image_count/good_count)*y_dim)), "E: "+str(round_to_n(ent))+"\nK: "+str(round_to_n(cnt))+"\nN: "+str(round_to_n(len(uniq_kmer)))+"\nF: "+str(round_to_n(cnt/max_cnt))+"\nS: "+str(round_to_n(sim))+"\nL:"+"".join(pos_list)+"\nC:"+str(clv), fill=(0), font=font)
		except ValueError:
			print(["external vars (x_dim, w, y_dim,h, image_count,good_count, ent, cnt, cnt/max_cnt, sim):", x_dim, w, y_dim,h,image_count, good_count, ent,cnt,cnt/max_cnt, sim])
			raise ValueError
	else:
		im.paste(img, (int((x_pos*x_dim)-w/2), int(y_pos*y_dim)) )

		draw.text((int((x_pos*x_dim)-w/2), int(y_pos*y_dim)), "E: "+str(round_to_n(ent))+"\nK: "+str(round_to_n(cnt))+"\nN: "+str(round_to_n(len(uniq_kmer)))+"\nF: "+str(round_to_n(cnt/max_cnt))+"\nS: "+str(round_to_n(sim))+"\nL:"+"".join(pos_list)+"\nC:"+str(clv), fill=(0), font=font) 
	image_count += 1

im.save(working_dir+"/"+cmd_args.output_root+".png", "PNG")

	
  #path = os.path.expanduser(file)
  #img = Image.open(path)
  #img.thumbnail((400, 400), Image.ANTIALIAS)
  #x = index // 2 * 400
  #y = index % 2 * 400
  #w, h = img.size
  #print('pos {0},{1} size {2},{3}'.format(x, y, w, h))
  #result.paste(img, (x, y, x + w, y + h))

#result.save(os.path.expanduser('~/image.jpg'))

import numpy as np
import copy
import argparse 
import sys

def match(str_):
	tmp = 0
	tmp+=str_.count(",")
	tmp+=str_.count(".")
	return tmp
def mismatch(str_):
	tmp=0
	tmp+=str_.count("A")
	tmp+=str_.count("C")
	tmp+=str_.count("G")
	tmp+=str_.count("T")
	tmp+=str_.count("N")
	tmp+=str_.count("a")
	tmp+=str_.count("c")
	tmp+=str_.count("g")
	tmp+=str_.count("t")
	tmp+=str_.count("n")
	return tmp

if __name__ == "__main__":
	num_cells = 0
	out_path = ""
	in_path = ""
	ms = 0
	nmc = 0

	ap = argparse.ArgumentParser()
	ap.add_argument("-names","--cell names", required=False, help="file containing the cell names")
	ap.add_argument("-n","--number of cells",required=False, help="the number of cells ")
	# ap.add_argument("-lst","--list of the cell names",required=False,
		# help="a file containing the names of the cells at each row, each row must be of the form \'cell CELL_ID\', for example \'cell 200\' denotes the cell with the ID 200")
	ap.add_argument("-out","--output file name",required=False, help="path of the summary mpileup file")
	ap.add_argument("-in","--input mpileup file",required=True, help="path to the input file")
	ap.add_argument("-ms","--minimum coverage",required=False, help="The minimum number of reads required to support the alternative, default value 3")
	ap.add_argument("-nmc","--minimum cells",required=False,
		help="The number of cells having a minimum number of alternative reads needed to consider a genomic location a mutation candidate locus, default value 2")
	args = vars(ap.parse_args())

	if args['number of cells']==None and args['cell names']==None:
		print("Please either provide the list of cell names or enter the number of cells \nUsage: python -in <input mpileup file> -n <number of cells> | -names <txt file containing the cell names>")
		sys.exit()
	elif args['number of cells']!=None:
		num_cells = int(args['number of cells'])
		cell_names = open(out_path+'cellNames.txt', 'w')
		for cell_indx in range(num_cells):
			cell_names.write('cell'+str(cell_indx)+'\n')
		cell_names.close()
	elif args['cell names']!=None:
		path_cell_names = args['cell names']
		num_cells = 0
		with open(path_cell_names, 'r') as cell_names:
			for line in cell_names:
				if len(line.strip())!=0:
					num_cells+=1


	if args['input mpileup file']!=None:
		in_path = args['input mpileup file']
	else:
		print('Please enter the path to the mpileup file \n Usage: python loci_filter.py -in <path to the mpileup file>')
		sys.exit()

	if args['output file name']!=None:
		out_path = args['output file name']
	else:
		out_path = "./candidate_loci.mpileup"
	if args['minimum coverage']!=None:
		ms = int(args['minimum coverage'])
	else:
		ms = 3
	if args['minimum cells']!=None:
		nmc = int(args['minimum cells'])
	else:
		nmc = 2

	print("# cells "+str(num_cells))

	ed_pileup = open(out_path, "w")

	print("Parse the mpileup file")	
	with open(in_path, "r") as infile:
		for line in infile:
			tmp_arr = []
			nmc_count = 0
			currecnt_pos = int(line.strip().split('\t')[1])
			split = line.strip().split('\t')[3:]
			for i in range(num_cells):
				if int(split[3*i])==0:
						# flag = True
					tmp_arr.append([0,0])
				else:
						# print split[3*i+1]
					tmp_arr.append([match(split[3*i+1]),mismatch(split[3*i+1])])
					if mismatch(split[3*i+1])>=ms:
						nmc_count+=1

			if nmc_count>=nmc:
				ed_pileup.write(line)
	ed_pileup.close()

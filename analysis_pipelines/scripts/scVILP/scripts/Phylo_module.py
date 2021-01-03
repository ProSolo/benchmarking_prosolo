import numpy as np
import copy
import sys
import os
from ete3 import Tree
import time 
import subprocess

def duplicates(lst, item):
	return [i for i, x in enumerate(lst) if x == item]

def gen_Newick(genotype, PerfectPhy_path, out_dir_path, names_):
	###########################################################################################
	###################### Remove the columns with only one state #############################
	###########################################################################################
	R = genotype
	R = np.array(R)
	l = R.shape[1]
	removal_R = []
	num_of_char = 2 ## This is a binary case
	for j in range(l):	
		result_count = len(set(R[:,j]))
		if result_count<num_of_char:
			removal_R.append(j)

	print(str(len(set(removal_R)))+" characters have been removed from the result")
	R = np.delete(R, removal_R, 1)
	#### The number of columns has changed
	n = len(R)
	l_R_removal = len(R[0])
	############################ Second Post Processing Step ###############################
	#### Build the input matrix of PerfectPhy without columns containing only 0s or 1s #####
	########################################################################################
	if not PerfectPhy_path.endswith("/"):
		PerfectPhy_path+="/"
	result_name = PerfectPhy_path+"src/tree_result.txt"
	if not out_dir_path.endswith("/"):
		out_dir_path+="/"
	result_newick = out_dir_path+"phylogeny.nex"
	cmd = 'rm '+result_name
	os.system(cmd)
	cmd ='rm '+result_newick
	os.system(cmd)	
	tree_result = open(result_name, "w")
	tree_result.write("0"+"\n")
	tree_result.write(str(n)+"\t"+str(l_R_removal)+"\n")


	####################### Third Post Processing Step ###############################
	################## Find the duplicated taxa in matrix/matrices ###################
	##################################################################################
	R_taxa = []
	for i in range(n):
		temp2 = ""
		for j in range(l_R_removal):
			temp2+=str(int(R[i][j]))
			
		R_taxa.append(temp2)

	R_taxa=np.array(R_taxa)

	##########################################################
	####### Prepare the input file/files of PerfectPhy #######
	##########################################################

	for i in range(n):

		for j in range(l_R_removal):
			if j==l_R_removal-1 and i!=n-1:
				tree_result.write(str(int(R[i][j]))+"\n")
			elif j==l_R_removal-1 and i==n-1:
				tree_result.write(str(int(R[i][j])))
			else:
				tree_result.write(str(int(R[i][j]))+"\t")
	tree_result.close()
	###################### Fourth Post Processing Step ############################
	######### Construct the perfect phylogeny using PerfectPhy ####################
	###############################################################################


	p_method = subprocess.Popen(PerfectPhy_path+'src/perfectphy -f '+result_name+' -newick', shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	data_method = ""
	for line in p_method.stdout.readlines():
		data_method=line
	retval = p_method.wait()
	data_method = data_method.replace(" ","")+";"

	p2_method = subprocess.Popen(PerfectPhy_path+'src/perfectphy -f '+result_name+' -unique', shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

	print("uniqueness of perfect phylogeny of the Inferred Tree:")
	for line in p2_method.stdout.readlines():
		print(line)
	##################### Build a Tree object from Newick(s) #######################
	tree_method = Tree(data_method, format=8)

	for leaf in tree_method:
		if leaf.name.replace("'","") in R_taxa:
			tmp_tuple = tuple(duplicates(R_taxa,leaf.name.replace("'","")))
			leaf.add_child(name=str(tmp_tuple[0]+1))
			for new_node in tmp_tuple[1:]:
				leaf.add_child(name=str(new_node+1))

	for leaf in tree_method:
		leaf.name = names_[int(leaf.name)-1]
	######## Name the internal nodes #############

	extra_R = n+1
	ntax = 0
	for node in tree_method.traverse("postorder"):
		ntax+=1
		if not node.is_leaf():
			node.name = "internal_"+str(extra_R)
			extra_R+=1

	data_method_write = tree_method.write(format=8)


	############################################################################
	#################### Generate the nexus format file ########################
	############################################################################	

	trees = open(result_newick, "w")
	trees.write("#NEXUS"+"\n")
	trees.write("begin taxa;"+"\n")
	trees.write("dimensions ntax="+str(ntax)+";"+"\n")
	trees.write("taxlabels"+"\n")
	for node in tree_method.traverse("postorder"):
		trees.write(node.name+"\n")
	trees.write(";"+"\n")
	trees.write("end;"+"\n")
	trees.write("Begin trees;"+"\n")	
	trees.write("Tree 'PAUP_Inferred' = [&U] ")
	trees.write(data_method_write+"\n")
	trees.write("END;")
	trees.close()
	return data_method_write
#!/usr/bin/env python
import numpy as np 
import argparse
from gurobipy import *
import copy
import sys
import os
import matplotlib.pyplot as plt
# plt.switch_backend('agg')
import seaborn as sns 
# sns.set_style("whitegrid")
import math
import scipy
import subprocess
from scipy.stats import *
from scipy.special import comb
from scipy.special import gammaln
from ete3 import Tree
import time 
from decimal import Decimal
import scipy.cluster.hierarchy as hier
import scipy.spatial.distance as dist
import parser
import VCF
import Phylo_module


def Binomial_pmf(k,n,p):
	''' calculates the pmf of binomial distribution '''
	k_decimal = Decimal(k)
	n_decimal = Decimal(n)
	p_decimal = Decimal(p)
	tmp = Decimal(gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1))+Decimal(k_decimal*p_decimal.ln()+(n_decimal-k_decimal)*Decimal(1-p_decimal).ln())
	return tmp.exp()

def check_InCompatibility(Matrix):

	'''Count the number of character pairs which violate
	 the infinite-sites assumption '''
	num_incompatibles = 0
	pairs = []
	Matrix = np.array(Matrix)
	n = Matrix.shape[0]
	l = Matrix.shape[1]
	B_prime = [[[0 for k in range(4)] for j in range(l)] for i in range(l)]
	for p in range(l):
		q = p+1
		while q<l:
			count01=0
			count10=0
			count11=0
			for cell in range(n):
				if count01+count10+count11==3:
					break  
				# if Matrix[cell][p]==0 and Matrix[cell][q]==0:
				# 	B_prime[p][q][0]=1
				if Matrix[cell][p]==0 and Matrix[cell][q]==1:
					B_prime[p][q][1]=1
				elif Matrix[cell][p]==1 and Matrix[cell][q]==0:
					B_prime[p][q][2]=1
				elif Matrix[cell][p]==1 and Matrix[cell][q]==1:
					B_prime[p][q][3]=1
			q+=1
	for p in range(l):
		q=p+1
		while q<l:
			s = sum(B_prime[p][q])
			if s==3:
				num_incompatibles+=1
				pairs.append((p,q))
			q+=1
	print(pairs)
	return num_incompatibles

def optimize(read_count_mat, fp, fn, missing_data_thr, K_vios, mu0, mu1):
    
    #########################################################################################################################
    ############ The arguments include the error rates, the read count matrix, and the threshold for missing data ###########
    #########################################################################################################################
    fp_decimal = Decimal(fp)
    fn_decimal = Decimal(fn)
    n = read_count_mat.shape[0]
    l = read_count_mat.shape[1]
    R = [[0 for i in range(l)] for j in range(n)]
    missing_data_threshold = missing_data_thr
    ######################################################################
    ######################### Build the model ############################
    ######################################################################
    model = Model("model")
    B = {}
    Y = []
    V = {}
    ########################################
    ### Add the variables to the model #####
    ########################################
    obj = LinExpr()
    vios = LinExpr()
    print("Add variables to the model")
    for i in range(n):
    	Y.append([])
    	for j in range(l):
    		Y[i].append(model.addVar(vtype=GRB.BINARY, name="Y[%d,%d]" % (i,j)))
    for p in range(l):
    	q=p+1
    	while q<l:
    		V["("+str(p)+","+str(q)+")"]=model.addVar(vtype=GRB.BINARY)
    		vios+=V["("+str(p)+","+str(q)+")"]
    		for k in range(3):
    			B["("+str(p)+","+str(q)+","+str(k+1)+")"]=model.addVar(vtype=GRB.BINARY, name="B["+str(p)+","+str(q)+","+str(k+1)+"]")
    		q+=1
    model.update()
    ######################################
    ### Add constraints to the model #####
    ######################################
    print("Add constraints to the model")
    for p in range(l):
    	q=p+1
    	while q<l:
    		
    		model.addConstr(V["("+str(p)+","+str(q)+")"]>=B["("+str(p)+","+str(q)+","+str(1)+")"]+B["("+str(p)+","+str(q)+","+str(2)+")"]+B["("+str(p)+","+str(q)+","+str(3)+")"]-2)
    		for taxon in range(n):
    			####### The constraints which control the B variables #######
    			model.addConstr(B["("+str(p)+","+str(q)+","+str(1)+")"]>=Y[taxon][q]-Y[taxon][p])
    			model.addConstr(B["("+str(p)+","+str(q)+","+str(2)+")"]>=Y[taxon][p]-Y[taxon][q])
    			model.addConstr(B["("+str(p)+","+str(q)+","+str(3)+")"]>=Y[taxon][p]+Y[taxon][q]-1)
    		q=q+1
    model.addConstr(vios<=K_vios)
    # mu0=1e-3
    # mu1=0.5
    #################################################################
    ################ Build the objective function ###################
    #################################################################
    print("Build the objective function")
    for i in range(n):
    	for j in range(l):
    		############# This line accounts for the missing data ##############
    		if read_count_mat[i][j][0]+read_count_mat[i][j][1]>=missing_data_threshold: 
    			r = int(read_count_mat[i][j][0])
    			v = int(read_count_mat[i][j][1])
    			AA = Binomial_pmf(v,r+v,mu0)
    			BB = Binomial_pmf(v,r+v,mu1)
    			#obj -= (Y[i][j])*np.float128(Decimal((fn_decimal/2)*AA+(1-fn_decimal/2)*BB).ln())
    			obj -= (Y[i][j])*np.float128(Decimal((fn_decimal)*AA+(1-fn_decimal)*BB).ln())
    			obj -= (1-Y[i][j])*np.float128(Decimal((1-fp_decimal)*AA+(fp_decimal)*BB).ln())
    		else:
    			pass
    model.update()
    ##################################################
    ########## Assign the objective function #########
    ##################################################
    print("Assign the objective function")
    model.setObjective(obj, GRB.MINIMIZE)
    #####################################################
    ######## Set the parameters of the model ############
    #####################################################
    #model.Params.timeLimit = 255000
    # model.Params.method=3
    #model.Params.Threads = 31
    #model.Params.ConcurrentMIP = 2
    # model.Params.MIPGap=0.01
    ########################################################
    ######### Optimize the model and report it #############
    ########################################################
    print("Optimize the model")
    model.optimize()
    print('IsMIP: %d' % model.IsMIP)
    if model.status == GRB.Status.INFEASIBLE:
    	print("The model is infeasible")
    print("Solved with MIPFocus: %d" % model.Params.MIPFocus)
    print("The noisy model has been optimized")
    print('Obj: %g' % model.objVal)
    #########################################################
    ##### Save the final array given by the ILP solver ######
    #########################################################
    for i in range(n):
    	for j in range(l):
    		R[i][j] = int(Y[i][j].x)
    gc.collect()
    return R

if __name__=="__main__":
	#########################################################################
	# Specify the path to the PerfectPhy directory here
	PerfectPhy_path_ = "./PerfectPhy"
	#########################################################################
	K_ = 0
	fn_given = 0.1
	fp_given = 1e-08
	missing_data_threshold = 10
	out_path = "./"
	data_path = ""
	cell_names_path = ""
	ap = argparse.ArgumentParser()
	ap.add_argument("-names","--cell names", required=True, help="file containing the cell names")
	ap.add_argument("-out","--output file name",required=False, help="path to the output directory")
	ap.add_argument("-in","--input mpileup file",required=True, help="path to the input file")
	ap.add_argument("-mdthr","--missing data threshold",required=False, help="minimum coverage for each ref-var pair, default value 10")
	ap.add_argument("-fp","--false positive rate",required=False, help="false positive error rate, default value 1e-08")
	ap.add_argument("-fn","--false negative rate",required=False, help="false negative error rate, default value 0.1")
	ap.add_argument("-vio","--maximum number of violations",required=False, help="maximum number of violations of infinite-sites assumption, default value 0")
	args = vars(ap.parse_args())

	if args['cell names']!=None:
		cell_names_path = args['cell names']
	else:
		print("Please enter the path to the cell names \nUsage: python scVILP_main.py -in <path to the mpileup file> -names <path to the list of cell names>")
		sys.exit()

	if args['output file name']!=None:
		out_path = args['output file name']
	if not out_path.endswith("/"):
		out_path+="/"
	if args['input mpileup file']!=None:
		data_path = args['input mpileup file']
	else:
		print("Please enter the path to the mpileup file\nUsage: python scVILP_main.py -in <path to the mpileup file> -names <path to the list of cell names>")
		sys.exit()
	if args['missing data threshold']!=None:
		missing_data_threshold = float(args['missing data threshold'])
	if args['false positive rate']!=None:
		fp_given = float(args['false positive rate'])
	if args['false negative rate']!=None:
		fn_given = float(args['false negative rate'])
	if args['maximum number of violations']!=None:
		K_ = int(args['maximum number of violations'])


	##############################################################################
	########################### Parse the mpileup file ###########################
	(read_counts, alts, refs, chroms, positions, names, depths) = parser.Parse(cell_names_path, data_path)
	n=read_counts.shape[0]
	l=read_counts.shape[1]

	print("# of taxa: %d" % n)
	print("# of mutations: %d" % l)

	print("false positive rate given: %f" %fp_given)
	print("false negative rate given: %f" %fn_given)



	mat_ = optimize(read_count_mat=read_counts,fp=fp_given,fn=fn_given,missing_data_thr=missing_data_threshold, K_vios=K_, mu0=1e-3, mu1=0.5)

	#############################################################################
	#################### Generate heatmap of the genotypes ######################
	#############################################################################
	mat_ = np.array(mat_)
	tmp_array = copy.copy(mat_)
	R1 = tmp_array.T
	distMatrix = dist.pdist(R1)
	distSquareMatrix = dist.squareform(distMatrix)
	linkageMatrix = hier.linkage(distMatrix,method='ward')
	dendro = hier.dendrogram(linkageMatrix)
	leaves1 = dendro['leaves']
	transformedData = R1[leaves1,:]


	R2=tmp_array
	distMatrix = dist.pdist(R2)
	distSquareMatrix = dist.squareform(distMatrix)
	linkageMatrix = hier.linkage(distMatrix,method='ward')
	dendro = hier.dendrogram(linkageMatrix)
	leaves2 = dendro['leaves']
	transformedData = transformedData[:,leaves2]
	##### leaves1 for the mutations sites
	##### leaves2 for the taxa
	fig_ = plt.figure(figsize=(6,6))
	ax_ = fig_.add_subplot(111)
	cax_ = ax_.matshow(transformedData,cmap='Blues',aspect="auto")
	ax_.set_ylabel('Cells')
	ax_.set_xlabel('Genomic Positions')


	# fig_.colorbar(cax_)
	fig_.savefig(out_path+"hier_clust_heatmap.png", dpi=1200)

	###########################################################################
	######################## Generate the VCF output ##########################
	###########################################################################
	VCF.gen_VCF(out_dir=out_path, genotype_mat=mat_, read_count_mat_=read_counts, chrs=chroms, posits=positions, alt_counts=alts, rfs=refs, ids=names, dps=depths)

	###################################################################################
	######################## Generate Perfect Phylogeny Newick ########################
	###################################################################################

	if K_==0:
		Phylo_module.gen_Newick(genotype=mat_, PerfectPhy_path=PerfectPhy_path_, out_dir_path=out_path, names_=names)



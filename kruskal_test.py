#This function Kruskal-Wallis H-test
from __future__ import division
import os
import pandas as pd
from scipy import stats
import numpy
import warnings
warnings.filterwarnings("ignore")
cdp=os.getcwd()
while 1:
	os.system('clear')
	print "At this stage, we perform the statistical Kruskal-Wallis H-test on normalized Metaphlan, KO, EC number, and Reaction data."
	print "that stored in the all_results folder."
	print "To get started, you first have to label the class of each sample in the class_label.txt file in the Read_files folder."
	print "The number of classes can be in the range [2:10]. For each sample, enter a separate row. for example:"
	print "class1"
	print "class2"
	print "class1"
	print "class3"
	print
	print "Enter g1 for running the Kruskal-Wallis H-test on the default annotated data."
	print "Enter g2 for running the Kruskal-Wallis H-test on the user extracted data."
	print "Enter 0 for back"
	xx=raw_input()
	#xx=int(c)
	if xx=='g1':
		import os
		os.system('clear')
		print "Please enter the p-value to filter the output. This value can be in the interval [0:1]"
		s=raw_input()
		pv=float(s)
		while pv<0 or pv>1:
			print "Select the p-value in renge of [0:1]="
			s=raw_input()
			pv=float(s)
		file_p=cdp+'/Read_files/class_label.txt'
		fcl=open(file_p,"r")
		cl_la=[]
		cl=fcl.readline().rstrip()
		while cl!="":
			cl_la.append(cl)
			cl=fcl.readline().rstrip()
		fcl.close()
		grps = pd.unique(cl_la)
		if len(grps)>10 or len(grps)<2:
			print "The number of classes can be in the range [2:10]"
			print "Please press the key to continue"	
			ccc=raw_input()
			continue
		nc=len(grps)
		print
		print "Please wait..."
		print
		##############################################Kruskal-Wallis H-test test metaphlan
		file_p=cdp+'/files/taxa_row.txt'
		fmt=open(file_p,"r")
		mt=fmt.readline().rstrip()
		taxa_row=[]
		while mt!="":
			taxa_row.append(mt)
			mt=fmt.readline().rstrip()
		fmt.close()


		file_p=cdp+'/all_results/metaphlan_taxa.txt'
		fmt=open(file_p,"r")
		mt=fmt.readline().rstrip()
		mt=fmt.readline().rstrip()
		meta_taxa1=[]
		while mt!="":
			meta_taxa1.append(mt)
			mt=fmt.readline().rstrip()
		meta_taxa=[]
		for ii in range(0,len(taxa_row)):
			x=int(taxa_row[ii])
			meta_taxa.append(meta_taxa1[x-1])
		ng=len(meta_taxa)
		fmt.close()
		file_p1=cdp+'/all_results/normal_matrix_metaphlan.txt'
		file_p=cdp+'/all_results/kruskal_test_metaphlan.txt'
		fag=open(file_p,"w")
		if nc==2:
			anv_mat=numpy.zeros((ng,3),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				F, p = stats.kruskal(y0,y1)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=p
		elif nc==3:
			anv_mat=numpy.zeros((ng,4),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				F, p = stats.kruskal(y0,y1,y2)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=p
		elif nc==4:
			anv_mat=numpy.zeros((ng,5),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				F, p = stats.kruskal(y0,y1,y2,y3)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=p
		elif nc==5:
			anv_mat=numpy.zeros((ng,6),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				F, p = stats.kruskal(y0,y1,y2,y3,y4)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=p
		elif nc==6:
			anv_mat=numpy.zeros((ng,7),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=p

		elif nc==7:
			anv_mat=numpy.zeros((ng,8),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=p
		elif nc==8:
			anv_mat=numpy.zeros((ng,9),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			a7 = [i for i, x in enumerate(cl_la) if x == grps[7]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				y7 = [x[index] for index in a7]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6,y7)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=sum(y7)/len(y7)
				anv_mat[i,8]=p
		elif nc==9:
			anv_mat=numpy.zeros((ng,10),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			a7 = [i for i, x in enumerate(cl_la) if x == grps[7]]
			a8 = [i for i, x in enumerate(cl_la) if x == grps[8]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				y7 = [x[index] for index in a7]
				y8 = [x[index] for index in a8]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6,y7,y8)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=sum(y7)/len(y7)
				anv_mat[i,8]=sum(y8)/len(y8)
				anv_mat[i,9]=p
		elif nc==10:
			anv_mat=numpy.zeros((ng,11),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			a7 = [i for i, x in enumerate(cl_la) if x == grps[7]]
			a8 = [i for i, x in enumerate(cl_la) if x == grps[8]]
			a9 = [i for i, x in enumerate(cl_la) if x == grps[9]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				y7 = [x[index] for index in a7]
				y8 = [x[index] for index in a8]
				y9 = [x[index] for index in a9]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6,y7,y8,y9)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=sum(y7)/len(y7)
				anv_mat[i,8]=sum(y8)/len(y8)
				anv_mat[i,9]=sum(y9)/len(y9)
				anv_mat[i,10]=p

		fag.write('%-35s'%'sample')	
		fag.write('\t')
		for i in range(0,len(grps)):
			ss='Avg in '+grps[i]
			fag.write('%-20s'%ss)	
			fag.write('\t')
		fag.write('p-value')
		fag.write('\n')
		c=numpy.argsort(anv_mat[:,nc])
		i=0
		while anv_mat[c[i],nc]<=pv:
			text=meta_taxa[c[i]]
			fag.write('%-35s'%text)	
			fag.write('\t')
			for j in range(0,nc+1):
				fag.write('%-20s'%str(anv_mat[c[i]][j]))	
				fag.write('\t')
			fag.write('\n')
			i+=1
		fag.close()
		##############################################Kruskal-Wallis H-test test gene
		file_p=cdp+'/files/gene_row.txt'
		fmt=open(file_p,"r")
		mt=fmt.readline().rstrip()
		gene_row=[]
		while mt!="":
			gene_row.append(mt)
			mt=fmt.readline().rstrip()
		fmt.close()
		ng=len(gene_row)

		file_p1=cdp+'/all_results/normal_matrix_gene.txt'
		file_p=cdp+'/all_results/kruskal_test_gene.txt'
		fag=open(file_p,"w")
		if nc==2:
			anv_mat=numpy.zeros((ng,3),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				F, p = stats.kruskal(y0,y1)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=p
		elif nc==3:
			anv_mat=numpy.zeros((ng,4),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				F, p = stats.kruskal(y0,y1,y2)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=p
		elif nc==4:
			anv_mat=numpy.zeros((ng,5),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				F, p = stats.kruskal(y0,y1,y2,y3)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=p
		elif nc==5:
			anv_mat=numpy.zeros((ng,6),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				F, p = stats.kruskal(y0,y1,y2,y3,y4)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=p
		elif nc==6:
			anv_mat=numpy.zeros((ng,7),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=p

		elif nc==7:
			anv_mat=numpy.zeros((ng,8),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=p
		elif nc==8:
			anv_mat=numpy.zeros((ng,9),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			a7 = [i for i, x in enumerate(cl_la) if x == grps[7]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				y7 = [x[index] for index in a7]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6,y7)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=sum(y7)/len(y7)
				anv_mat[i,8]=p
		elif nc==9:
			anv_mat=numpy.zeros((ng,10),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			a7 = [i for i, x in enumerate(cl_la) if x == grps[7]]
			a8 = [i for i, x in enumerate(cl_la) if x == grps[8]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				y7 = [x[index] for index in a7]
				y8 = [x[index] for index in a8]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6,y7,y8)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=sum(y7)/len(y7)
				anv_mat[i,8]=sum(y8)/len(y8)
				anv_mat[i,9]=p
		elif nc==10:
			anv_mat=numpy.zeros((ng,11),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			a7 = [i for i, x in enumerate(cl_la) if x == grps[7]]
			a8 = [i for i, x in enumerate(cl_la) if x == grps[8]]
			a9 = [i for i, x in enumerate(cl_la) if x == grps[9]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				y7 = [x[index] for index in a7]
				y8 = [x[index] for index in a8]
				y9 = [x[index] for index in a9]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6,y7,y8,y9)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=sum(y7)/len(y7)
				anv_mat[i,8]=sum(y8)/len(y8)
				anv_mat[i,9]=sum(y9)/len(y9)
				anv_mat[i,10]=p

		fag.write('%-14s'%'sample')	
		fag.write('\t')
		for i in range(0,len(grps)):
			ss='Avg in '+grps[i]
			fag.write('%-20s'%ss)	
			fag.write('\t')
		fag.write('p-value')
		fag.write('\n')
		c=numpy.argsort(anv_mat[:,nc])
		i=0
		while anv_mat[c[i],nc]<=pv:
			text='gene'+str(gene_row[c[i]])
			fag.write('%-14s'%text)	
			fag.write('\t')
			for j in range(0,nc+1):
				fag.write('%-20s'%str(anv_mat[c[i]][j]))	
				fag.write('\t')
			fag.write('\n')
			i+=1
		fag.close()
		file_p=cdp+'/files/catalog_name.txt'
		fcn=open(file_p,"r")
		cn=fcn.readline().rstrip()
		fcn.close()
		dd=numpy.zeros((i),dtype=int)
		for k in range(0,i):
			dd[k]=c[k]+1
		dd.sort()
		file_p=cdp+'/all_results/selected_genes.txt'
		fsg=open(file_p,"w")
		file_p=cdp+'/'+cn
		fgc=open(file_p,"r")
		gc=fgc.readline().rstrip('\n')
		for k in range(0,len(dd)):
			while gc[0]!='>' or int(dd[k])!=int(gc[5:]):
				gc=fgc.readline().rstrip('\n')
			fsg.write(gc)	
			fsg.write('\n')
			gc=fgc.readline().rstrip('\n')
			while gc!="" and gc[0]!='>':
				fsg.write(gc)	
				fsg.write('\n')
				gc=fgc.readline().rstrip()
		fgc.close()
		fsg.close()
		##############################################Kruskal-Wallis H-test test KO
		file_p=cdp+'/kegg_annotation/ko.txt'
		fmt=open(file_p,"r")
		mt=fmt.readline().rstrip()
		ko=[]
		while mt!="":
			ko.append(mt)
			mt=fmt.readline().rstrip()
		ng=len(ko)
		fmt.close()
		file_p1=cdp+'/all_results/normal_matrix_ko.txt'
		file_p=cdp+'/all_results/kruskal_test_ko.txt'
		fag=open(file_p,"w")
		if nc==2:
			anv_mat=numpy.zeros((ng,3),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				F, p = stats.kruskal(y0,y1)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=p
		elif nc==3:
			anv_mat=numpy.zeros((ng,4),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				F, p = stats.kruskal(y0,y1,y2)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=p
		elif nc==4:
			anv_mat=numpy.zeros((ng,5),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				F, p = stats.kruskal(y0,y1,y2,y3)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=p
		elif nc==5:
			anv_mat=numpy.zeros((ng,6),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				F, p = stats.kruskal(y0,y1,y2,y3,y4)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=p
		elif nc==6:
			anv_mat=numpy.zeros((ng,7),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=p

		elif nc==7:
			anv_mat=numpy.zeros((ng,8),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=p
		elif nc==8:
			anv_mat=numpy.zeros((ng,9),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			a7 = [i for i, x in enumerate(cl_la) if x == grps[7]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				y7 = [x[index] for index in a7]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6,y7)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=sum(y7)/len(y7)
				anv_mat[i,8]=p
		elif nc==9:
			anv_mat=numpy.zeros((ng,10),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			a7 = [i for i, x in enumerate(cl_la) if x == grps[7]]
			a8 = [i for i, x in enumerate(cl_la) if x == grps[8]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				y7 = [x[index] for index in a7]
				y8 = [x[index] for index in a8]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6,y7,y8)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=sum(y7)/len(y7)
				anv_mat[i,8]=sum(y8)/len(y8)
				anv_mat[i,9]=p
		elif nc==10:
			anv_mat=numpy.zeros((ng,11),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			a7 = [i for i, x in enumerate(cl_la) if x == grps[7]]
			a8 = [i for i, x in enumerate(cl_la) if x == grps[8]]
			a9 = [i for i, x in enumerate(cl_la) if x == grps[9]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				y7 = [x[index] for index in a7]
				y8 = [x[index] for index in a8]
				y9 = [x[index] for index in a9]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6,y7,y8,y9)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=sum(y7)/len(y7)
				anv_mat[i,8]=sum(y8)/len(y8)
				anv_mat[i,9]=sum(y9)/len(y9)
				anv_mat[i,10]=p

		fag.write('%-14s'%'sample')	
		fag.write('\t')
		for i in range(0,len(grps)):
			ss='Avg in '+grps[i]
			fag.write('%-20s'%ss)	
			fag.write('\t')
		fag.write('p-value')
		fag.write('\n')
		c=numpy.argsort(anv_mat[:,nc])
		i=0
		while anv_mat[c[i],nc]<=pv:
			text=ko[c[i]]
			fag.write('%-14s'%text)	
			fag.write('\t')
			for j in range(0,nc+1):
				fag.write('%-20s'%str(anv_mat[c[i]][j]))	
				fag.write('\t')
			fag.write('\n')
			i+=1
		fag.close()
		##############################################Kruskal-Wallis H-test test EC number
		file_p=cdp+'/kegg_annotation/def_ec.txt'
		fmt=open(file_p,"r")
		mt=fmt.readline().rstrip()
		ec=[]
		while mt!="":
			ec.append(mt)
			mt=fmt.readline().rstrip()
		ng=len(ec)
		fmt.close()
		file_p1=cdp+'/all_results/normal_matrix_ec.txt'
		file_p=cdp+'/all_results/kruskal_test_ec.txt'
		fag=open(file_p,"w")
		if nc==2:
			anv_mat=numpy.zeros((ng,3),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				F, p = stats.kruskal(y0,y1)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=p
		elif nc==3:
			anv_mat=numpy.zeros((ng,4),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				F, p = stats.kruskal(y0,y1,y2)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=p
		elif nc==4:
			anv_mat=numpy.zeros((ng,5),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				F, p = stats.kruskal(y0,y1,y2,y3)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=p
		elif nc==5:
			anv_mat=numpy.zeros((ng,6),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				F, p = stats.kruskal(y0,y1,y2,y3,y4)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=p
		elif nc==6:
			anv_mat=numpy.zeros((ng,7),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=p

		elif nc==7:
			anv_mat=numpy.zeros((ng,8),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=p
		elif nc==8:
			anv_mat=numpy.zeros((ng,9),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			a7 = [i for i, x in enumerate(cl_la) if x == grps[7]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				y7 = [x[index] for index in a7]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6,y7)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=sum(y7)/len(y7)
				anv_mat[i,8]=p
		elif nc==9:
			anv_mat=numpy.zeros((ng,10),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			a7 = [i for i, x in enumerate(cl_la) if x == grps[7]]
			a8 = [i for i, x in enumerate(cl_la) if x == grps[8]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				y7 = [x[index] for index in a7]
				y8 = [x[index] for index in a8]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6,y7,y8)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=sum(y7)/len(y7)
				anv_mat[i,8]=sum(y8)/len(y8)
				anv_mat[i,9]=p
		elif nc==10:
			anv_mat=numpy.zeros((ng,11),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			a7 = [i for i, x in enumerate(cl_la) if x == grps[7]]
			a8 = [i for i, x in enumerate(cl_la) if x == grps[8]]
			a9 = [i for i, x in enumerate(cl_la) if x == grps[9]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				y7 = [x[index] for index in a7]
				y8 = [x[index] for index in a8]
				y9 = [x[index] for index in a9]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6,y7,y8,y9)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=sum(y7)/len(y7)
				anv_mat[i,8]=sum(y8)/len(y8)
				anv_mat[i,9]=sum(y9)/len(y9)
				anv_mat[i,10]=p

		fag.write('%-14s'%'sample')	
		fag.write('\t')
		for i in range(0,len(grps)):
			ss='Avg in '+grps[i]
			fag.write('%-20s'%ss)	
			fag.write('\t')
		fag.write('p-value')
		fag.write('\n')
		c=numpy.argsort(anv_mat[:,nc])
		i=0
		while anv_mat[c[i],nc]<=pv:
			text=ec[c[i]]
			fag.write('%-14s'%text)	
			fag.write('\t')
			for j in range(0,nc+1):
				fag.write('%-20s'%str(anv_mat[c[i]][j]))	
				fag.write('\t')
			fag.write('\n')
			i+=1
		fag.close()
		##############################################Kruskal-Wallis H-test test reaction number
		file_p=cdp+'/kegg_annotation/def_re.txt'
		fmt=open(file_p,"r")
		mt=fmt.readline().rstrip()
		re=[]
		while mt!="":
			re.append(mt)
			mt=fmt.readline().rstrip()
		ng=len(re)
		fmt.close()
		file_p1=cdp+'/all_results/normal_matrix_re.txt'
		file_p=cdp+'/all_results/kruskal_test_re.txt'
		fag=open(file_p,"w")
		if nc==2:
			anv_mat=numpy.zeros((ng,3),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				F, p = stats.kruskal(y0,y1)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=p
		elif nc==3:
			anv_mat=numpy.zeros((ng,4),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				F, p = stats.kruskal(y0,y1,y2)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=p
		elif nc==4:
			anv_mat=numpy.zeros((ng,5),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				F, p = stats.kruskal(y0,y1,y2,y3)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=p
		elif nc==5:
			anv_mat=numpy.zeros((ng,6),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				F, p = stats.kruskal(y0,y1,y2,y3,y4)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=p
		elif nc==6:
			anv_mat=numpy.zeros((ng,7),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=p

		elif nc==7:
			anv_mat=numpy.zeros((ng,8),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=p
		elif nc==8:
			anv_mat=numpy.zeros((ng,9),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			a7 = [i for i, x in enumerate(cl_la) if x == grps[7]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				y7 = [x[index] for index in a7]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6,y7)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=sum(y7)/len(y7)
				anv_mat[i,8]=p
		elif nc==9:
			anv_mat=numpy.zeros((ng,10),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			a7 = [i for i, x in enumerate(cl_la) if x == grps[7]]
			a8 = [i for i, x in enumerate(cl_la) if x == grps[8]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				y7 = [x[index] for index in a7]
				y8 = [x[index] for index in a8]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6,y7,y8)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=sum(y7)/len(y7)
				anv_mat[i,8]=sum(y8)/len(y8)
				anv_mat[i,9]=p
		elif nc==10:
			anv_mat=numpy.zeros((ng,11),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			a7 = [i for i, x in enumerate(cl_la) if x == grps[7]]
			a8 = [i for i, x in enumerate(cl_la) if x == grps[8]]
			a9 = [i for i, x in enumerate(cl_la) if x == grps[9]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				y7 = [x[index] for index in a7]
				y8 = [x[index] for index in a8]
				y9 = [x[index] for index in a9]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6,y7,y8,y9)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=sum(y7)/len(y7)
				anv_mat[i,8]=sum(y8)/len(y8)
				anv_mat[i,9]=sum(y9)/len(y9)
				anv_mat[i,10]=p

		fag.write('%-14s'%'sample')	
		fag.write('\t')
		for i in range(0,len(grps)):
			ss='Avg in '+grps[i]
			fag.write('%-20s'%ss)	
			fag.write('\t')
		fag.write('p-value')
		fag.write('\n')
		c=numpy.argsort(anv_mat[:,nc])
		i=0
		while anv_mat[c[i],nc]<=pv:
			text=re[c[i]]
			fag.write('%-14s'%text)	
			fag.write('\t')
			for j in range(0,nc+1):
				fag.write('%-20s'%str(anv_mat[c[i]][j]))	
				fag.write('\t')
			fag.write('\n')
			i+=1
		fag.close()
		######################
		rr=[]
		for k in range(0,i):
			rr.append(re[c[k]])
		rr.sort()
		file_p=cdp+'/all_results/selected_reactions.txt'
		fsg=open(file_p,"w")
		file_p=cdp+'/kegg_annotation/def_re_eq.txt'
		fgc=open(file_p,"r")
		gc=fgc.readline().rstrip('\n')
		for k in range(0,len(rr)):
			while gc!="" and rr[k]!=gc[9:15]:
				gc=fgc.readline().rstrip('\n')
			fsg.write(gc)	
			fsg.write('\n')
			gc=fgc.readline().rstrip()
		fgc.close()
		fsg.close()
		###################
		print ""
		print "The ANOVA test was performed on all the normal data extracted from the samples and the results are stored"
		print "in the text files kruskal_test_metaphlan.txt, kruskal_test_gene.txt, Kruskal-Wallis H-test_test_ko.txt, Kruskal-Wallis H-test_test_ec.txt, and"
		print "kruskal_test_re.txt in the all_results folder, respectively."
		print ""
		print ""
		print "Also, the sequences of p-value filtered genes are stored in a text file named selected_genes.txt in folder all_results."
		print ""
		print "Finally, all filtered reactions with p-value and their definition and equation are stored in a file named selected_reactions.txt in folder all_results."
		print ""
		print "Please press the key to continue"
		c=raw_input()
	elif xx=='g2':
		import os
		os.system('clear')
		print "Please enter the p-value to filter the output. This value can be in the interval [0:1]"
		s=raw_input()
		pv=float(s)
		while pv<0 or pv>1:
			print "Select the p-value in renge of [0:1]="
			s=raw_input()
			pv=float(s)
		file_p=cdp+'/Read_files/class_label.txt'
		fcl=open(file_p,"r")
		cl_la=[]
		cl=fcl.readline().rstrip()
		while cl!="":
			cl_la.append(cl)
			cl=fcl.readline().rstrip()
		fcl.close()
		grps = pd.unique(cl_la)
		if len(grps)>10 or len(grps)<2:
			print "The number of classes can be in the range [2:10]"
			print "Please press the key to continue"	
			ccc=raw_input()
			continue
		nc=len(grps)
		print
		print "Please wait..."
		print
		##############################################Kruskal-Wallis H-test test metaphlan
		file_p=cdp+'/all_results/metaphlan_taxa.txt'
		fmt=open(file_p,"r")
		mt=fmt.readline().rstrip()
		mt=fmt.readline().rstrip()
		meta_taxa=[]
		while mt!="":
			meta_taxa.append(mt)
			mt=fmt.readline().rstrip()
		ng=len(meta_taxa)
		fmt.close()
		file_p1=cdp+'/all_results/normal_matrix_metaphlan.txt'
		file_p=cdp+'/all_results/kruskal_test_metaphlan.txt'
		fag=open(file_p,"w")
		if nc==2:
			anv_mat=numpy.zeros((ng,3),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				F, p = stats.kruskal(y0,y1)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=p
		elif nc==3:
			anv_mat=numpy.zeros((ng,4),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				F, p = stats.kruskal(y0,y1,y2)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=p
		elif nc==4:
			anv_mat=numpy.zeros((ng,5),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				F, p = stats.kruskal(y0,y1,y2,y3)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=p
		elif nc==5:
			anv_mat=numpy.zeros((ng,6),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				F, p = stats.kruskal(y0,y1,y2,y3,y4)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=p
		elif nc==6:
			anv_mat=numpy.zeros((ng,7),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=p

		elif nc==7:
			anv_mat=numpy.zeros((ng,8),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=p
		elif nc==8:
			anv_mat=numpy.zeros((ng,9),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			a7 = [i for i, x in enumerate(cl_la) if x == grps[7]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				y7 = [x[index] for index in a7]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6,y7)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=sum(y7)/len(y7)
				anv_mat[i,8]=p
		elif nc==9:
			anv_mat=numpy.zeros((ng,10),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			a7 = [i for i, x in enumerate(cl_la) if x == grps[7]]
			a8 = [i for i, x in enumerate(cl_la) if x == grps[8]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				y7 = [x[index] for index in a7]
				y8 = [x[index] for index in a8]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6,y7,y8)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=sum(y7)/len(y7)
				anv_mat[i,8]=sum(y8)/len(y8)
				anv_mat[i,9]=p
		elif nc==10:
			anv_mat=numpy.zeros((ng,11),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			a7 = [i for i, x in enumerate(cl_la) if x == grps[7]]
			a8 = [i for i, x in enumerate(cl_la) if x == grps[8]]
			a9 = [i for i, x in enumerate(cl_la) if x == grps[9]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				y7 = [x[index] for index in a7]
				y8 = [x[index] for index in a8]
				y9 = [x[index] for index in a9]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6,y7,y8,y9)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=sum(y7)/len(y7)
				anv_mat[i,8]=sum(y8)/len(y8)
				anv_mat[i,9]=sum(y9)/len(y9)
				anv_mat[i,10]=p

		fag.write('%-35s'%'sample')	
		fag.write('\t')
		for i in range(0,len(grps)):
			ss='Avg in '+grps[i]
			fag.write('%-20s'%ss)	
			fag.write('\t')
		fag.write('p-value')
		fag.write('\n')
		c=numpy.argsort(anv_mat[:,nc])
		i=0
		while anv_mat[c[i],nc]<=pv:
			text=meta_taxa[c[i]]
			fag.write('%-35s'%text)	
			fag.write('\t')
			for j in range(0,nc+1):
				fag.write('%-20s'%str(anv_mat[c[i]][j]))	
				fag.write('\t')
			fag.write('\n')
			i+=1
		fag.close()
		##############################################Kruskal-Wallis H-test test gene
		file_p=cdp+'/files/catalog_length.txt'
		fgl=open(file_p,"r")
		gl=fgl.readline().rstrip()
		ng=int(gl)
		fgl.close()
		file_p1=cdp+'/all_results/normal_matrix_gene.txt'
		file_p=cdp+'/all_results/kruskal_test_gene.txt'
		fag=open(file_p,"w")
		if nc==2:
			anv_mat=numpy.zeros((ng,3),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				F, p = stats.kruskal(y0,y1)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=p
		elif nc==3:
			anv_mat=numpy.zeros((ng,4),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				F, p = stats.kruskal(y0,y1,y2)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=p
		elif nc==4:
			anv_mat=numpy.zeros((ng,5),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				F, p = stats.kruskal(y0,y1,y2,y3)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=p
		elif nc==5:
			anv_mat=numpy.zeros((ng,6),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				F, p = stats.kruskal(y0,y1,y2,y3,y4)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=p
		elif nc==6:
			anv_mat=numpy.zeros((ng,7),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=p

		elif nc==7:
			anv_mat=numpy.zeros((ng,8),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=p
		elif nc==8:
			anv_mat=numpy.zeros((ng,9),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			a7 = [i for i, x in enumerate(cl_la) if x == grps[7]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				y7 = [x[index] for index in a7]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6,y7)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=sum(y7)/len(y7)
				anv_mat[i,8]=p
		elif nc==9:
			anv_mat=numpy.zeros((ng,10),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			a7 = [i for i, x in enumerate(cl_la) if x == grps[7]]
			a8 = [i for i, x in enumerate(cl_la) if x == grps[8]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				y7 = [x[index] for index in a7]
				y8 = [x[index] for index in a8]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6,y7,y8)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=sum(y7)/len(y7)
				anv_mat[i,8]=sum(y8)/len(y8)
				anv_mat[i,9]=p
		elif nc==10:
			anv_mat=numpy.zeros((ng,11),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			a7 = [i for i, x in enumerate(cl_la) if x == grps[7]]
			a8 = [i for i, x in enumerate(cl_la) if x == grps[8]]
			a9 = [i for i, x in enumerate(cl_la) if x == grps[9]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				y7 = [x[index] for index in a7]
				y8 = [x[index] for index in a8]
				y9 = [x[index] for index in a9]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6,y7,y8,y9)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=sum(y7)/len(y7)
				anv_mat[i,8]=sum(y8)/len(y8)
				anv_mat[i,9]=sum(y9)/len(y9)
				anv_mat[i,10]=p

		fag.write('%-14s'%'sample')	
		fag.write('\t')
		for i in range(0,len(grps)):
			ss='Avg in '+grps[i]
			fag.write('%-20s'%ss)	
			fag.write('\t')
		fag.write('p-value')
		fag.write('\n')
		c=numpy.argsort(anv_mat[:,nc])
		i=0
		while anv_mat[c[i],nc]<=pv:
			text='gene'+str(c[i]+1)
			fag.write('%-14s'%text)	
			fag.write('\t')
			for j in range(0,nc+1):
				fag.write('%-20s'%str(anv_mat[c[i]][j]))	
				fag.write('\t')
			fag.write('\n')
			i+=1
		fag.close()
		file_p=cdp+'/files/catalog_name.txt'
		fcn=open(file_p,"r")
		cn=fcn.readline().rstrip()
		fcn.close()
		dd=numpy.zeros((i),dtype=int)
		for k in range(0,i):
			dd[k]=c[k]+1
		dd.sort()
		file_p=cdp+'/all_results/selected_genes.txt'
		fsg=open(file_p,"w")
		file_p=cdp+'/'+cn
		fgc=open(file_p,"r")
		gc=fgc.readline().rstrip('\n')
		for k in range(0,len(dd)):
			while gc[0]!='>' or int(dd[k])!=int(gc[5:]):
				gc=fgc.readline().rstrip('\n')
			fsg.write(gc)	
			fsg.write('\n')
			gc=fgc.readline().rstrip('\n')
			while gc!="" and gc[0]!='>':
				fsg.write(gc)	
				fsg.write('\n')
				gc=fgc.readline().rstrip()
		fgc.close()
		fsg.close()
		##############################################Kruskal-Wallis H-test test KO
		file_p=cdp+'/kegg_annotation/ko.txt'
		fmt=open(file_p,"r")
		mt=fmt.readline().rstrip()
		ko=[]
		while mt!="":
			ko.append(mt)
			mt=fmt.readline().rstrip()
		ng=len(ko)
		fmt.close()
		file_p1=cdp+'/all_results/normal_matrix_ko.txt'
		file_p=cdp+'/all_results/kruskal_test_ko.txt'
		fag=open(file_p,"w")
		if nc==2:
			anv_mat=numpy.zeros((ng,3),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				F, p = stats.kruskal(y0,y1)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=p
		elif nc==3:
			anv_mat=numpy.zeros((ng,4),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				F, p = stats.kruskal(y0,y1,y2)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=p
		elif nc==4:
			anv_mat=numpy.zeros((ng,5),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				F, p = stats.kruskal(y0,y1,y2,y3)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=p
		elif nc==5:
			anv_mat=numpy.zeros((ng,6),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				F, p = stats.kruskal(y0,y1,y2,y3,y4)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=p
		elif nc==6:
			anv_mat=numpy.zeros((ng,7),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=p

		elif nc==7:
			anv_mat=numpy.zeros((ng,8),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=p
		elif nc==8:
			anv_mat=numpy.zeros((ng,9),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			a7 = [i for i, x in enumerate(cl_la) if x == grps[7]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				y7 = [x[index] for index in a7]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6,y7)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=sum(y7)/len(y7)
				anv_mat[i,8]=p
		elif nc==9:
			anv_mat=numpy.zeros((ng,10),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			a7 = [i for i, x in enumerate(cl_la) if x == grps[7]]
			a8 = [i for i, x in enumerate(cl_la) if x == grps[8]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				y7 = [x[index] for index in a7]
				y8 = [x[index] for index in a8]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6,y7,y8)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=sum(y7)/len(y7)
				anv_mat[i,8]=sum(y8)/len(y8)
				anv_mat[i,9]=p
		elif nc==10:
			anv_mat=numpy.zeros((ng,11),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			a7 = [i for i, x in enumerate(cl_la) if x == grps[7]]
			a8 = [i for i, x in enumerate(cl_la) if x == grps[8]]
			a9 = [i for i, x in enumerate(cl_la) if x == grps[9]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				y7 = [x[index] for index in a7]
				y8 = [x[index] for index in a8]
				y9 = [x[index] for index in a9]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6,y7,y8,y9)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=sum(y7)/len(y7)
				anv_mat[i,8]=sum(y8)/len(y8)
				anv_mat[i,9]=sum(y9)/len(y9)
				anv_mat[i,10]=p

		fag.write('%-14s'%'sample')	
		fag.write('\t')
		for i in range(0,len(grps)):
			ss='Avg in '+grps[i]
			fag.write('%-20s'%ss)	
			fag.write('\t')
		fag.write('p-value')
		fag.write('\n')
		c=numpy.argsort(anv_mat[:,nc])
		i=0
		while anv_mat[c[i],nc]<=pv:
			text=ko[c[i]]
			fag.write('%-14s'%text)	
			fag.write('\t')
			for j in range(0,nc+1):
				fag.write('%-20s'%str(anv_mat[c[i]][j]))	
				fag.write('\t')
			fag.write('\n')
			i+=1
		fag.close()
		###############################Kruskal-Wallis H-test test EC number
		file_p=cdp+'/kegg_annotation/ec.txt'
		fmt=open(file_p,"r")
		mt=fmt.readline().rstrip()
		ec=[]
		while mt!="":
			ec.append(mt)
			mt=fmt.readline().rstrip()
		ng=len(ec)
		fmt.close()
		file_p1=cdp+'/all_results/normal_matrix_ec.txt'
		file_p=cdp+'/all_results/kruskal_test_ec.txt'
		fag=open(file_p,"w")
		if nc==2:
			anv_mat=numpy.zeros((ng,3),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				F, p = stats.kruskal(y0,y1)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=p
		elif nc==3:
			anv_mat=numpy.zeros((ng,4),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				F, p = stats.kruskal(y0,y1,y2)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=p
		elif nc==4:
			anv_mat=numpy.zeros((ng,5),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				F, p = stats.kruskal(y0,y1,y2,y3)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=p
		elif nc==5:
			anv_mat=numpy.zeros((ng,6),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				F, p = stats.kruskal(y0,y1,y2,y3,y4)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=p
		elif nc==6:
			anv_mat=numpy.zeros((ng,7),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=p

		elif nc==7:
			anv_mat=numpy.zeros((ng,8),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=p
		elif nc==8:
			anv_mat=numpy.zeros((ng,9),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			a7 = [i for i, x in enumerate(cl_la) if x == grps[7]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				y7 = [x[index] for index in a7]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6,y7)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=sum(y7)/len(y7)
				anv_mat[i,8]=p
		elif nc==9:
			anv_mat=numpy.zeros((ng,10),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			a7 = [i for i, x in enumerate(cl_la) if x == grps[7]]
			a8 = [i for i, x in enumerate(cl_la) if x == grps[8]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				y7 = [x[index] for index in a7]
				y8 = [x[index] for index in a8]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6,y7,y8)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=sum(y7)/len(y7)
				anv_mat[i,8]=sum(y8)/len(y8)
				anv_mat[i,9]=p
		elif nc==10:
			anv_mat=numpy.zeros((ng,11),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			a7 = [i for i, x in enumerate(cl_la) if x == grps[7]]
			a8 = [i for i, x in enumerate(cl_la) if x == grps[8]]
			a9 = [i for i, x in enumerate(cl_la) if x == grps[9]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				y7 = [x[index] for index in a7]
				y8 = [x[index] for index in a8]
				y9 = [x[index] for index in a9]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6,y7,y8,y9)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=sum(y7)/len(y7)
				anv_mat[i,8]=sum(y8)/len(y8)
				anv_mat[i,9]=sum(y9)/len(y9)
				anv_mat[i,10]=p

		fag.write('%-14s'%'sample')	
		fag.write('\t')
		for i in range(0,len(grps)):
			ss='Avg in '+grps[i]
			fag.write('%-20s'%ss)	
			fag.write('\t')
		fag.write('p-value')
		fag.write('\n')
		c=numpy.argsort(anv_mat[:,nc])
		i=0
		while anv_mat[c[i],nc]<=pv:
			text=ec[c[i]]
			fag.write('%-14s'%text)	
			fag.write('\t')
			for j in range(0,nc+1):
				fag.write('%-20s'%str(anv_mat[c[i]][j]))	
				fag.write('\t')
			fag.write('\n')
			i+=1
		fag.close()
		##############################################Kruskal-Wallis H-test test reaction number
		file_p=cdp+'/kegg_annotation/re.txt'
		fmt=open(file_p,"r")
		mt=fmt.readline().rstrip()
		re=[]
		while mt!="":
			re.append(mt)
			mt=fmt.readline().rstrip()
		ng=len(re)
		fmt.close()
		file_p1=cdp+'/all_results/normal_matrix_re.txt'
		file_p=cdp+'/all_results/kruskal_test_re.txt'
		fag=open(file_p,"w")
		if nc==2:
			anv_mat=numpy.zeros((ng,3),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				F, p = stats.kruskal(y0,y1)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=p
		elif nc==3:
			anv_mat=numpy.zeros((ng,4),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				F, p = stats.kruskal(y0,y1,y2)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=p
		elif nc==4:
			anv_mat=numpy.zeros((ng,5),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				F, p = stats.kruskal(y0,y1,y2,y3)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=p
		elif nc==5:
			anv_mat=numpy.zeros((ng,6),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				F, p = stats.kruskal(y0,y1,y2,y3,y4)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=p
		elif nc==6:
			anv_mat=numpy.zeros((ng,7),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=p

		elif nc==7:
			anv_mat=numpy.zeros((ng,8),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=p
		elif nc==8:
			anv_mat=numpy.zeros((ng,9),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			a7 = [i for i, x in enumerate(cl_la) if x == grps[7]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				y7 = [x[index] for index in a7]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6,y7)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=sum(y7)/len(y7)
				anv_mat[i,8]=p
		elif nc==9:
			anv_mat=numpy.zeros((ng,10),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			a7 = [i for i, x in enumerate(cl_la) if x == grps[7]]
			a8 = [i for i, x in enumerate(cl_la) if x == grps[8]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				y7 = [x[index] for index in a7]
				y8 = [x[index] for index in a8]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6,y7,y8)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=sum(y7)/len(y7)
				anv_mat[i,8]=sum(y8)/len(y8)
				anv_mat[i,9]=p
		elif nc==10:
			anv_mat=numpy.zeros((ng,11),dtype=float)
			g_mat=numpy.loadtxt(file_p1, dtype=float)
			a0 = [i for i, x in enumerate(cl_la) if x == grps[0]]
			a1 = [i for i, x in enumerate(cl_la) if x == grps[1]]
			a2 = [i for i, x in enumerate(cl_la) if x == grps[2]]
			a3 = [i for i, x in enumerate(cl_la) if x == grps[3]]
			a4 = [i for i, x in enumerate(cl_la) if x == grps[4]]
			a5 = [i for i, x in enumerate(cl_la) if x == grps[5]]
			a6 = [i for i, x in enumerate(cl_la) if x == grps[6]]
			a7 = [i for i, x in enumerate(cl_la) if x == grps[7]]
			a8 = [i for i, x in enumerate(cl_la) if x == grps[8]]
			a9 = [i for i, x in enumerate(cl_la) if x == grps[9]]
			for i in range(0,ng):
				x = g_mat[i,]
				y0 = [x[index] for index in a0]
				y1 = [x[index] for index in a1]
				y2 = [x[index] for index in a2]
				y3 = [x[index] for index in a3]
				y4 = [x[index] for index in a4]
				y5 = [x[index] for index in a5]
				y6 = [x[index] for index in a6]
				y7 = [x[index] for index in a7]
				y8 = [x[index] for index in a8]
				y9 = [x[index] for index in a9]
				F, p = stats.kruskal(y0,y1,y2,y3,y4,y5,y6,y7,y8,y9)
				anv_mat[i,0]=sum(y0)/len(y0)
				anv_mat[i,1]=sum(y1)/len(y1)
				anv_mat[i,2]=sum(y2)/len(y2)
				anv_mat[i,3]=sum(y3)/len(y3)
				anv_mat[i,4]=sum(y4)/len(y4)
				anv_mat[i,5]=sum(y5)/len(y5)
				anv_mat[i,6]=sum(y6)/len(y6)
				anv_mat[i,7]=sum(y7)/len(y7)
				anv_mat[i,8]=sum(y8)/len(y8)
				anv_mat[i,9]=sum(y9)/len(y9)
				anv_mat[i,10]=p

		fag.write('%-14s'%'sample')	
		fag.write('\t')
		for i in range(0,len(grps)):
			ss='Avg in '+grps[i]
			fag.write('%-20s'%ss)	
			fag.write('\t')
		fag.write('p-value')
		fag.write('\n')
		c=numpy.argsort(anv_mat[:,nc])
		i=0
		while anv_mat[c[i],nc]<=pv:
			text=re[c[i]]
			fag.write('%-14s'%text)	
			fag.write('\t')
			for j in range(0,nc+1):
				fag.write('%-20s'%str(anv_mat[c[i]][j]))	
				fag.write('\t')
			fag.write('\n')
			i+=1
		fag.close()
		######################
		rr=[]
		for k in range(0,i):
			rr.append(re[c[k]])
		rr.sort()
		file_p=cdp+'/all_results/selected_reactions.txt'
		fsg=open(file_p,"w")
		file_p=cdp+'/kegg_annotation/re_eq.txt'
		fgc=open(file_p,"r")
		gc=fgc.readline().rstrip('\n')
		for k in range(0,len(rr)):
			while gc!="" and rr[k]!=gc[9:15]:
				gc=fgc.readline().rstrip('\n')
			fsg.write(gc)	
			fsg.write('\n')
			gc=fgc.readline().rstrip()
		fgc.close()
		fsg.close()
		###################
		print ""
		print "The ANOVA test was performed on all the normal data extracted from the samples and the results are stored"
		print "in the text files kruskal_test_metaphlan.txt, kruskal_test_gene.txt, kruskal_test_ko.txt, kruskal_test_ec.txt, and"
		print "kruskal_test_re.txt in the all_results folder, respectively."
		print ""
		print ""
		print "Also, the sequences of p-value filtered genes are stored in a text file named selected_genes.txt in folder all_results"
		print ""
		print "Finally, all filtered reactions with p-value and their definition and equation are stored in a file named selected_reactions.txt in folder all_results."
		print ""
		print "Please press the key to continue"
		c=raw_input()

	elif xx=='0':
		break


#This function Extract samples information in normalized bacteria, gene, KO, EC number and reaction matrix
from __future__ import division
import os
import numpy
import re
os.system('clear')
print "This function Extract samples information in normalized bacteria, gene, KO, EC number and reaction matrix"
print ""
cdp=os.getcwd()
################################
#this code extract number of samples
file_p=cdp+'/sample_file_names.txt'
fsn=open(file_p,"r")
snc=0
sn=fsn.readline().rstrip()
while sn!="":
	snc+=1
	sn=fsn.readline().rstrip()
fsn.close()
###############################
#this code extract number of taxa in samples
file_p=cdp+'/all_results/metaphlan_taxa.txt'
fmpt=open(file_p,"r")
mpc=0
mp=fmpt.readline().rstrip()
mp=fmpt.readline().rstrip()
while mp!="":
	mpc+=1
	mp=fmpt.readline().rstrip()
fmpt.close()
###############################
#this code extract number reads in samples
file_p=cdp+'/files/sample_reads.txt'
fsl=open(file_p,"r")
sa_le=[]
sl=fsl.readline().rstrip('\n')
while sl!="":
	sa_le.append(sl)
	sl=fsl.readline().rstrip()
fsl.close()
###############################
#this code extract length of genes in catalog
file_p=cdp+'/files/gene_length.txt'
fgl=open(file_p,"r")
ge_le=[]
gl=fgl.readline().rstrip()
while gl!="":
	ge_le.append(gl)
	gl=fgl.readline().rstrip()
gn=len(ge_le)
fgl.close()
###############################
#this code extract type of samples "paired end or single end"
file_p=cdp+'/files/catalog_align.txt'
fca=open(file_p,"r")
ca=fca.readline().rstrip()
ca=fca.readline().rstrip()
###############################
#this code open a file for normalized value of metaphlan outputs
file_p=cdp+'/all_results/normal_matrix_metaphlan.txt'
fmp=open(file_p,"w")

##############################################
#this code extract nornalized data from metaphlan outputs
if ca=='p':
	me_ph_mat=numpy.zeros((mpc,int(snc/2)),dtype=float)
else:
	me_ph_mat=numpy.zeros((mpc,snc),dtype=float)
file_p=cdp+'/all_results/total_metaphlan_results.txt'
fmp=open(file_p,"r")
mp=fmp.readline().rstrip()
mp=fmp.readline().rstrip()
mp=fmp.readline().rstrip()
mp=fmp.readline().rstrip()
i=0
while mp!="":
	j=0
	if ca=='p':
		k=0
		h=[m.start() for m in re.finditer('\t', mp)]
		#print(h)
		while k<len(h):
			#print(mp[h[k]+1:h[k+1]])
			a=float(mp[h[k]+1:h[k+1]])
			if (k+3)<len(h):
				#print(mp[h[k+1]+1:h[k+2]])
				b=float(mp[h[k+1]+1:h[k+2]])
			else:
				b=float(mp[h[k+1]+1:len(mp)])
				#print(mp[h[k+1]+1:len(mp)])
			me_ph_mat[i][j]=(a+b)/2
			j+=1
			k+=2
	else:
		k=0
		h=[m.start() for m in re.finditer('\t', mp)]
		#print(h)
		while k<len(h):
			if (k+1)<len(h):
				a=float(mp[h[k]+1:h[k+1]])
				#print(a)
			else:
				a=float(mp[h[k]+1:len(mp)])
				#print(a)
			me_ph_mat[i][j]=a
			j+=1
			k+=1
	i+=1
	mp=fmp.readline().rstrip()
fmp.close()
##############################calculate normalized value using Metagenomeseq1 with R package
file_p=cdp+'/all_results/normal_matrix_metaphlan.txt'
fmp=open(file_p,"w")
si=numpy.shape(me_ph_mat)
i=0
while i<si[0]:
	j=0
	while j<si[1]:
		fmp.write(str(me_ph_mat[i][j]))
		if j<si[1]:
			fmp.write('\t')	
		j+=1
	fmp.write('\n')	
	i+=1
fmp.close()
del me_ph_mat

file_p=cdp+'/files/param.txt'
fmp=open(file_p,"w")
if ca=='p':
	fmp.write(str(snc/2))
else:
	fmp.write(str(snc))
fmp.write('\n')	
print "Please enter minimum total count for each taxon. This value is percentage of all taxa, for example 1.2"
s=raw_input()
pv=float(s)
while pv<0:
	print "Please enter a positive floating point number"
	s=raw_input()
	pv=float(s)
fmp.write(str(pv))
fmp.close()

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 
from rpy2.robjects.packages import importr
metagenomeSeq = importr("metagenomeSeq")
edgeR=importr("edgeR")
from rpy2 import robjects
robjects.r('''
	s=getwd()
	n <- scan(paste0(s,'/files/param.txt'))
	n <- matrix(n, ncol = 1, byrow = TRUE)
	tot_tir <- scan(paste0(s,'/all_results/normal_matrix_metaphlan.txt'))
	tot_tir <- matrix(tot_tir, ncol = n[1], byrow = TRUE)
	dge <- DGEList(counts=tot_tir)
	keep <- filterByExpr(dge,min.count = 0, min.total.count = n[2])
	dge <- dge[keep, , keep.lib.sizes=FALSE]
	obj = newMRexperiment(dge)
	matx=cumNorm(obj, p = cumNormStatFast(obj))
	norm_mat_tot=MRcounts(matx, norm = TRUE, log = FALSE)
	row_name <- rownames(norm_mat_tot)
	write.table(norm_mat_tot, file=paste0(s,'/all_results/normal_matrix_metaphlan.txt'), sep="\t", col.names = F, row.names = F)
	write.table(as.integer(row_name), file=paste0(s,'/files/taxa_row.txt'), sep="\t", col.names = F, row.names = F)
	''')

###############################
#this code open a file for normalized value of metaphlan outputs
#file_p=cdp+'/all_results/normal_matrix_metaphlan.txt'
#fmp=open(file_p,"w")
#si=numpy.shape(me_ph_mat)
#i=0
#while i<si[0]:
	#j=0
	#while j<si[1]:
	#	fmp.write('%20s'%str(me_ph_mat[i][j]))
	#	if j<si[1]-1:
	#		fmp.write('\t')	
	#	j+=1
	#fmp.write('\n')	
	#i+=1
#fmp.close()
#del me_ph_mat
print "The normal values of the bacteria were extracted for samples and stored in the /all_results/normal_matrix_metaphlan.txt file."
print "In this file, the rows represent the bacteria and the columns represent the samples."
print ""
##############################################
#this code extract nornalized data from mosaik outputs
print ""
print ""
print "In this step, the abundance of genes is reading from the MOSAIK output files, which may take a long time."
print "please wait ..."
if ca=='p':
	mos_mat=numpy.zeros((gn,int(snc/2)),dtype=int)
else:
	mos_mat=numpy.zeros((gn,snc),dtype=int)
file_p=cdp+'/sample_file_names.txt'
fsn=open(file_p,"r")
sn=fsn.readline().rstrip()
j=0
while sn!="":
	ssss='The '+sn+'.sam file is now read'
	print ssss
	file_p=cdp+'/mosaik_outputs/'+sn+'.sam'
	fsam=open(file_p,"r")
	sam=fsam.readline().rstrip()
	while sam[0]=='@':
		sam=fsam.readline().rstrip()
	while sam!="":
		h=[m.start() for m in re.finditer('\t', sam)]
		a=int(sam[h[1]+5:h[2]])
		#print(a)
		mos_mat[a-1][j]+=1
		sam=fsam.readline().rstrip()
	j+=1
	sn=fsn.readline().rstrip()
	if ca=='p':
		sn=fsn.readline().rstrip()
fsn.close()
###############################
#this code open a file for normalized value of mosaik outputs
file_p=cdp+'/all_results/normal_matrix_gene.txt'
fge=open(file_p,"w")
si=numpy.shape(mos_mat)
i=0
while i<si[0]:
	j=0
	while j<si[1]:
		fge.write(str(mos_mat[i][j]))
		if j<si[1]:
			fge.write('\t')	
		j+=1
	fge.write('\n')	
	i+=1
fge.close()
del mos_mat

file_p=cdp+'/files/param.txt'
fmp=open(file_p,"w")
if ca=='p':
	fmp.write(str(snc/2))
else:
	fmp.write(str(snc))
fmp.write('\n')	
print "Please enter minimum total count of mapped reads per gene in total samples. For example 5"
s=raw_input()
pv=int(s)
while pv<0:
	print "Please enter a positive integer number"
	s=raw_input()
	pv=int(s)
fmp.write(str(pv))
fmp.close()

robjects.r('''
	s=getwd()
	n <- scan(paste0(s,'/files/param.txt'))
	n <- matrix(n, ncol = 1, byrow = TRUE)
	tot_tir <- scan(paste0(s,'/all_results/normal_matrix_gene.txt'))
	tot_tir <- matrix(tot_tir, ncol = n[1], byrow = TRUE)
	glr <- scan(paste0(s,'/files/gene_length.txt'))
	glr <- matrix(glr, ncol = 1, byrow = TRUE)
	dge <- DGEList(counts=tot_tir)
	keep <- filterByExpr(dge,min.count = 0, min.total.count = n[2])
	dge <- dge[keep, , keep.lib.sizes=FALSE]
	row_name <- as.integer(rownames(dge))
	glr=glr[row_name]
	tot_tir=tot_tir[row_name,]
	for (i in 1:ncol(tot_tir))
	{
		tot_tir[,i]=tot_tir[,i]/glr
	}
	obj = newMRexperiment(tot_tir)
	matx=cumNorm(obj, p = cumNormStatFast(obj))
	norm_mat_tot=MRcounts(matx, norm = TRUE, log = FALSE)
	write.table(norm_mat_tot, file=paste0(s,'/all_results/normal_matrix_gene.txt'), sep="\t", col.names = F, row.names = F)
	write.table(row_name, file=paste0(s,'/files/gene_row.txt'), sep="\t", col.names = F, row.names = F)
	''')


print "The normal values of the genes were extracted for samples and stored in the /all_results/normal_matrix_gene.txt file."
print "In this file, the rows represent the gene and the columns represent the samples."
print ""
##############################################
#this code extract nornalized data from KO optained from KEGG Database
file_p=cdp+'/kegg_annotation/ko.txt'
fko=open(file_p,"r")
lko=[]
ko=fko.readline().rstrip()
while ko!="":
	lko.append(ko)
	ko=fko.readline().rstrip()
fko.close()
koc=len(lko)
####################
file_p=cdp+'/kegg_annotation/gene_ko.txt'
fgeko=open(file_p,"r")
lge=[]
lgeko=[]
geko=fgeko.readline().rstrip()
while geko!="":
	h=geko.find('\t')
	lge.append(int(geko[4:h]))
	lgeko.append(geko[h+1:len(geko)])
	geko=fgeko.readline().rstrip()
fgeko.close()
####################
file_p=cdp+'/all_results/normal_matrix_gene.txt'
mos_mat=numpy.loadtxt(file_p, dtype=float)
file_p=cdp+'/files/gene_row.txt'
gene_row=numpy.loadtxt(file_p, dtype=float)
####################
if ca=='p':
	ko_mat=numpy.zeros((koc,int(snc/2)),dtype=float)
	sncc=int(snc/2)
else:
	ko_mat=numpy.zeros((koc,snc),dtype=float)
	sncc=snc

for i in range(0,koc):
	h=[index for index, value in enumerate(lgeko) if value == lko[i]]
	gene=[]
	for j in range(0,len(h)):
		gene.append(lge[h[j]])
	#glen=0
	#for j in range(0,len(gene)):
	#	glen+=int(ge_le[gene[j]-1])
	for j in range(0,sncc):
		sumg=0
		for k in range(0,len(gene)):
			hh = numpy.where(gene_row == gene[k])
			##print(hh[0])
			if hh[0].size>0:
				hhh=hh[0][0]
				sumg+=mos_mat[hhh][j]
		ko_mat[i][j]=sumg
###############################
#this code open a file for normalized value of KO
file_p=cdp+'/all_results/normal_matrix_ko.txt'
fko=open(file_p,"w")
si=numpy.shape(ko_mat)
i=0
while i<si[0]:
	j=0
	while j<si[1]:
		fko.write('%20s'%str(ko_mat[i][j]))
		if j<si[1]-1:
			fko.write('\t')	
		j+=1
	fko.write('\n')	
	i+=1
fko.close()
del lko
del lge 
del lgeko
del ko_mat
print "The normal values of the KOs were extracted for samples and stored in the /all_results/normal_matrix_ko.txt file."
print "In this file, the rows represent the KO and the columns represent the samples."
print ""
##############################################
#this code extract nornalized data from EC number optained from KEGG Database
file_p=cdp+'/kegg_annotation/def_ec.txt'
fec=open(file_p,"r")
lec=[]
ec=fec.readline().rstrip()
while ec!="":
	lec.append(ec)
	ec=fec.readline().rstrip()
fec.close()
ecc=len(lec)
####################
file_p=cdp+'/kegg_annotation/gene_ec.txt'
fgeec=open(file_p,"r")
lge=[]
lgeec=[]
geec=fgeec.readline().rstrip()
while geec!="":
	h=geec.find('\t')
	lge.append(int(geec[4:h]))
	lgeec.append(geec[h+1:len(geec)])
	geec=fgeec.readline().rstrip()
fgeec.close()
####################
if ca=='p':
	ec_mat=numpy.zeros((ecc,int(snc/2)),dtype=float)
	sncc=int(snc/2)
else:
	ec_mat=numpy.zeros((ecc,snc),dtype=float)
	sncc=snc

for i in range(0,ecc):
	h=[index for index, value in enumerate(lgeec) if value == lec[i]]
	gene=[]
	for j in range(0,len(h)):
		gene.append(lge[h[j]])
	for j in range(0,sncc):
		sumg=0
		for k in range(0,len(gene)):
			hh = numpy.where(gene_row == gene[k])
			if hh[0].size>0:
				hhh=hh[0][0]
				sumg+=mos_mat[hhh][j]
		ec_mat[i][j]=sumg
###############################
#this code open a file for normalized value of EC number
file_p=cdp+'/all_results/normal_matrix_ec.txt'
fec=open(file_p,"w")
si=numpy.shape(ec_mat)
i=0
while i<si[0]:
	j=0
	while j<si[1]:
		fec.write('%20s'%str(ec_mat[i][j]))
		if j<si[1]-1:
			fec.write('\t')	
		j+=1
	fec.write('\n')	
	i+=1
fec.close()
del lec
del lge 
del lgeec
del ec_mat
print "The normal values of the EC numbers were extracted for samples and stored in the /all_results/normal_matrix_ec.txt file."
print "In this file, the rows represent the EC number and the columns represent the samples."
print ""
##############################################
#this code extract nornalized data from reactions optained from KEGG Database
file_p=cdp+'/kegg_annotation/def_re.txt'
fre=open(file_p,"r")
lre=[]
re1=fre.readline().rstrip()
while re1!="":
	lre.append(re1)
	re1=fre.readline().rstrip()
fre.close()
rec=len(lre)
#print rec
####################
file_p=cdp+'/kegg_annotation/gene_re.txt'
fgere=open(file_p,"r")
lge=[]
lgere=[]
gere=fgere.readline().rstrip()
while gere!="":
	h=gere.find('\t')
	lge.append(int(gere[4:h]))
	lgere.append(gere[h+1:len(gere)])
	gere=fgere.readline().rstrip()
fgere.close()
####################
if ca=='p':
	re_mat=numpy.zeros((rec,int(snc/2)),dtype=float)
	sncc=int(snc/2)
else:
	re_mat=numpy.zeros((rec,snc),dtype=float)
	sncc=snc

for i in range(0,rec):
	h=[index for index, value in enumerate(lgere) if value == lre[i]]
	gene=[]
	for j in range(0,len(h)):
		gene.append(lge[h[j]])
	for j in range(0,sncc):
		sumg=0
		for k in range(0,len(gene)):
			hh = numpy.where(gene_row == gene[k])
			if hh[0].size>0:
				hhh=hh[0][0]
				sumg+=mos_mat[hhh][j]
		re_mat[i][j]=sumg
###############################
#this code open a file for normalized value of reactions
file_p=cdp+'/all_results/normal_matrix_re.txt'
fre=open(file_p,"w")
si=numpy.shape(re_mat)
i=0
while i<si[0]:
	j=0
	while j<si[1]:
		fre.write('%20s'%str(re_mat[i][j]))
		if j<si[1]-1:
			fre.write('\t')	
		j+=1
	fre.write('\n')	
	i+=1
fre.close()
del lre
del lge 
del lgere
del re_mat
print "The normal values of reactions were extracted for samples and stored in the /all_results/normal_matrix_re.txt file."
print "In this file, the rows represent the reactions and the columns represent the samples."
print ""

print "Please press any key to continue"

c=raw_input()












#This function extract annotated information that has already been extracted from the KEGG database
import os
os.system('clear')
cdp=os.getcwd()
print "KEGG annotation" 
print 
print "We extracted all the annotated information related to the KOs and EC numbers and reactions to the date 2019/6/30"
print "from the KEGG database and stored in the 'kegg_annotation' folder in the text files that start with the def prefix."
print
print "At this step, we extract all the annotated information associated with the sequences step by step."
print "Do you want to continue [y/n]:"
s=raw_input()
if s=='y':
	print "In this step, extract KO-related EC numbers from local datasets"
	print
	print
	file_p=cdp+'/kegg_annotation/gene_ko.txt'
	fgeko=open(file_p,"r")
	file_p=cdp+'/kegg_annotation/def_ko_ec.txt'
	fkoec=open(file_p,"r")
	file_p=cdp+'/kegg_annotation/gene_ec.txt'
	fgeec=open(file_p,"w")
	file_p=cdp+'/files/catalog_length.txt'
	fca=open(file_p,"r")
	clen=int(fca.readline().rstrip('\n'))
	fca.close()
	import re
	import numpy
	ge_ko_mat=numpy.zeros((clen,10), dtype=numpy.uint16)
	lko=[]
	lec=[]
	koec=fkoec.readline().rstrip('\n')
	while koec !="":
		h1=koec.find('\t')
		lko.append(koec[0:h1])
		lec.append(koec[h1+1:len(koec)])
		koec=fkoec.readline().rstrip('\n')
	lkoi=numpy.zeros((len(lko)),dtype=numpy.uint16)
	j=0
	while j<len(lko):
		s=lko[j];
		lkoi[j]=int(s[1:len(s)])
		j=j+1
	geko=fgeko.readline().rstrip('\n')
	while geko !="":
		h1=geko.find('\t')
		ge=int(geko[4:h1])
		ko=int(geko[h1+2:len(geko)])
		j=0
		while ge_ko_mat[ge-1][j]>0:
			j=j+1
		ge_ko_mat[ge-1][j]=ko
		geko=fgeko.readline().rstrip('\n')
		i=0;
	while i<clen:
		ec=[]
		j=0;
		while ge_ko_mat[i][j]>0:
			h1=numpy.where(lkoi==ge_ko_mat[i][j])
			h2=h1[0]
			#print(h2)
			if len(h2)>0 and lec[h2[0]] not in ec: 
				ec.append(lec[h2[0]])
			j=j+1	
		if len(ec)>0:
			ec.sort()
			j=0
			while j<len(ec):
				fgeec.write('gene'+str(i+1)+'\t')
				fgeec.write(ec[j]+'\n')
				j=j+1
		i+=1
	del ge_ko_mat
	del lko
	del lec
	del lkoi			
	fgeko.close()
	fkoec.close()
	fgeec.close()
	print "EC numbers related information was successfully extracted."
	print "The results are visible in the kegg_annotation folder in the gene_ec.txt and def_ko_ec.txt and def_ec.txt files."
	print "Please press the key to continue"	
	c=raw_input()
	####################################################
	print "In this step, extract EC-related reactions from local datasets"
	print
	print
	file_p=cdp+'/kegg_annotation/gene_ec.txt'
	fgeec=open(file_p,"r")
	file_p=cdp+'/kegg_annotation/def_ec_re.txt'
	fecre=open(file_p,"r")
	file_p=cdp+'/kegg_annotation/gene_re.txt'
	fgere=open(file_p,"w")
	file_p=cdp+'/files/catalog_length.txt'
	fca=open(file_p,"r")
	clen=int(fca.readline().rstrip('\n'))
	fca.close()
	import re
	import numpy
	ge_ec_mat=numpy.zeros((clen,10), dtype=numpy.uint16)
	lec=[]
	lre=[]
	ecre=fecre.readline().rstrip('\n')
	while ecre !="":
		h1=ecre.find('\t')
		lec.append(ecre[0:h1])
		lre.append(ecre[h1+1:len(ecre)])
		ecre=fecre.readline().rstrip('\n')
	#lkoi=numpy.zeros((len(lko)),dtype=numpy.uint16)
	#j=0
	#while j<len(lko):
	#	s=lko[j];
	#	lkoi[j]=int(s[1:len(s)])
	#	j=j+1
	geec=fgeec.readline().rstrip('\n')
	while geec !="":
		h1=geec.find('\t')
		ge=int(geec[4:h1])
		ec=geec[h1+1:len(geec)-1]
		if ec in lec:
			h2=[index for index, value in enumerate(lec) if value == ec]
			#print h2
			#h2=lec.index(ec)
			#h2=lec.find(ec)
			if len(h2)>0:
				j=0
				f=0
				y=0
				while y<len(h2):
					while ge_ec_mat[ge-1][j]>0:
						if ge_ec_mat[ge-1][j]-1==h2[y]:
							f=1
						j=j+1
					if f==0:
						ge_ec_mat[ge-1][j]=h2[y]+1
					y+=1
		geec=fgeec.readline().rstrip('\n')
	#print ge_ec_mat
	i=0;
	while i<clen:
		re1=[]
		j=0;
		while ge_ec_mat[i][j]>0:
			re1.append(lre[ge_ec_mat[i][j]-1])
			j=j+1	
		if len(re1)>0:
			re1.sort()
			j=0
			while j<len(re1):
				fgere.write('gene'+str(i+1)+'\t')
				fgere.write(re1[j]+'\n')
				j=j+1
		i+=1
	del ge_ec_mat
	del lec
	del lre
	fgeec.close()
	fecre.close()
	fgere.close()
	print "Reactions related information was successfully extracted."
	print "The results are visible in the kegg_annotation folder in the gene_re.txt and def_ec_re.txt and def_re.txt files."
	print "Please press the key to continue"
	c=raw_input()
###############################
	print 
	print
	print "Finally, The reaction-related equations are visible in the kegg_annotation folder in the def_re_eq.txt file."
	print "Please press the key to continue"
	c=raw_input()
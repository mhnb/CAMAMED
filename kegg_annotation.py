#This function extract annotated information from KEGG databases
def ko_read(text):
	url='http://www.genome.jp/dbget-bin/www_bfind_sub?mode=bfind&max_hit=1000&locale=en&serv=kegg&dbkey=orthology&keywords='+text
	#print(url)
	f = urllib.urlopen(url)
	myfile = f.read()
	#print(myfile)
	h1=myfile.find('style="margin-left:2em"> ')
	h2=myfile.find('[EC:')
	h3=h2
	if h3>-1:
		while myfile[h3]!=']':
			h3=h3+1
		st1=myfile[h2+4:h3]
	else:
		st1=""
	return st1
######################################end of ko_read() function
def ce_read(text):
	url='http://www.genome.jp/dbget-bin/www_bget?ec:'+text
	#print(url)
	f = urllib.urlopen(url)
	myfile = f.read()
	h1=myfile.find('Reaction(IUBMB)')
	h2=myfile.find('Reaction(KEGG)')
	st=myfile[h1:h2]
	return st
####################################end of ec_read() function
while 1:
	import os
	os.system('clear')
	print "KEGG annotation" 
        print 
	print "We extracted all the annotated information related to the KOs and EC numbers and reactions to the date 2019/6/30"
	print "from the KEGG database and stored in the 'kegg_annotation' folder in the text files that start with the def prefix."
	print "But if you want to extract this information yourself, you can use these options."
	print
	print "****These steps may take a long time****"
	print
	print
	print "Enter e41 to extract KO-related EC numbers from the KEGG database."
	print "Enter e42 to extract EC-related reactions from the KEGG database."
	print "Enter e43 to extract reaction-related equation from the KEGG database."
	print "Enter 0 for back"
	x=raw_input()
	#x=int(c)
	if x=='e41':
		import os
		os.system('clear')
		print "This step, extract KO-related EC numbers from the KEGG database."
		print "To do this, have a file named ko.txt in the kegg_annotation folder. This file contains ko obtained for the sequence of genes from KEGG databases."	
		print ""
		print "Extract this KO-related EC numbers"
		import urllib
		cdp=os.getcwd()
		file_p=cdp+'/kegg_annotation/ko.txt'
		fko=open(file_p,"r")
		file_p=cdp+'/kegg_annotation/ko_ec.txt'
		fkoec=open(file_p,"w")
		file_p=cdp+'/kegg_annotation/ec.txt'
		fec=open(file_p,"w")
		lec=[]
		while 1:
			text=fko.readline().rstrip('\n')
			if text=="":
				break
			else:
				print(text)	
				st1=ko_read(text)	
				#print(st1)
				if st1!="" and (' ' in st1)==False: 
					if ('-' in st1)==False:
						fkoec.write(text[0:6])
						fkoec.write('\t')	
						fkoec.write(st1)
						fkoec.write('\n')
						if st1 not in lec: 
							lec.append(st1)				
				elif st1!="" and (' ' in st1)==True:
					while st1!="":
						h1=st1.find(' ')
						if h1!=-1:
							st2=st1[0:h1]
						else:
							st2=st1
							st1=""
						#print(st2)
						st1=st1[h1+1:len(st1)]
						#print(st1)
						#input('ver=')
						if ('-' in st2)==False:
							fkoec.write(text[0:6])
							fkoec.write('\t')
							fkoec.write(st2)
							fkoec.write('\n')
							if st2 not in lec: 
								lec.append(st2)
		fko.close()
		fkoec.close()
		lec.sort()
		i=0
		while i<len(lec):
			fec.write(lec[i])
			fec.write('\n')
			i+=1
		fec.close()
		#####################################################
		file_p=cdp+'/kegg_annotation/gene_ko.txt'
		fgeko=open(file_p,"r")
		file_p=cdp+'/kegg_annotation/ko_ec.txt'
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

		print "This information was successfully extracted."
		print "The results are visible in the kegg_annotation folder in the gene_ec.txt and ko_ec.txt and ec.txt files"
		print "Please press the key to continue"	
		c=raw_input()
	elif x=='e42':
		import os
		os.system('clear')
		print "This step, optain EC-related reactions from the KEGG database."
		print "To do this, have a file named ec.txt in the kegg_annotation folder. This file contains EC numbers extracted for the KO numbers from KEGG databases."	
		print ""
		print "Extract this EC numbers-related reactions"
		import urllib
		cdp=os.getcwd()
		file_p=cdp+'/kegg_annotation/ec.txt'
		fec=open(file_p,"r")
		file_p=cdp+'/kegg_annotation/ec_re.txt'
		fecre=open(file_p,"w")
		file_p=cdp+'/kegg_annotation/re.txt'
		fre=open(file_p,"w")
		lre=[]
		while 1:
			text=fec.readline().rstrip('\n')
			if text=="":
				break
			else:
				print(text)
				st=ce_read(text)		
				h3=st.find('href="/dbget-bin/www_bget?rn:')
				while h3>-1:
		        		fecre.write(text)
					fecre.write('\t')	
					#print(st[h3+29:h3+35])
					fecre.write(st[h3+29:h3+35])
					fecre.write('\n')
					if st[h3+29:h3+35] not in lre: 
						lre.append(st[h3+29:h3+35])
					st=st[h3+35:len(st)]	
					h3=st.find('href="/dbget-bin/www_bget?rn:')
		fec.close()
		fecre.close()
		lre.sort()
		i=0
		while i<len(lre):
			fre.write(lre[i])
			fre.write('\n')
			i+=1
		fre.close()
		#####################################################
		file_p=cdp+'/kegg_annotation/gene_ec.txt'
		fgeec=open(file_p,"r")
		file_p=cdp+'/kegg_annotation/ec_re.txt'
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
			ec=geec[h1+1:len(geec)]
			if ec in lec:
				h2=[index for index, value in enumerate(lec) if value == ec]
				#h2=lec.index(ec)
				#h2=lec.find(ec)
				#print(h2)
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

		print "This information was successfully extracted."
		print "The results are visible in the kegg_annotation folder in the gene_re.txt and ec_re.txt and re.txt files"
		print "Please press the key to continue"
		c=raw_input()
	elif x=='e43':
		import os
		os.system('clear')
		print "This step, extract reaction-related equations from the KEGG database."
		print "To do this, have a file named re.txt in the kegg_annotation folder. This file contains reaction numbers extracted from KEGG databases."	
		print ""
		print "Extract this reactions-related equation"
		import urllib
		cdp=os.getcwd()
		file_p=cdp+'/kegg_annotation/re.txt'
		fre=open(file_p,"r")
		file_p=cdp+'/kegg_annotation/re_eq.txt'
		freeq=open(file_p,"w")
		while 1:
			text=fre.readline().rstrip('\n')
			if text=="":
				break
			else:
				print(text)		
				url='http://www.genome.jp/dbget-bin/www_bget?rn:'+text
				#print(url)
				f = urllib.urlopen(url)
				myfile = f.read()
				#print(myfile)
				h1=myfile.find('Definition')
				h1=h1+236
				h2=h1
			        while myfile[h2+1]!='<':
					h2=h2+1;
				freeq.write('Reaction:'+text+'\t')
				s=myfile[h1:h2+1].replace("&lt;=&gt;", "<=>")
				freeq.write('Definition:'+s+'\t')
				h1=myfile.find('Equation')
				h1=h1+175
				import re
				h3=[m.start() for m in re.finditer('<br>', myfile)]
				j=0
				while h3[j]<h1:
					j=j+1
				h2=h3[j]-1
				st=myfile[h1:h2].replace("&lt;=&gt;", "<=>")
				#s=st.replace('<a href="/dbget-bin/www_bget?cpd:', '')
				#print(st)
				freeq.write('Equation:')
				h1=[m.start() for m in re.finditer('<a href="/dbget-bin/www_bget?',st)]
				#print(h1)
				h3=st.find('<=>')
				j=0
				while h1[j]<h3:
					if j==0:
						ss=st[h1[j]+41:h1[j]+47]+'+'
					else:
						ss=ss+st[h1[j]+41:h1[j]+47]+'+'
					j=j+1	
				ss=ss+"<=>"
				while j<len(h1)-1:
					ss=ss+st[h1[j]+41:h1[j]+47]+'+'
					j=j+1
				ss=ss+st[h1[j]+41:h1[j]+47]
				s=ss.replace("+<=>", " <=> ")
				#print(s)
				freeq.write(s)
				freeq.write('\n')
		fre.close()
		freeq.close()
		print "This information was successfully extracted."
		print "The results are visible in the kegg_annotation folder in the re_eq.txt file."
		print "Please press the key to continue"
		c=raw_input()
	elif x=='0':
		break

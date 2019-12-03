#this function is a preprocessing step in pipeline
import os
cdp=os.getcwd()
while 1:
	import os
	os.system('clear')
	print "This function is a preprocessing step in pipeline"
	print
	print "Enter a1 for running SAR toolkit to convert SRA files to Fastq or Fasta files.  (**optional**)"
	print "Enter a2 to receive information on sequences and gene catalogs."
	print "Enter a3 for running CD-HIT on the gene catalog to remove redundant genes. (**optional**)"
	print "Enter a4 preprocessing gene catalog."
	print "Enter 0 for Back"
	x=raw_input()
	#x=int(c)
	if x=='a1':
		import subprocess
		subprocess.call("./sra_samlpe.sh")
	elif x=='a2':
		os.system('clear')
		print "You must first select the gene catalog. Gene sequences must have Fastaasta format."
		print "First copy the gene catalog into the metagenomic_data_analysis folder and then enter the filename." 
		print "For example, gene_catalog.fa:"
		s=raw_input()
		filename, file_extension = os.path.splitext(s)
		while file_extension!='.fa' and file_extension!='.fasta':
			print "The selected file type is not correct"
			print "select .fa or .fasta files"
			s=raw_input()
			filename, file_extension = os.path.splitext(s)
		file_p=cdp+'/files/def_catalog_name.txt'
		cot_na=open(file_p,"w")
		cot_na.write(s)		
		cot_na.close()
		os.system('clear')
		print "At this step, you must first enter the name of the files for the samples in the 'sample_file_names.txt' text file as below."
		print "Then copy the sample file according to its name in the 'Read_files' folder."
		print "single end fastq files"
		print "file1.fastq"
		print "file2.fastq"
		print "file3.fastq"
		print
		print "single end fasta files"
		print "file1.fasta"
		print "file2.fasta"
		print "file3.fasta"
		print
		print "paired end fasta files"
		print "file1_1.fasta"
		print "file1_2.fasta"
		print "file2_1.fasta"
		print "file2_2.fasta"
		print
		print "paired end fastq files"
		print "file1_1.fastq"
		print "file1_2.fastq"
		print "file2_1.fastq"
		print "file2_2.fastq"
		file_p=cdp+'/files/catalog_align.txt'
		cot_al=open(file_p,"w")
		print "select type of read files[fastq/fasta]= "
		q=raw_input()
		#cot_al.write('q=')
		while q!='fastq' and q!='fasta':
			print "The selected file type is not correct"
			print "select type of read files[fastq/fasta]= "
			q=raw_input()
		if q=='fastq':
			cot_al.write('q')
		elif q=='fasta':
			cot_al.write('a')
		cot_al.write('\n')
		print "select type of read files- paired end or single end [p/s]= "
		p=raw_input()
		while p!='p' and p!='s':
			print "The selected read type is not correct"
			print "select type of read files paired-end or single-end [p/s]= "
			p=raw_input()
		cot_al.write(p)
		cot_al.write('\n')
		print "select number of cores= "
		co=raw_input()
		cot_al.write(co)
		cot_al.write('\n')
		if p=='p':
			print "For paired end data, enter the median fragment length (insert size an integer number) or ignor enter -1 ="
			s=raw_input()
			xx=int(s)
			while xx<-1:
				print "Insert size should be an integer greater than or equal to -1, please re-enter"
				s=raw_input()
				xx=int(s)
			#cot_al.write('s=')
			cot_al.write(str(xx))
			cot_al.write('\n')
		cot_al.close()
		print "This information was successfully saved"
		print "Please click a key to continue"
		s=raw_input()
	elif x=='a3':
		print "You can use this option if you think your gene catalog has redundant sequences."
		print "Do you want to continue [y/n]:"
		s=raw_input()
		while s!='y' and s!='n':
			print "You can use this option if you think your gene catalog has redundant sequences."
			print "Do you want to continue [y/n]:"
			s=raw_input()
		if s=='y':
			print "Select the sequence identity threshold. This value can be in the range of 0.8 to 1 and its default value is 0.9="
			s=raw_input()
			c=float(s)
			#print c
			#print(int())
			while c<0.8 or c>1:
				print "Select the sequence identity threshold in renge of [0.8:1]="
				s=raw_input()
				c=float(s)
			file_p=cdp+'/files/def_catalog_name.txt'
			mf3=open(file_p,"r")
			ca_na=mf3.readline()
			mf3.close()
			bashcommand='cdhit-est -i '+ca_na+' -o cd_hit_gene_catalog -c '+s+' -T 0 -M 0'
			import subprocess
			#subprocess.call(bashcommand)
			subprocess.check_call(bashcommand, shell=True)
			file_p=cdp+'/files/def_catalog_name.txt'
			mf3=open(file_p,"w")
			mf3.write('cd_hit_gene_catalog')
			mf3.close()
			print "CD-HIT runs on the main gene catalog."
			print "After the duplicate genes have been deleted, the new gene catalog is stored with the name cd_hit_gene_catalog."
			print "Also, clustered genes are stored in a file named cd_hit_gene_catalog.clstr, and the genes of the head cluster are marked with *."
			print "Please click a key to continue"
			s=raw_input()
	elif x=='a4':
		print "At this point, the names of the genes are deleted from the gene catalog and stored in the gene_name.txt file."
		print "and for the genes the gene1, gene2 and ... are respectively selected."
		print
		print "Please wait..."
		file_p=cdp+'/files/def_catalog_name.txt'
		mf3=open(file_p,"r")
		ca_na=mf3.readline()
		#print ca_na
		mf3.close()
		mf1=open(ca_na,"r")
		#print mf1
		file_p=cdp+'/files/gene_length.txt'
		gene=open(file_p,"w")
		i=0
		s=0;
		text=mf1.readline()
		while 1:
			if text=="":
				break
			else:
				text=mf1.readline()
				l=0
				while text!="" and text[0]!='>':
					l=l+len(text)
					text=mf1.readline()
				gene.write(str(l))
				gene.write('\n')
				i=i+1
				s=s+l;
		s=s/i
		mf1.close()
		gene.close()
		file_p=cdp+'/files/catalog_length.txt'
		ca_le=open(file_p,"w")
		ca_le.write(str(i))
		ca_le.close()
		#mf1=open("catalog_name.txt","r")
		#ca_na=mf1.readline()
		#print ca_na
		#mf1.close()
		mainf=open(ca_na,"r")
		mainwf=open("main_gene_catalog.fa","w")
		file_p=cdp+'/files/gene_name.txt'
		gene=open(file_p,"w")
		i=1
		j=0
		while 1:
			#if i%1000==0:
			#	print(i)
			text=mainf.readline()
			if text=="":
				break
			if text[0]=='>':
				gene.write('>gene'+str(i)+'\n')
				gene.write(text)
				mainwf.write('>gene'+str(i)+'\n')
				i=i+1
			else:
				ix=0
				l=len(text)
				while ix < l-1:
					if text[ix]=='a':
						t=text[0:ix]+'A'+text[ix+1:l]
						text=t
					elif text[ix]=='t':
						t=text[0:ix]+'T'+text[ix+1:l]
						text=t
					elif text[ix]=='c':
						t=text[0:ix]+'C'+text[ix+1:l]
						text=t
					elif text[ix]=='g':
						t=text[0:ix]+'G'+text[ix+1:l]
						text=t
					elif text[ix]=='u':
						t=text[0:ix]+'T'+text[ix+1:l]
						text=t
					elif text[ix]=='U':
						t=text[0:ix]+'T'+text[ix+1:l]
						text=t
					elif text[ix]!='A' and text[ix]!='T' and text[ix]!='C' and text[ix]!='G':
						#print text[ix]
						t=text[0:ix]+'N'+text[ix+1:l]
						text=t
						j=j+1	
					ix=ix+1
				mainwf.write(text)
		mainwf.close()
		gene.close()	  
		mainf.close()
		file_p=cdp+'/files/catalog_name.txt'
		mf2=open(file_p,"w")
		mf2.write('main_gene_catalog.fa')
		mf2.close()
		#import os
		#os.remove(ca_na)
		#import os
		#os.rename('IGCnew.fa',ca_na)
		print
		print "The preprocessing state has been completed on gene catalog and last gene catalog name is main_gene_catalog.fa" 
		print ("There are "+str(i-1)+" genes in this catalog, The length of genes is "+str(s)+"bp on average")
		print "Please click a key to continue"
		s=raw_input()
		
	elif x=='0':
		break
#!/usr/bin/env python
#this function is a preprocessing step in pipeline
import sys
import os
cdp=os.getcwd()
#################################################
def help_fun():
	print('''
	This function is a preprocessing step in pipeline

	Command: ./camamed_pre_processing.py method [options]

	      method:
		  -h	Shows help related to this function
		  -sra  Running SAR toolkit to convert SRA files to Fastq or Fasta files.  (**optional**)
			Copy SRA samples in the /CAMAMED/sra_files folder.
			Copy the sample names in the /CAMAMED/sra_files/sra_file_names.txt.
			for example:
				file1.sra
				file2.sra
				file3.sra
			After executing SRA-Toolkit, Fastq or Fasta files are automatically copied to folder /Read_files.
			For Example: ./camamed_pre_processing.py -sra		

		  -gec  Get information on sequences and gene catalogs.
			options:
			   -gc (gene_catalog): Gene catalog sequences must have Fasta format.
				       	  First copy the gene catalog into the /CAMAMED folder.
			   -rff (read_files_format)[fastq/fasta]:
					  Copy the sample files in the /CAMAMED/Read_files folder.
					  Enter the file name of the samples in the /CAMAMED/Read_files/sample_file_names.txt
						For single end fastq files:
						    file1.fastq
						    file2.fastq
						    file3.fastq
						For paired end fastq files:
						    file1_1.fastq
						    file1_2.fastq
						    file2_1.fastq
						    file2_2.fastq
			   -rt (read_type): paired end or single end [p/s]	
			   -is (insert_size): For paired end sequences (An integer number) or ignor enter -1 (default=-1).		
			For Example: ./camamed_pre_processing.py -gec -gc gene_catalog.fa -rff fastq -rt p
				     ./camamed_pre_processing.py -gec -gc gene_catalog.fa -rff fastq -rt p -is 350

		  -cdh  Running CD-HIT on the gene catalog to remove redundant genes. (**optional**)
			You can use this option if you think your gene catalog has redundant sequences.
			options:
			   -sit (sequence_identity_threshold): This value can be in the range of 0.8 to 1 (default=0.9).
			For Example: ./camamed_pre_processing.py -cdh
				     ./camamed_pre_processing.py -cdh -sit 0.95

		  -pgc  Peprocessing gene catalog.
			At this point, the names of the genes are deleted from the gene catalog and stored in the /CAMAMED/files/gene_name.txt
			and for the genes the gene1, gene2 and ... are respectively selected.
			For Example: ./camamed_pre_processing.py -pgc
	''')
#################################################
if not os.path.exists('files'):
	os.mkdir('files')
if len(sys.argv)<2:
	help_fun()
elif str(sys.argv[1])=='-h':
	help_fun()	
elif str(sys.argv[1])=='-sra':
	import subprocess
	subprocess.call("./sra_samlpe.sh")
elif str(sys.argv[1])=='-gec':
	if len(sys.argv)<8:
		print("Not enough input")
	else:
		sf=0
		sk=0
		ss=sys.argv
		try:
			i = ss.index('-gc')
		except ValueError:
			i=-1
			print("syntax error")
		if i+1<len(sys.argv) and i!=-1:
			sk+=1
			s=str(sys.argv[i+1])
			filename, file_extension = os.path.splitext(s)
			if file_extension!='.fa' and file_extension!='.fasta':
				print("The selected file type is not correct")
				print("select .fa or .fasta files")
				sf=1
			else:
				file_p=cdp+'/files/def_catalog_name.txt'
				cot_na=open(file_p,"w")
				cot_na.write(s)		
				cot_na.close()
		try:
			i = ss.index('-rff')
		except ValueError:
			i=-1
			print("syntax error")
		if i+1<len(sys.argv) and i!=-1:
			sk+=1			
			q=str(sys.argv[i+1])
			if q!='fastq' and q!='fasta':
				print("The selected file type is not correct")
				print("select type of read files[fastq/fasta]")
				sf=1
		try:
			i = ss.index('-rt')
		except ValueError:
			i=-1
			print("syntax error")
		p='z'
		if i+1<len(sys.argv) and i!=-1:
			sk+=1
			p=str(sys.argv[i+1])
			if p!='p' and p!='s':
				print("The selected read type is not correct")
				print("select type of read files paired-end or single-end [p/s]")
				sf=1
		xx=-1
		if p=='p' and len(sys.argv)>8:
			try:
				i = ss.index('-is')
			except ValueError:
				i=-1
				print("syntax error")
				sk-=1
			if i+1<len(sys.argv) and i!=-1:
				xx=int(sys.argv[i+1])
				if xx<1:
					print("Insert size should be an integer greater than zero")
					sf=1
		if sf==0 and sk==3:
			file_p=cdp+'/files/catalog_align.txt'
			cot_al=open(file_p,"w")
			if q=='fastq':
				cot_al.write('q')
			elif q=='fasta':
				cot_al.write('a')
			cot_al.write('\n')
			cot_al.write(p)
			cot_al.write('\n')
			cot_al.write(str(xx))
			cot_al.write('\n')
			cot_al.close()
			print("This information was successfully saved")
elif str(sys.argv[1])=='-cdh':
	if os.path.isfile('cd_hit_gene_catalog'):
		os.remove('cd_hit_gene_catalog')
	if os.path.isfile('cd_hit_gene_catalog.clstr'):
		os.remove('cd_hit_gene_catalog.clstr')
	if len(sys.argv)<2:
		print("Not enough input")
	else:
		c=0.7
		if len(sys.argv)>2:
			ss=sys.argv
			try:
				i = ss.index('-sit')
			except ValueError:
				i=-1
				print("syntax error")
			if i+1<len(sys.argv) and i!=-1:
				c=float(sys.argv[i+1])
			if c<0.8 or c>1:
				print("Select the sequence identity threshold in renge of [0.8:1]")
		else:
			c=0.9
		if c>=0.8 and c<=1:
			file_p=cdp+'/files/def_catalog_name.txt'
			mf3=open(file_p,"r")
			ca_na=mf3.readline()
			mf3.close()
			bashcommand='cdhit-est -i '+ca_na+' -o cd_hit_gene_catalog -c '+str(c)+' -T 0 -M 0'
			import subprocess
			subprocess.check_call(bashcommand, shell=True)
			file_p=cdp+'/files/def_catalog_name.txt'
			mf3=open(file_p,"w")
			mf3.write('cd_hit_gene_catalog')
			mf3.close()
			print("CD-HIT runs on the main gene catalog.")
			print("After the duplicate genes have been deleted, the new gene catalog is stored with the name cd_hit_gene_catalog in /CAMAMED.")
			print("Also, clustered genes are stored in a file named cd_hit_gene_catalog.clstr, and the genes of the head cluster are marked with *.")
elif str(sys.argv[1])=='-pgc':
	print("Please wait...")
	file_p=cdp+'/files/def_catalog_name.txt'
	mf3=open(file_p,"r")
	ca_na=mf3.readline()
	mf3.close()
	mf1=open(ca_na,"r")
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
	s=round(s/i)
	mf1.close()
	gene.close()
	file_p=cdp+'/files/catalog_length.txt'
	ca_le=open(file_p,"w")
	ca_le.write(str(i))
	ca_le.close()
	mainf=open(ca_na,"r")
	mainwf=open("main_gene_catalog.fa","w")
	file_p=cdp+'/files/gene_name.txt'
	gene=open(file_p,"w")
	i=1
	j=0
	while 1:
		text=mainf.readline()
		if text=="":
			break
		if text[0]=='>':
			if (i%100000)==0:
				print('Processed up to gene'+str(i))
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
	print()
	print("The preprocessing state has been completed on gene catalog and last gene catalog name is main_gene_catalog.fa")
	print("There are "+str(i-1)+" genes in this catalog, The length of genes is "+str(s)+"bp on average")
else: 
	help_fun()
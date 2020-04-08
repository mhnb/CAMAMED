#!/usr/bin/env python
#This function executes fastqc quality control and SeqKit program on sample sequences
import sys
import os
cdp=os.getcwd()
#################################################
def help_fun():
	print('''
	This function executes fastqc quality control and SeqKit tools on sample sequences.
	At this state, if the sequences have Fastq format, they will be quality controlled and their statistical information extracted.
	But if they have Fasta format, only their statistical information obtained.
	Before running this step, the sample data should be in the /CAMAMED/Read_files folder
	and the file names should be written in the text file /CAMAMED/Read_files/sample_file_names.txt, respectively.

	Command: ./camamed_quality_control.py method

	      method:
		  -h	Shows help related to this function.
		  -fqc  Execut FastQC quality control only for Fastq files.
			For Example: ./camamed_quality_control.py -fqc		

		  -sek  Execut SeqKit to extract information from sample files.
			For Example: ./camamed_quality_control.py -sek

		  -eti  Extract total information from SeqKit outputs.
			For Example: ./camamed_quality_control.py -eti
	''')
#################################################
if not os.path.exists('fastqc_output'):
	os.mkdir('fastqc_output')
if not os.path.exists('seqkit_output'):
	os.mkdir('seqkit_output')
if not os.path.exists('all_results'):
	os.mkdir('all_results')
file_p=cdp+'/files/catalog_align.txt'
fid=open(file_p,"r")
data_ty=fid.readline().rstrip('\n')
fid.close()
if len(sys.argv)<2:
	help_fun()
elif len(sys.argv)>2:
	print("syntax error")
	help_fun()
elif str(sys.argv[1])=='-h':
	help_fun()
elif str(sys.argv[1])=='-fqc':
	if data_ty=='a':
		print("Because the data format is Fasta you can only select -sek and -eti methods")
	else:
		os.system('clear')
		print("At this state FastQC software runs on samples")
		import subprocess
		subprocess.call("./fastqc_samlpe.sh")
		print()
		print("The FastQC software was run on samples and its results are visible in the /CAMAMED/fastqc_output folder.")
elif str(sys.argv[1])=='-sek':
	os.system('clear')
	print("In this step, SeqKit software runs on samples and statistical information of the sample files is extracted.")
	import subprocess
	subprocess.call("./seqkit_samlpe.sh")
	print()
	print("The SeqKit software was run on samples and its results are visible in the /CAMAMED/seqkit_output folder.")
elif str(sys.argv[1])=='-eti':
	cdp=os.getcwd()
	file_p=cdp+'/Read_files/'+'sample_file_names.txt'
	sam_na=open(file_p,"r")
	file_p=cdp+'/all_results/'+'total_sample_info.txt'
	ffq=open(file_p,"w")
	file_p=cdp+'/files/sample_reads.txt'
	fsr=open(file_p,"w")
	ffq.write('Total sample information Results\n')
	ffq.write('file              format  type  num_seqs    sum_len  min_len  avg_len  max_len\n')
	while 1:
		san=sam_na.readline().rstrip()
		if san=="":
			break
		else:
			file_p=cdp+'/seqkit_output/'+san
			fkit=open(file_p,"r")
			lkit=fkit.readline().rstrip()
			lkit=fkit.readline().rstrip()
			st=cdp+'/Read_files/'
			l=len(st)
			lkit2=lkit[l:]
			lkit3=lkit2.replace(",", "")
			a=[int(s) for s in lkit3.split() if s.isdigit()]
			fsr.write(str(a[0]))
			fsr.write('\n')
			ffq.write(lkit2)
			ffq.write('\n')
			fkit.close()
	print()
	print("This information was successfully extracted.")
	print("The results are visible in the /CAMAMED/all_results/total_sample_info.txt file")
	ffq.close()	
	fsr.close()
	sam_na.close()	
else:
	help_fun()

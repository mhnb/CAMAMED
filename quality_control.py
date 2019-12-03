#This function executes fastqc quality control program on sample data
import os
cdp=os.getcwd()
file_p=cdp+'/files/catalog_align.txt'
fid=open(file_p,"r")
data_ty=fid.readline().rstrip('\n')
fid.close()
while 1:
	import os
	os.system('clear')
	print "At this state, if the sequences have Fastq format, they will be quality controlled and their statistical information extracted."
	print "But if they have Fasta format, only their statistical information obtained."
	print "Before running this step, the sample data should be in the 'Read_files' folder" 
	print "and the file names should be written in the text file 'sample_file_names.txt', respectively."
	print
	print "Enter b1 to execut FastQC quality control only for Fastq files."
	print "Enter b2 to execut SeqKit to extract information from sample files."
	print "Enter b3 to extract total information from SeqKit outputs."
	print "Enter 0 to back"
	x=raw_input()
	#x=int(c)
	while data_ty=='a' and x=='b1':
		print "Because the data format is Fasta you can only select options b2 and b3"
		print
		print "Enter b1 for executes FastQC quality control only for Fastq files."
		print "Enter b2 for executes SeqKit to extract information from sample files."
		print "Enter b3 for extract total information from SeqKit outputs."
		print "Enter 0 for back"
		x=raw_input()
	if x=='b1':
		import os
		os.system('clear')
		print "At this state FastQC software runs on samples"
		import subprocess
		subprocess.call("./fastqc_samlpe.sh")
		print
		print "The FastQC software was run on samples and its results are visible in the fastqc_output folder."
		print "Please press the key to continue"					
		c=raw_input()
	elif x=='b2':
		import os
		os.system('clear')
		print "In this step, SeqKit software runs on samples and statistical information of the sample files is extracted."
		import subprocess
		subprocess.call("./seqkit_samlpe.sh")
		print
		print "The SeqKit software was run on samples and its results are visible in the seqkit_output folder."
		print "Please press the key to continue"					
		c=raw_input()
	elif x=='b3':
		sam_na=open("sample_file_names.txt","r")
		cdp=os.getcwd()
		file_p=cdp+'/all_results/'+'total_sample_info.txt'
		ffq=open(file_p,"w")
		file_p=cdp+'/files/sample_reads.txt'
		fsr=open(file_p,"w")
		ffq.write('Total sample information Results\n')
		ffq.write('file              format  type  num_seqs    sum_len  min_len  avg_len  max_len\n')
		while 1:
			san=sam_na.readline().rstrip()
			#filename, file_extension = os.path.splitext(san)
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
				#os.remove(file_p)
		print
		print "This information was successfully extracted."
		print "The results are visible in the all_results folder in the total_sample_info.txt files"
		print "Please press the key to continue"					
		ffq.close()	
		fsr.close()
		sam_na.close()	
		c=raw_input()
	elif x=='0':
		break

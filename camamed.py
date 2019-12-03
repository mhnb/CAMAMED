#This program is an open source program used for analyzing metagenome data
#This program takes sample reads files as well as the reference catalog of the genes and performs analyzes on them. 
#Meanwhile, the program uses the KEGG database for data analysis.
import subprocess
subprocess.check_call('chmod +x MosaikAligner', shell=True)
subprocess.check_call('chmod +x MosaikBuild', shell=True)
subprocess.check_call('chmod +x MosaikJump', shell=True)
subprocess.check_call('chmod +x MosaikText', shell=True)
subprocess.check_call('chmod +x seqkit', shell=True)
subprocess.check_call('chmod +x ./FastQC/fastqc', shell=True)
subprocess.check_call('chmod +x fastqc_samlpe.sh', shell=True)
subprocess.check_call('chmod +x metaphlan_samlpe.sh', shell=True)
subprocess.check_call('chmod +x mosaik_build_ref.sh', shell=True)
subprocess.check_call('chmod +x mosaik_read_aligner.sh', shell=True)
subprocess.check_call('chmod +x seqkit_samlpe.sh', shell=True)
subprocess.check_call('chmod +x sra_samlpe.sh', shell=True)


while 1:
	import os
	os.system('clear')
	print "This program is an open source program used for analyzing metagenome data"
	print "This program takes sample reads files as well as the reference catalog of the genes and performs analyzes on them. "
	print "Meanwhile, the program uses the kegg database for data analysis." 
	print
	print "Enter a for Preprocessing"   
	print "Enter b for FastQC quality control" 
	print "Enter c for Extracting bacterial frequency using MetaPhlan"   
	print "Enter d for Mapping Reads to reference gene catalog using MOSAIK."  
	print "Enter e for Extracting samples information using the user's gene catalog and KEGG database."
	print "Enter f for Extracting samples information from normalized bacteria, gene, KO, EC number and reaction matrix."
	print "Enter g for performing statistical Kruskal-Wallis test on normalized data."  
	print "Enter 0 for Exit"  
	x=raw_input()
	#x=int(c)
	if x=='a':
		execfile("pre_processing.py")
	elif x=='b':
		execfile("quality_control.py")
	elif x=='c':
		execfile("metaphlan_profiling.py")
	elif x=='d':
		execfile("mapping_mosaik.py")
	elif x=='e':
		execfile("annotate_user_genes.py")
	elif x=='f':
		execfile("data_normalization.py")
	elif x=='g':
		execfile("kruskal_test.py")
	elif x=='0':
		break
	else: 
		print "This is an incorrect input"






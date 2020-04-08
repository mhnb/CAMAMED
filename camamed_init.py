#!/usr/bin/env python
#This program is an open source program used for analyzing metagenome data
#This program takes sample reads files as well as the reference catalog of the genes and performs analyzes on them. 
#Meanwhile, the program uses the KEGG database for data analysis.
import subprocess
import os
import sys


subprocess.check_call('tar -xvf FastQC.tar.gz', shell=True)
os.system('clear')
subprocess.check_call('chmod +x MosaikAligner', shell=True)
subprocess.check_call('chmod +x MosaikBuild', shell=True)
subprocess.check_call('chmod +x MosaikJump', shell=True)
subprocess.check_call('chmod +x seqkit', shell=True)
subprocess.check_call('chmod +x ./FastQC/fastqc', shell=True)
subprocess.check_call('chmod +x fastqc_samlpe.sh', shell=True)
subprocess.check_call('chmod +x metaphlan_samlpe.sh', shell=True)
subprocess.check_call('chmod +x mosaik_build_ref.sh', shell=True)
subprocess.check_call('chmod +x mosaik_read_aligner.sh', shell=True)
subprocess.check_call('chmod +x seqkit_samlpe.sh', shell=True)
subprocess.check_call('chmod +x sra_samlpe.sh', shell=True)
subprocess.check_call('chmod +x camamed_pre_processing.py', shell=True)
subprocess.check_call('chmod +x camamed_quality_control.py', shell=True)
subprocess.check_call('chmod +x camamed_metaphlan_profiling.py', shell=True)
subprocess.check_call('chmod +x camamed_mapping_mosaik.py', shell=True)
subprocess.check_call('chmod +x camamed_kegg_annotation.py', shell=True)
subprocess.check_call('chmod +x camamed_data_normalization.py', shell=True)
subprocess.check_call('chmod +x camamed_statistical_test.py', shell=True)

exist = subprocess.call('command -v fastqc >> /dev/null', shell=True)
if exist !=0:
	print("installing FastQC")
	subprocess.check_call('sudo apt install fastqc', shell=True)

exist = subprocess.call('command -v fastq-dump >> /dev/null', shell=True)
if exist !=0:
	print("installing sra-toolkit")
	subprocess.check_call('sudo apt-get install sra-toolkit', shell=True)

exist = subprocess.call('command -v samtools >> /dev/null', shell=True)
if exist !=0:
	print("installing samtools")	
	subprocess.check_call('sudo apt-get install -y samtools', shell=True)


exist = subprocess.call('command -v cdhit-est --version >> /dev/null', shell=True)
if exist !=0:
	print("installing cdhit-est")
	subprocess.check_call('sudo apt-get install cd-hit', shell=True)


if sys.version_info.major==3:
	import importlib.util
	try:
		import pip
	except ImportError:
		print("installing pip3")
		subprocess.check_call('sudo apt install python3-pip', shell=True)
	exist = importlib.util.find_spec('pandas')
	if exist is None:
		print("installing pandas")
		subprocess.check_call('sudo pip3 install pandas', shell=True)
	exist = importlib.util.find_spec('numpy')
	if exist is None:
		print("installing numpy")
		subprocess.check_call('sudo pip3 install numpy', shell=True)
	exist = importlib.util.find_spec('scipy')
	if exist is None:
		print("installing scipy")
		subprocess.check_call('sudo pip3 install scipy', shell=True)
	exist = importlib.util.find_spec('Bio')
	if exist is None:
		print("installing biopython")
		subprocess.check_call('sudo pip3 install biopython --upgrade', shell=True)
else:
	try:
		import pip
	except ImportError:
		print("installing pip")
		subprocess.check_call('sudo apt update', shell=True)
		subprocess.check_call('sudo apt install python-pip', shell=True)
	try:
		import pandas
	except ImportError:
		print("installing pandas")
		subprocess.check_call('sudo pip install pandas', shell=True)
	try:
		import numpy
	except ImportError:
		print("installing numpy")
		subprocess.check_call('sudo pip install numpy', shell=True)
	try:
		import scipy
	except ImportError:
		print("installing scipy")
		subprocess.check_call('sudo pip install scipy', shell=True)
	try:
		import Bio
	except ImportError:
		print("installing biopython")
		subprocess.check_call('sudo pip install biopython --upgrade', shell=True)

exist = subprocess.call('command -v metaphlan2 >> /dev/null', shell=True)
if exist !=0:
	print("installing metaphlan2")
	subprocess.check_call('sudo apt install metaphlan2', shell=True)

print
print ("This program is an open source software used for analyzing metagenome data")
print ("This program takes sample reads files as well as the reference catalog of the genes and performs analyzes on them.")
print ("Meanwhile, the program uses the kegg database for data analysis.")
print
print ("***********Initial setup is done and the software is ready to run***********")
print
print
print
print

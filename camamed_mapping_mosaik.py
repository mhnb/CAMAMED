#!/usr/bin/env python 
#This function maps the reads to the genes catalog using the MOSAIK software
import sys
import os
cdp=os.getcwd()
#################################################
def help_fun():
	print('''
	This function maps the reads to the genes catalog using the MOSAIK software.
	To mapping the read sequences to the gene catalog, the following two steps should be implemented.
	In the first step, the gene catalog should be converted into an acceptable form for the MOSAIK software.
	In the second step, the reads are mapped to the gene catalog. 

	Command: ./camamed_mapping_mosaik.py method [options]

	      method:
		  -h	Shows help related to this function
		  -cag  Creating an acceptable form of gene catalog.
			options:
			   -hw (hash word): Select the length of the hash word that can be in range 4 to 32 (default=15).
			For Example: ./camamed_mapping_mosaik.py -cag		
				     ./camamed_mapping_mosaik.py -cag -hw 17

		  -mrc  Mapping reads to the gene catalog using MOSAIK software.
			options:
			   -c (core_number): Select number of cores for MOSAIK execution (default=1).
			    The number of cores must be a positive integer.
			For Example: ./camamed_mapping_mosaik.py -mrc		
				     ./camamed_mapping_mosaik.py -mrc -c 5

	''')
#################################################
if not os.path.exists('mosaik_outputs'):
	os.mkdir('mosaik_outputs')
if len(sys.argv)<2:
	help_fun()
elif str(sys.argv[1])=='-h':
	help_fun()
elif str(sys.argv[1])=='-cag':
	h=15
	sk=0
	if len(sys.argv)>2:
		ss=sys.argv
		try:
			i = ss.index('-hw')
		except ValueError:
			i=-1
			print("syntax error")
			sk=1
		if i+1<len(sys.argv) and i!=-1:
			h=int(sys.argv[i+1])
			if h<4 or h>32:
				sk=1
				print("You must select the length of the hash word from 4 to 32")
	if sk==0:
		os.system('clear')
		import subprocess
		cmd1 = './mosaik_build_ref.sh '+str(h)
		subprocess.call(cmd1, shell=True)
elif str(sys.argv[1])=='-mrc':
	co=1
	sk=0
	if len(sys.argv)>2:
		ss=sys.argv
		try:
			i = ss.index('-c')
		except ValueError:
			i=-1
			print("syntax error")
			sk=1
		if i+1<len(sys.argv) and i!=-1:
			co=int(sys.argv[i+1])
			if co<1:
				sk=1
				print("The number of cores must be a positive integer.")
	if sk==0:
		os.system('clear')
		import subprocess
		cmd1 = './mosaik_read_aligner.sh '+str(co)
		subprocess.call(cmd1, shell=True)		
else: 
	help_fun()
#This function maps the reads to the genes catalog using the MOSAIK software
import os
cdp=os.getcwd()
while 1:
	os.system('clear')
	print "To mapping the read sequences to the gene catalog, the following two steps should be implemented"
	print "In the first step, the gene catalog should be converted into an acceptable form for the MOSAIK program."
	print "In the second step, the reads are mapped to the gene catalog. Refer to the help file for this step"
	print
	print "Enter d1 for creating an acceptable form of gene catalog."
	print "Enter d2 for mapping reads to gene catalog."
	print "Enter 0 for back"
	x=raw_input()
	#x=int(c)
	if x=='d1':
		import os
		os.system('clear')
		print "At this point, you must specify the length of the hash word that can be in range 4 to 32 and the default value is 15"
		h=raw_input()
		hh=int(h)
		if hh<4 or hh>32:
			print "You must select the length of the hash word from 4 to 32"
			h=raw_input()
			print
			print
		else:
			file_p=cdp+'/files/catalog_hash.txt'
			cot_ha=open(file_p,"w")
			cot_ha.write('-hs ')
			cot_ha.write(h)
			cot_ha.write('\n')
			cot_ha.close()
			file_p=cdp+'/files/catalog_name.txt'
			cot_na=open(file_p,"r")
			text=cot_na.readline()
			cot_na.close()
			print("Your default gene catalog is '"+text+"' and the length of the hash word is "+h+". Do you agree[y/n]=")
			h=raw_input()
			if h=='y':
				import subprocess
				subprocess.call("./mosaik_build_ref.sh")
			print
			print "Please press any key to continue"
			c=raw_input()
	elif x=='d2':
		import os
		os.system('clear')
		file_p=cdp+'/files/catalog_align.txt'
		cot_al=open(file_p,"r")
		p=cot_al.readline().rstrip()
		if p=='q':
			q='fastq'
		else:
			q='fasta'
		p=cot_al.readline().rstrip()
		co=cot_al.readline().rstrip()
		s=cot_al.readline().rstrip()
		cot_al.close()
		if p=='p':
			pp="paired end "
		else:
			pp="single end "
		if p=='p':
			print (pp+q+" file type and insert size is="+s+" and number of cores="+co+", Do you agree[y/n]=")
			h=raw_input()
		elif p=='s':
			print (pp+q+" file type"+" and number of cores="+co+", Do you agree[y/n]=")
			h=raw_input()
		if h=='y':
			import subprocess
			subprocess.call("./mosaik_read_aligner.sh")
	elif x=='0':
		break
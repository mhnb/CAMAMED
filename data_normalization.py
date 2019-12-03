#This function Extract samples information in normalized bacteria, gene, KO, EC number and reaction matrix


while 1:
	import os
	os.system('clear')
	print "This function Extract samples information in normalized bacteria, gene, KO, EC number and reaction matrix"
	print
	print "Enter f1 for normalizing data using local information previously extracted from the KEGG database."   
	print "Enter f2 for normalizing data using information extracted by the user from the KEGG database." 
	print "Enter 0 for back"  
	x=raw_input()
	#x=int(c)
	if x=='f1':
		execfile("data_normalization_local.py")
	elif x=='f2':
		execfile("data_normalization_kegg.py")
	elif x=='0':
		break
	else: 
		print "This is an incorrect input"






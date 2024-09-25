import sys

# This script is to try and generate a network list to figure out the distribution of (1) 1E->1P. (2) 1E->MP, (3) ME->1P, (4) ME->MP. 
# Need to remember that the grouping is for each E-P pair but uses the entire E-P interactions for the selecting the group (not merging all E-P interactions into a single E group)
# (1) E only interacts with 1P and that P only interacts with that E
# (2) E interacts with multiple P but for the E-P pair it is the only E that interacts with that P 
# (3) E exclusive to P but P interacts with multiple E (inverse of 2)
# (4) E interacts with multiple P and that P interacts with multiple E
dicte = {}
dictp = {}
dictdat = {} 
# First assumption: no header line. Second assumption: enhancer has no name (just bed position for columns 1-3), gene in column 4, ABC interaction score (not required) in column 5, and cluster name is in column 6 (Cell Type Label) and sorted so that clusters are all together
clust = ""
for line in open(sys.argv[1], 'r'):
	linet = line.strip().split("\t")
	lint = linet[3]
	lintc = linet[0] + ":" + linet[1] + ":" + linet[2]
	# Simple check to see if it is the start of the program and properly set the cluster checker
	if clust == "":
		clust = linet[5]
	# If the cluster changes then need to print the C# data for the cluster since we want it separate per cluster then reset the dictionaries
	if linet[5] != clust:
		for en, pr in dicte.items():
			if(len(pr) == 1):
				if(len(dictp[pr[0]]) == 1):
					print(dictdat[en + "_" + pr[0]] + "\t1")
				else:
					print(dictdat[en + "_" + pr[0]] + "\t3")
			else:
				for prom in pr:
					if(len(dictp[prom]) == 1):
						print(dictdat[en + "_" + prom] + "\t2")
					else:
						print(dictdat[en + "_" + prom] + "\t4")	
		clust = linet[5]
		dicte = {}
		dictp = {}
		dictdat = {} 		
	# Now we do the actual separating data to figure out later what C# they belong to
	if(lintc not in dicte.keys()):
		dicte[lintc] = []
	if(lint not in dicte[lintc]):
		dicte[lintc].append(lint)

	if(lint not in dictp.keys()):
		dictp[lint] = []
	if(lintc not in dictp[lint]):
		dictp[lint].append(lintc)
	# Save the entire line to a dictionary so we can print it later instead of iterating the file again and appending the number that way
	lintrol = lintc + "_" + lint
	dictdat[lintrol] = "\t".join(linet)

# Once the loop is done we still need to handle the final cluster that was read in 
for en, pr in dicte.items():
	if(len(pr) == 1):
		if(len(dictp[pr[0]]) == 1):
			print(dictdat[en + "_" + pr[0]] + "\t1")
		else:
			print(dictdat[en + "_" + pr[0]] + "\t3")
	else:
		for prom in pr:
			if(len(dictp[prom]) == 1):
				print(dictdat[en + "_" + prom] + "\t2")
			else:
				print(dictdat[en + "_" + prom] + "\t4")	



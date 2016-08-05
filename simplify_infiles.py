#simply_infiles.py
#This script re-writes the output files from the BLATs to only include the headers,
#first hit (which is sometimes a self-hit), and the next hit's worth of hits

#The three files that are parsed are output files from BLATs of novel contigs versus:
#1) unmasked novel contigs
#2) masked novel contigs
#3) canFam3.1+chrUn

###############################################################################
#Parsing unmasked novel v unmasked novel BLAT results:
infile = '/home/jmkidd/kidd-lab/ampend-projects/Novel_Sequence_Analysis/RedundantNovelContigs/results/BLAT_to_UnmaskedNovelContigs/Total_UnmaskedNovelContig_BLATResults.psl'
#For testing:
#infile = '/home/jmkidd/kidd-lab/ampend-projects/Novel_Sequence_Analysis/RedundantNovelContigs/results/BLAT_to_UnmaskedNovelContigs/test.psl'
inFile = open(infile, 'r')

#Generating outfiles
outfile = '/home/jmkidd/kidd-lab/ampend-projects/Novel_Sequence_Analysis/RedundantNovelContigs/results/BLAT_to_UnmaskedNovelContigs/Total_UnmaskedNovelContig_BLATResults.psl_simplified'
#For testing:
#outfile = '/home/jmkidd/kidd-lab/ampend-projects/Novel_Sequence_Analysis/RedundantNovelContigs/results/BLAT_to_UnmaskedNovelContigs/test.psl_simplified'
outFile = open(outfile, 'w')

print 'Processing novel v unmasked...\n'

count = 1 

subIDList = {}

for lineFull in inFile:
	line = lineFull.rstrip()
	line = line.split()
	if '#' in line:
		outFile.write(lineFull)
		continue

	ID = line[0]
	ID = ID.split('_') #splits ID and length by _
	queryID = ID[0]
	qLength = int(ID[1])
	subID = line[1] #Subject contig ID

	count += 1

	if queryID == subID: #if self-self hit, ignore
		outFile.write(lineFull)
		continue
	
	if queryID not in subIDList:
		subIDList[queryID] = []
		subIDList[queryID].append(subID)
		outFile.write(lineFull)
		continue
	else:
		if subID in subIDList[queryID]:
			outFile.write(lineFull)
			continue
		else:
			continue		
print 'Done processing novel v unmasked...\n'

###############################################################################
#Parsing unmasked novel v masked novel BLAT results:
infile = '/home/jmkidd/kidd-lab/ampend-projects/Novel_Sequence_Analysis/RedundantNovelContigs/results/BLAT_to_MaskedNovelContigs/Total_MaskedNovelContig_BLATResults.psl'
#For testing
#infile = '/home/jmkidd/kidd-lab/ampend-projects/Novel_Sequence_Analysis/RedundantNovelContigs/results/BLAT_to_MaskedNovelContigs/test.psl'
inFile = open(infile, 'r')

#Generating outfile
outfile = '/home/jmkidd/kidd-lab/ampend-projects/Novel_Sequence_Analysis/RedundantNovelContigs/results/BLAT_to_MaskedNovelContigs/Total_MaskedNovelContig_BLATResults.psl_simplified'
#for testing:
#outfile = '/home/jmkidd/kidd-lab/ampend-projects/Novel_Sequence_Analysis/RedundantNovelContigs/results/BLAT_to_MaskedNovelContigs/test.psl_simplified'
outFile = open(outfile, 'w')

print 'Processing novel v masked...\n'

count = 1 

subIDList = {}

for lineFull in inFile:
	line = lineFull.rstrip()
	line = line.split()
	if '#' in line:
		outFile.write(lineFull)
		continue

	ID = line[0]
	ID = ID.split('_') #splits ID and length by _
	queryID = ID[0]
	qLength = int(ID[1])
	subID = line[1] #Subject contig ID

	count += 1

	if queryID == subID: #if self-self hit, ignore
		outFile.write(lineFull)
		continue
	
	if queryID not in subIDList:
		subIDList[queryID] = []
		subIDList[queryID].append(subID)
		outFile.write(lineFull)
		continue
	else:
		if subID in subIDList[queryID]:
			outFile.write(lineFull)
			continue
		else:
			continue		
print 'Done processing novel v masked...\n'




###############################################################################
#Parsing unmasked novel v canFam3+chrUn BLAT results:
infile = '/home/jmkidd/kidd-lab/ampend-projects/Novel_Sequence_Analysis/RedundantNovelContigs/results/BLAT_to_canFam3withchrUn/Total_canFamNovelContig_BLATResults.psl'
#for testing
#infile = '/home/jmkidd/kidd-lab/ampend-projects/Novel_Sequence_Analysis/RedundantNovelContigs/results/BLAT_to_canFam3withchrUn/test.psl'
inFile = open(infile, 'r')

#Generating output file:
outfile = '/home/jmkidd/kidd-lab/ampend-projects/Novel_Sequence_Analysis/RedundantNovelContigs/results/BLAT_to_canFam3withchrUn/Total_canFamNovelContig_BLATResults.psl_simplified'
#for testing
#outfile = '/home/jmkidd/kidd-lab/ampend-projects/Novel_Sequence_Analysis/RedundantNovelContigs/results/BLAT_to_canFam3withchrUn/test.psl_simplified'
outFile = open(outfile, 'w')

print 'Processing novel v canFam3.1+chrUn...\n'

count = 1 

subIDList = {}
countList = {}

for lineFull in inFile:
	line = lineFull.rstrip()
	line = line.split()
	if '#' in line:
		outFile.write(lineFull)
		continue

	ID = line[0]
	ID = ID.split('_') #splits ID and length by _
	queryID = ID[0]
	qLength = int(ID[1])
	subID = line[1] #Subject contig ID

	count += 1
	
	if queryID == subID: #if self-self hit, write as first hit
		outFile.write(lineFull)
		continue
	
	if queryID not in subIDList:
		subIDList[queryID] = []
		subIDList[queryID].append(subID)
		
		countList[queryID] = 1 #Add one hit for the contig
		outFile.write(lineFull)
		continue
	else:
		if countList[queryID] <= 15: #If hit count is greater than 15, these are not likely strong hits. Only want top 15 hits
			if subID in subIDList[queryID]:
				countList[queryID] += 1 #Add one to hit count
				outFile.write(lineFull)
				continue
		else:
			continue
		
print 'Done processing novel v masked...\n'
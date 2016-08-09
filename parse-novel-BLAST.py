# parse-novel-BLAST.py
# Parses BLAT outputs of the novel contigs against
#	themselves and canFam3.1+chrUn

import operator
import sys
import genutils

###############################################################################
#Example Output:
"""
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
zoey-scaffold-5220	zoey-scaffold-5220	100.00	543	0	0	1	543	1	543	2.1e-313	1070.0
zoey-scaffold-5220	wolf-scaffold-6825	91.89	74	6	0	417	490	317	390	7.3e-28	122.0
zoey-scaffold-5220	wolf-scaffold-6825	100.00	26	0	0	417	442	33	89.7e-07	51.0

queryID = line[0]
subID = line[1]
perID = int(line[2])
alnLength = int(line[3])
mismatch = int(line[4])
gaps = int(line[5])
qStart = line[6]
qEnd = line[7]
sStart = line[8]
sEnd = line[9]
eval = 
bitScore = 

"""
###############################################################################
def writingUnmaskedFasta(seqs):
	#Writing unmasked chrNovel Fasta file
	finalfile = '/home/jmkidd/kidd-lab/ampend-projects/Novel_Sequence_Analysis/RedundantNovelContigs/Final_chrNovel_Fasta/chrNovel_NonRedundant.fa'
	finalFile = open(finalfile, 'w')

	finalFile.write('>chrNovel\n')

	#Opening the original unmasked novel contig fasta sequence
	with open('/home/jmkidd/kidd-lab/ampend-projects/Novel_Sequence_Analysis/NovelSequence/novel.v2.fa') as fasta:
		prev_seq = []
		for line in fasta:
			if line.startswith('>'):
				line = line.rstrip()
				line = line.split()	
				ID = line[0]			
				shortID = ID.replace('>','')
				if shortID in passList:
					seqs.append("".join(prev_seq))
					prev_seq = []
			else:
				if shortID in passList:
					prev_seq.append(line.rstrip())

	# This appends the lines after the last ">"
	seqs.append("".join(prev_seq))
	#Writing out
	finalFile.write(Nstring.join(seqs))
###############################################################################
def writingMaskedFasta(seqs):
	#Writing unmasked chrNovel 
	finalfile = '/home/jmkidd/kidd-lab/ampend-projects/Novel_Sequence_Analysis/RedundantNovelContigs/Final_chrNovel_Fasta/chrNovel_NonRedundant.fa.masked'
	finalFile = open(finalfile, 'w')
	
	finalFile.write('>chrNovel\n')
	
	#Opening the original unmasked novel contig fasta sequence
	with open('/home/jmkidd/kidd-lab/ampend-projects/Novel_Sequence_Analysis/NovelSequence/novel.v2.fa') as fasta:
		prev_seq = []
		for line in fasta:
			if line.startswith('>'):
				line = line.rstrip()
				line = line.split()	
				ID = line[0]			
				shortID = ID.replace('>','')
				if shortID in passList:
					seqs.append("".join(prev_seq))
					prev_seq = []
			else:
				if shortID in passList:
					prev_seq.append(line.rstrip())

	# This appends the lines after the last ">"
	seqs.append("".join(prev_seq))
	#Writing out
	finalFile.write(Nstring.join(seqs))
###############################################################################
def makeGapDict(gapList):
	#HARD CODED: This gap file contains the gaps of all chromosomes and then for the merged chrUn
	gapFile = open('/home/ampend/kidd-lab/ampend-projects/Novel_Sequence_Analysis/RedundantNovelContigs/canFam3.1.withUnk.gaps.bed','r')
	
	for line in gapFile:
		line = line.rstrip()
		line = line.split()
		
		chromGap = line[0]
		startGap = int(line[1])
		endGap = int(line[2])
		
		if chromGap not in gapList:
			gapList[chromGap] = []
			gapList[chromGap].append((startGap, endGap))
		else:
			gapList[chromGap].append((startGap, endGap))
	return gapList
###############################################################################
def findGapOverlap(queryID, sStart, sEnd, chromNum, gapList, dict):
	#Add chr to chromosome Number for processing gaps
	if chromNum == 'Un':
		longChr = 'chrUn'
	else:
		longChr = 'chr' + chromNum

	StartIDX = 0 #Setting start to 0 
	EndIDX = len(gapList[longChr]) - 1 #Setting end index length-1
	MidIDX = 0
	gapPresent = False
	#print '%s\t%s\t%s\t%s\t%s' % (queryID,longChr,StartIDX,EndIDX,MidIDX)

	while True:  	
		MidIDX = int(EndIDX + StartIDX) / 2
		start_mid = int(gapList[longChr][MidIDX][0])
		end_mid = int(gapList[longChr][MidIDX][1])
		#THOSE ARE NARROWING SEARCH FOR GAPS
		if sEnd < start_mid:
			EndIDX = int(MidIDX)
			#print '%s\t%s\t%s\t%s\t%s' % (queryID,longChr,StartIDX,EndIDX,MidIDX)
			if EndIDX - StartIDX <= 1: #no gap overlap possible
				break
			else:
				continue
		if sStart > end_mid:
			StartIDX = int(MidIDX)
			#print '%s\t%s\t%s\t%s\t%s' % (queryID,longChr,StartIDX,EndIDX,MidIDX)
			if EndIDX - StartIDX <= 1: #no gap overlap possible
				break
			else:
				continue
		#THOSE BELOW ARE OVERLAPPING WITH GAPS
		if sStart < start_mid and sEnd > end_mid:
			dict[queryID][18] = True
			gapPresent = True
			#print '%s\t%s\t%s\t%s\t%s' % (queryID,chr,StartIDX,EndIDX,MidIDX)
			break
		if sStart == start_mid or sStart == end_mid:
			dict[queryID][18] = True
			gapPresent = True
			#print '%s\t%s\t%s\t%s\t%s' % (queryID,chr,StartIDX,EndIDX,MidIDX)
			break
		if sEnd == start_mid or sEnd == end_mid:
			dict[queryID][18] = True
			gapPresent = True
			#print '%s\t%s\t%s\t%s\t%s' % (queryID,chr,StartIDX,EndIDX,MidIDX)
			break			
		#This means there's an overlap if all conditions fail
		if EndIDX - StartIDX <= 1:
			#no gap overlap possible
			break
	return dict[queryID][18], gapPresent
	
"""#SAVING OLD ROUTINE
def findGapOverlap(queryID, sStart, sEnd, chromNum, gapList, dict):
	#Add chr to chromosome Number for processing gaps
	if chromNum == 'Un':
		longChr = 'chrUn'
	else:
		longChr = 'chr' + chromNum
	
	StartIDX = 0 #Setting start to 0 
	EndIDX = len(gapList[longChr]) - 1 #Setting end index length-1
	MidIDX = 0
	gapPresent = False
	#print '%s\t%s\t%s\t%s\t%s' % (queryID,longChr,StartIDX,EndIDX,MidIDX)

	while True:  	
		MidIDX = int(EndIDX + StartIDX) / 2
		start_mid = int(gapList[longChr][MidIDX][0])
		end_mid = int(gapList[longChr][MidIDX][1])
		#THOSE ARE NARROWING SEARCH FOR GAPS
		if sEnd < start_mid:
			EndIDX = int(MidIDX)
			#print '%s\t%s\t%s\t%s\t%s' % (queryID,longChr,StartIDX,EndIDX,MidIDX)
			if EndIDX - StartIDX <= 1: #no gap overlap possible
				break
			else:
				continue
		if sStart > end_mid:
			StartIDX = int(MidIDX)
			#print '%s\t%s\t%s\t%s\t%s' % (queryID,longChr,StartIDX,EndIDX,MidIDX)
			if EndIDX - StartIDX <= 1: #no gap overlap possible
				break
			else:
				continue
		#THOSE BELOW ARE OVERLAPPING WITH GAPS
		if sStart < start_mid and sEnd > end_mid:
			dict[queryID][18] = True
			gapPresent = True
			#print '%s\t%s\t%s\t%s\t%s' % (queryID,chr,StartIDX,EndIDX,MidIDX)
			break
		if sStart == start_mid or sStart == end_mid:
			dict[queryID][18] = True
			gapPresent = True
			#print '%s\t%s\t%s\t%s\t%s' % (queryID,chr,StartIDX,EndIDX,MidIDX)
			break
		if sEnd == start_mid or sEnd == end_mid:
			dict[queryID][18] = True
			gapPresent = True
			#print '%s\t%s\t%s\t%s\t%s' % (queryID,chr,StartIDX,EndIDX,MidIDX)
			break			
		#This means there's an overlap if all conditions fail
		if EndIDX - StartIDX <= 1:
			#no gap overlap possible
			break
	return dict[queryID][18], gapPresent"""
###############################################################################
###############################################################################
#Generating output file that totals all the results per contig
outfile = '/home/jmkidd/kidd-lab/ampend-projects/Novel_Sequence_Analysis/RedundantNovelContigs/results/BLATResultsTable_AllContigs.txt'
outFile = open(outfile, 'w')

###############################################################################
#Parsing unmasked novel v unmasked novel BLAT results:
infile = '/home/jmkidd/kidd-lab/ampend-projects/Novel_Sequence_Analysis/RedundantNovelContigs/results/BLAT_to_UnmaskedNovelContigs/Total_UnmaskedNovelContig_BLATResults.psl_simplified'
#For testing:
#infile = '/home/jmkidd/kidd-lab/ampend-projects/Novel_Sequence_Analysis/RedundantNovelContigs/results/BLAT_to_UnmaskedNovelContigs/test.psl_simplified'
inFile = open(infile, 'r')

count = 1
dict = {} #Opening up dictionary for storing contig data
LOG_EVERY_N = 10000 #For counting processed contigs during analysis

print 'Processing novel v unmasked...\n'

for lineFull in inFile:
	if (count % LOG_EVERY_N) == 0:
		print 'Processed novel v unmasked hits #:',count
	
	line = lineFull.rstrip()
	line = line.split()
	
	if '#' in line: #skipping header lines
		continue 
	count += 1

	#The 1st column of blat hit looks like the following: #wolf-scaffold-1_435, where the _435 = length of contig #We want to keep track of the length of the contig for parsing
	ID = line[0]
	ID = ID.split('_') #splits ID and length by _ character
	queryID = ID[0]
	qLength = int(ID[1])
	subID = line[1] #Subject contig ID
	
	#Open up dictionary for queryID
	if queryID not in dict.keys():######DICTIONARY STRUCTURE KEY: ######
		#0=ID, 1=AllPass (default = True), 
		#2=BLAT1 Pass (default = True), 3=BLAT1 Contig Hit, 4=BLAT1 %ID
		#5=BLAT2 Pass (default = True), 6=BLAT2 Contig Hit, 7=BLAT2 %ID
		#8=BLAT3 Pass (default = True), 9=BLAT3 canFam3 Top Hit, 10=BLAT3 %ID, 11=Start Coord, 12=End Coord, 13=canFam3 alignment length (bp), 
		#14 = BLAT3 Pass Adding Multiple hits (default = True), 15 = Total Sequenced UNaligned to CanFam3+chrUn, 16 = Proportion of contig UNaligned
		#17 = Contig Length 
		#18 = Gap Present (default = False), 19=Gap Counts (default=0)
		#AllPass is default = True, because without having a hit to parse in BLAT3 (to canFam3+chrUn)
		#	there wont be a chance for the contig to go through the for loop assigning all pass versus fail (see below)
		#BLAT3 is default = True because if there's no hit, then the contig will not receive an outfile to be parsed
		#	and there needs to be a result there
		#Setting default values for dictionary
		dict[queryID] = [queryID,True,True,0,0,True,0,0,True,0,0,0,0,0,True,0,0,0,False,0]
		dict[queryID][17] = qLength

	#BLAT results		
	perID = float(line[2])
	alnLength = int(line[3])
	propAligned = float(alnLength)/qLength
	
	maxHit = 0.8 * qLength #max allowed hit is > 80% of query hit
	
	if queryID != subID and dict[queryID][2] is True:
		if perID > 75.0 and alnLength > maxHit: #If it hits with >75% identity AND the alignment length is >80% of total length
			dict[queryID][2] = False #This contig FAILS BLAT against unmasked novel contigs
			dict[queryID][3] = subID #Hit ID added
			dict[queryID][4] = propAligned #Proportion of contig aligned with percent ID > 75%
			dict[queryID][1] = False #Therefore FAILS OVERALL to be used
		else:
			dict[queryID][2] = True #Passes, but below info will still be added
			dict[queryID][3] = subID #Hit ID added
			if perID > 75.0:
				dict[queryID][4] = propAligned #Proportion of contig aligned with percent ID > 75%

print 'Done processing novel v unmasked...\n'

###############################################################################
#Parsing unmasked novel v masked novel BLAT results:
infile = '/home/jmkidd/kidd-lab/ampend-projects/Novel_Sequence_Analysis/RedundantNovelContigs/results/BLAT_to_MaskedNovelContigs/Total_MaskedNovelContig_BLATResults.psl_simplified'
#For testing
#infile = '/home/jmkidd/kidd-lab/ampend-projects/Novel_Sequence_Analysis/RedundantNovelContigs/results/BLAT_to_MaskedNovelContigs/test.psl_simplified'

inFile = open(infile, 'r')

print 'Processing novel v masked...\n'
count = 1 

for lineFull in inFile:
	line = lineFull.rstrip()
	line = line.split()
	if '#' in line: #skipping headers
		continue 

	ID = line[0]
	ID = ID.split('_') #splits ID and length by _
	queryID = ID[0]
	qLength = int(ID[1])
	subID = line[1] #Subject contig ID

	count += 1

	#BLAT results
	perID = float(line[2])
	alnLength = int(line[3])
	propAligned = float(alnLength)/qLength
		
	maxHit = 0.8 * qLength #max allowed hit is > 80% of query hit

	if queryID != subID and dict[queryID][5] is True:
		if perID > 75.0 and alnLength > maxHit: #If it hits with >75% identity AND the alignment length is >80% of total length
			dict[queryID][5] = False #This contig FAILS BLAT against unmasked novel contigs
			dict[queryID][6] = subID #Hit ID added
			dict[queryID][7] = propAligned #Proportion of contig aligned with percent ID > 75%
			dict[queryID][1] = False #Therefore FAILS OVERALL to be used
		else:
			dict[queryID][5] = True #Passes, but below info will still be added
			dict[queryID][6] = subID #Hit ID added
			if perID > 75.0:
				dict[queryID][7] = propAligned #Proportion of contig aligned with percent ID > 75%

	if (count % LOG_EVERY_N) == 0:
		print 'Processed novel v masked hits #:',count
print 'Done processing novel v masked...\n'

###############################################################################
#Steps for finding overlaps:	
	#Want to take qStart and qEnd, if qStart less than valueStart, replace
		#else keep valueStart
	#Then want to take qEnd is less than valueEnd, replace.
		#Else, keep valueStart
	#Do this for each value in dictionary until there are no more checks to make
	#When done, count the amount of nucleotides that align to canFam3 by adding
	#	up all start,end in dictionary
###############################################################################
###############################################################################
#Generating dictionary of gap coordinates
gapList = {}
makeGapDict(gapList)


###############################################################################
#Parsing unmasked novel v canFam3+chrUn BLAT results:
infile = '/home/jmkidd/kidd-lab/ampend-projects/Novel_Sequence_Analysis/RedundantNovelContigs/results/BLAT_to_canFam3withchrUn/Total_canFamNovelContig_BLATResults.psl_simplified'
#for testing
#infile = '/home/jmkidd/kidd-lab/ampend-projects/Novel_Sequence_Analysis/RedundantNovelContigs/results/BLAT_to_canFam3withchrUn/test.psl'

inFile = open(infile, 'r')
prevList = []
skipList = []
hitList = []
hits = {}
hitCount = 0
chromList = []

count = 1
check = 0
print 'Processing novel v canFam3+chrUn...\n'

for lineFull in inFile:
	line = lineFull.rstrip()
	line = line.split()
	if '#' in line: #Skipping header
		continue 
			
	#The 1st column of blat hit looks like the following: wolf-scaffold-1_435, where the _435 = length of contig. We want to keep track of the length of the contig for parsing
	ID = line[0]
	ID = ID.split('_') #splits ID and length by '_'
	queryID = ID[0]
	qLength = int(ID[1])
	subID = line[1] #canFam chromosome
	chromNum = subID.split('chr')
	chromNum = chromNum[1]
	dict[queryID][17] = qLength
		
	count += 1 #hit count being processed overall
	
	if queryID == subID: #if self-self hit, ignore
		dict[queryID][8] = True
		continue
	if queryID not in hitList: #If not, add to processed list
		if perID > 0.75:
			hits[queryID] = {}
			hitList.append(queryID)
			gapCount = 0 #Set gap counts = 0 for a new contig ID
			prevID = 0
			prevStart = 0
			prevEnd = 0
			prevChrom = 0
		else:
			continue
	if queryID not in hits:
		continue #breaks loop if it didn't meet the first requirement
			
	#BLAT results
	perID = float(line[2])
	alnLength = int(line[3])
	qStart = line[6]
	qEnd = line[7]
	
	#Determining hit start/end
	if int(line[8]) > int(line[9]): #If BLAT results are in (-) orientation:
		sStart = int(line[9])
		sEnd = int(line[8])
	else:
		sStart = int(line[8])
		sEnd = int(line[9])
	
	#Defining maximum alignment length
	maxHit = 0.8 * qLength #max allowed hit is > 80% of query hit

	#########################################	
	#Processing THIS hit
	if dict[queryID][8] is True:
		if dict[queryID][10] < 1: #Checking to make sure that only the top hit goes into the outFile, default for dict[10] = 0
			dict[queryID][9] = subID
			dict[queryID][10] = perID
			dict[queryID][11] = sStart
			dict[queryID][12] = sEnd
			dict[queryID][13] = alnLength
			#Saving this hit as the one to compare for all additional 
			prevHitID = subID
			prevHitStart = int(sStart)
			prevHitEnd = int(sEnd)
			prevHitChrom = subID
			prevHitqStart = qStart
			prevHitqEnd = qEnd
			if chromNum not in hits[queryID].keys():
				hits[queryID][chromNum]=[]
			else:
				continue
		
		#Processing whether or not the novel contig aligns to canFam3 where a gap is present
		if chromNum in hits[queryID].keys():
			if 'M' not in chromNum: #chrM does not have gaps
				temp = findGapOverlap(queryID, sStart, sEnd, chromNum, gapList, dict)
				gapPresent = temp[1]
				if gapPresent is True:
					dict[queryID][19] += 1
					gapCount += 1

	#Comparing THIS hit to PREVIOUS hit
	if queryID in hitList: #Has this hit been processed already?
		if prevHitChrom == subID: #If chromosomes are ok, keep going
			if queryID not in skipList: #If the ID is not in skipList, keep going
				if perID > 75.0 or dict[queryID][18] is True: #If percent identity > 80%, add as a hit OR If there is a gap, add as a hit
					hits[queryID][chromNum].append((qStart, qEnd)) #append start, end to list by chromosome
				#if dict[queryID][18] is True: #If there is a gap, add as a hit
				#	hits[queryID][chromNum].append((qStart, qEnd)) #append start, end to list by chromosome
				else:
					del hits[queryID]
		else:
			skipList.append(queryID)
	
	if (count % LOG_EVERY_N) == 0:
		print 'Processed canFam hits #:',count

	#if count > 50:
	#	break
		
#############################################################
#Reminder of dictionary structure:
		#0=ID, 1=AllPass (default = True), 
		#2=BLAT1 Pass (default = True), 3=BLAT1 Contig Hit, 4=BLAT1 %ID
		#5=BLAT2 Pass (default = True), 6=BLAT2 Contig Hit, 7=BLAT2 %ID
		#8=BLAT3 Pass (default = True), 9=BLAT3 canFam3 Top Hit, 10=BLAT3 %ID, 11=Start Coord, 12=End Coord, 13=canFam3 alignment length (bp), 
		#14 = BLAT3 Pass Adding Multiple hits (default = True), 15 = Total Sequence UNaligned to CanFam3+chrUn, 16 = Proportion of contig UNaligned
		#17 = Contig Length 
		#18 = Gap Present (default = False), 19=Gap Counts (default=0)
#############################################################

print '\nMerging coordinates...' 

for key in dict.keys(): #for contig (queryID) in dictionary
	if key in hits: #for contig (queryID) in hits (only those that hit to canFam significantly)
		#Merging the coordinates per novel contig to calculate total length aligned to canFam3
		for chr in hits[key].keys(): #because we have chromosome as key for hits[queryID]
			hits[key][chr] = [(int(x), int(y)) for x, y in hits[key][chr]] #saves the coordinates x= start, y=end of BLAT alignments
			hits[key][chr].sort(key=operator.itemgetter(0)) #This sorts the coordinates added
			new_list = [hits[key][chr][0]] #have to create new list for next steps

			for start, end in hits[key][chr][1:]: #for x, y (start, end) in coordinates per contig
				last_start, last_end = new_list[-1]
				if start > last_end: #if the start comes after the previous end, then append the pair of coordinates to the new list
					new_list.append((start, end))
				else: #if not, then the coordinate pair is now (last_start, and the maximum of last_end and end)
					new_list[-1] = last_start, max(last_end, end)
			
			#print(hits[key][chr]) #--> old, without merging of overlapping coordinates
			#print(new_list) #--> new list, with merging of overlapping coordinates
			
			#Writing all left and right coordinates to a list, per contig
			left=[m[0] for m in new_list]
			right=[m[1] for m in new_list]
			length_each=[abs(int(m[1])-int(m[0])+1) for m in new_list] #calculates the distance between the left and right coordinate along the list
			alignedTotal = sum(length_each) #adds up all lengths between coordinates to become the total amount spanned by the contig on canFam3
			contigLength = int(dict[key][17])
			unalignedTotal = contigLength - alignedTotal
			dict[key][15]  = unalignedTotal
			
			#Calculating proportion of contig that aligns to CanFam3+chrUn
			propAligned = float(alignedTotal)/contigLength
			dict[key][16] = propAligned #Saving propAligned to dictionary in index #16
			
			#NEW SECTION
			#GAPS
			gappresent = dict[key][18] #gappresent can be true or false, default = false
			if gappresent is False: #If gaps present = False
				if propAligned > .8: #If greater than 75% of the contig is unaligned to canFam3+chrUn 
					dict[key][1] = False #Then, it fails to pass and is considered a rendundant contig
					dict[key][8] = False #Then it also fails canFam3+chrUn BLAT
			if gappresent is True: #If gaps present = True
				#if dict[key][1] is True: #If this contig passed the other two BLAT conditions
				dict[key][1] = True #It should be considered as pass for the all conditions regardless of canFam3 pass
		"""
		#Now choosing between the two scaffolds that hit to one another, choosing the contig that's longest
		if dict[key][2] is False: #If contig fails self BLAT against unmasked contigs
			contigHitID = dict[key][3] #ContigID of hit
			#print 'contig: %s\thitID: %s' % (key, contigHitID)
			if contigHitID in hits:
				if dict[contigHitID][8] is False and dict[key][8] is True:
					dict[key][1] = True
					dict[contigHitID][1] = False
				if dict[contigHitID][8] is True and dict[key][8] is False: #if ContigID of hit did not pass canFam BLAT, stays False
					dict[key][1] = False
					dict[contigHitID][1] = True
				if dict[contigHitID][17] > dict[key][17]: #If the hit is longer than the current contig being analyzed, keep the longer of the two
					dict[contigHitID][1] = True
				else:
					dict[key][1] = True
		if dict[key][5] is False: #If contig fails self BLAT against the masked contigs
			if dict[key][1] is not False: #If it's already false, then skip because it failed for a different reason
				contigHitID = dict[key][3] #ContigID of hit
				if contigHitID in hits:
					if dict[contigHitID][17] > dict[key][17]: #If the hit is longer than the current contig being analyzed, keep the longer of the two
						dict[contigHitID][1] = True
					else:
						dict[key][1] = True
		
		#Checking percent ID of hit against canFam3+chrUn to determine whether the contig passes overall
		perID = dict[key][10] #% ID of contig BLAT hit to top canFam3.1+chrUn hit
		if perID > 75.0:
			dict[key][8] = False
			if gappresent is not True:
				dict[key][1] = False
			if gappresent is True:
				dict[key][1] = True
		"""
	else:
		continue
#############################################################
#writing final bedfile	
bedfile = '/home/ampend/kidd-lab/ampend-projects/Novel_Sequence_Analysis/RedundantNovelContigs/results/canFam3_novelcontig_hitsOnly.bed'
bedFile = open(bedfile, 'w')
gapsFile = '/home/ampend/kidd-lab/ampend-projects/Novel_Sequence_Analysis/RedundantNovelContigs/canFam3.1.withUnk.gaps.bed'
bedOutfile = '/home/ampend/kidd-lab/ampend-projects/Novel_Sequence_Analysis/RedundantNovelContigs/results/Intersect_canFam3Hits_AgainstGaps.txt'

for key in dict.keys():
	if key in hits: #alignment to genome dictionary 
		chrom = dict[key][9]
		start = dict[key][11]
		bed_Start = start - 2
		end = dict[key][12]
		bed_end = end + 2
		ID = dict[key][0]
		
		bedFile.write('%s\t%i\t%i\t%s\n' % (chrom,start,end,ID))

print 'running bedtools window against gaps within 2bp: '		
cmd = 'bedtools window -w 5 -a %s -b %s > %s' % (bedfile, gapsFile, bedOutfile)

print cmd
#genutils.runCMD(cmd)
print 'Done with intersection'
print 'Results written to: ', bedOutfile


#Have to choose be
for key in dict.keys():
	gappresent = dict[key][18]
	if key in hits:
		#Now choosing between the two scaffolds that hit to one another, choosing the contig that's longest
		#	or the one that passes the canFam3 blat
		if dict[key][2] is False: #If contig fails self BLAT against unmasked contigs
			contigHitID = dict[key][3] #ContigID of hit
			#print 'contig: %s\thitID: %s' % (key, contigHitID)
			if contigHitID in hits: #ContighitID=contigA, key = contigB
				if dict[contigHitID][8] is False and dict[key][8] is True: #if contigA fails canFam3 blat but contigB passes
					dict[key][1] = True #then contigB passes overall
					dict[contigHitID][1] = False #then contigA fails overall
					continue
				if dict[contigHitID][8] is True and dict[key][8] is False: #if contigA passes canFam3 blat but contigB fails
					dict[key][1] = False #then contigB fails overall
					dict[contigHitID][1] = True #then contigA passes overall
					continue
				if dict[contigHitID][8] is True and dict[key][8] is True: #if both contigA and contigB pass canFam BLAT
					if dict[contigHitID][17] > dict[key][17]: #If the contigB is longer than the contigA, keep the longer of the two
						dict[contigHitID][1] = True #longest passes overall
						dict[key][1] = False #shortest fails overall
						continue
					if dict[key][17] > dict[contigHitID][17]: #If the contigA is longer than the contigB, keep the longer of the two
						dict[contigHitID][1] = False #shortest fails overall
						dict[key][1] = True #longest passes overall
						continue
				if dict[contigHitID][8] is False and dict[key][8] is False: #if both contigA and contigB fail canFam BLAT
						if dict[contigHitID][18] is False and dict[key][18] is False: #if both contigA and contigB dont span gaps
							dict[contigHitID][1] = False #then both fail overall
							dict[key][1] = False #then both fail overall
							continue
				#else:
				#	dict[key][1] = True
		if dict[key][5] is False: #If contig fails self BLAT against the masked contigs
			contigHitID = dict[key][3] #ContigID of hit
			if contigHitID in hits: #ContighitID=contigA, key = contigB
				if dict[contigHitID][8] is False and dict[key][8] is True: #if contigA fails canFam3 blat but contigB passes
					dict[key][1] = True #then contigB passes overall
					dict[contigHitID][1] = False #then contigA fails overall
					continue
				if dict[contigHitID][8] is True and dict[key][8] is False: #if contigA passes canFam3 blat but contigB fails
					dict[key][1] = False #then contigB fails overall
					dict[contigHitID][1] = True #then contigA passes overall
					continue
				if dict[contigHitID][8] is True and dict[key][8] is True: #if both contigA and contigB pass canFam BLAT
					if dict[contigHitID][17] > dict[key][17]: #If the contigB is longer than the contigA, keep the longer of the two
						dict[contigHitID][1] = True #longest passes overall
						dict[key][1] = False #shortest fails overall
						continue
					if dict[key][17] > dict[contigHitID][17]: #If the contigA is longer than the contigB, keep the longer of the two
						dict[contigHitID][1] = False #shortest fails overall
						dict[key][1] = True #longest passes overall
						continue
				if dict[contigHitID][8] is False and dict[key][8] is False: #if both contigA and contigB fail canFam BLAT
						if dict[contigHitID][18] is False and dict[key][18] is False: #if both contigA and contigB dont span gaps
							dict[contigHitID][1] = False #then both fail overall
							dict[key][1] = False #then both fail overall
							continue
		"""if dict[key][1] is True: #If it's already false, then skip because it failed for a different reason
				contigHitID = dict[key][3] #ContigID of hit
				if contigHitID in hits:
					if dict[contigHitID][17] > dict[key][17]: #If the hit is longer than the current contig being analyzed, keep the longer of the two
						dict[contigHitID][1] = True
					else:
						dict[key][1] = True
			else: 
				continue
		"""
		#Checking percent ID of hit against canFam3+chrUn to determine whether the contig passes overall
		perID = dict[key][10] #% ID of contig BLAT hit to top canFam3.1+chrUn hit
		propAligned = dict[key][16] #calculated above, this is the proportion of the novel contig that aligns to canFam3 (sum of all relevant hits to canFam3)
		if perID > 75.0 and propAligned > 0.8:
			dict[key][8] = False
			if gappresent is False:
				dict[key][1] = False
	#GAP CHECK
	if dict[key][18] is True: #if there is a gap being spanned by a contig, it automatically passes overall
		dict[key][1] = True #turns overall pass to true

#############################################################
#############################################################

#Print out table where:
#Column 1 = queryID = Scaffold ID
#Column 2 = Passes all criteria (True in all BLATs)
#Column 3 = True if it has no significant hit in novel vs unmasked novel
#Column 4 = True if it has no significant hit in novel vs masked novel
#Column 5 = True if it has no significant hit in novel vs canFam3.1+chrUn

#Writing out header
outFile.write('ScaffoldID\tPassAll\tPassUnmaskedBLAT\tUnmaskedContigHit\tProportionAligningToUnmaskedContigGT75PercentID\t')
outFile.write('PassMaskedBLAT(T/F)\tMaskedContigHit\tProportionAligningToMaskedContigGT75PercentID\t')
outFile.write('PassCanFamBLAT(T/F)\tCanFamHitChr\tCanFamHitPercentIdentity\tCanFamHitStart\tCanFamHitEnd\tCanFamAlignmentLength\t')
outFile.write('CanFamBLATMultipleHits(T/F)\tProportionAlignedToCanFam\tTotalAlignmentOfMultipleHitsToCanFam\t')
outFile.write('NovelContigLength\t')
outFile.write('GapPresent(T/F)\tGapsSpanned\n')

#Writing out results
for keys in dict.keys():
	outFile.write("\t".join(map(str,dict[keys])) + '\n')

#############################################################

print 'Processed %i contigs' % (count)
print 'DONE PARSING!\n\n'

"""
#############################################################
#############################################################
#Now moving around the fasta files that pass
outFile.close()
resultsfile = '/home/jmkidd/kidd-lab/ampend-projects/Novel_Sequence_Analysis/RedundantNovelContigs/results/BLATResultsTable_AllContigs.txt'
resultsFile = open(resultsfile, 'r')

fastaDir = '/home/jmkidd/kidd-lab/ampend-projects/Novel_Sequence_Analysis/RedundantNovelContigs/fasta/'

print 'Reading in results from:', resultsfile

passCount = 0
failCount = 0
#more detailed fails
failUMCount = 0
failMCount = 0
failBothContigCount = 0
failCanFamCount = 0
failAll = 0

passList = []

for line in resultsFile:
	line = line.rstrip()
	line = line.split()

	#Summary Stat calculations of fail versus pass
	#All-out fail	
	if line[1] == 'False':
		failCount += 1 
		#Different types of failure
		if line[2] == 'False' and line[5] == 'True':
			failUMCount += 1
		if line[2] == 'True' and line[5] == 'False':
			failMCount += 1
		if line[2] == 'False' and line[5] == 'False':
			failBothContigCount += 1	
		if line[2] == 'True' and line[5] == 'True' and line[8] == 'False':
			failCanFamCount += 1
		if line[2] == 'False' and line[5] == 'False' and line[8] == 'False':
			failAll += 1
		continue #Continue if any are false, remaining will have FASTAs processed below.
		
	#Contigs that pass BLAT thresholds
	if line[1] == 'True':
		passCount += 1 

	contigID = line[0]
	
	passList.append(contigID)		

	#if passCount > 10:
		#break


print '\n'
print 'Contig Stats:'
print '%i contigs pass all BLATs\n' % (passCount)		
print '%i contigs fail altogether' % (failCount)
print '%i contigs fail all three BLATs' % (failAll)
print '%i contigs fail when BLATed against unmasked novel contigs' % (failUMCount)
print '%i contigs fail when BLATed against masked novel contigs' % (failMCount)
print '%i contigs fail when BLATed against unmasked AND masked novel contigs' % (failBothContigCount)
print '%i contigs fail when BLATed against canFam3.1 + chrUn\n' % (failCanFamCount)

print 'Math check:'
sumFail = failUMCount + failMCount + failBothContigCount + failCanFamCount + failAll
print 'Failed contigs:',sumFail
sumAll = sumFail + passCount
print 'Total analyzed contigs',sumAll	

###############################################################
#This goes through line-by-line of the FASTA file of the novel contigs and:
#1) Asks if the novel contig ID is in the "pass list" generated above
#2) If so, it takes the fasta sequence and adds it to the seq list
#3) It then concatenates each fasta with user-defined length of N's (here, 1000)

#For adding 100 'N's between FASTA sequences to generate chrNovel
Nstring = 1000 * 'N'

#Writing unmasked chrNovel 
seqs = []
#writingUnmaskedFasta(seqs)

#Writing masked chrNovel 
seqs = []
#writingMaskedFasta(seqs)
"""



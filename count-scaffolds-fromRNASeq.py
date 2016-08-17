# count-scaffolds-fromRNASeq.py
# Tracks expression of each novel contig from
#	Cufflinks/CuffMerge output

import re
import sys 
from optparse import  OptionParser

###############################################################################
USAGE = """
count-scaffolds-fromRNASeq.py
			 	--input < File with filenames needed to go into analysis > 
				--fasta < Novel contig fasta file >
input == File with filenames needed to go into analysis
fasta == FASTA file that has the novel contig sequences

"""

parser = OptionParser(USAGE)
parser.add_option('--input',dest='input', help = 'File with filenames needed to go into analysis')
parser.add_option('--fasta',dest='fasta', help = 'FASTA file that has the novel contig sequences')
(options, args) = parser.parse_args()

parser = OptionParser(USAGE)
if options.input is None:
	parser.error('input parsed SV file name not given')
if options.fasta is None:
	parser.error('FASTA file that has the novel contig sequences')
############################################################################
def parse_infile(inFile,contigDict):
	lineCount = 0
	for line in inFile:
		lineCount += 1 
		if lineCount == 1: #skips header line
			continue
		line = line.rstrip()
		line = line.split()
		
		if 'scaffold' in line[6]: #ensures that this is a novel contig ID and not a gene or cufflinks ID
			scaffID = line[6]
			tmp = scaffID.split(':')
			scaffID = tmp[0]
		else: #ignores "genes" that aren't novel contigs
			continue
			
		if scaffID in contigDict:
			if runID in contigDict[scaffID]:
				continue
			else:
				contigDict[scaffID].append(runID)
############################################################################
def write_output(contigDict,runArray,outFile):
	header = []
	header.append('ContigID')
	
	yesTotalCount = 0
	for r in runArray:
		header.append(r)
	outFile.write('\t'.join(map(str,header)) + '\n')

	for key in contigDict:
		line = []
		line.append(key)
		yesCount = 0
		for r in runArray:
			tmp = r.split('_')
			runID = tmp[0]
			if runID in contigDict[key]:
				line.append('Yes')
				yesCount += 1
			else:
				line.append('No')
		if yesCount >= 3:
			yesTotalCount += 1
		outFile.write('\t'.join(map(str,line)) + '\n')
		
	print 'The results written to %s have a total of %i contigs with at expression detected in at least 3 libraries' % (outFile,yesTotalCount)
############################################################################
#############################################################################

#Dictionary for saving ALL data
dict = {}

##############################
# Saving IDs in contig IDs in
#	dictionary
##############################
contigFasta = open(options.fasta, 'r')
print 'Reading in fasta sequences from the novel contig fasta file: \n', options.fasta

contigDict = {}
count = 0

for line in contigFasta:
	line = line.rstrip()

	if '>' in line:
		count += 1
		m = re.match(r">(.+)", line)
		contigID = m.group(1)	
	contigDict[contigID] = []
print 'Read in %s contig IDs into dictionary\n\n' % len(contigDict)

##############################
# Defining paths for finding
#	cufflinks results
##############################
novelDir = '/home/ampend/kidd-lab/ampend-projects/Novel_Sequence_Analysis/rna-seq/results/canFam3.1-Novel/'
canFamDir = '/home/ampend/kidd-lab/ampend-projects/Novel_Sequence_Analysis/rna-seq/results/canFam3.1-withUn/'
canFamNovelDir = '/home/ampend/kidd-lab/ampend-projects/Novel_Sequence_Analysis/rna-seq/results/canFam3.1-Novel-Unknown/'

##############################
# Defining output files
##############################
dir = '/home/ampend/kidd-lab/ampend-projects/Novel_Sequence_Analysis/rna-seq/results/ContigCounts/'
novelOut = dir + 'Novel_Cufflinks_ContigCounts.txt'
canFamNovelOut = dir + 'CanFam+Novel_Cufflinks_ContigCounts.txt'

##############################
# Defining infiles
##############################
sampleInfo = open(options.input, 'r')
print 'Reading in sample information from file: \n', options.input

##############################
# Parsing NOVEL 
##############################
count = 0 #keeping track of line number
runArray = [] #keeping track of runIDs that have been parsed (SRR...)
outFile = open(novelOut, 'w')
	
for line in sampleInfo:
	count += 1
	if count ==1: #skips header line
		continue 
	line = line.rstrip()
	line = line.split()
	tissue = line[5]
	runID = line[4]

	if runID not in runArray:
		runArray.append(runID + '_' + tissue)
	inFile = open(novelDir + runID + '/genes.fpkm_tracking', 'r')
	parse_infile(inFile,contigDict)
write_output(contigDict,runArray,outFile)
outFile.close()
sampleInfo.close()

"""##############################
# Parsing CANFAM + NOVEL
##############################
sampleInfo = open(options.input, 'r')
count = 0 #keeping track of line number
runArray = [] #keeping track of runIDs that have been parsed (SRR...)
outFile = open(canFamNovelOut, 'w')

for line in sampleInfo:
	count += 1
	if count ==1: #skips header line
		continue 
	line = line.rstrip()
	line = line.split()
	tissue = line[5]
	runID = line[4]

	if runID not in runArray:
		runArray.append(runID + '_' + tissue)
	inFile = open(canFamNovelDir + runID + '/genes.fpkm_tracking', 'r')
	#outFile = open(canFamNovelOut, 'w')
	parse_infile(inFile,contigDict)
write_output(contigDict,runArray,outFile)
outFile.close()

print '\n\nDONE!!!!\n\n'
"""

##########################################################################################
##########################################################################################
##########################################################################################

##############################
# Setting up new dictionary
#	for parsing
##############################

for key in contigDict:
	dict[key] = [False,0,0,0,False,0,0,0,False,0,0,0]
##############################
##0 = BLAT Hit (True/False)
#1 = chrom, 2 = start, 3 = end, 
##3 = BLAST Hit (True/False)
#5 = hit (genbank) , 6 = long description of hit , 7 = evalue of hit, 
##8 = Expressed in RNA-Seq (True/False)
#9 = Number of Libraries Expressed


##############################
# Defining TOTAL outfile and 
#	relevant in files
##############################
#File in which all important information will be written that includes:
#	BLAT results, RNA-Seq results, and alignment results
totalFile = open('/home/ampend/kidd-lab/ampend-projects/Novel_Sequence_Analysis/rna-seq/results/TotalParsedContigs_RNASeq_BLAT_AlignmentData.txt', 'w')

#BLAST output file of masked contigs versus non-redundant database
blastFile = open('/home/ampend/kidd-lab/ampend-projects/Novel_Sequence_Analysis/blastx_Novel_To_NR/results/Parsed_MASKED_BLASTx_toNR.txt', 'r')

#BLAT output file that was used to identify redundant contigs 
#SOME MAY NOT BE IN THIS FINAL FILE:
blatFile = open('/home/ampend/kidd-lab/ampend-projects/Novel_Sequence_Analysis/RedundantNovelContigs/results/BLATResultsTable_AllContigs.txt', 'r')


##############################
# Parsing BLAT outputs
##############################
##0 = BLAT Hit (True/False)
#1 = chrom, 2 = start, 3 = end, 
##3 = BLAST Hit (True/False)
#5 = hit (genbank) , 6 = long description of hit , 7 = evalue of hit, 
##8 = Expressed in RNA-Seq (True/False)
#9 = Number of Libraries Expressed

print '\n\nParsing BLAT output:'
lineCount = 0
contigCount = 0
for line in blatFile:
	lineCount += 1
	if lineCount == 1: #skips header
		continue
	
	line = line.rstrip()
	line = line.split()
	
	contig=line[0]
	chrom = line[9]
	start = line[11]
	end = line[12]
	if contig in dict:
		if 'chr' in chrom: #If the contig aligns to a 'CHR', then it becomes True for a possible alignment
			dict[contig][0] = True
		else: #If the contig does not have chr in chr column, then it becomes False for a possible alignment
			dict[contig][0] = False #(default)
		dict[contig][1] = chrom
		dict[contig][2] = start
		dict[contig][3] = end
		contigCount += 1
	else:
		continue
print 'Found BLAT coordinates (or empty coordinates) for %i contigs' % (contigCount)

##############################
# Parsing BLAST outputs
##############################
##0 = BLAT Hit (True/False)
#1 = chrom, 2 = start, 3 = end, 
##4 = BLAST Hit (True/False)
#5 = hit (genbank) , 6 = long description of hit , 7 = evalue of hit, 
##8 = Expressed in RNA-Seq (True/False)
#9 = Number of Libraries Expressed

print '\nParsing BLAST output:'
lineCount = 0
contigCount = 0

for line in blastFile:
	lineCount += 1
	if lineCount == 1: #skips header
		continue

	#line = line.rstrip()
	line = line.split('\t')
	
	contig = line[0]
	hit = line[1]
	evalue = line[11]
	longID = line[13]

	if contig in dict:
		dict[contig][4] = True
		dict[contig][5] = hit
		dict[contig][6] = longID
		dict[contig][7] = evalue
		contigCount += 1
	else:
		continue
print 'Found BLAST information for %i contigs' % (contigCount)

##############################
# Parsing RNASeq outputs
##############################
##0 = BLAT Hit (True/False)
#1 = chrom, 2 = start, 3 = end, 
##4 = BLAST Hit (True/False)
#5 = hit (genbank) , 6 = long description of hit , 7 = evalue of hit, 
##8 = Expressed in RNA-Seq (True/False)
#9 = Number of Libraries Expressed
inFile = open(novelOut, 'r')

print '\nParsing RNA-Seq table output:'
lineCount = 0
contigCount = 0

for line in inFile:
	lineCount += 1
	if lineCount == 1: #skips header
		continue

	line = line.rstrip()
	line = line.split()
	yesCount = 0
	
	contig = line[0]
	
	for k in range(1,len(line)):
		if 'Yes' in line[k]:
			yesCount += 1

	if yesCount >= 3: #If a contig is expressed in 3 or more libraries, I classify this as expressed
		dict[contig][8] = True
	dict[contig][9] = yesCount


##############################
# Writing new dictionary out
#	to outfile
##############################

header = 'ContigID\t' + 'BLAT_Alignment(T/F)\t' + 'Chrom\t' + 'Start\t' + 'End\t'
header = header + 'BLAST_Hit(T/F)\t' + 'HitID\t' + 'LongID\t' + 'Evalue\t'
header = header + 'Expressed(T/F)\t' + 'LibrariesExpressed\t'
totalFile.write(header)

for key in dict:
	totalFile.write('%s\t' % (key))
	totalFile.write('\t'.join(map(str,dict[key])) + '\n')

print '\n\nDONE!!!\n\n\n'





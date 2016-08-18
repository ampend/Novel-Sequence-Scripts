# parse_gtf_file.py
# Tracks expression of each novel contig from
#	Cufflinks/CuffMerge output

import re
import sys 
from optparse import  OptionParser

###############################################################################
USAGE = """
python parse_gtf_file.py
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

#Dictionary for saving ALL data
dict = {}

##############################
# Saving IDs in contig IDs in
#	dictionary
##############################
contigFasta = open(options.fasta, 'r')
print 'Reading in fasta sequences from the novel contig fasta file: \n', options.fasta

contigDict = {}
runDict = {}
count = 0

for line in contigFasta:
	line = line.rstrip()

	if '>' in line:
		count += 1
		m = re.match(r">(.+)", line)
		contigID = m.group(1)
	#Sets default exon counts for all contigs = 0. Will overwrite to greater 
	#	number if its seen higher after parsing transcripts.gtf file
	contigDict[contigID] = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
print 'Read in %s contig IDs into dictionary\n\n' % len(contigDict)

##############################
# Saving run IDs into list
##############################
sampleInfo = open(options.input, 'r')
print 'Reading in sample information from file: \n', options.input
runArray = []
count = 0
for line in sampleInfo:
	count += 1
	if count ==1: #skips header line
		continue 
	line = line.rstrip()
	line = line.split()
	runID = line[4]
	if runID not in runArray:
		runArray.append(runID)	
##############################
# Parsing NOVEL GTF file
##############################
#Defining path for finding cufflinks results
novelDir = '/home/ampend/kidd-lab/ampend-projects/Novel_Sequence_Analysis/rna-seq/results/canFam3.1-Novel/'

keyList = []
exonList = []

keyCount = 0

#Running through each runID (SRR ID)
for i in range(len(runArray)):
	runID = runArray[i] #keeps track of position, should be same position in the contigDict also
	#open GTF file
	gtfFile = open(novelDir + runID + '/transcripts.gtf', 'r')
	contigArray = [] #clears contig Array with each run

	#Processing GTF file:
	for line in gtfFile: #goes through line by line of transcripts.gtf file
		line = line.rstrip()
		line = line.split()
	
		if 'exon' not in line[2]: #only want to keep track of exons. Each entry has one 'transcript' line, and 1 or more exon lines (each exon line pertains to individual exons)
			continue
		
		contig=line[0] #contig ID is in first column
		if contig in contigArray: #If the contig ID has been seen before, then add 1 to exon counts
			exonCount += 1
			contigDict[contig][i] = exonCount #overwrites the exon count in a given index position (matches position of runID) for exon #
		else: #if contigID hasn't been seen before, sets exon counts to 1
			exonCount = 1
			contigArray.append(contig)	
			contigDict[contig][i] = exonCount #changes exon count in dictionary to 1

#Opening outfile to print table to
outFile = open('/home/ampend/kidd-lab/ampend-projects/Novel_Sequence_Analysis/rna-seq/results/ExonCounts/ExonCounts_NovelCufflinks.txt','w')

#Writing Header
outFile.write('ContigID\t')
outFile.write('\t'.join(map(str,runArray)) + '\n')
#Writing out exon counts
for key in contigDict:
	outFile.write('%s\t' % (key))
	outFile.write('\t'.join(map(str,contigDict[key])) + '\n')







########################################################
# 2016-08-04
# write-final-fasta.py
# A. Pendleton
# Finds contigs that pass parsing requirements and
#	writes their FASTA sequences to new files
########################################################

from optparse import  OptionParser
import sys

###############################################################################
USAGE = """
python write-final-fasta.py		
				--directory <Directory to write output files>
				--infile <Novel Contig FASTA>

directory == Directory to write output files
infile == Infile that contains the basic probe stats -- UNCORRECTED

"""

parser = OptionParser(USAGE)
parser.add_option('--infile',dest='infile', help = 'Infile that contains the basic probe stats -- UNCORRECTED')
parser.add_option('--directory',dest='directory', help = 'Directory to write output files')

(options, args) = parser.parse_args()

if options.infile is None:
    parser.error('probe stats infile not given')
if options.directory is None:
    parser.error('output directory not given')
###############################################################################
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
	with open('/home/jmkidd/kidd-lab/ampend-projects/Novel_Sequence_Analysis/NovelSequence/novel.v2.fa.masked') as fasta:
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
infile = options.infile 
inFile = open(infile, 'r')
print '\nReading in parsed data from file: ', infile 

outdir = options.directory
print '\nWriting output files to the directory: ', outdir



#____INPUT FILES_____#
#fastafile - unmasked novel contig fasta file
fastafile = '/home/jmkidd/kidd-lab/ampend-projects/Novel_Sequence_Analysis/NovelSequence/novel.v2.fa'
#fastaFile = open(fastafile, 'r')
print 'Unmasked FASTA sequences from: ', fastafile
#masked novel contig fasta file
maskedfastafile = '/home/jmkidd/kidd-lab/ampend-projects/Novel_Sequence_Analysis/NovelSequence/novel.v2.fa.masked'
#maskedfastaFile = open(maskedfastafile, 'r')
print 'Masked FASTA sequences from: ', maskedfastafile

#____OUTPUT FILES____#
#out_unmasked = open(outdir + 'chrNovel_NonRedundant.fa', 'w')
#out_masked = open(outdir + 'chrNovel_NonRedundant.fa.masked', 'w')

#Defining variables
count = 0 #for keeping track of total contig counts
passList = []

for line in inFile:
	line = line.rstrip()
	line = line.split()
	
	if 'ScaffoldID' in line[0]: #skips header
		continue
	
	count += 1 # number of contigs that were in table
		
	contigID = line[0]
	contigPass = line[1]
	if contigPass is False:
		continue
	else:
		passList.append(contigID)

print '\nNumber of contigs: ', count
print '\nNumber of PASS contigs: ', len(passList)

#PARSING UNMASKED CONTIG FASTA
lineCount = 0
Nstring = 1000 * 'N'


###############################################################
#This goes through line-by-line of the FASTA file of the novel contigs and:
#1) Asks if the novel contig ID is in the "pass list" generated above
#2) If so, it takes the fasta sequence and adds it to the seq list
#3) It then concatenates each fasta with user-defined length of N's (here, 1000)

#For adding 100 'N's between FASTA sequences to generate chrNovel
Nstring = 1000 * 'N'

#Writing unmasked chrNovel 
print 'Writing unmasked fasta sequence...\n'
seqs = []
writingUnmaskedFasta(seqs)

#Writing masked chrNovel 
print 'Writing unmasked fasta sequence...\n'
seqs = []
writingMaskedFasta(seqs)


print 'DONE!!!!\n\n'

"""for line in fastaFile:
	line = line.rstrip()
	line = line.split()
	lineCount += 1 
	
	if '>' in line[0]:
		if lineCount > 1 and contigID in contigPassList:
			'\n'.join(str(s) for s in seq)
			out_unmasked.write('>%s\n%s\n' % (previousID, seq))
		ID = line[0]
		ID = ID.split('>') #splits ID from > symbol (i.e. >zoey-scaffold-686 3009)
		contigID = ID[1]
		previousID = contigID 
		seq = [] #clears sequence list
	else:
		seq.append(line) #appends sequence to line		
	
	if lineCount >1000:
		break
"""		
	
	
	
		














































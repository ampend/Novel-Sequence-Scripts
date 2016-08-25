import sys
from optparse import  OptionParser
import genutils

###############################################################################
#changes based on PS_template_April2010 value
def display_edit_setup(line):
     if '%%DocumentFonts: Helvetica-Bold Helvetica' in line:
         line = '%%DocumentFonts: Arial-Bold Arial\n'
     
     if '%%BoundingBox:' in line:
         line = '%%BoundingBox: 0 0 612 792\n'
     if '/repwidth' in line:
         line = '/repwidth 0.25 cm def\n'
     if '/linkwidth' in line:
         line = '/linkwidth 0 def\n'
     if '/gapwidth' in line:
         line = '/gapwidth 0 def\n'
     if '/leftmargin' in line:
         line = '/leftmargin 0.5 cm def\n'
     if '/bottommargin' in line:
         line = '/bottommargin 2 cm def\n'
     if '/pageheight' in line:
         line = '/pageheight 26 cm def\n'

     if '/Helvetica-Bold findfont' in line:
         line = '/Arial-Bold findfont 7 scalefont setfont\n\n'
          
     if '/Helvetica findfont' in line:
         line = '/Arial findfont 7 scalefont setfont\n'
     if '0 pageheight 2 cm sub moveto' in line:
         line = '0 pageheight 1.5 cm sub moveto\n'
     if '/graphicmargin 17.5' in line:
         line = '/graphicmargin 17.5 cm def\n'
         
     if '(Longest Sequence ' in line:
         line = '\n'    
     if '( Threshold Score ' in line:
         line = '\n'    

     if 'curveto' in line:
         line= line.replace('curveto','lineto')
     return line
###############################################################################
def read_rm_file(rmFileName):
    rmLines = []
    inFile = open(rmFileName,'r')
    for line in inFile:
        if line == '\n':
            continue
        line = line.rstrip()
        line = line.split()
        if line[0] == 'There':
            return []
        if line[0] == 'SW':
            continue
        if line[0] == 'score':
            continue
        rmLines.append(line)
    inFile.close()
    return rmLines
###############################################################################
def repeat_class_to_name(r):
    if 'SINE' in r:
        return 'SINE'
    if 'ARTEFACT' in r:
        return 'ARTEFACT'
    if 'DNA' in r:
        return 'DNA'
    if 'LINE' in r:
        return 'LINE'
    if 'Low_complexity' in r:
        return 'Low_complexity'
    if 'LTR' in r:
        return 'LTR'
    if 'Other' in r:
        return 'Other'
    if 'rRNA' in r:
        return 'rRNA'
    if 'scRNA' in r:
        return 'scRNA'
    if 'snRNA' in r:
        return 'snRNA'
    if 'srpRNA' in r:
        return 'srpRNA'
    if 'tRNA' in r:
        return 'tRNA'
    if 'RNA' in r:
        return 'RNA'
    if 'Satellite' in r:
        return 'Satellite'
    if 'Simple_repeat' in r:
        return 'Simple_repeat'
    if 'Unknown' in r:
        return 'Unknown'
    if 'Retroposon' in r:
        return 'SINE'    
    print 'repeat class unknown for',r
    return 'Unknown'
###############################################################################
def sort_and_merge_repeats(repeats):
    repeats.sort()
    newRep = []
    for r in repeats:
        s = r[0]
        e = r[1]
        orientation = r[2]
        repClass = r[3]
        if len(newRep) == 0:
            newRep.append(r)
        else:
            lr = newRep[-1]
            ls = lr[0]
            le = lr[1]
            lorient = lr[2]
            lrepClass = lr[3]
            # overlap, need to extend
            if le > s and lrepClass == repClass and lorient == orientation:
                n = [ls,e,orientation,repClass]
                newRep[-1] = n
            else:
                newRep.append(r)
    return newRep    
###############################################################################
def process_repeat_file(fn):
    repeatLines = read_rm_file(fn)
    repeats = []
    for R in repeatLines:
        s = int(R[5])
        e = int(R[6])
        orientation = R[8]
        
        repClass = R[10]
        if orientation == 'C':
            orientation = '-'
        repClass = repeat_class_to_name(repClass)
        reps = [s,e,orientation,repClass]
        repeats.append(reps)
    repeats = sort_and_merge_repeats(repeats)
    return repeats
###############################################################################
def process_gap_file(fn):
    gapLines = []
    inFile = open(fn,'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        s = int(line[1])
        e = int(line[2])
        gapLines.append([s,e])
    inFile.close()
    return gapLines
###############################################################################
def process_exon_file(fn):
    exonLines = []
    inFile = open(fn,'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        s = int(line[1])
        e = int(line[2])
        name = line[3]
        exonLines.append([name,s,e])
        exonLines.sort()
    inFile.close()
    return exonLines  
###############################################################################
def process_primer_file(fn):
    primerLines = []
    inFile = open(fn,'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        s = int(line[1])
        e = int(line[2])
        dir = line[3]
        primerLines.append([s,e,dir])
    inFile.close()
    return primerLines
###############################################################################
USAGE = """
python annotate-miropeats-3seqs-AP.py    
--miroin <in miropeat ps file>  
--topRM <RM of top> 
--bottomRM <RM of bottom>
--topName <name of top> 
--bottomName <name of bottom>
--blat <blat output of canFam vs fosmid, blast9 output format>
"""

parser = OptionParser(USAGE)
parser.add_option('--miroin',dest='miropeatsInput',help='input file of miropeats ps file')
parser.add_option('--topRM',dest='topRM',help='repeat mask out file for top sequence')
parser.add_option('--bottomRM',dest='bottomRM',help='repeat mask out file for bottom sequence')
parser.add_option('--topName',dest='topName',help='name of top seqence')
parser.add_option('--bottomName',dest='bottomName',help='name of bottom seqence')
parser.add_option('--blat',dest='blat',help='blat of canFam versus fosmid results, blast9 output format')

(options,args)=parser.parse_args()
if options.miropeatsInput is None:
    parser.error('miropeats input file name not given')
if options.topRM is None:
    parser.error('top RM input file name not given')
if options.bottomRM is None:
    parser.error('bottom RM input file name not given')
if options.topName is None:
    parser.error('top seq name not given')
if options.bottomName is None:
    parser.error('bottom seq name not given')
if options.blat is None:
    parser.error('blat output not given, blast9 output format')
    
###############################################################################
bottomName = options.bottomName
##### Info files #####
#PSTemplate = '/home/jmkidd/kidd-lab/jmkidd-projects/people-projects/jwilds-projects/RetroSeq-HGDP/align/scripts/PS_template_July2008.txt' 
PSTemplate = '/home/jmkidd/kidd-lab/jmkidd-projects/people-projects/jwilds-projects/RetroSeq-HGDP/align/scripts/PS_template_July2008_withlines.txt'

#### setup the repeat to color dictionary#####
color = {}
color['Other'] = 'Black'
color['Simple Repeat'] = 'DkGray'
color['Low Complexity'] = 'LtGray'
color['DNA'] = 'Pink'
color['LTR'] = 'Orange'
color['LINE'] = 'Green'
color['SINE'] = 'Purple'
color['GAP'] = 'Black' #AP: Changed
color['EXON'] = 'Blue' #AP: Changed

#height of ref and fos sequence, in "Y" coordinates
#fos is what is on the bottom
#ref_line=1.2
#fos_line=0.4
#####bottomLine = 0.3
bottomLine = .6 #0.3
#middleLine = 0.8
####topLine = 1.0
topLine = 1.3
# for direction...
#arrowEndBp = 1000
arrowEndBp = 100

###############################################################################

miropeatsInFile = options.miropeatsInput
topRepeatMaskFile = options.topRM
bottomRepeatMaskFile = options.bottomRM
miropeatsOutFile = miropeatsInFile + '.annotated.ps'

# these are for coloring breakpoints... we will skip this for now...
chrm_breaks = {}
chrm_breaks['not do'] = 1
clone_breaks = {}
clone_breaks['not do'] = 1

inFile = open(miropeatsInFile,'r')
outFile = open(miropeatsOutFile,'w')

# read in to tagends function, then add in the template
while True:
    line = inFile.readline()
    if line[0:8] == '/tagends':
        break
    else:
        line = display_edit_setup(line)
        outFile.write(line)
#here, line is the begin of the tagends function
# here is a good place to print out other info

outFile.write('/Arial-Bold findfont 7 scalefont setfont\n')
#outFile.write('0 pageheight 1.5 cm sub moveto\n')
#outFile.write('0 pageheight 2.0 cm sub moveto\n')

print 'Adding PS template info.....'
inTEMP = open(PSTemplate,'r')
for t in inTEMP:
    outFile.write(t)
inTEMP.close()
print 'added the template info!'
# to put us back...
outFile.write('/Arial findfont 8 scalefont setfont\n')    

#################################################
bedfile = open('/home/jmkidd/kidd-lab/ampend-projects/fosmids/Processing/miropeats/Total_FosmidcanFamCoordinates.bed', 'r')

for d in bedfile:
	d = d.rstrip()
	d = d.split()
	if options.bottomName in d[3]:
		#Print header line that shows the name of the fosmid and the coordinates
		outFile.write('(%s -- %s:%s-%s) show' % (options.bottomName,d[0],d[1],d[2]))

#################################################
#Identifies correction factor of start of Zoey fosmid vs canFam3
blatFile = open(options.blat,'r')

for b in blatFile:
	b = b.rstrip()
	b = b.split()
	if '#' in b[0]:
		continue
	if 'Zoey' in b[0] and 'canFam' in b[1]:#hard coded for my purposes, fosmid in column 1 and the canFam in column 2
		offset = int(b[8])
		print 'canFam/fosmid offset equals:', offset
		break

#################################################
# print out rest of info down to where coords will go
outFile.write(line)
started_pc = 0
more  = 1
while more == 1:
    line = inFile.readline()
    
    if 'printcontig' in line: #AP: Reads through the coordinate alignment file and assigns
    							# bottom or top line to 1 or 0	
        line = line.rstrip()
        started_pc = 1
        f = line.split()
        # figure out the line
        if options.topName in f[1]:
            f[0] = topLine
        elif options.bottomName in f[1]:
            f[0] = bottomLine
        else:
            print line
            print 'where am i?'
            sys.exit()

        if options.topName in f[6]:
            f[5] = topLine
        elif options.bottomName in f[6]:
            f[5] = bottomLine
        else:
            print line
            print 'where am i?'
            sys.exit()
        
        p_cmd = 'printcontig'

        # The below is if we are going to color potential chromosome breakpoints that we have called based on the
        # miropeats output directly
        #chrom first
        if 'canFam' in f[1]: #AP: I Changed this from 'chr' for my formatting
            if f[2] in chrm_breaks and f[7] in clone_breaks:
                p_cmd = 'printcontig_left'
            if f[3] in chrm_breaks and f[8] in clone_breaks:
                p_cmd = 'printcontig_right'
            if f[3] in chrm_breaks and f[8] in clone_breaks and f[2] in chrm_breaks and f[7] in clone_breaks:
                p_cmd = 'printcontig_both'
        else:
            if f[7] in chrm_breaks and f[2] in clone_breaks:
                p_cmd = 'printcontig_left'
            if f[8] in chrm_breaks and f[3] in clone_breaks:
                p_cmd = 'printcontig_right'
            if f[3] in chrm_breaks and f[8] in clone_breaks and f[2] in chrm_breaks and f[7] in clone_breaks:
                p_cmd = 'printcontig_both'
        
        print p_cmd
        f[10] = p_cmd
        f = [str(j) for j in f]
        line = ' '.join(f) + '\n'
    if (('printcontig' in line) is False) and (started_pc == 1):  # means that we are done with the printcontig commands
        more = 0
    if more == 1:
        outFile.write(line)

#################################################
#Identifies gaps in the canFam3 alignment from bedtools intersect output
exonFile = open('/home/jmkidd/kidd-lab/ampend-projects/fosmids/Processing/miropeats/EXONS_Intersect_withFosmids.txt','r') #HARD CODED!!!
out1 = open('/home/jmkidd/kidd-lab/ampend-projects/fosmids/Processing/miropeats/%s/exons_canFam.exons'%(bottomName), 'w')
out2 = open('/home/jmkidd/kidd-lab/ampend-projects/fosmids/Processing/miropeats/%s/exons_Zoey.exons'%(bottomName), 'w')


for e in exonFile:
	e = e.rstrip()
	e = e.split()
	#Just in case the wrong contig ID was input that misses the 'canFam' in the contig ID
	#	this step will add the canFam into the name
	if 'canFam' not in e[3]:
		contigID = 'canFam_' + e[3]
	#if the contig ID is in the bedfile
	if options.topName in contigID:
		exonStart = int(e[5]) 
		exonEnd = int(e[6]) 
		exonName = e[7]
		canFamStart = int(e[1]) 
		canFamEnd = int(e[2])
		#canFam coordinates of exon
		start = exonStart - canFamStart + 1
		end = exonEnd - canFamStart + 1 
		out1.write('%s\t%i\t%i\t%s\n' % (options.topName, start, end, exonName)) 
		#fosmid coordinates of exon
		start = exonStart - canFamStart - offset + 1
		end = exonEnd - canFamStart - offset + 1
		if start > 0 and end < canFamEnd: 
			out2.write('%s\t%i\t%i\t%s\n' % (options.bottomName, start, end, exonName))
out1.close()
out2.close()		

#################################################
#Identifies aCGH probes in the canFam3 alignment from bedtools intersect output
probeFile = open('/home/jmkidd/kidd-lab/ampend-projects/fosmids/Processing/miropeats/PROBES_Intersect_withFosmids.txt','r') #HARD CODED!!!
out1 = open('/home/jmkidd/kidd-lab/ampend-projects/fosmids/Processing/miropeats/%s/probes_canFam.probes'%(bottomName), 'w')
out2 = open('/home/jmkidd/kidd-lab/ampend-projects/fosmids/Processing/miropeats/%s/probes_Zoey.probes'%(bottomName), 'w')

for e in probeFile:
	e = e.rstrip()
	e = e.split()
	#Just in case the wrong contig ID was input that misses the 'canFam' in the contig ID
	#	this step will add the canFam into the name
	if 'canFam' not in e[3]:
		contigID = 'canFam_' + e[3]
	#if the contig ID is in the bedfile
	if options.topName in contigID:
		#print options.topName
		probeStart = int(e[5]) 
		probeEnd = int(e[6]) 
		probeName = e[7]
		canFamStart = int(e[1]) 
		canFamEnd = int(e[2])
		#canFam coordinates of probe
		#start = probeStart - canFamStart + 7000 + 1 
		#end = probeEnd - canFamStart + 7000 + 1 
		start = probeStart - canFamStart +  1
		end = probeEnd - canFamStart +  1
		out1.write('%s\t%i\t%i\t%s\n' % (options.topName, start, end, probeName)) 
	 	print '%s\t%i\t%i\t%s\t' % (options.topName, start, end, probeName)
		#fosmid coordinates of probe
		start = probeStart - canFamStart + offset + 1
		end = probeEnd - canFamStart + offset + 1
		if start > 0 and end < canFamEnd: 
			out2.write('%s\t%i\t%i\t%s\n' % (options.bottomName, start, end, probeName))
		 	print '%s\t%i\t%i\t%s\n' % (options.bottomName, start, end, probeName)
out1.close()
out2.close()
        
#OK now time for repeats info                       
print 'Annotating Reference (top) Sequence...'
print 'Adding repeats...'
###############################################################################
# put reading in RM *out file, merging, etc into function below
repeats = process_repeat_file(topRepeatMaskFile)
print 'did sort and the population'
print 'adding repeats'
n = 0
ypos = topLine + 0.03
for R in repeats:
    rs = R[0]
    re = R[1]
    rd = R[2]
    rt = R[3]
    over = 1 # has to be on this chrm at this point, so just force it
    # may change later if we want to draw subset of coords
    if over == 1:
        if rt in color:
            c = color[rt]
        else:
            c = color['Other']
        eLen = re - rs  + 1
        outFile.write('%s\n' % c)
        if rd == '+':
            if eLen < arrowEndBp:
                outFile.write('%f %i %i RIGHT_end\n' % (ypos,rs,re))
                n += 1
            else:
                ar = re - arrowEndBp
                outFile.write('%f %i %i exon\n' % (ypos,rs,ar))
                outFile.write('%f %i %i RIGHT_end\n' % (ypos,ar,re))
        elif rd == '-':
            if eLen < arrowEndBp:
                outFile.write('%f %i %i LEFT_end\n' % (ypos,rs,re))
                n += 1
            else:
                ar = rs + arrowEndBp
                outFile.write('%f %i %i LEFT_end\n' % (ypos,rs,ar))
                outFile.write('%f %i %i exon\n' % (ypos,ar,re))
        else:
            print 'repeat error!  what dir???'
            print R
            sys.exit()
outFile.write('Black\n')
tpos = ypos + 0.02
outFile.write('/Arial findfont 9 scalefont setfont\n')
outFile.write('%f (Repeats) printname_right\n' % (tpos))


###############################################################################

# would do DUPMASK here, if we had it

###############################################################################
# print out the repeats key, on the left
m = topLine + 5*0.05
ys = m + 0.1 # print key 2 lines above htis
outFile.write('/yh drop %f mul def\n' % (ys)) # transform to right coordinates
# print key at yh
outFile.write('0 yh moveto\n')
outFile.write('Black\n')
outFile.write('(Other     ) 100 string cvs show\n')
outFile.write('DkGray\n')
outFile.write('(Simple Repeat     ) 100 string cvs show\n')
outFile.write('LtGray\n\n')
outFile.write('(Low Complexity     ) 100 string cvs show\n')
outFile.write('Pink\n')
outFile.write('(DNA     ) 100 string cvs show\n')
outFile.write('Orange\n')
outFile.write('(LTR     ) 100 string cvs show\n')
outFile.write('Green\n')
outFile.write('(LINE     ) 100 string cvs show\n')
outFile.write('Purple\n')
outFile.write('(SINE     ) 100 string cvs show\n')
outFile.write('Black\n')
###############################################################################

#Now annotating the middle allele
"""
#print 'Now annotating the (middle)'
#repeats = process_repeat_file(middleRepeatMaskFile)
#print 'Did sort'
#print 'Adding repeats....'
#print 'adding repeats'
#n = 0
#ypos = middleLine+0.05
#ypos = middleLine+0.03
#for R in repeats:
    #rs = R[0]
    re = R[1]
    rd = R[2]
    rt = R[3]
    over = 1 # has to be on this chrm at this point, so just force it
    # may change later if we want to draw subset of coords
    if over == 1:
        if rt in color:
            c = color[rt]
        else:
            c = color['Other']
        eLen = re - rs  + 1
        outFile.write('%s\n' % c)
        if rd == '+':
            if eLen < arrowEndBp:
                outFile.write('%f %i %i RIGHT_end\n' % (ypos,rs,re))
                n += 1
            else:
                ar = re - arrowEndBp
                outFile.write('%f %i %i exon\n' % (ypos,rs,ar))
                outFile.write('%f %i %i RIGHT_end\n' % (ypos,ar,re))
        elif rd == '-':
            if eLen < arrowEndBp:
                outFile.write('%f %i %i LEFT_end\n' % (ypos,rs,re))
                n += 1
            else:
                ar = rs + arrowEndBp
                outFile.write('%f %i %i LEFT_end\n' % (ypos,rs,ar))
                outFile.write('%f %i %i exon\n' % (ypos,ar,re))
        else:
            print 'repeat error!  what dir???'
            print R
            sys.exit()
outFile.write('Black\n')
tpos = ypos + 0.02
outFile.write('/Arial findfont 9 scalefont setfont\n')
outFile.write('%f (Repeats) printname_right\n' % (tpos))
"""
######## GAPS ##########

print 'Now annotating gaps in the canFam contig...\n'
topGaps = '/home/jmkidd/kidd-lab/ampend-projects/fosmids/Processing/miropeats/%s/gaps_canFam.gaps'%(bottomName)
#topGaps = 'gaps_canFam.gaps'

gaps = process_gap_file(topGaps)
n = 0
ypos = topLine + 0.09
gap_pos = topLine + 0.09
for R in gaps:
    rs = R[0]
    re = R[1]
    over = 1 # has to be on this chrm at this point, so just force it
    if over == 1:
        c = 'Black'
        eLen = re - rs  + 1
        outFile.write('%s\n' % c)
        outFile.write('%f %i %i exon\n' % (ypos,rs,re))
        n += 1
outFile.write('Black\n')
tpos = ypos + 0.02
outFile.write('%f (Gaps) printname_right\n' % (tpos))

#########################
######## EXONS ##########
#########################

exonList = []

print 'Now annotating exons in the canFam contig...\n'

#topExons = 'exons_canFam.exons'
topExons = '/home/jmkidd/kidd-lab/ampend-projects/fosmids/Processing/miropeats/%s/exons_canFam.exons'%(bottomName) 

exonLines = []
#exons = process_gap_file(topExons)
exons = process_exon_file(topExons)

n = 0
#ypos = topLine + 0.15
#ypos = topLine + .05
ypos = gap_pos
exon_pos = gap_pos + 0.03
exons.sort()

for R in exons:
	rs = R[1]
	re = R[2]
	rname = R[0]
	if rname in exonList:
		print R
		over = 1 # has to be on this chrm at this point, so just force it
		if over == 1:
			c = 'Blue'
			eLen = re - rs  + 1
			outFile.write('%s\n' % c)
			outFile.write('%f %i %i exon\n' % (ypos,rs,re))
			n += 1
			outFile.write('Blue\n')
			tpos = ypos + 0.02
			outFile.write('%f (%s) printname_right\n' % (tpos, rname))
	else:
		exonList.append(rname)
		ypos = ypos + 0.03
		print R
		over = 1 # has to be on this chrm at this point, so just force it
		if over == 1:
			c = 'Blue'
			eLen = re - rs  + 1
			outFile.write('%s\n' % c)
			outFile.write('%f %i %i exon\n' % (ypos,rs,re))
			n += 1		
			outFile.write('Blue\n')
			tpos = ypos + 0.02
			outFile.write('%f (%s) printname_right\n' % (tpos, rname))
exonPos = ypos

##########################
######## PROBES ##########
##########################

print 'Now annotating aCGH probes in the canFam contig...\n'
#topProbes = 'probes_canFam.exons'
topProbes = '/home/jmkidd/kidd-lab/ampend-projects/fosmids/Processing/miropeats/%s/probes_canFam.probes'%(bottomName)

probes = process_gap_file(topProbes)
n = 0
#ypos = topLine + 0.09
probe_pos = exonPos + 0.03
ypos  = exonPos + 0.03
for R in probes:
    rs = R[0]
    re = R[1]
    over = 1 # has to be on this chrm at this point, so just force it
    if over == 1:
        c = 'Red'
        eLen = re - rs  + 1
        outFile.write('%s\n' % c)
        outFile.write('%f %i %i exon\n' % (ypos,rs,re))
        n += 1
outFile.write('Red\n')
tpos = ypos + 0.02
outFile.write('%f (Probes) printname_right\n' % (tpos))

###############################################################################
#Now annotating the bottom allele
###############################################################################

print 'Now annotating the (bottom)'
repeats = process_repeat_file(bottomRepeatMaskFile)
print 'Did sort'
print 'Adding repeats....'
print 'adding repeats'
n = 0
ypos = bottomLine-0.15
for R in repeats:
    rs = R[0]
    re = R[1]
    rd = R[2]
    rt = R[3]
    over = 1 # has to be on this chrm at this point, so just force it
    # may change later if we want to draw subset of coords
    if over == 1:
        if rt in color:
            c = color[rt]
        else:
            c = color['Other']
        eLen = re - rs  + 1
        outFile.write('%s\n' % c)
        if rd == '+':
            if eLen < arrowEndBp:
                outFile.write('%f %i %i RIGHT_end\n' % (ypos,rs,re))
                n += 1
            else:
                ar = re - arrowEndBp
                outFile.write('%f %i %i exon\n' % (ypos,rs,ar))
                outFile.write('%f %i %i RIGHT_end\n' % (ypos,ar,re))
        elif rd == '-':
            if eLen < arrowEndBp:
                outFile.write('%f %i %i LEFT_end\n' % (ypos,rs,re))
                n += 1
            else:
                ar = rs + arrowEndBp
                outFile.write('%f %i %i LEFT_end\n' % (ypos,rs,ar))
                outFile.write('%f %i %i exon\n' % (ypos,ar,re))
        else:
            print 'repeat error!  what dir???'
            print R
            sys.exit()
outFile.write('Black\n')
tpos = ypos + 0.02
outFile.write('/Arial findfont 9 scalefont setfont\n')
outFile.write('%f (Repeats) printname_right\n' % (tpos))

########################
######## GAPS ##########
########################

print 'Now annotating gaps in fosmid...\n'
#bottomGaps = 'gaps_Zoey.gaps'
bottomGaps = '/home/jmkidd/kidd-lab/ampend-projects/fosmids/Processing/miropeats/%s/gaps_Zoey.gaps'%(bottomName)

gaps = process_gap_file(bottomGaps)

n = 0
print 'GAP start ypos',ypos
ypos = bottomLine - 0.18

for R in gaps:
    rs = R[0]
    re = R[1]
    over = 1 # has to be on this chrm at this point, so just force it
    if over == 1:
        c = 'Black'
        eLen = re - rs  + 1
        outFile.write('%s\n' % c)
        outFile.write('%f %i %i exon\n' % (ypos,rs,re))
        n += 1
outFile.write('Black\n')
tpos = ypos + 0.02
print 'GAP end ypos',ypos
gap_pos = ypos
outFile.write('%f (Gaps) printname_right\n' % (tpos))

#################################
######## NOVEL CONTIGS ##########
#################################

print 'Now annotating novel contigs in fosmid...\n'

#novelContigFasta = '~/kidd-lab/ampend-projects/Novel_Sequence_Analysis/NovelSequence/novel.v2.fa.masked'
#New non-reundant Fasta
novelContigFasta = '~/kidd-lab/ampend-projects/Novel_Sequence_Analysis/RedundantNovelContigs/Final_chrNovel_Fasta/novelContigs_NonRedundant.fa.masked'
bottomRM = options.bottomRM
masked_bottomRM = bottomRM.replace(".out",".masked")
contigfile = 'BLAT_novelContigs_vs_fosmid.blat'
 
#cmd = 'blat -fine -minMatch=1 -minScore=10 -out=blast9 %s %s %s' % (novelContigFasta,masked_bottomRM,contigfile)
cmd = 'blat -noHead %s %s %s' % (novelContigFasta,masked_bottomRM,contigfile)
print cmd
genutils.runCMD(cmd) 
n = 0
#ypos = bottomLine - 0.25
exon_pos = gap_pos - 0.03
ypos = exon_pos

contigFile = open('BLAT_novelContigs_vs_fosmid.blat','r')

contigList = []

print 'NOVEL CONTIG start ypos', ypos 

for b in contigFile:
	b = b.rstrip()
	b = b.split()
	if b[0].isdigit() is False:
		continue
  	contigLen = int(b[14])
  	match = int(b[0])
	prop = float(match)/float(contigLen) 	
	ID = b[13]
	ns = int(b[11])
	ne = int(b[12])
	if prop > 0.9 or int(b[15]) == 0 and prop > 0.5 or int(b[16])==contigLen and prop > 0.5 :
		tpos=ypos
		print '\nNovel contig with >90% alignment to fosmid:'
		if ID in contigList:
			append.contigList(ID)
			over = 1 # has to be on this chrm at this point, so just force it
			if over == 1:
				c = 'Purple'
				eLen = ne - ns  + 1
				outFile.write('%s\n' % c)
				outFile.write('%f %i %i exon\n' % (ypos,ns,ne))
				print '%f %i %i exon\n' % (ypos,ns,ne)
				n += 1
				outFile.write('Purple\n')
				tpos = ypos + 0.015
				outFile.write('/Arial findfont 7 scalefont setfont\n')
				outFile.write('%f (%s) printname_right\n' % (tpos, ID))
				print '%f (%s) printname_right\n' % (ypos, ID)
				print '%s\t%i\t%i' % (ID, ne, ns)
				
		else:
			contigList.append(ID)
			ypos = ypos - 0.03
			print b[13]
			over = 1 # has to be on this chrm at this point, so just force it
			if over == 1:
				c = 'Purple'
				eLen = ne - ns  + 1
				outFile.write('%s\n' % c)
				outFile.write('%f %i %i exon\n' % (ypos,ns,ne))
				print '%f %i %i exon\n' % (ypos,ns,ne)
				n += 1		
				outFile.write('Purple\n')
				tpos = ypos + 0.015
				outFile.write('/Arial findfont 7 scalefont setfont\n')
				outFile.write('%f (%s) printname_right\n' % (tpos, ID))
				print '%f (%s) printname_right\n' % (ypos, ID)
				print '%s\t%i\t%i' % (ID, ne, ns)
print 'NOVEL CONTIG end ypos', ypos 
				
print 'Done adding novel contigs\n'
	
outFile.write('Black\n')
tpos = ypos + 0.02



"""
print 'Now annotating novel contigs in fosmid...\n'

novelContigFasta = '~/kidd-lab/ampend-projects/Novel_Sequence_Analysis/NovelSequence/novel.v2.fa.masked'
bottomRM = options.bottomRM
masked_bottomRM = bottomRM.replace(".out",".masked")
contigfile = '/home/jmkidd/kidd-lab/ampend-projects/fosmids/Processing/miropeats/%s/BLAT_novelContigs_vs_fosmid.blat'%(bottomName)

cmd = 'blat -noHead %s %s %s' % (novelContigFasta,masked_bottomRM,contigfile)
print cmd
genutils.runCMD(cmd) 
n = 0
print 'bottomLine =', bottomLine
print 'CONTIG start ypos=', ypos
ypos = bottomLine - 0.30
print 'CONTIG end ypos=', ypos

contigFile = open(contigfile,'r')
contigList = []
for b in contigFile:
	b = b.rstrip()
	b = b.split()
	if b[0].isdigit() is False:
		continue
  	contigLen = int(b[14])
  	match = int(b[0])
	prop = float(match)/float(contigLen) 	
	ID = b[13]
	ns = int(b[11])
	ne = int(b[12])
	if prop > 0.9 or prop > 0.5 and int(b[15]) == 0 or prop > 0.5 and int(b[16]) == contigLen:
		print 'Novel contigs with >90% alignment to fosmid:\n'
		if ID in contigList:
			c = 'Purple'
			eLen = ne - ns  + 1
			outFile.write('%s\n' % c)
			outFile.write('%f %i %i exon\n' % (ypos,ns,ne))
			n += 1
			outFile.write('Purple\n')
			tpos = ypos + 0.02
			outFile.write('%f (%s) printname_right\n' % (tpos, ID))
			print '%s\t%i\t%i' % (ID, ne, ns)
			print tpos
		else:
			contigList.append(ID)
			c = 'Purple'
			eLen = ne - ns  + 1
			outFile.write('%s\n' % c)
			outFile.write('%f %i %i exon\n' % (ypos,ns,ne))
			n += 1
			outFile.write('Purple\n')
			tpos = ypos + 0.02
			outFile.write('%f (%s) printname_right\n' % (tpos, ID))
			print '%s\t%i\t%i' % (ID, ne, ns)
			print tpos
		ypos = tpos

print 'Done adding novel contigs\n'
	
outFile.write('Black\n')
tpos = ypos + 0.02
#outFile.write('%f (Novel Contigs) printname_right\n' % (tpos))
"""
####################################################

# Dupmask would go here if we were doing it for the fosmid


######ADD IN THE BREAK POINTS#####################

outFile.write(str(line)) # last line we didn't do yet, should be blank
outFile.write('Black\n')

# set scale bar at the appropiate position

outFile.write('/y1 drop %f mul def\n' % bottomLine)
outFile.write('/ybot y1 -0.3 cm add def\n')
outFile.write('/ybot2 y1 -0.5 cm add def\n')

while True:
    l = inFile.readline()
    if l == '':
        break
    # replace hard coded positions with pos relative to the fos_line
    l = l.replace('-0.3 cm','ybot')
    l = l.replace('-.3 cm','ybot')

    l = l.replace('-0.5 cm','ybot2')
    l = l.replace('-.5 cm','ybot2')
    outFile.write(l)


inFile.close()
outFile.close()

print 'DONE with annotation adding!'
# ok, now I think we may be done....


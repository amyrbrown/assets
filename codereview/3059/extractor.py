################################################
# Base count extractor                         #
# A simple script using the pySAM library      #
# to access SAMTools to obtain high quality    #
# base counts from an aligned BAM file         #
#                                              #
# david.eyre@ndm.ox.ac.uk                      #
# 04 March 2013                                #
################################################


### options specified here

#numeric list of positions of variable sites in reference genome, numbered from zero
POSLIST = [0, 1, 1234]

#path to bam files for extraction of base counts
BAMPATH = ''

#list of file names for bam files, without .bam suffix, these will be recycled for output names
FILES = ['','']

#output path
OUTPATH = ''

### end of options





import pysam, os

def Bases_At_Pos(samfile, pos, chromname, minbasequal, minmapqual):
    'Return a string of the bases at that position.'
    position = 0
    coverage = 0
    bases = ""
    for pileupcolumn in samfile.pileup(reference=chromname, start=pos-1, end=pos):
        if ((pileupcolumn.pos+1) >= pos and (pileupcolumn.pos+1) <= pos):
            position = int(pileupcolumn.pos+1)
            coverage = int(pileupcolumn.n)
            for pileupread in pileupcolumn.pileups:
                if (pileupread.indel == 0 and pileupread.is_del == 0 and \
                (ord(pileupread.alignment.qual[pileupread.qpos])-33) >= minbasequal and \
                float(pileupread.alignment.mapq) >= minmapqual):
                    bases += pileupread.alignment.seq[pileupread.qpos]
    return position, coverage, bases

POSLIST.sort()

#loop over each of the bam files
for file in FILES:
	inbampath = BAMPATH+file+'.bam'
	if not os.path.exists(inbampath): 
		#skip this file
		continue
	
	#open bam file with pySAM library
	samfile = pysam.Samfile(inbampath, 'rb')
	
	#open output file and write headers
	with open (OUTPATH+file+'_bc.txt', 'w') as f:
		f.write('site\tA\tC\tG\tT\n')
		#for each position of interest write base counts
		for pos in POSLIST:
			hq_position, hq_coverage, hq_bases = Bases_At_Pos(samfile, int(pos), ref, 30.0, 30.0)
			hq_basecounts = [ len([x for x in hq_bases if x == 'A']), len([x for x in hq_bases if x == 'C']), len([x for x in hq_bases if x == 'G']), len([x for x in hq_bases if x == 'T']) ]
			if(len([i for i in hq_basecounts if i>0])>=1):
				f.write('%s\t%s\n'%(hq_position, '\t'.join([str(i) for i in hq_basecounts])))
		f.close()
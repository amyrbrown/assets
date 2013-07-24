"""
This file is part of DEAP.

DEAP is free software: you can redistribute it and/or modify it under 
the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

DEAP is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with DEAP.  If not, see <http://www.gnu.org/licenses/>.

This file was written by Winston Haynes. If you have any questions, 
please contact me at: winston.haynes@seattlechildresn.org
"""

import gc
import re
import random
import stats
import time
import os
import sys
from sets import Set
from pathwayStructures import *
import rpy2.robjects as robjects
from deapHelpers import *
import threading
import Queue
import logging
import copy
import string

gc.enable()
###Read in user arguments
if len(sys.argv)==6:
	pwd=sys.argv[1]
	subdir=sys.argv[2]
	datadir=sys.argv[3]
	edgedir=sys.argv[4]
	isPaired=string.lower(sys.argv[5])
else:
	print 'Incorrect number of arguments. aborting.\nProper format is: python DEAP.py [working directory] [analysis name] [data directory] [edge file absolute path] [is data paired? y/N]'
	sys.exit()

if isPaired=='Y' or isPaired=='y':
	setIsPaired(True)
else:
	setIsPaired(False)

###Get all expression files
expressionFiles=[datadir.split('/')[-1]]
datadir=re.sub(expressionFiles[0],'',datadir)

###Check for existence of all output directories. If they don't exist, make them.
outputdir=os.path.join(pwd,'output/')
if not os.path.exists(outputdir):
	os.makedirs(outputdir)
logdir=os.path.join(outputdir,subdir,'logs/')
if not os.path.exists(logdir):
	os.makedirs(logdir)
outdir=os.path.join(outputdir,subdir,'results/')
if not os.path.exists(outdir):
	os.makedirs(outdir)
nulldir=os.path.join(outputdir,subdir,'null/')
if not os.path.exists(nulldir):
	os.makedirs(nulldir)

###Set up logging
beginTime=re.sub(r'\.','',str(time.time()))
logging.basicConfig(filename=os.path.join(logdir,subdir+'.log'),level=logging.INFO, format='[%(levelname)s] (%(threadName)-10s) %(message)s',)

###Iterate over all expression files
for expressionData in expressionFiles:	
	###Open output file, write header line, define function for future writes
	filename=expressionData.split('.')[0]
	outfile=open(os.path.join(outdir,filename)+'.txt','w')
	outfile.write('PathwayName\t')
	outfile.write('||AbsValue\tpVal\tnullMean\tnullStDev\t||PathwaySubset\n')
	def writeResult(result):
		xmlfile=result[0]
		statvals=result[1]
		outfile.write(xmlfile+'\t')
		outfile.write(str(statvals[2].curval)+'\t'+str(statvals[2].qval)+'\t'+str(statvals[2].mean)+'\t'+str(statvals[2].stdev)+'\t')
		outfile.write(str(statvals[2].pathSubset)+'\n')
		outfile.flush()

	###Load expression data
	expdict = createExpressionDict(os.path.join(pwd,datadir,expressionData))

	###Perform data manipulations. In this case, data rotation
	timeb=time.time()
	logging.info('rotating data')
	rotatedData=buildRotatedData(expdict)
	logging.info('data rotated in '+str(time.time()-timeb)+' secs.')
	
	### Perform transforms on data based on user input

	expdict=meanExpDict(expdict)
	rotatedData=meanDictList(rotatedData)

	###Define threadable class for performing all operations on each xmlfile
	class xmlFileBuilder(threading.Thread):
		###Initialize with no result
		def __init__(self,xmlfile):
			self.xmlfile=xmlfile
			self.result=None
			threading.Thread.__init__(self,name=xmlfile)

		###Called after thread completion to get the result
		def get_result(self):
			return [self.xmlfile,self.result]

		###The thread processing
		def run(self):
			timea=time.time()
			logging.info('Building for file: '+self.xmlfile)
			protedges=loadEdgesFromFile(os.path.join(edgedir,self.xmlfile))
			logging.info('Graph built for '+self.xmlfile+' in '+str(time.time()-timea)+' secs')
			timea=time.time()
			scores=calculateScores(expdict,protedges)
			logging.info('Scores for '+self.xmlfile+' calculated in '+str(time.time()-timea)+' secs:')
			for score in scores:
				logging.info('\tmin= '+str(score[0])+'\tmax= '+str(score[1]))
			nullfile=os.path.join(nulldir,filename+'_'+self.xmlfile)
			timea=time.time()
			nullvalues=buildNull(rotatedData,protedges,nullfile)
			logging.info('null distribution built in '+str(time.time()-timea)+' secs')
			self.result=calculateStats(scores,nullvalues)
			logging.info('Stats calculated')
			logging.info('Build completed')

	###Begin a thread for each file and add it to the queue
	def producer(q,files):
		for xmlfile in files:
			thread=xmlFileBuilder(xmlfile)
			thread.start()
			q.put(thread,True)

	###As files finish in the queue, pull them into the list
	finished=[]
	def consumer(q,numFiles):
		while numFiles > len(finished):
			thread=q.get(True)
			thread.join()
			finished.append(thread.get_result())
	
	###Grab list of all graph files, assuming all .edg files in directory are supposed to be parsed
	xmlfiles=[]
	for i in os.listdir(os.path.join(pwd,edgedir)):
		if re.search('\.edg',i):
			xmlfiles.append(i)
	logging.info(xmlfiles)	

	###Begin multithreading of xml file parsing
	timec=time.time()
	try:
		ou=open('a.log','w')
		q=Queue.Queue()
		prod_thread=threading.Thread(target=producer, args=(q,xmlfiles))
		cons_thread=threading.Thread(target=consumer, args=(q,len(xmlfiles)))
		prod_thread.start()
		cons_thread.start()
		prod_thread.join()
		cons_thread.join()
	except Exception as inst:
		ou.write(type(inst))
		ou.write(inst.args)
		ou.write(inst)
		sys.exit()
	logging.info('all graphs built in '+str(time.time()-timec)+' secs.')

	###After all threads are completed, sort and write out the results
	for piece in sorted(finished,key=lambda result:result[1][2].qval):
		writeResult(piece)
	outfile.close()
	del rotatedData
	del expdict
	del finished
	gc.collect()

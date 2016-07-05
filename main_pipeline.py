#Python libraries
import sys
import os
import re
import json
import collections as cl
from multiprocessing import Pool
import copy

#Pipeline libraries
lib_dir = sys.path[0] + '/lib/'
sys.path.append(lib_dir)
import raw_processor as raw
import annIO
import annotation as ann

usage = "Requires an input file and control file"
argv = sys.argv[1:]

if len(argv) == 2:
	inputFile = argv.pop(0)
	configFile = argv.pop(0)
	print '\n'
else:
	sys.exit(usage)

#
#Regex precompiles
#

hfinder = re.compile('>(\S+)')					#Identifies headers
gfinder = re.compile('^([^\|]+)\|')				#Identifies genome tags
pfinder = re.compile('^\S+\|([^\|]+)\|')				#Identifies plasmid tag
pfinder2 = re.compile('^([^\|]+)\|.+plasmid_([^_]+)')		#Identifies genome tag + plasmid from fastas
p_nameRE = re.compile('plasmid_([^_]+)')			#Finds plasmid name

#Used for pulling info from blast hits
idREG = re.compile('\s+Identities = \d+/(\d+)\s+\((\d+)')
specREG = re.compile('\s+\[[^\]]+\]\s*$')

#
#Objects
#

#Container for information about a particular database (HMMer, BLAST, LAST)
class Database(object):
	__slots__ = ('name','path','type','multi','TRUST_SCORE','NOISE_SCORE','MIN_LEN',
				'MAX_LEN','FAM','param','source','override')

	def __init__(self,DATA,name):
		self.name = str(name)
		self.path = DATA['path']
		self.source = str(DATA['source'].lower())
		self.multi = DATA['allow_multiple']
		self.param = str(DATA['param'])
		self.type = str(DATA['type'].lower())
		self.override = DATA['override']

		self.TRUST_SCORE = {}
		if DATA['score_trust']:
			with open(DATA['score_trust'],'r') as f:
				for ln in f:
					ln = ln.replace('\n','')
					ln = ln.split('\t')
					self.TRUST_SCORE[ln[0]] = float(ln[1])

		self.NOISE_SCORE = {}
		if DATA['score_noise']:
			if DATA['score_noise'] == 'GA':
				self.param += ' --cut_ga'
			elif DATA['score_noise'] == 'NC':
				self.param += ' --cut_nc'
			elif DATA['score_noise'] == 'TC':
				self.param += ' --cut_tc'
			else:
				with open(DATA['score_noise'],'r') as f:
					for ln in f:
						ln = ln.replace('\n','')
						ln = ln.split('\t')
						self.NOISE_SCORE[ln[0]] = float(ln[1])

		self.MIN_LEN = {}
		if DATA['min_len']:
			with open(DATA['min_len'],'r') as f:
				for ln in f:
					ln = ln.replace('\n','')
					ln = ln.split('\t')
					self.MIN_LEN[ln[0]] = float(ln[1])

		self.MAX_LEN = {}
		if DATA['max_len']:
			with open(DATA['max_len'],'r') as f:
				for ln in f:
					ln = ln.replace('\n','')
					ln = ln.split('\t')
					self.MAX_LEN[ln[0]] = float(ln[1])

		self.FAM = {}
		if DATA['fam_names']:
			with open(DATA['fam_names'],'r') as f:
				for ln in f:
					ln = ln.replace('\n','')
					ln = ln.split('\t')
					self.FAM[ln[0]] = ln[1]

#Temporary container for sorting out overlapping rRNA hits
class rna(object):
	__slots__ = ('contig','start','end','id','type','strand','length','data')
	def __init__(self):
		self.contig = ''
		self.start = 0.0
		self.end = 0.0
		self.id = 0.0
		self.type = ''
		self.strand = ''
		self.length = 0
		self.data = ''



#		Processing Config File 		#
DB = cl.OrderedDict()
ctrl = open(configFile,'r')
try:
	CTRL = json.load(ctrl,object_pairs_hook=cl.OrderedDict)
except ValueError as err:
	annIO.json_error(err,ctrl,"control")
ctrl.close()

#Check for sensible parameter inputs
PARAM = CTRL['param'].copy()

#Build database objects
for i in CTRL:
	if i != 'param':
		DB[str(i)] = Database(CTRL[i],i)


#		Build directories		#
cmdLst = []
cmdLst.append("mkdir -p {0}/chromosomes".format(PARAM['outputFolder']))
cmdLst.append("mkdir -p {0}/plasmids".format(PARAM['outputFolder']))
cmdLst.append("mkdir -p {0}/json_features".format(PARAM['outputFolder']))
cmdLst.append("mkdir -p {0}/json_genome".format(PARAM['outputFolder']))

for cmd in cmdLst:
	os.system(cmd)

#		Read existing genome JSON files		#
GENOME = {}

GCHANGE = {} #Keep track of which genomes have changes to write new JSONs later
GENOME = annIO.importGenomeData('{0}/{1}'.format(PARAM['outputFolder'],'json_genome'))

#		Processing new inputs		#
inF = open(inputFile,'r')
try:
	INPUT = json.load(inF)
except ValueError as err:
	annIO.json_error(err,inF,"input")
inF.close()

newG = []
names = []
NEED_ANTIFAM = {}
ADD = {}
for g in INPUT:
	if not GENOME.has_key(g):
		newG.append(g)

if newG:
	if INPUT[g]['inputType'].lower() == 'fasta':
		ADD = raw.fasta_extract(INPUT,newG,PARAM)
	elif INPUT[g]['inputType'].lower() == 'coords':
		ADD = raw.coord_extract(INPUT,newG,PARAM)
	elif INPUT[g]['inputType'].lower() == 'genbank':
		ADD = raw.gbk_extract(INPUT,newG,PARAM)
if ADD:
	for i in ADD:
		NEED_ANTIFAM[i] = ''
		GENOME[i] = ADD[i].copy()
	del ADD


LOCUS = {}	#Locus tag list
for g in GENOME:
	for rep in GENOME[g]:
		if not LOCUS.has_key(g):
			LOCUS[g] = GENOME[g][rep].tag

#Cleaning up temporary folders and files
cmdLst.append("rm -f *.prodigal")
cmdLst.append("rm -f junk.out")

for cmd in cmdLst:
	#os.system(cmd)
	pass


#			DNA Sequence Gathering		#


#Genome() object lists for chromosomes and plasmid sequences
CHROM = {}
PLASMID = {}

genome_list = []

os.system('mkdir -p {0}'.format(PARAM['outputFolder']))
os.system('mkdir -p {0}/db_data'.format(PARAM['outputFolder']))

print 'Processing chromosomes and plasmids...'

#Compiling contig data
contigNameRE = re.compile('_([^_]+$)')
plasmidContigRE = re.compile('(plasmid_\S+)')

for g in GENOME:
	for rep in GENOME[g]:
		totalLength = 0
		with open(GENOME[g][rep].fastaPath, 'r') as f:
			for ln in f:
				if '>' in ln:
					ln = ln.translate(None, '>\n')
					contigName = ln.split('_')[-1]

					contigLength = 0
					GENOME[g][rep].contigs[contigName] = [0,0,0] #Length, Number, Start
					GENOME[g][rep].contigs[contigName][1] = len(GENOME[g][rep].contigs)
					GENOME[g][rep].contigs[contigName][2] = totalLength + 1
				else:
					length = len(ln.replace('\n','')) 
					totalLength += length


# 		Feature Data Gathering		#
print 'Processing feature data files...'

#Builds feature lists from JSON files
for g in GENOME:
	for rep in GENOME[g]:

		f = '{0}/json_features/{2}_{1}.json'.format(PARAM['outputFolder'],rep,GENOME[g][rep].tag)
		try:
			inF = open(f,'r')
			inF.close()
			genes = annIO.importFeatureData(f,'all',PARAM)
			for gene in genes:
				GENOME[g][rep].feat[gene.id] = gene
			del genes
		except IOError as (errno, strerror):
			print "I/O error({0}): {1}".format(errno, strerror)
			pass

#		Screen new genomes with antifam		#
for g in NEED_ANTIFAM:
	GENOME[g], DEBUG = raw.runAntifam(GENOME[g],PARAM)

#		Functional Annotations, & ORF Curration			#
for g in GENOME:
	for rep in GENOME[g]:
		#Set absolute start site
		for feat in GENOME[g][rep].feat:
			if not GENOME[g][rep].feat[feat].absStart > 0:
				contig = GENOME[g][rep].feat[feat].contig
				start = GENOME[g][rep].feat[feat].start
				GENOME[g][rep].feat[feat].absStart = start + GENOME[g][rep].contigs[contig][2] - 1
				GCHANGE[g] = ''
 
#Build annotation dicts for all proteins where needed
for g in GENOME:
	for rep in GENOME[g]:
		for prot in GENOME[g][rep].feat:
			if GENOME[g][rep].feat[prot].type == 'cds':
				for db in DB:
					if DB[db].type == 'cds':
						if not GENOME[g][rep].feat[prot].ann.has_key(DB[db].name):
							GENOME[g][rep].feat[prot].ann[DB[db].name] = ""
							GCHANGE[g] = ''

						if not GENOME[g][rep].feat[prot].score.has_key(DB[db].name):
							GENOME[g][rep].feat[prot].score[DB[db].name] = 0.0
							GCHANGE[g] = ''

						if not GENOME[g][rep].feat[prot].eval.has_key(DB[db].name):
							GENOME[g][rep].feat[prot].eval[DB[db].name] = 1
							GCHANGE[g] = ''

						if not GENOME[g][rep].feat[prot].putative.has_key(DB[db].name):
							GENOME[g][rep].feat[prot].putative[DB[db].name] = False
							GCHANGE[g] = ''

PROTLEN = {}
for g in GENOME:
	for rep in GENOME[g]:
		for prot in GENOME[g][rep].feat:
			if GENOME[g][rep].feat[prot].type == 'cds':
				PROTLEN[prot] = GENOME[g][rep].feat[prot].length

#		Export JSON files of annotated features		#
for g in GCHANGE:
	if GCHANGE.has_key(g):
		for rep in GENOME[g]:
			outFile = '{0}/json_features/{2}_{1}.json'.format(PARAM['outputFolder'],rep,GENOME[g][rep].tag)
			debug = annIO.exportJsonData(GENOME[g][rep].feat,outFile,'id')
			if not debug:
				print "Failed to export JSON data"

#
#Running proteins against DBs
#
for db in DB:
	NEED = cl.OrderedDict()
	if DB[db].type == 'cds':
	
		#HMMer seems to be unstable when every core on a computer is running its own HMMer instance
		if DB[db].source == 'hmmer':
			fNum = PARAM['cores'] / 2
		else:
			fNum = PARAM['cores']

		#Identify protein coding genes that lack annotations from this database to extract sequences
		for g in GENOME:
			for rep in GENOME[g]:
				for prot in GENOME[g][rep].feat:
					if GENOME[g][rep].feat[prot].type == 'cds':
						if  DB[db].override == True or not GENOME[g][rep].dbsDone.has_key(DB[db].name):
							GCHANGE[g] = ''
							if not NEED.has_key(GENOME[g][rep].fastaPath):
								NEED[GENOME[g][rep].fastaPath] = cl.deque()
							seq = annIO.SeqFetch()
							seq.locus = GENOME[g][rep].feat[prot].id
							seq.start = GENOME[g][rep].feat[prot].start - 1
							seq.end = GENOME[g][rep].feat[prot].end
							seq.strand = GENOME[g][rep].feat[prot].strand
							seq.contig = GENOME[g][rep].feat[prot].contig

							NEED[GENOME[g][rep].fastaPath].append(seq)

		#Extracting sequences, minimize memory usage by not storing bulk sequence data in ram
		if NEED:
			DEBUG = annIO.extractSeq(NEED,PARAM['outputFolder'],'prot',fNum)

			fastaNames = os.listdir('{0}/tempFastas'.format(PARAM['outputFolder']))
			fastas = []
			for f in fastaNames:
				fastas.append('{0}/tempFastas/{1}'.format(PARAM['outputFolder'],f))

			#Run algorithm and extract annotation data
			PROT_ANN, debug = ann.db_run(DB[db],PARAM,fastas,LEN=PROTLEN)

			#Set annotation in master data structure
			if debug == 0:
				dbName = DB[db].name
				for sub in PROT_ANN:
					for g in GENOME:
						for rep in GENOME[g]:
							for prot in GENOME[g][rep].feat:
								if PROT_ANN[sub].has_key(prot):
									GENOME[g][rep].feat[prot].ann[dbName] = PROT_ANN[sub][prot].ann[dbName]
									GENOME[g][rep].feat[prot].eval[dbName] = PROT_ANN[sub][prot].eval[dbName]
									GENOME[g][rep].feat[prot].score[dbName] = PROT_ANN[sub][prot].score[dbName]
									GENOME[g][rep].feat[prot].putative[dbName] = PROT_ANN[sub][prot].putative[dbName]
									GCHANGE[g] = ""

									if PROT_ANN[sub][prot].orphan == 0:
										GENOME[g][rep].feat[prot].orphan = 0
				del PROT_ANN
			else:
				print "Error running and parsing {0} database, exiting".format(DB[db].name)
				sys.exit()

			#Mark genome as run against this db, and mark proteins with no hits
			for g in GENOME:
				for rep in GENOME[g]:
					if not GENOME[g][rep].dbsDone.has_key(DB[db].name):
						GENOME[g][rep].dbsDone[DB[db].name] = ''
						GCHANGE[g] = ""
					for prot in GENOME[g][rep].feat:
						if not GENOME[g][rep].feat[prot].ann[dbName]:
							GENOME[g][rep].feat[prot].ann[dbName] = "-"
							GCHANGE[g] = ""

	elif DB[db].type == "rna":
		fastas = []
		for g in GENOME:
			for rep in GENOME[g]:
				if not GENOME[g][rep].dbsDone.has_key(DB[db].name):
					GCHANGE[g] = ''
					fastas.append(GENOME[g][rep].fastaPath)
		if fastas:
			NEW_FEAT, debug = ann.db_run(DB[db], PARAM, fastas, featCount =GENOME[g][rep].rnaCount)
			for g in GENOME:
				for rep in NEW_FEAT:
					if GENOME[g].has_key(rep):
						GENOME[g][rep].dbsDone[DB[db].name] = ""
						for n in NEW_FEAT[rep]:
							contig = '.'.join(n.split('.')[:-1])
							GENOME[g][rep].rnaCount += 1
							newFeat = copy.copy(NEW_FEAT[g][n])
							locusFinal = '%07d' % (GENOME[g][rep].rnaCount)
							locusFinal = '{0}|{0}_r{1}'.format(GENOME[g][rep].tag,locusFinal)
							newFeat.id = locusFinal
							GENOME[g][rep].feat[locusFinal] = newFeat
							GCHANGE[g] = ""

	#ADD THIS IN
	elif DB[db].type == "misc":
		pass

#		Data Finalization		#
print 'Finalizing feature assignments and annotations...'
#		Export JSON files of annotated features		#
for g in GENOME:
	if GCHANGE.has_key(g):
		for rep in GENOME[g]:
			outFile = '{0}/json_features/{2}_{1}.json'.format(PARAM['outputFolder'],rep,GENOME[g][rep].tag)
			debug = annIO.exportJsonData(GENOME[g][rep].feat,outFile,'id')
			if not debug == 1:
				print "Failed to export JSON data"

#
#Set final protein annotations, prioritized based on config file order
#
if PARAM['newGeneNames'] == True:
	for g in GENOME:
		for rep in GENOME[g]:
			for prot in GENOME[g][rep].feat:
				for db in DB:
					if DB[db].type == 'cds' and GENOME[g][rep].feat[prot].type == 'cds':
						if GENOME[g][rep].feat[prot].ann[DB[db].name] != '-' and GENOME[g][rep].feat[prot].ann[DB[db].name] != '' and not GENOME[g][rep].feat[prot].function:

							#Check if it has a nice human-readable gene name
							if DB[db].multi == False:
								if DB[db].FAM.has_key(GENOME[g][rep].feat[prot].ann[DB[db].name]):
										name = DB[db].FAM[GENOME[g][rep].feat[prot].ann[DB[db].name]]
										if GENOME[g][rep].feat[prot].putative[DB[db].name] == True:
											GENOME[g][rep].feat[prot].function = 'putative {0}'.format(name)
										else:
											GENOME[g][rep].feat[prot].function = name

							#If it only has a database ID
							if not GENOME[g][rep].feat[prot].function:
								if DB[db].multi == False:
									name = GENOME[g][rep].feat[prot].ann[DB[db].name]
									if GENOME[g][rep].feat[prot].putative[DB[db].name]:
										GENOME[g][rep].feat[prot].function = 'putative {0} family protein'.format(name)
									else:
										GENOME[g][rep].feat[prot].function = '{0} family protein'.format(name)
								else:
									fams = GENOME[g][rep].feat[prot].ann[DB[db].name]
									fams = filter(None, fams)

									if len(fams) > 1:
										name = '-'.join(fams)
									else:
										name = ''.join(fams)
									GENOME[g][rep].feat[prot].function = '{1} {0} domain protein'.format(name,DB[db].name)

				if GENOME[g][rep].feat[prot].function == "":
					GENOME[g][rep].feat[prot].function = 'hypothetical protein'

				#Tag partial proteins
				if GENOME[g][rep].feat[prot].type == 'cds':
					if GENOME[g][rep].feat[prot].partial != '00':
						GENOME[g][rep].feat[prot].function += ' partial protein'
					if 'protein protein' in GENOME[g][rep].feat[prot].function:
						GENOME[g][rep].feat[prot].function = GENOME[g][rep].feat[prot].function.replace('protein protein','protein')
					if 'putative putative' in GENOME[g][rep].feat[prot].function:
						GENOME[g][rep].feat[prot].function = GENOME[g][rep].feat[prot].function.replace('putative putative','putative')
					if 'protein partial protein' in GENOME[g][rep].feat[prot].function:
						GENOME[g][rep].feat[prot].function = GENOME[g][rep].feat[prot].function.replace('protein partial protein','partial protein')

#		Print locus tag table		#
f = open('{0}/genome_locus.txt'.format(PARAM['outputFolder']),'w')
for i in sorted(LOCUS):
	f.write('{0}\t{1}\n'.format(LOCUS[i],i))
f.close()

#		Export JSON files of annotated features		#
for g in GENOME:
	if GCHANGE.has_key(g):

		outFile = '{0}/json_genome/{2}_{1}.json'.format(PARAM['outputFolder'],g,LOCUS[g])
		annIO.exportJsonData(GENOME[g],outFile,"name")

		for rep in GENOME[g]:
			outFile = '{0}/json_features/{2}_{1}.json'.format(PARAM['outputFolder'],rep,GENOME[g][rep].tag)
			debug = annIO.exportJsonData(GENOME[g][rep].feat,outFile,'id')
			if not debug == 1:
				print "Failed to export JSON data"
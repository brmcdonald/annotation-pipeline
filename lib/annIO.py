import os
import collections as cl
import re
import json
import sys
import os

lib_dir = sys.path[0] + '/lib/'
sys.path.append(lib_dir)
import raw_processor as raw

#Used for chromosomes and plasmids
class Replicon(object):
	__slots__ = ('name','rnaCount','proteinCount','partials','contigs',
				'contigLengths','feat','tag','length','parent',
				'plasmidName','fastaPath','genePath','protPath','rnaPath',
				'miscCount','rawFastaPath',"type","tempPath",
				"dbsDone","orgType")

	def __init__(self):
		self.name = ''
		self.rnaCount = 0
		self.proteinCount = 0
		self.miscCount = 0
		self.partials = 0
		self.contigs = cl.OrderedDict()						#{contig} = [length,number,startPos]
		self.contigLengths = cl.OrderedDict()		#List of contig lengths
		self.feat = cl.OrderedDict()								#Annotated list of features is collected here
		self.dbsDone = {}				#Set of databases this genome has been run against
		self.tag = ''								#Locus tag used
		self.length = 0								#Total length
		self.parent = ''							#Parent genome, used for plasmids
		self.plasmidName = ''						#Plasmid name
		self.rawFastaPath = ''							#Path to original source DNA fasta file
		self.fastaPath = ''							#Path to DNA fasta file
		self.genePath = ''							#Path to DNA protein coding sequences
		self.protPath = ''							#Path to amino acid protein sequences
		self.rnaPath = ''							#path to RNA coding sequences
		self.type = ""				#Is it a chromsome or plasmid?
		self.tempPath = ""			#attribute for storing paths to temporary files
		self.orgType = ""		#Bacteria, Archea, Eukaryote, Virus, Mixed

#Class for protein coding genes
class Protein(object):
	__slots__ = ('contig','strand','tag','locusNum','function','plasmidName','plasmidTag',
				'start','end','length','external','partial','absStart','ann','score','aaSeq',
				'nucSeq','geneTag','putative','orphan','type','coords','id',
				'eval','failedAntifam','genome','replicon')

	def __init__(self):
		self.type = 'cds'
		self.contig = ""			#Name of contig it's found on
		self.genome = ''		#name of genome it's from
		self.replicon = ""		#name of replicon it's found on
		self.strand = ""		#Forward or reverse (1,-1)
		self.tag = ""				#genome tag
		self.locusNum = ""			#Locus number
		self.function = ""			#Primary gene function annotation, used to name gene
		self.plasmidName = ""		#Plasmid name, if there is one
		self.plasmidTag = ""
		self.id = ""		#Final locus id to use for gene

		self.start = 0			#Start on parent contig (note always 5 to 3 prime on genome seq)
		self.end = 0				#End on parent contig
		self.length = 0	#Gene length
		self.external = 0			#Whether the feature comes from an external source ie GBK
		self.partial = ""
		self.failedAntifam = False			#Whether the feature has failed antifam curration

		#Start relative to the start of the parent fasta file, in case there are multiple contigs
		self.absStart = 0

		self.ann = {}			#Annotations for each database (string)
		self.score = {}			#Bitscores for the top hit in each db (float)
		self.coords = {}		#For DBs with multiple hits per gene, coords of each hit
		self.eval = {}			#Evalue for hits
		self.aaSeq = ""			#Amino acid sequence
		self.nucSeq = ""		#Nucleotide sequence
		self.geneTag = ""		#Short name for a gene, ie dnaA or metB
		self.putative = {}		#How sure are we of each db hit? True/False
		self.orphan = 1			#Used to track whether a gene has a decent annotation

#Class for RNAs or noncoding regions
class Feature(object):
	__slots__ = ('contig','strand','name','tag','locusNum','function',
				'plasmidName','plasmidTag','start','end','length','external',
				'partial','abs_start','type','anti','amino','geneTag','id','nucSeq',
				'type','score','eval','full','putative')

	def __init__(self,feat_type):
		self.type = feat_type
		self.contig = ""			#Name of contig
		self.strand = 0				#Forward or reverse (1,-1)
		self.name = ""				#Name of gene/tRNA in the source fasta or output file
		self.tag = ""				#Locus tag
		self.locusNum = ""			#Locus number
		self.ann = ""			#Final gene annotation
		self.plasmidName = ""		#Plasmid name, if there is one
		self.plasmidTag = ""
		self.id = ""		#final locus id

		self.start = 0				#Start on parent contig (note always 5 to 3 prime on genome seq)
		self.end = 0				#End on parent contig
		self.length = 0				#Gene length
		self.external = 0			#Whether the feature comes from an external source ie GBK
		self.partial = '00'
		self.geneTag = ""			#Short name for a gene, ie dnaA or metB
		self.nucSeq = ""
		self.score = 0.0
		self.eval = ""
		self.full = ""				#Full output line from raw output file
		self.putative = False

		#Start relative to the start of the parent fasta file, in case there are multiple contigs
		self.abs_start = 0

		#For tRNAs only
		self.anti = ""		#Anticodon
		self.amino = ""		#Amino acid for tRNA

#Temporary container for extracting sequence data from DNA fasta files
class SeqFetch(object):
	__slots__ = ('locus','start','end','strand','contig')
	def __init__(self):
		self.locus = ''
		self.start = 0
		self.end = 0
		self.strand = 0
		self.contig = ''

def importFeatureData(f,featType,PARAM):
	"""extract feature data from json data files"""
	data = open(f,'r')
	try:
		DATA = json.load(data,object_pairs_hook=cl.OrderedDict)
	except ValueError as err:
		json_error(err,data,"feature data")

	data.close()
	genes = []

	for g in DATA:
		genomeTag = DATA[g]['tag']
		if DATA[g]['type'] == 'cds' and (featType == 'all' or featType == 'cds'):
			gene = Protein()

			for attr in DATA[g]:
				if attr == 'function' and PARAM['newGeneNames'] == False:
					if not DATA.has_key('function'):
						gene.function = 'Hypothetical protein'
					elif DATA[g]['function']:
						gene.final_ann = DATA[g]['function']
						gene.geneTag = genomeTag
						gene.external = 1
						gene.orphan = 0
					else:
						gene.final_ann = 'Hypothetical protein'
				elif type(DATA[g][attr]) == unicode:
					gene.__setattr__(attr,str(DATA[g][attr]))
				elif type(DATA[g][attr]) == cl.OrderedDict:
					gene.__setattr__(attr,_decode_dict(DATA[g][attr]))
				elif type(DATA[g][attr]) == list:
					gene.__setattr__(attr,_decode_list(DATA[g][attr]))
				else:
					gene.__setattr__(attr,DATA[g][attr])

		elif featType == "RNA":
			pass
		genes.append(gene)
	return genes

def importGenomeData(folder):
	"""extract genome data from JSON files"""
	jsons = os.listdir(folder)
	GENOME = {}
	ALL = {}
	for f in jsons:
		data = open('{0}/{1}'.format(folder,f),'r')
		try:
			DATA = json.load(data,object_pairs_hook=cl.OrderedDict)
		except ValueError as err:
			json_error(err,data,"genome data")

		data.close()

		for rep in DATA:
			replicon = Replicon()
			for attr in DATA[rep]:
				if type(DATA[rep][attr]) == unicode:
					replicon.__setattr__(attr,str(DATA[rep][attr]))
				elif type(DATA[rep][attr]) == cl.OrderedDict:
					replicon.__setattr__(attr,_decode_dict(DATA[rep][attr]))
				elif type(DATA[rep][attr]) == list:
					replicon.__setattr__(attr,_decode_list(DATA[rep][attr]))
				else:
					replicon.__setattr__(attr,DATA[rep][attr])

			if not GENOME.has_key(replicon.parent):
				GENOME[replicon.parent] = cl.OrderedDict()
			GENOME[replicon.parent][replicon.name] = replicon
	return GENOME

def exportJsonData(FEAT,outFile,name,full=False,source='Std'):
	"""Export data from pipeline objects to JSON format."""
	out = open(outFile,'w')
	NEVER = {'contigs':'','contigLengths':''}	#Avoid writing contig sequences in JSON, too much text
	FULL = {"feat":"","miscFeat":"","rnaFeat":""}	#Only write feature data for feature JSON files
	WRITE = cl.OrderedDict()

	inputType = type(FEAT)

	if inputType == cl.deque:
		for f in FEAT:
			WRITE[getattr(f,name)] = cl.OrderedDict()
			for att in dir(f):
				if full == False:
					if not att.startswith("__") and not NEVER.has_key(att) and not FULL.has_key(att):
						WRITE[getattr(f,name)][att] = getattr(f,att)
				else:
					if not att.startswith("__") and not NEVER.has_key(att):
						WRITE[getattr(f,name)][att] = getattr(f,att)

	elif	inputType == cl.OrderedDict and source == 'Std':
		for f in FEAT:
			WRITE[getattr(FEAT[f],name)] = cl.OrderedDict()
			for att in dir(FEAT[f]):
				if full == False:
					if not att.startswith("__") and not NEVER.has_key(att) and not FULL.has_key(att):
						WRITE[getattr(FEAT[f],name)][att] = getattr(FEAT[f],att)
				else:
					if not att.startswith("__") and not NEVER.has_key(att):
						WRITE[getattr(FEAT[f],name)][att] = getattr(FEAT[f],att)

	elif	inputType == cl.OrderedDict and source == 'genome':
		for f in FEAT:
			WRITE[getattr(FEAT[f],name)] = cl.OrderedDict()
			for att in dir(FEAT[f]):
				if full == False:
					if not att.startswith("__") and not NEVER.has_key(att) and not FULL.has_key(att):
						WRITE[getattr(FEAT[f],name)][att] = getattr(FEAT[f],att)
				else:
					if not att.startswith("__") and not NEVER.has_key(att):
						WRITE[getattr(FEAT[f],name)][att] = getattr(FEAT[f],att)
	else:
		f = FEAT
		WRITE[getattr(f,name)] = cl.OrderedDict()
		for att in dir(f):
			if full == False:
				if not att.startswith("__") and not NEVER.has_key(att) and not FULL.has_key(att):
					WRITE[getattr(f,name)][att] = getattr(f,att)
			else:
				if not att.startswith("__") and not NEVER.has_key(att):
					WRITE[getattr(f,name)][att] = getattr(f,att)
	json.dump(WRITE,out,indent=4)
	out.close()
	return 1

#Decoding JSON unicode to strings
def _decode_list(data):
	rv = []
	for item in data:
		if isinstance(item, unicode):
			item = item.encode('utf-8')
		elif isinstance(item, list):
			item = _decode_list(item)
		elif isinstance(item, dict):
			item = _decode_dict(item)
		rv.append(item)
	return rv

def _decode_dict(data):
	rv = cl.OrderedDict()
	for key, value in data.iteritems():
		if isinstance(key, unicode):
			key = key.encode('utf-8')
		if isinstance(value, unicode):
			value = value.encode('utf-8')
		elif isinstance(value, list):
			value = _decode_list(value)
		elif isinstance(value,cl.OrderedDict):
			value = _decode_dict(value)
		rv[key] = value
	return rv


def extractSeq(GENE,outFolder,seqType,fNum):
	"""Extract DNA or Protein sequences based on coordinate"""
	os.system('rm -r {0}/tempFastas > /dev/null'.format(outFolder))
	os.system('mkdir -p {0}/tempFastas'.format(outFolder))
	SEQ = cl.OrderedDict()
	
	for f  in GENE:
		CONTIG = {}
		with open(f,'r') as data:
			for ln in data:
				if '>' in ln:
					ln = ln.translate(None,'\n>')
					ln = ln.split('_')[-1]
					contig = ln
					CONTIG[contig] = cl.deque()
				else:
					CONTIG[contig].append(ln.replace('\n',''))
		for c in CONTIG:
			CONTIG[c] = ''.join(CONTIG[c])
		for seq in GENE[f]:
			if CONTIG.has_key(seq.contig):
				SEQ[seq.locus] = CONTIG[seq.contig][seq.start:seq.end]
			else:
				print "Cant find contig {0} in replicon {1}".format(seq.contig,f)
			if seq.strand == '-1':
				SEQ[seq.locus] = raw.revcomp(SEQ[seq.locus])
				pass
			if seqType == 'prot':
				SEQ[seq.locus] = raw.translate(SEQ[seq.locus])
				pass


	seqEachFile = len(SEQ) / fNum + 1

	seqCount = 0
	out = open('{0}/tempFastas/seq{1}.fasta'.format(outFolder,fNum),'w')
	for s in SEQ:
		out.write('>{0}\n{1}\n'.format(s,SEQ[s]))
		seqCount += 1
		if seqCount == seqEachFile:
			out.close()
			fNum = fNum - 1
			out = open('{0}/tempFastas/seq{1}.fasta'.format(outFolder,fNum),'w')
			seqCount = 0
	if out:
		out.close()
	return CONTIG

def json_error(err,f,fileType):
	"""Processing JSON import errors, suggesting common JSON syntax errors"""
	err = str(err)
	if "Expecting property name" in err:
		m = re.search('line (\d+)',err)
		if m:
			badLn = m.group(1)
		else:
			badLn = '?'
		print "Error in {2} JSON file: {0} on line {1}, unnecessary comma?".format(f,badLn,fileType)
	elif "Expecting , delimiter" in err:
		m = re.search('line (\d+)',err)
		if m:
			badLn = m.group(1)
		else:
			badLn = '?'
		print "Error in {2} JSON file: {0} on line {1}, missing a comma?".format(f,badLn,fileType)
	elif "No JSON object could be decoded" in err:
		print "Error in {1} JSON file: {0}. missing value for a variable?".format(f,fileType)
	else:
		print "Error in parsing {2} JSON file {0}: {1}".format(f,badLn,fileType)
	sys.exit()
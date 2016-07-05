import re
import os
import sys
import json
import collections as cl
from multiprocessing import Pool

antifamPath = sys.path[0] + '/antifam/'
lib_dir = sys.path[0] + '/lib/'
sys.path.append(lib_dir)
import annIO
import annotation as ann


#Dictionaries and lists of things to clean from headers
#Should be run in order: Remove, Regex remove, Replace
hRemove = [
	"_genomic_scaffold_super", "_chr_", "subsp._", "sp._", "sv._"
	"str._", "serovar_", "sv_", "(", ")", "\'", ".scaffold", "chromosome", "sequence", "sp_",
	"complete_genome", "strain_", "str_", "sp.", "str._",'_DNA'
]
hRegex = [
	"_gcontig_\d+","\,[^\,]+$","\s*NZ_\S+","\s*sq\S+","\s*contig\S+","\s*plasmid.+",
	"\s+\d\d\d\d\d\d+$","\S+\.Contig\d+","_main","_super","SM8_Contig[^\,]+","\s*cont\d+\.\S+",
	"\s*cont\d+","_\:[^\:]+$","NODE_\S+","_scaffold_\d+_Cont\d+","_genomic_scaffold_[^\,]*",
	"Scaffold[^\,]*","_Contig[^\,]+","\.Reformatted\S+","\.454AllContigs\S*","_supercont\S* ",
	"\.fna\s*\d*"
]
hReplace = {
	"Tu_":"Tu", "NBRC_":"NBRC", "NRRL_":"NRRL", "NRRLB_":"NRRLB", "ATCC_":"ATCC",
	"IFM_":"IFM", "DSM_":"DSM", "Tu_":"Tu", "sp_":"", "NCTC_":"NCTC","NCPPB_":"NCPPB",
	"ATCC_":"ATCC", "IFM_":"IFM", "DSM_":"DSM", "/":"-", "_super_genomic_scaffold":"",
	"JCM_":"JCM", ".":"-", "PAMC_":"PAMC", "str_":""
}

#DNA->protein translation table
CODONS = {
	"ACC":"T","ATG":"M","ACA":"T","ACG":"T","ATC":"I","AAC":"N","ATA":"I","AGG":"R",
	"CCT":"P","CTC":"L","AGC":"S","AAG":"K","AGA":"R","CAT":"H","AAT":"N","ATT":"I",
	"CTG":"L","CTA":"L","ACT":"T","CAC":"H","AAA":"K","CAA":"Q","AGT":"S","CCA":"P",
	"CCG":"P","CCC":"P","CTT":"L","TAT":"Y","GGT":"G","TGT":"C","CGA":"R","CAG":"Q",
	"CGC":"R","GAT":"D","CGG":"R","TTT":"F","TGC":"C","GGG":"G","TGA":"*","GGA":"G",
	"TGG":"W","GGC":"G","TAC":"Y","TTC":"F","TCG":"S","TAG":"*","TTG":"L","CGT":"R",
	"GAA":"E","TCA":"S","GCA":"A","GTA":"V","GCC":"A","GTC":"V","GCG":"A","GTG":"V",
	"GAG":"E","GTT":"V","GCT":"A","TTA":"L","GAC":"D","TCC":"S","TAA":"*","TCT":"S",
	'NNN':'x'}

#multithreadable command runner
def cmd_runner(cmd):
	os.system(cmd)

def translate(seq):
	"""Translates nucleotide to amino acid sequence"""
	global CODONS

	aaSeq = cl.deque()
	codons = [seq[i:i+3] for i in range(0, len(seq), 3)]

	for i in codons:
		i = i.upper()
		if len(i) != 3:
			return ''.join(aaSeq)
		elif not CODONS.has_key(i):
			aaSeq.append('X')
		else:
			aaSeq.append(CODONS[i])
	return ''.join(aaSeq)

def revcomp (seq):
	"Generates reverse complement"

	if not isinstance(seq,str):
		raise TypeError("revcomp() requires a string input")
	f_lst = ("G", "C", "A", "T", "g", "c", "a", "t")
	r_lst = ("C", "G", "T", "A", "c", "g", "t", "a")
	r_seq = []
	seq = list(seq.strip())
	for i in reversed(seq):
		for n in range(len(f_lst)):
			if f_lst[n] is i:
				r_seq.append(r_lst[n])
				break
	return ''.join(r_seq)


def coord_extract(GENOME,newGenomes,PARAM):
	"""Extracts proteins from fasta files based on provided gene coordinates. Acceptable coordinate file formats:
	gm: metagenemark orfs output
	gff: gff format, ie IMG metagnomes"""

	def metagenemark_parser(f):
		"""Parses metagenemark format coordinate file"""
		GENE = cl.OrderedDict()
		record = 0
		with open(f,'r') as data:
			for ln in data:
				ln = ln.replace('\n','')

				if 'FASTA definition line:' in ln:
					contig = ln.split('line: ')[-1]
				elif 'Predicted proteins:' in ln or 'Model information' in ln:
					record = 0
				elif 'LeftEnd    RightEnd' in ln:
					record = 1

				elif record == 1 and '#' not in ln and ln:
					ln = ln.split()
					try:
						geneNum = ln[0]
						strand = ln[1]
						leftEnd = ln[2]
						rightEnd = ln[3]
					except IndexError:
						print "coud not read line, exiting:\n{0}\n".format(ln)
						sys.exit()

					name = contig + '_' + geneNum
					GENE[name] = annIO.Protein()
					GENE[name].contig = contig
					GENE[name].id = name
					GENE[name].start = int(leftEnd.replace('<',''))
					GENE[name].end = int(rightEnd.replace('>',''))
					GENE[name].length = GENE[name].end - GENE[name].start
					GENE[name].strand = strand
					GENE[name].type = 'cds'

					if '<' in leftEnd:
						GENE[name].partial = '10'
					if '>' in rightEnd:
						if GENE[name].partial == '10':
							GENE[name].partial = '11'
						else:
							GENE[name].partial = '01'
		return GENE

	def gff_parser(f):
		"""Parses gff format coordinate file"""
		GENE = cl.OrderedDict()

		with open(f,'r') as data:
			for ln in data:
				ln = ln.replace('\n','')
				ln = ln.split('\t')
				featType = ln[2]
				
				names = ln[-1].split(';')
				for n in names:
					if 'ID=' in n:
						name = n.replace('ID=','')
					if 'locus_tag=' in n:
						name = n.replace('locus_tag=','')

				if featType == 'CDS':
					GENE[name] = annIO.Protein()
					GENE[name].type = 'cds'
				else:
					GENE[name] = annIO.Feature()
					GENE[name].type = featType


				GENE[name].contig = ln[0]
				GENE[name].id = name
				GENE[name].start = int(ln[3])
				GENE[name].end = int(ln[4])
				GENE[name].length = GENE[name].end - GENE[name].start

				if '-' in ln[6]:
					GENE[name].strand = '-1'
				else:
					GENE[name].strand = '1'
		return GENE

	#Generating output folders
	cmdLst = []
	cmdLst.append("mkdir -p {0}/chromosomes".format(PARAM['outputFolder']))
	for cmd in cmdLst:
		os.system(cmd)

	REPLICON = cl.OrderedDict()
	CONTIG_LOOKUP = {}

	for gName in newGenomes:
		for idx, replicon in enumerate(GENOME["chromosomeFastas"]):
			gObj = annIO.Replicon()
			gObj.name = gName
			numFastas =  len(GENOME["chromosomeFastas"])
			numNames =  len(GENOME['chromosomeNames'])
			if numFastas == numNames and len(str(GENOME["chromosomeNames"][idx])) > 0:
				gObj.name = gObj.name + "_" + GENOME["chromosomeNames"][idx]
			gObj.tag = GENOME['locusTag']
			gObj.type = 'metagenome'
			gObj.parent = gName
			gObj.fastaPath = os.path.abspath('{0}/chromosomes/{1}.fna'.format(PARAM['outputFolder'],gObj.name))

			outChrom = open(gObj.fastaPath,'w')
			contigCount = 0
			totalLength = 0
			with open('{0}/{1}'.format(GENOME["inputFolder"],replicon),'r') as data:
				for ln in data:
					if '>' in ln:
						ln = ln.replace('\n','')
						ln = ln.replace('\r','')
						contigCount += 1
						contigNameFinal = ''

						if ">gi|" in ln:
							data = ln.split("|")
							contigName = data[3]
							contigName = re.sub("\.\S+$","",contig_name)
							contigNameFinal = contig_name.replace("_",".")
						if not contigNameFinal:
							contigNameFinal = 'contig.' + str(contigCount)
						contigName = ln.replace('>','')

						gObj.contigs[contigNameFinal] = [0,contigCount,totalLength + 1]
						CONTIG_LOOKUP[contigName] = contigNameFinal
						outChrom.write('>{2}|{0}_{1}\n'.format(gName,contigNameFinal,gObj.tag))
					else:
						outChrom.write(ln)
						gObj.contigs[contigNameFinal][0] += len(ln) - 1
						totalLength += len(ln) - 1

			REPLICON[gObj.name] = gObj
			outChrom.close()

	#		Process Feature Data		#
	for rep in REPLICON:
		if GENOME['coordType'] == 'gm':
			for idx, thing in enumerate(GENOME["chromosomeFastas"]):
				FEAT = metagenemark_parser("{0}/{1}".format(GENOME["inputFolder"],GENOME["coordFiles"][idx]))
				for coord_feature in FEAT:
					REPLICON[rep].feat[coord_feature] = FEAT[coord_feature]
		elif GENOME['coordType'] == 'gff':
			for idx, thing in enumerate(GENOME["chromosomeFastas"]):
				FEAT = gff_parser("{0}/{1}".format(GENOME["inputFolder"],GENOME["coordFiles"][idx]))
				for coord_feature in FEAT:
					REPLICON[rep].feat[coord_feature] = FEAT[coord_feature]

	protCount = 0
	rnaCount = 0
	TAG = {}
	for rep in REPLICON:
		for feat in REPLICON[rep].feat:
			REPLICON[rep].feat[feat].failedAntifam=False
			REPLICON[rep].feat[feat].genome = gName
			REPLICON[rep].feat[feat].replicon = rep
			REPLICON[rep].feat[feat].contig = CONTIG_LOOKUP[REPLICON[rep].feat[feat].contig]

			if PARAM['newLocusTags'] == False:
				REPLICON[rep].feat[feat].locusNum = rawLocus.split('_')[-1]
			else:
				if REPLICON[rep].feat[feat].type == 'cds':
					protCount += 1
					locus_id_final = '%07d' % (protCount,)
					REPLICON[rep].feat[feat].locusNum = locus_id_final
					if REPLICON[rep].type == 'plasmid':
						data = (REPLICON[g].tag,REPLICON[rep].plasmidName,locus_id_final)
						REPLICON[rep].feat[feat].id = '{0}|{1}_{2}'.format(*data)
					else:
						REPLICON[rep].feat[feat].id = '{0}|{0}_{1}'.format(REPLICON[rep].tag,locus_id_final)
				elif 'RNA' in REPLICON[rep].feat[feat].type:
						rnaCount += 1
						locus_id_final = '%07d' % (rnaCount,)
						REPLICON[rep].feat[feat].locusNum = locus_id_final
						if REPLICON[rep].type == 'plasmid':
							data = (REPLICON[g].tag,REPLICON[rep].plasmidName,locus_id_final)
							REPLICON[rep].feat[feat].id = '{0}|{1}_{2}'.format(*data)
						else:
							REPLICON[rep].feat[feat].id = '{0}|{0}_{1}'.format(REPLICON[rep].tag,locus_id_final)
				
	#		Write replicon and feature JSON files		#
	GENOME = {}
	for i in REPLICON:
		REPLICON[i].tempPath = ""
		parent = REPLICON[i].parent
		if not GENOME.has_key(parent):
			GENOME[parent] = cl.OrderedDict()
			TAG[parent] = REPLICON[i].tag
		GENOME[parent][REPLICON[i].name] = REPLICON[i]

		outDF = "{0}/json_features/{2}_{1}.json".format(PARAM['outputFolder'],REPLICON[i].name,TAG[parent])
		annIO.exportJsonData(GENOME[parent][REPLICON[i].name].feat,outDF,'id')

	for i in GENOME:
		outFile = '{0}/json_genome/{2}_{1}.json'.format(PARAM['outputFolder'],i,TAG[i])
		annIO.exportJsonData(GENOME[i],outFile,"parent")
	return GENOME


def fasta_extract(GENOME,newGenomes,PARAM):
	"""Uses prodigal to identify protein coding regions from fasta files"""
	global hRemove
	global hRegex
	global hReplace

	#Cleaners for prodigal headers
	contigReg = re.compile(">(\S+)_\d+_*\d*")
	aReg = re.compile("contig\S+")
	bReg = re.compile("scaffold\S+")
	cReg = re.compile("^>([^\|]+)\|(\S+)")			#Find genome tag
	dReg = re.compile("plasmid_([^_]+)")			#Find plasmid name
	eReg = re.compile("plasmid_\S+")				#Remove excess plasmid name stuff
	debugData = cl.deque()

	def prodigal_head_cleaner(i):
		"""Returns Genome, Length, partial (0/1), gene, start, end, forward/reverse, genome tag, gene num, 
			plasmdid name,contig"""
		DATA = {}

		if ">" in i:
			data = i.split(" # ")
			data2 = data[4].split(";")
			DATA['rawName'] = data[0].replace('>','')
			DATA['contig'] = data[0].split('_')[-2]
			DATA['start'] = int(data[1])
			DATA['end'] = int(data[2])
			DATA['length'] = (int(data[2]) - int(data[1])) / 3
			DATA['strand'] = data[3]

			if "partial=00" in i:
				DATA['partial'] = "00"
			elif "partial=10" in i:
				DATA['partial'] = "10"
			elif "partial=01" in i:
				DATA['partial'] = "01"
			elif "partial=11" in i:
				DATA['partial'] = "11"

			if "ID=" in data2[0]:
				DATA['id'] =  data2[0].replace("ID=","")

			return DATA

	plasmidFlag = 0		#Whether we need to deal with plasmid sequences or not
	errFlag = 0
	plasmid_count = 0

	#Error file if needed
	ERR = open("genome_processing.err", 'w')
	REPLICON = cl.OrderedDict()
	SEQ_LEN = {}

	print "Reading raw sequence files..."
	for newG  in newGenomes:
		for idx, replicon in enumerate(GENOME[newG]["chromosomeFastas"]):
			gObj = annIO.Replicon()
			gObj.name = newG
			numFastas =  len(GENOME[newG]["chromosomeFastas"])
			numNames =  len(GENOME[newG]['chromosomeNames'])
			if numFastas == numNames and len(str(GENOME[newG]["chromosomeNames"][idx])) > 0:
				gObj.name = gObj.name + "_" + GENOME[newG]["chromosomeNames"][idx]
			gObj.tag = GENOME[newG]['locusTag']
			gObj.type = 'chromosome'
			gObj.parent = newG
			gObj.fastaPath = os.path.abspath('{0}/chromosomes/{1}.fna'.format(PARAM['outputFolder'],gObj.name))
			gObj.tempPath = os.path.abspath('{0}/prodigal_temp/{1}.faa'.format(PARAM['outputFolder'],gObj.name))

			with open('{0}/{1}'.format(GENOME[newG]["inputFolder"],replicon),'r') as data:
				contigName = ''
				for ln in data:
					ln = ln.replace('\n','')
					ln = ln.replace('\r','')
					if '>' in ln:
						contigName = ln.replace('>','')
						gObj.contigs[contigName] = cl.deque()
					elif contigName:
						gObj.contigs[contigName].append(ln)
			REPLICON[gObj.name] = gObj

			outChrom = open(gObj.fastaPath,'w')
			contigCount = 0

			for contig in REPLICON[gObj.name].contigs:
				contigName = contig.replace(' ','_')
				contigNameFinal = ''
				if ">gi|" in contigName:
					data = ln.split("|")
					contigName = data[3]
					contigName = re.sub("\.\S+$","",contig_name)
					contigNameFinal = contig_name.replace("_",".")

				if not contigNameFinal:
					contigCount += 1
					contigNameFinal = 'contig.' + str(contigCount)

				outChrom.write('>{2}|{0}_{1}\n'.format(newG,contigNameFinal,REPLICON[gObj.name].tag))
				outChrom.write('\n'.join(REPLICON[gObj.name].contigs[contig]))
				outChrom.write('\n')
			outChrom.close()

		for idx, replicon in enumerate(GENOME[newG]["plasmidFastas"]):
			if replicon:
				gObj = annIO.Replicon()
				gObj.name = newG + "_plasmid_" + GENOME[newG]['plasmidNames'][idx]
				gObj.tag = GENOME[newG]['locusTag']
				gObj.type = 'plasmid'
				gObj.plasmidName = GENOME[newG]['plasmidNames'][idx]
				gObj.parent = newG
				gObj.fastaPath = os.path.abspath('{0}/plasmids/{1}.fna'.format(PARAM['outputFolder'],gObj.name))
				gObj.tempPath = os.path.abspath('{0}/prodigal_temp/{1}.faa'.format(PARAM['outputFolder'],gObj.name))


				with open('{0}/{1}'.format(GENOME[newG]["inputFolder"],GENOME[newG]['plasmidFastas'][idx]),'r') as data:
					contigName = ''
					for ln in data:
						ln = ln.replace('\n','')
						ln = ln.replace('\r','')
						if '>' in ln:
							contigName = ln.replace('>','')
							gObj.contigs[contigName] = cl.deque()
						elif contigName:
							gObj.contigs[contigName].append(ln)
				REPLICON[gObj.name] = gObj

				contigCount = 0
				for contig in REPLICON[gObj.name].contigs:
					contigNameRaw = contig.replace(' ','_')
					contigNameFinal = ''
					if ">gi|" in contigNameRaw:
						data = ln.split("|")
						contigNameRaw = data[3]
						contigNameRaw = re.sub("\.\S+$","",contig_nameRaw)
						contigNameFinal = contig_nameRaw.replace("_",".")

					if not contigNameFinal:
						contigCount += 1
						contigNameFinal = 'contig.' + str(contigCount)

					outP = open(gObj.fastaPath,'w')
					data = (GENOME[newG]['locusTag'],REPLICON[gObj.name].plasmidName,REPLICON[gObj.name].name,contigNameFinal)
					outP.write('>{0}|{1}|{2}_{3}\n'.format(*data))
					outP.write('\n'.join(REPLICON[gObj.name].contigs[contig]))
					outP.write('\n')
				outP.close()

	#
	#PRODIGAL RUNNING
	#

	#Calculate length of each sequence
	for i in REPLICON:
		with open(REPLICON[i].fastaPath,'r') as f:
			seq = 0
			for ln in f:
				if ">" not in ln:
					seq += len(ln.replace('\n',''))
			REPLICON[i].length = seq

	#Build prodigal directories: raw outputs, clean outputs, partial gene/protein clean outputs
	cmdList = []
	folderLst = []

	tempPath = "{0}/prodigal_temp".format(PARAM['outputFolder'])
	os.system("mkdir -p {0}".format(tempPath))

	#Run through genome files and set up prodigal commands
	cmdList = []
	trainerLst = []
	defaultFlag = 0

	for g in REPLICON:
		found = 0
		if REPLICON[g].length >= 100000:
			data = (PARAM['prodigal_path'],REPLICON[g].tempPath,REPLICON[g].fastaPath)
			cmdList.append("{0} -m -q -a  {1}  -i {2} -o /dev/null".format(*data))
			found = 1
		else:
			#If not, check for the parent genome and use that to train prodigal
			parent = REPLICON[g].parent
			if REPLICON.has_key(parent) and REPLICON[parent].length >= 100000:
				data = (REPLICON[parent].fastaPath,parent,PARAM['prodigal_path'],tempPath)
				trainerLst.append("{2} -m -q -t {3}/{1}_training.prodigal -i {0} -o /dev/null".format(*data))

				data = (PARAM['prodigal_path'],REPLICON[g].tempPath,REPLICON[g].fastaPath,tempPath,parent)
				cmdList.append("{0} -m -q  -t {3}/{4}_training.prodigal -a {1}  -i {2} -o /dev/null".format(*data))

				found = 1
			#If still none, use the default training sequence
			if found == 0:
				data = (PARAM['prodigal_path'],REPLICON[g].tempPath,REPLICON[g].fastaPath,tempPath)
				cmdList.append("{0} -q -m -t {3}/def_training.prodigal -a {1}  -i {2} -o /dev/null".format(*data))
				defaultFlag = 1

	#Run prodigal training file building and gene prediction on multiple threads
	cores = PARAM['cores']
	if defaultFlag == 1:
		data = (PARAM['prodigal_path'],REPLICON[PARAM['defaultGenome']].tempPath,REPLICON[g].fastaPath,tempPath)
		trainerLst.append("{0} -m -q -t {3}/def_training.prodigal -i {0} -o junk.out".format(*data))

	print "Running prodigal on {0} cores...".format(cores)
	p=Pool(cores)
	for cmd in trainerLst:
		p.apply_async(cmd_runner, args=(cmd,))
		pass
	p.close()
	p.join()

	p=Pool(cores)
	for cmd in cmdList:
		p.apply_async(cmd_runner, args=(cmd,))
		pass
	p.close()
	p.join()

	#
	#PRODIGAL HEADER CLEANING
	#
	COUNT = {}

	for i in REPLICON:
		COUNT[REPLICON[i].parent] = 0

	print "Processing prodigal outputs..."
	for g in REPLICON:
		outDF = "{0}/json_features/{2}_{1}.json".format(PARAM['outputFolder'],REPLICON[g].name,REPLICON[g].tag)
		protData = cl.deque()

		#JSON data tables
		with open(REPLICON[g].tempPath, 'r') as f:
			for ln in f:
				if '>' in ln:
					c = prodigal_head_cleaner(ln)
					if c:
						protObj = annIO.Protein()
						protObj.id = c['rawName']
						protObj.length=c['length']
						protObj.partial=c['partial']
						protObj.start=c['start']
						protObj.end=c['end']
						protObj.strand=c['strand']
						protObj.contig=c['contig']
						protObj.type='cds'
						protObj.genome = REPLICON[g].parent
						protObj.replicon = g

						if PARAM['newLocusTags'] == False:
							protObj.locusNum = rawLocus.split('_')[-1]
						else:
							COUNT[protObj.genome] += 1
							locus_id_final = '%07d' % (COUNT[protObj.genome],)
							protObj.locusNum = locus_id_final
							if REPLICON[g].type == 'plasmid':
								data = (REPLICON[g].tag,REPLICON[g].plasmidName,locus_id_final)
								protObj.id = '{0}|{1}_{2}'.format(*data)
							else:
								protObj.id = '{0}|{0}_{1}'.format(REPLICON[g].tag,locus_id_final)

						if REPLICON[g].type == 'plasmid':
							protObj.plasmidName= REPLICON[g].name
							protObj.plasmidTag= REPLICON[g].plasmidName
						protData.append(protObj)
		REPLICON[g].proteinCount = len(protData)
		annIO.exportJsonData(protData,outDF,'id')

	#		Writing JSON data tables for replicons		#
	GENOME = cl.OrderedDict()
	TAG = {}

	for i in REPLICON:
		REPLICON[i].tempPath = ""
		parent = REPLICON[i].parent
		if not GENOME.has_key(parent):
			GENOME[parent] = cl.OrderedDict()
			TAG[parent] = REPLICON[i].tag
		GENOME[parent][i] = REPLICON[i]
	for i in GENOME:
		outFile = '{0}/json_genome/{2}_{1}.json'.format(PARAM['outputFolder'],i,TAG[i])
		annIO.exportJsonData(GENOME[i],outFile,"name")
	print "Successfuly extracted data from new fasta genomes\n"
	return GENOME

def extract_gene_seq(CONTIG,GENE):
	"""Extracts DNA sequence for a gene from contigs, and translates to protein seq"""

	missing = cl.deque()
	for g in GENE:
		if CONTIG.has_key(GENE[g].contig):
			start = GENE[g].start - 1
			end = GENE[g].end
			seq = CONTIG[GENE[g].contig][start:end]

			if GENE[g].strand == '-':
				seq = revcomp(seq)
			GENE[g].nuc = seq

			if GENE[g].coding == 'cds':
				GENE[g].aa = translate(seq)
			GENE[g].length = len(GENE[g].aa.replace('\n',''))
		else:
			missing.append(GENE[g].contig)
			pass
	if missing:
		missing = set(missing)
		for i in missing:
			print 'Missing contig in fasta file: {0}'.format(i)

	#Adding newlines to sequence strings for printing
	for g in GENE:
		GENE[g].nuc = break_lines(GENE[g].nuc,60)
		GENE[g].aa = break_lines(GENE[g].aa,60)

	return GENE

def fasta2dict(f):
	"""Generates a dictionary from a fasta file: dict[name] = sequence"""
	DATA = cl.OrderedDict()
	if os.path.isdir(f):
		return
	with open(f,'r') as fasta:
		for ln in fasta:
			ln = ln.replace('\n','')
			if '>' in ln:
				name = ln.split()[0].replace('>','')
				DATA[name] = cl.deque()
			else:
				DATA[name].append(ln)
	for i in DATA:
		DATA[i] = ''.join(DATA[i])
	return DATA

def break_lines(seq,num):
	"""Adds newlines to a long sequence at set intervals"""

	seqLen = len(seq)
	seqPrint = cl.deque()

	start = 0
	end = int(num)
	while end < seqLen:
		seqPrint.append("{0}\n".format(seq[start:end]))
		start += num
		end += num

	seqPrint.append("{0}\n".format(seq[start:]))
	return ''.join(seqPrint)

def runAntifam(GENOME,PARAM):
	if not PARAM['hmmscan_path']:
		print '\nNo hmmscan path, skipping Antifam curration\n'
		return GENOME

	NEED = {}
	coreNum = 2
	fNum = PARAM['cores'] / 2
	hmm_param = '--cut_tc --noali --cpu {0} --tblout '.format(coreNum)

	#Extracting all sequences
	orgType = ''
	for rep in GENOME:
		if not orgType:
			if GENOME[rep].orgType.lower() == 'bacteria':
				antifamDB = 'AntiFam_Bacteria.hmm'
			elif GENOME[rep].orgType.lower() == 'archaea':
				antifamDB = 'AntiFam_Archaea.hmm'
			elif GENOME[rep].orgType.lower() == 'eukaryote':
				antifamDB = 'AntiFam_Eukaryota.hmm'
			elif GENOME[rep].orgType.lower() == 'virus':
				antifamDB = 'AntiFam_Virus.hmm'
			else:
				antifamDB = 'AntiFam_All.hmm'

		for prot in GENOME[rep].feat:
			if GENOME[rep].feat[prot].type == 'cds':
				if not NEED.has_key(GENOME[rep].fastaPath):
					NEED[GENOME[rep].fastaPath] = cl.deque()
				seq = annIO.SeqFetch()
				seq.locus = GENOME[rep].feat[prot].id
				seq.start = GENOME[rep].feat[prot].start - 1
				seq.end = GENOME[rep].feat[prot].end
				seq.strand = GENOME[rep].feat[prot].strand
				seq.contig = GENOME[rep].feat[prot].contig
				NEED[GENOME[rep].fastaPath].append(seq)
	DEBUG = annIO.extractSeq(NEED,PARAM['outputFolder'],'prot',fNum)

	fastaNames = os.listdir('{0}/tempFastas'.format(PARAM['outputFolder']))
	fastas = []
	for f in fastaNames:
		fastas.append('{0}/tempFastas/{1}'.format(PARAM['outputFolder'],f))

	cmd = []
	for f in fastas:
		if os.path.getsize(f) > 0:
			name = str(f)
			name = name.split('/')[-1].split('.')[:-1]
			name = '.'.join(name)

			inputs = (PARAM['hmmscan_path'],PARAM['outputFolder'],name,hmm_param,f,antifamPath,antifamDB)
			cmd.append('{0} {3} {1}/db_data/{2}.temp.antifam.out {5}/{6} {4} > /dev/null'.format(*inputs))

	if fNum > 1:
		p=Pool(fNum)
		for i in cmd:
			p.apply_async(cmd_runner, args=(i,))
			pass #for debugging
		p.close()
		p.join()
	else:
		for c in cmd:
			os.system(c)

	#Output files
	outs = os.listdir('{0}/db_data'.format(PARAM['outputFolder']))

	count = 0
	for f in outs:
		if '.temp.{0}.out'.format("antifam") in f:
			count += 1
			PROT = ann.antifam_best_hit("{0}/db_data/{1}".format(PARAM['outputFolder'],f))
		for prot in PROT:
			if PROT[prot].ann['antiFam']:
				for contig in GENOME:
					if GENOME[contig].feat.has_key(prot):
						GENOME[contig].feat[prot].failedAntifam = True
	if count == 0:
		print "Failed to run antifam"
		sys.exit()
	else:
		return GENOME, PROT

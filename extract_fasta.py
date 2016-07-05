import sys
import os
import json
import collections as cl

usage = """\n

Generates fasta files for features from the dataset. Assumes standard genetic code.

<dataset folder>
<Feature type to extract: cds, rna, misc, or all>
<sequence type to output: aa or nuc>
<output folder name>\n\n"""
argv = sys.argv[1:]

if len(argv) == 4:
	inFolder= argv.pop(0)
	featType = argv.pop(0)
	seqType = argv.pop(0)
	outFolder = argv.pop(0)
else:
	sys.exit(usage)

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

def codon_table():
	"""Builds a lookup table for codon -> aa"""
	CODON_TABLE = {
		"ACC":"T","ATG":"M","ACA":"T","ACG":"T","ATC":"I","AAC":"N","ATA":"I","AGG":"R",
		"CCT":"P","CTC":"L","AGC":"S","AAG":"K","AGA":"R","CAT":"H","AAT":"N","ATT":"I",
		"CTG":"L","CTA":"L","ACT":"T","CAC":"H","AAA":"K","CAA":"Q","AGT":"S","CCA":"P",
		"CCG":"P","CCC":"P","CTT":"L","TAT":"Y","GGT":"G","TGT":"C","CGA":"R","CAG":"Q",
		"CGC":"R","GAT":"D","CGG":"R","TTT":"F","TGC":"C","GGG":"G","TGA":"*","GGA":"G",
		"TGG":"W","GGC":"G","TAC":"Y","TTC":"F","TCG":"S","TAG":"*","TTG":"L","CGT":"R",
		"GAA":"E","TCA":"S","GCA":"A","GTA":"V","GCC":"A","GTC":"V","GCG":"A","GTG":"V",
		"GAG":"E","GTT":"V","GCT":"A","TTA":"L","GAC":"D","TCC":"S","TAA":"*","TCT":"S",
		'NNN':'X'}
	return CODON_TABLE

def fasta2dict(f):
	"""Generates a dictionary from a fasta file: dict[name] = sequence"""
	import collections as cl

	DATA = {}
	if os.path.isdir(f):
		print "ERROR: cannot read a directory"
		return
	with open(f,'r') as fasta:
		for ln in fasta:
			ln = ln.replace('\n','')
			if '>' in ln:
				name = ln.split()[0].replace('>','')
				name = name.split('_')[-1]
				DATA[name] = cl.deque()
			else:
				DATA[name].append(ln)
	for i in DATA:
		DATA[i] = ''.join(DATA[i])
	return DATA

def translate(seq,CODON_TABLE):
	aaSeq = []
	ln = cl.deque(seq)
	while len(ln) > 2:
		codon = ''
		codon += ln.popleft()
		codon += ln.popleft()
		codon += ln.popleft()
		if CODON_TABLE.has_key(codon.upper()):
			aaSeq.append(CODON_TABLE[codon.upper()])
		else:
			aaSeq.append('X')
	return ''.join(aaSeq)


os.system('mkdir -p {0}'.format(outFolder))
os.system('mkdir -p {0}/replicons'.format(outFolder))

if seqType == 'aa':
	CODON = codon_table()
	ext = 'faa'
else:
	ext = 'fna'

FASTAS = {}
chroms = os.listdir('{0}/chromosomes'.format(inFolder))
for i in chroms:
	FASTAS[i.replace('.fna','')] = '{0}/chromosomes/{1}'.format(inFolder,i)
plasmids = os.listdir('{0}/plasmids'.format(inFolder))
for i in plasmids:
	FASTAS[i.replace('.fna','')] = '{0}/plasmids/{1}'.format(inFolder,i)

jsons = os.listdir('{0}/json_features'.format(inFolder))
for j in jsons:
	name = j.replace('.json','')
	repName = '_'.join(name.split('_')[1:])
	seq = []
	if FASTAS.has_key(repName):
		RAW = fasta2dict(FASTAS[repName])
	else:
		print "Cant find fasta file for feature table {0}".format(j)

	with open('{0}/json_features/{1}'.format(inFolder,j),'r') as f:
		data = json.load(f,object_pairs_hook=cl.OrderedDict)
		for g in data:
			if not featType.lower() == 'all':
				if not data[g]['type'] == featType.lower():
					continue
			gName = '{0} {1}'.format(g,data[g]['function'])
			contig = data[g]['contig'].split('_')[-1]
			contigSeq = RAW[contig]
			geneSeq = contigSeq[data[g]['start'] - 1:data[g]['end']]
			geneSeq = ''.join(geneSeq)
			if data[g]['strand'] == "-1":
				geneSeq = revcomp(geneSeq)
			if seqType == 'aa':
				geneSeq = translate(geneSeq,CODON)
			seq.append((gName,geneSeq))

	outFile = open('{0}/replicons/{1}.{2}'.format(outFolder,name,ext),'w')
	for g in seq:
		outFile.write('>{0}\n{1}\n'.format(g[0],g[1]))
	outFile.close()

os.system('cat {0}/replicons/*.{1} > {0}/all_seqs.{1}'.format(outFolder,ext))

	
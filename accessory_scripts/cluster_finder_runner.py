import sys
import os
import json
import collections as cl
from multiprocessing import Pool

####################################################################################################
usage = """\n
Runs ClusterFinder on a set of genomes. Genes on the 5' and 3' ends of clusters
are removed if they don't contain a domain above the minimum probabilty cutoff.
Can also separate clusters that lack a gene with at least one antiSMASH domain
hit if desired.

Pfam name conversion table should look like this:
PF04951.8	Peptidase_M55
PF11663.3	Toxin_YhaV
PF08677.5	GP11
PF06362.6	DUF1067
etc.

Requires annotation folder with PFAM database annotations. AntiSMASH annotations
also required if these are used to sort clusters.


<Annotation folder>
<Genome list, or ALL>
<Pfam name conversion table>
<Min probability to retain edge genes>
<Min gene number for a cluster>
<antiSMASH annotated clusters only? (Y/N)>
<cores>
<output folder>\n\n"""
argv = sys.argv[1:]

if len(argv) == 8:
	inFolder = argv.pop(0)
	gList = argv.pop(0)
	pfamTbl = argv.pop(0)
	minScore = float(argv.pop(0))
	if minScore > 1:
		minScore = minScore / 100
	minGenes = int(argv.pop(0))
	smashFlag = argv.pop(0)
	cores = int(argv.pop(0))
	outFolder = argv.pop(0)
else:
	sys.exit(usage)

class Genome():
	def __init__(self):
		self.name = ''
		self.genes = cl.OrderedDict()
		self.clusters = cl.OrderedDict()

class Gene():
	def __init__(self):
		self.id = ''
		self.seq_status = 'draft'
		self.org_name = ''
		self.contig = ''
		self.org_id = ''
		self.locus_tag = ''
		self.start = ''
		self.end = ''
		self.strand = ''
		self.pfam_t_start = 0
		self.pfam_t_end = 0
		self.pfam_start = 0
		self.pfam_end = 0
		self.pfam_id = ''
		self.pfam_score = 0
		self.enzyme_id = ''
		self.antismash = 0		#Flag for presence of an antismash annotated domain

class Cluster(object):
	def __init__(self):
		self.name = ''

		#List indices for the following lists of genes below the probability cutoff
		self.badGenes = {} 

		#list of lists, internal lists are data lines for each domain in each prot
		self.dataLines = cl.OrderedDict()

		#list of gene names
		self.genes = cl.OrderedDict()
 		
		#list of lists, each internal list contains the data on all domains in a prot
		self.probs = cl.OrderedDict()		#Probabilities
		self.starts = cl.OrderedDict()		#Start coords
		self.ends = cl.OrderedDict()		#End coords
		self.domains = cl.OrderedDict()		#Pfam domain ids
		self.goodGenes = 0					#Count of good genes after currating
		self.noSmash = 0				#Flag false positives if antismash domains not present

	def currate(self,cutoff):
		"""Identify genes whose domains all have low probability scores"""
		for g in self.probs:
			maxScore = 0.0
			for p in self.probs[g]:
				if p > maxScore:
					maxScore = p
			if maxScore < cutoff:
				self.badGenes[g] = ''
		self.goodGenes = len(self.genes) - len(self.badGenes)

def fasta2dict(f):
	"""Generates a dictionary from a fasta file: dict[name] = sequence"""
	DATA = {}
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

def getHeader(f,out):
	"""Get fasta headers from a file"""
	o = open(out,'w')

	with open(f,'r') as data:
		for ln in data:
			if '>' in ln:
				o.write(ln)
	o.close()

def files2fasta(folder,output):
	DATA = {}
	files = os.listdir(folder)
	for f in files:
		name = f.split('.')[0]
		DATA[name] = []
		with open('{0}/{1}'.format(folder,f),'r') as things:
			for thing in things:
				thing = thing.translate(None,'>\n')
				DATA[name].append(thing)
	out = open(output,'w')
	for thing in DATA:
		out.write('>{0}\n{1}\n'.format(thing,','.join(DATA[thing])))
	out.close()


#Simple command multithreader
def multi_run(cmd):
	os.system(cmd)

#antiSMASH HMMer models for fatty acid biosynthesis false positives
FALSE_POSITIVE = {'bt1fas':'','ft1fas':'','t2fas':'','fabH':''}
####################################################################################################


#Generating output folders and reading basic data tables
os.system('mkdir -p {0}'.format(outFolder))
os.system('mkdir -p {0}/inputs'.format(outFolder))
os.system('mkdir -p {0}/outputs'.format(outFolder))
os.system('mkdir -p {0}/currated_clusters'.format(outFolder))
os.system('mkdir -p {0}/raw_outputs'.format(outFolder))
os.system('mkdir -p {0}/cluster_prots'.format(outFolder))
os.system('mkdir -p {0}/cluster_genes'.format(outFolder))
CF_PATH = '~/tools_and_software/ClusterFinder-master/ClusterFinder.py'

PFAM_TRANS = {}
with open(pfamTbl,'r') as f:
	for ln in f:
		ln = ln.strip().split('\t')
		PFAM_TRANS[ln[1]] = ln[0].split('.')[0]

GENOME = {}
if gList != 'ALL':
	with open(gList,'r') as f:
		for ln in f:
			ln = ln.replace('\n','')
			ln = ln.split()
			GENOME[ln[0]] = Genome()

#Getting coordinates of PFAM domain
COORDS = {}
coord_files = os.listdir('{0}/db_data/PFAM'.format(inFolder))
for f in coord_files:
	with open('{0}/db_data/PFAM/{1}'.format(inFolder,f),'r') as data:
		for lnD in data:
			ln = lnD.strip().split()
			if len(ln) != 5:
				print 'error in coord file:\t{0}'.format(f)
				print ln
				break
			gene = ln.pop(0)
			if OLD2NEW.has_key(gene):
				geneName = OLD2NEW[gene]
				if not COORDS.has_key(geneName):
					COORDS[geneName] = {'start': [], 'end': []}
				COORDS[geneName]['start'].append(ln[1])
				COORDS[geneName]['end'].append(ln[2])

#Parsing data from the annotation tables
anns = os.listdir('{0}/annotation_tables'.format(inFolder))
for f in anns:
	data = open('{0}/annotation_tables/{1}'.format(inFolder,f),'r')
	D = data.readlines()
	data.close()

	genome_name = f.split('.')[0]
	if 'plasmid' in genome_name:
		genome_name = genome_name.split('_plasmid')[0]

	header = D.pop(0)
	header = header.split()
	pfam_col = 0
	pfam_score_col = 0
	for idx, i in enumerate(header):
		if i == 'PFAM':
			pfam_col = idx
		elif i == 'PFAM_bitscore':
			pfam_score_col = idx
		elif i == 'antiSMASH':
			antismashCol = idx

	for gene in D:
		gene = gene.strip().split('\t')
		genome = gene[0].split('|')[0]
		if gList == 'ALL' and not GENOME.has_key(genome):
			GENOME[genome] = Genome()
		if not GENOME.has_key(genome):
			break

		if NEW2OLD.has_key(gene[0]) and gene[pfam_col] != '-':
			g = Gene()
			g.id = gene[0]
			g.org_name = genome_name
			g.contig = gene[7]
			g.org_id = genome
			g.locus_tag = gene[0]
			g.start = gene[2]
			g.end = gene[3]
			if gene[6] == '1':
				g.strand = '+'
			elif gene[6] == '-1':
				g.strand = '-'
			g.pfam_score = gene[pfam_score_col].split(',')
			g.pfam_id = gene[pfam_col].split(',')
			g.pfam_start =  COORDS[g.id]['start']
			g.pfam_end = COORDS[g.id]['end']
			if smashFlag == 'Y':
				if gene[antismashCol] != '-':
					smashHit = gene[antismashCol].split(',')
					for i in smashHit:
						if not FALSE_POSITIVE.has_key(i):
							g.antismash = 1
			else:
				g.antismash = 1

			GENOME[genome].genes[g.id] = g

#Generating ClusterFinder inputs
for genome in GENOME:
	out = open('{0}/inputs/{1}.cf_input'.format(outFolder,genome),'w')
	for gene in GENOME[genome].genes:
		g = GENOME[genome].genes[gene]
		data_all = (g.id,'Genome',g.org_name,g.contig,g.org_id,g.locus_tag,g.start,g.end,g.strand)

		for idx, pfam in enumerate(g.pfam_id):
			pfamTrans = PFAM_TRANS[g.pfam_id[idx].split('/')[0]]
			data_pfam = ('0','0',g.pfam_start[idx],g.pfam_end[idx],pfamTrans,g.pfam_score[idx],'na')

			out.write('\t'.join(data_all) + '\t' + '\t'.join(data_pfam) + '\n')
	out.close()

#Running Clusterfinder
cf_inputs = os.listdir('{0}/inputs'.format(outFolder))
cmd_run = []
for i in cf_inputs:
	name = i.split('.')[0]
	cmd_run.append('python {0} {1}/inputs/{2} {1}/outputs/{3}'.format(CF_PATH,outFolder,i,name))

p=Pool(cores)
for i in cmd_run:
	p.apply_async(multi_run, args=(i,))
	pass
p.close()
p.join()

os.system('mv {0}/inputs/*.out {0}/raw_outputs/'.format(outFolder))

#Ugly parser for uncurrated ClusterFinder outputs
rawClusters = os.listdir('{0}/outputs'.format(outFolder))
for f in rawClusters:
	with open(outFolder + '/outputs/' + f,'r') as data:
		oldc = ''
		oldg = ''
		g = ''
		c = ''

		for ln in data:
			ln = ln.replace('\n','')
			lnData = ln.split()
			if c:
				oldc = c
				oldg = g
			c = lnData[0]
			g = lnData[2]
			genome = g.split('|')[0]
			prob = float(lnData[-1])
			domain = lnData[-2]
			start = lnData[3]
			end = lnData[4]

			if not GENOME[genome].clusters.has_key(c):
				if oldc:
					GENOME[genome].clusters[oldc].probs[oldg] = probList
					GENOME[genome].clusters[oldc].domains[oldg] = domList
					GENOME[genome].clusters[oldc].starts[oldg] = startList
					GENOME[genome].clusters[oldc].ends[oldg] = endList
					GENOME[genome].clusters[oldc].dataLines[oldg] = lnList
				probList = []
				domList = []
				startList = []
				endList = []
				lnList = []
				GENOME[genome].clusters[c] = Cluster()
				GENOME[genome].clusters[c].name = c

			if not GENOME[genome].clusters[c].genes.has_key(g):
				GENOME[genome].clusters[c].genes[g] = ''
				if probList:
					GENOME[genome].clusters[c].probs[oldg] = probList
					GENOME[genome].clusters[c].domains[oldg] = domList
					GENOME[genome].clusters[c].starts[oldg] = startList
					GENOME[genome].clusters[c].ends[oldg] = endList
					GENOME[genome].clusters[c].dataLines[oldg] = lnList
				probList = []
				domList = []
				startList = []
				endList = []
				lnList = []
			probList.append(prob)
			domList.append(domain)
			startList.append(start)
			endList.append(end)
			lnList.append(ln)
		if probList:
			GENOME[genome].clusters[c].probs[g] = probList
			GENOME[genome].clusters[c].domains[g] = domList
			GENOME[genome].clusters[c].starts[g] = startList
			GENOME[genome].clusters[c].ends[g] = endList
			GENOME[genome].clusters[c].dataLines[g] = lnList

#Currating clusters
for genome in GENOME:
	for c in GENOME[genome].clusters:
		GENOME[genome].clusters[c].currate(minScore)

#Check for antismash domains if required
if smashFlag == 'Y':
	for genome in GENOME:
		for c in GENOME[genome].clusters:
			clean = 0
			for gene in GENOME[genome].clusters[c].genes:
				if GENOME[genome].genes[gene].antismash:
					clean = 1
			if not clean:
					GENOME[genome].clusters[c].noSmash = 1


#
#Printing clusters that passed all curration steps
#

#Printing currated output files
for g in GENOME:
	out = open(outFolder + '/currated_clusters/' + g + '.cf','w')
	for c in GENOME[g].clusters:
		if GENOME[g].clusters[c].goodGenes >= minGenes and GENOME[g].clusters[c].noSmash == 0:
			for gene in GENOME[g].clusters[c].dataLines:
				if not GENOME[g].clusters[c].badGenes.has_key(gene):
					out.write('\n'.join(GENOME[g].clusters[c].dataLines[gene]) + '\n')
	out.close()

#Printing cluster domain content
out = open('{0}/cluster_domains.lst'.format(outFolder),'w')
for g in GENOME:
	for c in GENOME[g].clusters:
		if GENOME[g].clusters[c].goodGenes >= minGenes and GENOME[g].clusters[c].noSmash == 0:
			out.write('>{0}|{1}\n'.format(g,c))
			goodDomains = cl.deque()
			for gene in GENOME[g].clusters[c].domains:
				if not GENOME[g].clusters[c].badGenes.has_key(gene):
					goodDomains.append(','.join(GENOME[g].clusters[c].domains[gene]))
			goodDomains = ','.join(goodDomains)
			out.write(goodDomains + '\n')
out.close()

#Printing cluster sequences
PROTS = fasta2dict('{0}/seqs/Omes/all_prots.faa'.format(inFolder))
PARTIAL = fasta2dict('{0}/seqs/Omes/all_partial.faa'.format(inFolder))
for g in GENOME:
	for c in GENOME[g].clusters:
		if GENOME[g].clusters[c].goodGenes >= minGenes and GENOME[g].clusters[c].noSmash == 0:
			out = open('{0}/cluster_prots/{1}.faa'.format(outFolder,c),'w')
			for gene in GENOME[g].clusters[c].genes:
				if PROTS.has_key(gene) and not GENOME[g].clusters[c].badGenes.has_key(gene):
					out.write('>{0}\n{1}\n'.format(gene,PROTS[gene]))
				elif PARTIAL.has_key(gene) and not GENOME[g].clusters[c].badGenes.has_key(gene):
					out.write('>{0}\n{1}\n'.format(gene,PARTIAL[gene]))
			out.close()
del PROTS
del PARTIAL

NUC = fasta2dict('{0}/seqs/Omes/all_genes.fna'.format(inFolder))
PARTIAL = fasta2dict('{0}/seqs/Omes/all_partial.fna'.format(inFolder))
for g in GENOME:
	for c in GENOME[g].clusters:
		if GENOME[g].clusters[c].goodGenes >= minGenes and GENOME[g].clusters[c].noSmash == 0:
			out = open('{0}/cluster_genes/{1}.fna'.format(outFolder,c),'w')
			for gene in GENOME[g].clusters[c].genes:
				if NUC.has_key(gene) and not GENOME[g].clusters[c].badGenes.has_key(gene):
					out.write('>{0}\n{1}\n'.format(gene,NUC[gene]))
				elif PARTIAL.has_key(gene) and not GENOME[g].clusters[c].badGenes.has_key(gene):
					out.write('>{0}\n{1}\n'.format(gene,PARTIAL[gene]))
			out.close()
del NUC
del PARTIAL

os.system('mkdir -p {0}/cluster_gene_lists'.format(outFolder))
clusterFastas = os.listdir('{0}/cluster_prots/'.format(outFolder))
for fasta in clusterFastas:
	name = fasta.split('.')[0]
	getHeader(outFolder + '/cluster_prots/' + fasta,'{0}/cluster_gene_lists/{1}.lst'.format(outFolder,name))
files2fasta('{0}/cluster_gene_lists'.format(outFolder),'{0}/genes_by_cluster.out'.format(outFolder))

#
#Printing clusters with no antismash domains if this was used to currate clusters
#

if smashFlag == 'Y':
	os.system('mkdir -p {0}/missing_antismash_clusters'.format(outFolder))
	os.system('mkdir -p {0}/missing_antismash_cluster_prots'.format(outFolder))
	os.system('mkdir -p {0}/missing_antismash_cluster_genes'.format(outFolder))

	#Printing currated output files
	for g in GENOME:
		out = open(outFolder + '/missing_antismash_clusters/' + g + '.cf','w')
		for c in GENOME[g].clusters:
			if GENOME[g].clusters[c].goodGenes >= minGenes and GENOME[g].clusters[c].noSmash == 1:
				for gene in GENOME[g].clusters[c].dataLines:
					if not GENOME[g].clusters[c].badGenes.has_key(gene):
						out.write('\n'.join(GENOME[g].clusters[c].dataLines[gene]) + '\n')
		out.close()

	#Printing cluster domain content
	out = open('{0}/missing_antismash_cluster_domains.lst'.format(outFolder),'w')
	for g in GENOME:
		for c in GENOME[g].clusters:
			if GENOME[g].clusters[c].goodGenes >= minGenes and GENOME[g].clusters[c].noSmash == 1:
				out.write('>{0}|{1}\n'.format(g,c))
				goodDomains = cl.deque()
				for gene in GENOME[g].clusters[c].domains:
					if not GENOME[g].clusters[c].badGenes.has_key(gene):
						goodDomains.append(','.join(GENOME[g].clusters[c].domains[gene]))
				goodDomains = ','.join(goodDomains)
				out.write(goodDomains + '\n')
	out.close()

	#Printing cluster sequences
	PROTS = fasta2dict('{0}/seqs/Omes/all_prots.faa'.format(inFolder))
	PARTIAL = fasta2dict('{0}/seqs/Omes/all_partial.faa'.format(inFolder))
	for g in GENOME:
		for c in GENOME[g].clusters:
			if GENOME[g].clusters[c].goodGenes >= minGenes and GENOME[g].clusters[c].noSmash == 1:
				out = open('{0}/missing_antismash_cluster_prots/{1}.faa'.format(outFolder,c),'w')
				for gene in GENOME[g].clusters[c].genes:
					if PROTS.has_key(gene) and not GENOME[g].clusters[c].badGenes.has_key(gene):
						out.write('>{0}\n{1}\n'.format(gene,PROTS[gene]))
					elif PARTIAL.has_key(gene) and not GENOME[g].clusters[c].badGenes.has_key(gene):
						out.write('>{0}\n{1}\n'.format(gene,PARTIAL[gene]))
				out.close()
	del PROTS
	del PARTIAL

	NUC = fasta2dict('{0}/seqs/Omes/all_genes.fna'.format(inFolder))
	PARTIAL = fasta2dict('{0}/seqs/Omes/all_partial.fna'.format(inFolder))
	for g in GENOME:
		for c in GENOME[g].clusters:
			if GENOME[g].clusters[c].goodGenes >= minGenes and GENOME[g].clusters[c].noSmash == 1:
				out = open('{0}/missing_antismash_cluster_genes/{1}.fna'.format(outFolder,c),'w')
				for gene in GENOME[g].clusters[c].genes:
					if NUC.has_key(gene) and not GENOME[g].clusters[c].badGenes.has_key(gene):
						out.write('>{0}\n{1}\n'.format(gene,NUC[gene]))
					elif PARTIAL.has_key(gene) and not GENOME[g].clusters[c].badGenes.has_key(gene):
						out.write('>{0}\n{1}\n'.format(gene,PARTIAL[gene]))
				out.close()
	del NUC
	del PARTIAL

os.system('mkdir -p {0}/missing_antismash_cluster_gene_lists'.format(outFolder))
clusterFastas = os.listdir('{0}/missing_antismash_cluster_prots/'.format(outFolder))
for fasta in clusterFastas:
	name = fasta.split('.')[0]
	getHeader(outFolder + '/missing_antismash_cluster_prots/' + fasta,'{0}/missing_antismash_cluster_gene_lists/{1}.lst'.format(outFolder,name))
files2fasta('{0}/missing_antismash_cluster_gene_lists'.format(outFolder),'{0}/missing_antismash_genes_by_cluster.out'.format(outFolder))

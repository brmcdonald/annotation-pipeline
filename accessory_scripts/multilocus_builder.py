import sys
import os
from multiprocessing import Pool
import collections as cl

usage = """\n
<annotation table folder>
<amino acid fasta file, or NA if not needed>
<nucleotide fasta file, or NA if not needed>
<List of gene families to align>
<Database name>
<Locus tag list of genomes to include, or ALL>
<Num cores to use>
<Bootstrap num (1 = no bootstraps)>
<(P)rotein, (N)ucleotide, or (C)odon alignment>
<(R)un tree or (O)utput alignment and RAxML command>
<output folder>\n\n"""
argv = sys.argv
junk = argv.pop(0)

if len(sys.argv) == 11:
	ann_folder = argv.pop(0)
	aaFile = argv.pop(0)
	nucFile = argv.pop(0)
	fam_lst = argv.pop(0)
	db_flag = argv.pop(0)
	tag_lst = argv.pop(0)
	cores = int(argv.pop(0))
	boot_num = int(argv.pop(0))
	seq_flag = argv.pop(0)
	run_flag = argv.pop(0)
	out_folder = argv.pop(0)
else:
	sys.exit(usage)

#
#Mafft path
#
mafft_path = "mafft --retree 2 --quiet --maxiterate 1000 --thread 2"

#
#RAxML Path
#
raxml = "/home/bradon/tools_and_software/RAxML-8.1.24/raxmlHPC-PTHREADS-SSE3"

#Multithreader
def cmd_runner(cmd):
	os.system(cmd)

class Gene():
	def __init__(self,name,fam,score,tag):
		self.name = name
		self.fam = fam
		self.score = float(score)
		self.aa = ""
		self.nuc = ""
		self.genome = tag

class Genome():
	def __init__(self,tag):
		self.tag = tag
		self.name = ""
		self.final = {}		#Dict of gene objects that will be used in the tree
		self.genes = {}		#Dict of gene objects

	#Identifies the gene with the highest bit score in each gene family
	def best_fam_member(self):
		for g in self.genes:
			if self.final.has_key(self.genes[g].fam):
				if self.genes[g].score > self.final[self.genes[g].fam].score:
					self.final[self.genes[g].fam] = self.genes[g]
			else:
				self.final[self.genes[g].fam] = self.genes[g]

def fasta2dict(f):
	"""Generates a dictionary from a fasta file: dict[name] = sequence"""
	DATA = {}
	if os.path.isdir(f):
		print "ERROR: cannot read a directory"
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

#Output name
outName = out_folder.split('/')[-1]

####################################################################################################
#######################################   Script Start   ###########################################
####################################################################################################

#Which columns of the annotation tables to get data from

os.system("mkdir -p {0}".format(out_folder))
if seq_flag == "P" or seq_flag == "C":
	os.system("mkdir -p {0}/faa".format(out_folder))
	os.system("mkdir -p {0}/faa_align".format(out_folder))
if seq_flag == "N" or seq_flag == "C":
	os.system("mkdir -p {0}/fna".format(out_folder))
	os.system("mkdir -p {0}/fna_align".format(out_folder))

GENOME = {}		#Dict of genome objects, which each contain gene objects
FAM = {}		#Dict to store which families we're aligning


fam_col = 0
partial_idx = 0
score_col = 0

ann_files = os.listdir(ann_folder)


#Record genomes that need to be in the tree
if tag_lst != "ALL":
	f = open(tag_lst,'r')
	f = f.readlines()
	for i in f:
		i = i.replace('\n','')
		if "\t" in i:
			i = i.split()
			x = Genome(i[0])
			x.name = i[1]
			GENOME[i[0]] = x
		else:
			i = i.replace('\n','')
			if i:
				x = Genome(i)
				GENOME[i] = x
else:
	for f in ann_files:
		f = f.replace(".ann","")
		genome = '_'.join(f.split('_')[1:])
		tag = f.split('_')[0]
		x = Genome(tag)
		GENOME[tag] = x

#Record gene families to build the tree with
f = open(fam_lst,'r')
f = f.readlines()
for i in f:
	i = i.split()
	FAM[i[0]] = ""

#Ugly hack to get list of annotation tables to read
read = []
for f in ann_files:
	if not 'plasmid' in f:
		genome = '.'.join(f.split('.')[:-1])
		fTag = genome.split('_')[0]
		found = 0
		for tag in GENOME:
			if tag_lst == 'ALL' or GENOME.has_key(fTag):
				found = 1
		if found == 1:
			read.append(f)
			found = 0

print "\nExtracting annotation data..."
#Pull bitscore data for each gene in the desired families from the annotation tables
for f in read:
	fi = open("{0}/{1}".format(ann_folder,f),'r')
	fi = fi.readlines()
	header = fi.pop(0).split('\t')

	#Figure out which columns the data for the desired database are stored in
	if fam_col == 0:
		for idx, h in enumerate(header):
			if h == db_flag:
				fam_col = idx
				score_col = idx+1
			if 'Partial' in h:
				partial_idx = idx

	for ln in fi:
		ln = ln.replace("\n","")
		data = ln.split("\t")

		#Complete genes only
		if data[fam_col] and int(data[partial_idx]) == 0:
			if FAM.has_key(data[fam_col]):
				tag = data[0].split('|')[0]
				gene = Gene(data[0],data[fam_col],data[score_col],tag)
				GENOME[tag].genes[data[0]] = gene

#Find best family member for each genome and extract sequences from fasta files
if seq_flag == "N" or seq_flag == 'C':
	if nucFile == 'NA':
		print "Need nucleotide sequences to make this tree"
		sys.exit()
	NUC_SEQ = fasta2dict(nucFile)

if seq_flag == "P" or seq_flag == 'C':
	if aaFile == 'NA':
		print "Need amino acid sequences to make this tree"
		sys.exit()
	AA_SEQ = fasta2dict(aaFile)

print 'Identifying orthologous genes in each genome...'
for genome in GENOME:
	record = 0
	GENOME[genome].best_fam_member()

	#Record nuc seq if needed
	if seq_flag == "N" or seq_flag == 'C':
		for fam in GENOME[genome].final:
			if NUC_SEQ.has_key(GENOME[genome].final[fam].name):
				GENOME[genome].final[fam].nuc = NUC_SEQ[GENOME[genome].final[fam].name]
			else:
				print "No nucleotide sequence for {0}".format(GENOME[genome].final[fam].name)

	#Record prot seq if needed
	if seq_flag == "P" or seq_flag == 'C':
		for fam in GENOME[genome].final:
			if AA_SEQ.has_key(GENOME[genome].final[fam].name):
				GENOME[genome].final[fam].aa = AA_SEQ[GENOME[genome].final[fam].name]
			else:
				print "No amino acid sequence for {0}".format(GENOME[genome].final[fam].name)



print "Extracting sequence data for best orthologs"
#Print fasta files for the best ortholog in each gene family
for fam in FAM:
	if seq_flag == "N":
		outN = open("{0}/fna/{1}.fna".format(out_folder,fam),'w')
		for genome in GENOME:
			if GENOME[genome].final.has_key(fam) and GENOME[genome].final[fam].nuc:
				g = GENOME[genome].final[fam].name
				outN.write(">{0}\n{1}\n".format(g,''.join(GENOME[genome].final[fam].nuc)))
		outN.close()

	elif seq_flag == "P":
		outP = open("{0}/faa/{1}.faa".format(out_folder,fam),'w')
		for genome in GENOME:
			if GENOME[genome].final.has_key(fam) and GENOME[genome].final[fam].aa:
				g = GENOME[genome].final[fam].name
				outP.write(">{0}\n{1}\n".format(g,''.join(GENOME[genome].final[fam].aa)))
		outP.close()

	elif seq_flag == "C":
		outN = open("{0}/fna/{1}.fna".format(out_folder,fam),'w')
		outP = open("{0}/faa/{1}.faa".format(out_folder,fam),'w')
		for genome in GENOME:
			if GENOME[genome].final.has_key(fam) and GENOME[genome].final[fam].nuc:
				g = GENOME[genome].final[fam].name
				outP.write(">{0}\n{1}\n".format(g,''.join(GENOME[genome].final[fam].aa)))
				outN.write(">{0}\n{1}\n".format(g,''.join(GENOME[genome].final[fam].nuc)))
		outP.close()
		outN.close()

print 'Running alignments...'
#Run mafft on each sequence
if seq_flag == "P" or seq_flag == "C":
	dirList = os.listdir("{0}/faa".format(out_folder))
	cmdLst = []
	for f in dirList:
		name = f
		name = name.replace(".fna",'')
		name = name.replace(".faa",'')
		name = name.replace(".fasta",'')
		name = name.replace("_clean.list", '')

		size = os.path.getsize('{0}/faa/{1}'.format(out_folder,f))
		if size > 0:
			cmd = "{0} {1}/faa/{3} > ./{1}/faa_align/{2}.mafft".format(mafft_path,out_folder,name,f)
			cmdLst.append(cmd)

elif seq_flag == "N":
	dirList = os.listdir("{0}/fna".format(out_folder))
	cmdLst = []
	for f in dirList:
		name = f
		name = name.replace(".fna",'')
		name = name.replace(".faa",'')
		name = name.replace(".fasta",'')
		name = name.replace("_clean.list", '')

		size = os.path.getsize('{0}/fna/{1}'.format(out_folder,f))
		if size > 0:
			cmd = "{0} {1}/fna/{3} > ./{1}/fna_align/{2}.mafft".format(mafft_path,out_folder,name,f)
			cmdLst.append(cmd)

mafft_cores = cores / 2
if mafft_cores < 1:
	mafft_cores = 1
p=Pool(mafft_cores)
for cmd in cmdLst:
	p.apply_async(cmd_runner, args=(cmd,))
	#pass
p.close()
p.join()

g_list = open("{0}/genomes.lst".format(out_folder),'w')
for i in GENOME:
	g_list.write("{0}\n".format(i))
g_list.close()

#Get rid of any gene families that didn't have alignments
if seq_flag == "P" or seq_flag == "C":
	aligns = os.listdir("{0}/faa_align".format(out_folder))
	for f in aligns:
		size = os.path.getsize('{0}/faa_align/{1}'.format(out_folder,f))
		if size == 0:
			os.system("rm -f {0}/faa_align/{1}".format(out_folder,f))
elif seq_flag == "N":
	aligns = os.listdir("{0}/fna_align".format(out_folder))
	for f in aligns:
		size = os.path.getsize('{0}/fna_align/{1}'.format(out_folder,f))
		if size == 0:
			os.system("rm -f {0}/fna_align/{1}".format(out_folder,f))

if seq_flag == "C":
	#Convert protein to nucleotide
	print 'Converting to codon alignments...'
	converter = "python ~/scripts/fasta_manipulation/aa2nuc_align.py"
	os.system("{0} {1}/faa_align {2} {1}/fna_align".format(converter,out_folder,nucFile))

if seq_flag == "N" or seq_flag == "C":
	#Combine nucleotide alignments
	combiner = "perl ~/scripts/Z_Random_Old/concatinate_by_genome.pl"
	os.system("{0} {1}/fna_align {1}/genomes.lst P N".format(combiner,out_folder))

if seq_flag == "P":
	#Combine protein alignments
	combiner = "perl ~/scripts/Z_Random_Old/concatinate_by_genome.pl"
	os.system("{0} {1}/faa_align {1}/genomes.lst P N".format(combiner,out_folder))

os.system("mv full_alignment.out {0}/{1}_align.phylip".format(out_folder,outName))
os.system('perl ~/scripts/fasta_manipulation/phylip2fasta.pl {0}/{1}_align.phylip > {0}/{1}_align.fasta'.format(out_folder,outName))

#Run RAxML
if run_flag == "R":
	if seq_flag == "N" or seq_flag == "C":
		if boot_num == 1:
			os.system("{0} -s {1}/{3}_align.phylip -n multilocus_{1} -T {2} -m GTRGAMMA -p 53715".format(raxml,out_folder,cores,outName))
		else:
			os.system("{0} -s {1}/{4}_align.phylip -n multilocus_{1} -T {2} -m GTRGAMMA  -p 53715 -x 3524627 -f a -N {3}".format(raxml,out_folder,cores,boot_num,outName))
	elif seq_flag == "P":
		if boot_num == 1:
			os.system("{0} -s {1}/{3}_align.phylip -n multilocus_{1} -T {2} -m PROTGAMMABLOSUM62 -p 53715".format(raxml,out_folder,cores,outName))
		else:
			os.system("{0} -s {1}/{4}_align.phylip -n multilocus_{1} -T {2} -m PROTGAMMABLOSUM62 -p 53715 -x 3524627 -f a -N {3}".format(raxml,out_folder,cores,boot_num,outName))

#Print RAxML input command
elif run_flag == "O":
	if seq_flag == "N" or seq_flag == "C":
		if boot_num == 1:
			print "{0} -s {1}/{3}_align.phylip -n multilocus_{1} -T {2} -m GTRGAMMA -p 53715 ".format(raxml,out_folder,cores,outName)
		else:
			print "{0} -s {1}/{4}_align.phylip -n multilocus_{1} -T {2} -m GTRGAMMA -p 53715 -x 3524627 -f a -N {3}".format(raxml,out_folder,cores,boot_num,outName)
	elif seq_flag == "P":
		if boot_num == 1:
			print "{0} -s {1}/{3}_align.phylip -n multilocus_{1} -T {2} -m PROTGAMMABLOSUM62".format(raxml,out_folder,cores,outName)
		else:
			print "{0} -s {1}/{4}_align.phylip -n multilocus_{1} -T {2} -m PROTGAMMABLOSUM62 -x 3524627 -f a -N {3}".format(raxml,out_folder,cores,boot_num,outName)

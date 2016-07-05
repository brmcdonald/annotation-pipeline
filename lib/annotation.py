import re
import os
import sys
import collections as cl
from multiprocessing import Pool

lib_dir = sys.path[0] + '/lib/'
antifamPath = sys.path[0]+'/antifam/'
sys.path.append(lib_dir)
import annIO

#Simple command multithreader
def multi_run(cmd):
	os.system(cmd)


def db_run(db,PARAM,subFastas,LEN={},featCount=0):
	"""Runs sequences against a database and parses according to any provided cutoffs"""

	debug = ''
	print 'Running and processing {0}...'.format(db.name)


	#Setting number of cores to use per file, and total files to run at once
	if db.source == 'hmmer':
		coreNum = 2
		fNum = PARAM['cores'] / 2

		if db.multi == True:
			hmmTbl = '--domtblout'
		else:
			hmmTbl = '--tblout'
	elif len(subFastas) < PARAM['cores']:
		coreNum = PARAM['cores'] / len(subFastas)
		fNum =  len(subFastas)
	else:
		fNum = PARAM['cores']
		coreNum = PARAM['cores']

	#Setting parameters for appropriate algorithm
	if db.source == 'hmmer':
		if '--cut_ga' in db.param or '--cut_tc' in db.param or '--cut_nc' in db.param:
			data = (PARAM['hmmscan_path'],db.param,hmmTbl,coreNum)
			hmm_param = '{0} {1} --noali --cpu {3} {2} '.format(*data)
		else:
			data = (PARAM['hmmscan_path'],db.param,hmmTbl,PARAM['hmmer_min_eval'],coreNum)
			hmm_param = '{0} {1} -E {3} --noali --cpu {4} {2} '.format(*data)
	elif db.source == 'blastp':
		data = (PARAM['blast_path'],db.param,PARAM['blast_min_eval'],coreNum)
		blastp_param = '{0} {1} -E {2} -F f -p blastp -v 1 -m 8 -b 1 -a {3}'.format(*data)
	elif db.source == 'blastn':
		data = (PARAM['blast_path'],db.param,PARAM['blast_min_eval'],coreNum)
		blastn_param = '{0} {1} -E {2} -F f -p blastn -v 1 -m 8 -b 1 -a {3}'.format(*data)
	elif db.source == 'infernal':
		if '--cut_ga' in db.param or '--cut_tc' in db.param or '--cut_nc' in db.param:
			data = (PARAM['cmscan_path'],db.param,coreNum)
			cm_param = '{0} {1} --noali --cpu {2}'.format(*data)
		else:
			data = (PARAM['cmscan_path'],db.param,PARAM['cmscan_min_eval'],coreNum)
			cm_param = '{0} {1} -E {2} --noali --cpu {3}'.format(*data)

	os.system('mkdir -p {0}/db_data/{1}'.format(PARAM['outputFolder'],db.name))

	#If it's just a table to read in, parse it and exit <---- NOT YET IMPLEMENTED
	if db.source == 'table':
		tbl_parse(db.path)
		return

	cmd = []
	for f in subFastas:
		if os.path.getsize(f) > 0:
			name = str(f)
			name = name.split('/')[-1].split('.')[:-1]
			name = '.'.join(name)

			#HMMer command for single best hits per protein, disgard full output and keep tbl
			if db.source == 'hmmer' and db.multi == False:
				inputs = (PARAM['outputFolder'],name, db.path,f,db.name,hmm_param)
				cmd.append('{5} {0}/db_data/{1}.temp.{4}.out {2} {3} > /dev/null'.format(*inputs))

			#HMMer for multiple hits (ie domain database), need full output to get hit coordinates
			elif db.source == 'hmmer' and db.multi == True:
				inputs = (PARAM['outputFolder'],name, db.path,f,db.name,hmm_param)
				cmd.append('{5} {0}/db_data/{1}.temp.{4}.out {2} {3} > {0}/db_data/{1}.full_{4}.out'.format(*inputs))

			#BLAST commands
			elif db.source == 'blastp':
				inputs = (blastp_param,db.path,PARAM['outFolder'],f,name,db.name)
				cmd.append('{0} {1} -d {1} -i {3} -o {2}/db_data/{4}.temp_{5}.out'.format(*inputs))
			elif db.source == 'blastn':
				inputs = (blastn_param,db.path,PARAM['outFolder'],f,name,db.name)
				cmd.append('{0} {1} -d {1} -i {3} -o {2}/db_data/{4}.temp_{5}.out'.format(*inputs))

			#Infernal
			elif db.source == 'infernal':
				inputs = (cm_param,PARAM['outputFolder'],db.name,name,f, db.path)
				cmd.append('{0} --tblout {1}/db_data/{3}.temp.{2}.out {5}  {4} > /dev/null'.format(*inputs))

			#LAST command
			elif db.source == 'last':
				inputs = (last,db.path,out_f,name,db.name,name)
				cmd.append('{0} -f 0 -o {2}/db_data/{3}.temp.{4}.out {1} {5}'.format(*inputs))
			else:
				print "Don't know what to do with database {0}, unrecongized source program {1}".format(db.name,db.source)

	if fNum > 1:
		p=Pool(fNum)
		for i in cmd:
			p.apply_async(multi_run, args=(i,))
			pass #for debugging
		p.close()
		p.join()
	else:
		for c in cmd:
			os.system(c)

	#Output files
	outs = os.listdir('{0}/db_data'.format(PARAM['outputFolder']))

	PROT = {}
	for f in outs:
		name = f.replace('.temp.{0}.out'.format(db.name),'')
		hmmerData = [PARAM,LEN,'{0}/db_data/{1}'.format(PARAM['outputFolder'],f),db,name,'prot']

		if '.temp.{0}.out'.format(db.name) in f:
			if db.source == 'hmmer' and db.multi == False:
				PROT[name] = hmmer_best_hit(*hmmerData)
			elif db.source == 'blastp':
				PROT[name] = blast_best_hit(*inputData.append('blastp'))
			elif db.source == 'last':
				PROT[name] = last_best_hit(*inputData)
			elif db.source == 'infernal':
				inputData = (PARAM,'{0}/db_data/{1}'.format(PARAM['outputFolder'],f),db,featCount)
				PROT[name] = infernal_currator(*inputData)

		elif '.full_{0}.out'.format(db.name) in f and '{0}.clean'.format(db.name) not in f:
			if db.source == 'hmmer' and db.multi == True:
				outFile = '{0}/db_data/{1}'.format(PARAM['outputFolder'],f)
				PROT[name] = hmmer_multi_hit(PARAM,LEN,outFile, db)

	#os.system('rm -f {0}/db_data/seq_*.temp.{1}*'.format(PARAM['outputFolder'],db.name))
	#os.system('rm -f {0}/db_data/seq_*.{1}*'.format(PARAM['outputFolder'],db.name))
	#os.system('rm -f {0}/db_data/seq_*.full_{1}*'.format(PARAM['outputFolder'],db.name))
	#os.system('rm -f {0}/temp.faa'.format(PARAM['outputFolder']))
	#os.system('rm -f -r {0}/tempFastas'.format(PARAM['outputFolder']))

	if PROT:
		return (PROT,0)
	else:
		return (inputData,1)

def antifam_best_hit(f):
	"""Simple parser for antifam hits"""
	PROT = {}
	with open(f, 'r') as f_data:
		for ln in f_data:
			if '#' not in ln:
				data = ln.split()
				cat = data[0]
				gene = data[2]
				score = float(data[8])
				evalue = data[7]

				if not PROT.has_key(gene):
					PROT[gene] = annIO.Protein()
					PROT[gene].ann["antiFam"] = ""
					PROT[gene].score["antiFam"] = 0
					PROT[gene].putative["antiFam"] = False
					PROT[gene].eval["antiFam"] = 1
				if score > PROT[gene].score["antiFam"]:
					PROT[gene].ann["antiFam"] = cat
					PROT[gene].score["antiFam"] = score
					PROT[gene].eval["antiFam"] = float(evalue)
					PROT[gene].orphan = 0
	return PROT


#Parsers for protein sequence comparison programms
def hmmer_best_hit(PARAM,LEN,f, db,name,seqType):
	"""Identifies best single HMMer hit for each gene in a tbl formatted output"""

	PROT = {}
	with open(f, 'r') as f_data:
		for ln in f_data:
			if '#' not in ln:
				data = ln.split()
				cat = data[0]
				gene = data[2]
				score = float(data[8])
				evalue = data[7]
				hit = ''
				prot_len = LEN[gene] / 3

				#If we have a noise cutoff, get rid of things that are noise
				if db.NOISE_SCORE.has_key(cat):
					if score < db.NOISE_SCORE[cat]:
						hit = 'noise'

					#If we have a trust cutoff, we can say putative or trusted for things above
					#the noise
					elif db.TRUST_SCORE.has_key(cat):
						if score < db.TRUST_SCORE[cat]:
							hit = 'putative'
						else:
							hit = 'trust'
					#If we only have a noise cutoff and it's above the noise, it's good
					#as far as we can tell
					else:
						hit = 'trust'

				#If we have no score cutoff info, just run with it
				else:
					hit = 'trust'

				#Check length requirments next
				if db.MIN_LEN.has_key(cat):
					if prot_len < db.MIN_LEN[cat]:
						hit = 'noise'
				if db.MAX_LEN.has_key(cat):
					if prot_len > db.MAX_LEN[cat]:
						hit = 'noise'

				#Assigning hit to protein
				if hit != 'noise':
					if not PROT.has_key(gene):
						PROT[gene] = annIO.Protein()
						PROT[gene].ann[db.name] = ""
						PROT[gene].score[db.name] = 0
						PROT[gene].putative[db.name] = False
						PROT[gene].eval[db.name] = 1
					if score > PROT[gene].score[db.name]:
						PROT[gene].ann[db.name] = cat
						PROT[gene].score[db.name] = score
						PROT[gene].eval[db.name] = float(evalue)
						PROT[gene].orphan = 0
						if hit == 'putative':
							PROT[gene].putative[db.name] = True
						else:
							PROT[gene].putative[db.name] = False
	return PROT

#
#TO DO: Modify this to use the tab delimited HMMer output, need to get around grouping hits by gene
#
def hmmer_multi_hit(PARAM,LEN,f, db):
	"""Identifies best non-overlapping HMMer hits, needs standard hmmer output format"""

	#Temporary container for sorting out multiple hits to the same protein
	class Hit(object):
		__slots__ = ('model','prot','length','start','end','score','type','eval')
		def __init__(self,name,query):
			self.model = name
			self.prot = query
			self.length = 0
			self.start = 0
			self.end = 0
			self.score = 0.0
			self.type = ''
			self.eval = 0.0

	#Path for the tab-delimited output file
	tabF = f.replace('full_{0}'.format(db.name),'temp.{0}'.format(db.name))
	PROT = {}	#Final protein annotations to return


	queryRE = re.compile('Query:\s+(\S+)\s+\[L=(\d+)\]')

	with open(f, 'r') as f_data:
		for ln in f_data:
			#Hits for each gene have to be compiled as we go
			if 'Query:' in ln:
				model = ''
				m = queryRE.search(ln)
				if m:
					q = m.group(1)

					genome = q.split('|')[0]
					HMM = {}		#Initialize a dict for all the hits to this gene
					count = 0
				else:
					print ln
			elif '>>' in ln:
				model = ln.split()[1].replace('.hmm','')

			#Found a decent hit
			elif '!' in ln and model:
				ln = ln.replace('\n', '')
				data = ln.split()

				count += 1
				evalue = float(data[5])
				if evalue <= PARAM["hmmer_min_eval"]:
					x = Hit(model,q)
					x.start = int(data[9])
					x.end = int(data[10])
					x.length = x.end - x.start
					x.score = float(data[2])
					x.eval = evalue
					HMM[count] = x


			#All hits for a gene compiled, now we can sort through them
			elif 'Internal pipeline statistics summary:' in ln and HMM:
				model = ''

				#Now check remaining hits against any score/length cutoffs, sort by start coord
				#(see comments in HMM best hit parser)
				for h in sorted(HMM.items(), key=lambda x: x[1].start):
					i = h[0]

					if db.NOISE_SCORE.has_key(HMM[i].model):
						if HMM[i].score < db.NOISE_SCORE[HMM[i].model]:
							HMM[i].type = 'noise_score'
						elif db.TRUST_SCORE.has_key(HMM[i].model):
							if HMM[i].score < db.TRUST_SCORE[HMM[i].model]:
								HMM[i].type = 'putative'
							else:
								HMM[i].type = 'trust'
						else:
							HMM[i].type = 'trust'
					else:
						HMM[i].type = 'trust'

					#Check length requirments next
					if db.MIN_LEN.has_key(HMM[i].model):
						if HMM[i].length < db.MIN_LEN[HMM[i].model]:
							HMM[i].type = 'noise_short'
					if db.MAX_LEN.has_key(HMM[i].model):
						if HMM[i].length > db.MAX_LEN[HMM[i].model]:
							HMM[i].type = 'noise_long'

				#Remove overlapping hits
				rm_lst = cl.deque()
				for i in HMM:
					if 'noise' not in HMM[i].type:
						for j in HMM:
							if 'noise' not in HMM[j].type:
								if i != j:
									iCoords = set(range(HMM[i].start, HMM[i].end))
									jCoords = set(range(HMM[j].start, HMM[j].end))

									if len(iCoords.intersection(jCoords)) > 5:

										#If two good hits overlap, pick the better scoring one
										if HMM[i].score < HMM[j].score:
											HMM[i].type = 'noise_overlap'
										elif HMM[i].score > HMM[j].score:
											HMM[j].type = 'noise_overlap'
										elif HMM[i].score == HMM[j].score:
											model2 = HMM[i].model.replace('.hmm','')
											HMM[j].type = 'dual_{0}'.format(model2)
											HMM[i].type = 'noise_overlap'

				#Final assignment for the hit
				for h in sorted(HMM.items(), key=lambda x: x[1].start):
					i = h[0]

					if not 'noise' in HMM[i].type:
						if db.type == 'cds':
							if not PROT.has_key(q):
								PROT[q] = annIO.Protein()
								PROT[q].coords[db.name] = {}
								PROT[q].coords[db.name]['starts'] = []
								PROT[q].coords[db.name]['ends'] = []
								PROT[q].ann[db.name] = []
								PROT[q].score[db.name] = []
								PROT[q].eval[db.name] = []


							model = HMM[i].model.replace('.hmm','')
							if 'dual' in HMM[i].type:
								second_match = HMM[i].type.split('_')[1]
								model = model + '/' + second_match
							PROT[q].ann[db.name].append(model)
							PROT[q].score[db.name].append(str(HMM[i].score))
							PROT[q].coords[db.name]['starts'].append(str(HMM[i].start))
							PROT[q].coords[db.name]['ends'].append(str(HMM[i].end))
							PROT[q].eval[db.name].append(HMM[i].eval)
							PROT[q].orphan = 0
							if HMM[i].type == 'putative':
								PROT[q].putative[db.name] = True
							else:
								PROT[q].putative[db.name] = False
	return PROT

#
#TO DO: Get blast and last parsers working
#
def blast_best_hit(PARAM,LEN,f, db, name,blast_type):
	"""Identifies best single BLAST hit for each gene in a tbl formatted output"""

	LN = {}
	with open('{0}/db_data/{1}'.format(PARAM['outputFolder'],f), 'r') as f_data:
		for ln in f_data:
			raw = ln
			ln = ln.split('\t')
			if not PROT.has_key(ln[0]):
				print 'Cannot find gene feature for {0} in blast'.format(ln[0])
				pass

			#Length of protein in AAs
			prot_len = PROT[ln[0]].length / 3
			#Length of hit in AAs
			perc_hit = (float(ln[3]) / prot_len) * 100
			perc_id = float(ln[2])

			if perc_id >= PARAM['blast_min_id'] and perc_hit > PARAM['blast_min_cov']:
				if float(ln[11]) > PROT[ln[0]].score[db.name]:
					PROT[gene].ann[db.name] = ln[1]
					PROT[gene].score[db.name] = score
					PROT[gene].orphan = 0

					if clean_flag == 0:
						LN[ln[0]] = raw

	if clean_flag == 0:
		out = open('{0}/db_data/{1}.{2}.clean'.format(PARAM['outputFolder'],name,db.name),'w')
		for gene in LN:
			out.write(LN[gene])
		out.close()
	return PROT

def last_best_hit(PARAM,LEN,f, db, name,blast_type):
	"""Identifies best single LAST hit for each gene in a tbl formatted output"""

	LN = {}
	with open('{0}/db_data/{1}'.format(PARAM['outputFolder'],f), 'r') as f_data:
		for ln in f_data:
			raw = ln
			ln = ln.split('\t')
			if not PROT.has_key(ln[0]):
				print 'Cannot find gene feature for {0} in blast'.format(ln[0])
				pass

			#Length of protein in AAs
			prot_len = PROT[ln[0]].len / 3
			#Length of hit in AAs
			perc_hit = (float(ln[3]) / prot_len) * 100
			perc_id = float(ln[2])

			if perc_id >= PARAM['blast_min_id'] and perc_hit > PARAM['blast_min_cov']:
				if float(ln[11]) > PROT[ln[0]].score[db.name]:
					PROT[gene].ann[db.name] = ln[1]
					PROT[gene].score[db.name] = score
					PROT[gene].orphan = 0

					if clean_flag == 0:
						LN[ln[0]] = raw

	if clean_flag == 0:
		out = open('{0}/db_data/{1}.{2}.clean'.format(PARAM['outputFolder'],name,db.name),'w')
		for gene in LN:
			out.write(LN[gene])
		out.close()
	return PROT


def infernal_currator(PARAM,f,db,featCount):
	"""Identifies best non-overlapping Infernal hits"""

	#Temporary container for sorting out multiple hits to the same region
	class Hit(object):
		__slots__ = ('model','length','start','end','score','type','eval','contig','strand','id')
		def __init__(self,name,contig,idTag):
			self.model = name
			self.id = idTag
			self.length = 0
			self.start = 0
			self.end = 0
			self.score = 0.0
			self.eval = 0.0
			self.contig = contig.split('_')[-1]
			self.strand = ""
			self.type = ''

	def currate_hits(HMM,db,featCount):
		"""Currate RNA hits on a contig, removing poor hits and overlaps"""
		model = ''

		#Now check remaining hits against any score/length cutoffs, sort by start coord
		#(see comments in HMM best hit parser)
		for h in sorted(HMM.items(), key=lambda x: x[1].start):
			i = h[0]

			if db.NOISE_SCORE.has_key(HMM[i].model):
				if HMM[i].score < db.NOISE_SCORE[HMM[i].model]:
					HMM[i].type = 'noise_score'
				elif db.TRUST_SCORE.has_key(HMM[i].model):
					if HMM[i].score < db.TRUST_SCORE[HMM[i].model]:
						HMM[i].type = 'putative'
					else:
						HMM[i].type = 'trust'
				else:
					HMM[i].type = 'trust'
			else:
				HMM[i].type = 'trust'

			#Check length requirments next
			if db.MIN_LEN.has_key(HMM[i].model):
				if HMM[i].length < db.MIN_LEN[HMM[i].model]:
					HMM[i].type = 'noise_short'
			if db.MAX_LEN.has_key(HMM[i].model):
				if HMM[i].length > db.MAX_LEN[HMM[i].model]:
					HMM[i].type = 'noise_long'

		#Remove overlapping hits
		rm_lst = cl.deque()
		for i in HMM:
			if 'noise' not in HMM[i].type:
				for j in HMM:
					if 'noise' not in HMM[j].type:
						if i != j:
							iCoords = set(range(HMM[i].start, HMM[i].end))
							jCoords = set(range(HMM[j].start, HMM[j].end))

							if len(iCoords.intersection(jCoords)) > 5:
								#If two good hits overlap, pick the better scoring one
								if HMM[i].score < HMM[j].score:
									HMM[i].type = 'noise_overlap'
								elif HMM[i].score > HMM[j].score:
									HMM[j].type = 'noise_overlap'
								elif HMM[i].score == HMM[j].score:
									model2 = HMM[i].model.replace('.hmm','')
									HMM[j].model = 'dual_{0}'.format(model2)
									HMM[i].type = 'noise_overlap'

		#Final assignment for the hit
		RNA = {}
		for h in sorted(HMM.items(), key=lambda x: x[1].start):
			i = h[0]
			if not 'noise' in HMM[i].type:
				model = HMM[i].model.replace('.hmm','')

				if 'LSU_rRNA_bacteria' in model or "LSU_rRNA_archaea" in model:
					model = 'rRNA-23s'
				elif 'SSU_rRNA_bacteria' in model or "SSU_rRNA_archaea" in model:
					model = 'rRNA-16s'
				elif '5S_' in model:
					model = 'rRNA-5s'
				elif ' 5_8S_rRNA' in model:
					model = 'rRNA-5.8s'
				elif 'LSU_rRNA_eukarya' in model:
					model = 'rRNA-28s'
				elif "SSU_rRNA_eukarya" in model:
					model = 'rRNA-18s'
				elif 'tRNA' in model:
					model = 'tRNA'
				else:
					model = ln[0]

				RNA[HMM[i].id] = annIO.Feature('rna')

				if 'dual' in HMM[i].type:
					second_match = HMM[i].type.split('_')[1]
					model = model + '/' + second_match
				RNA[HMM[i].id].ann = model
				RNA[HMM[i].id].score = float(HMM[i].score)
				RNA[HMM[i].id].start = HMM[i].start
				RNA[HMM[i].id].end = HMM[i].end
				RNA[HMM[i].id].eval = float(HMM[i].eval)
				RNA[HMM[i].id].contig = HMM[i].contig

				if HMM[i].type == 'putative':
					RNA[HMM[i].id].putative = True
				else:
					RNA[HMM[i].id].putative = False
		return RNA

	RNA_FINAL = {}
	with open(f, 'r') as f_data:
		HMM = {} 			#Hits for each contig have to be compiled as we go
		count = 0
		contig = ''
		for ln in f_data:
			if "#" not in ln:
				ln = ln.replace('\n', '')
				data = ln.split()

				#Finished looking at all the RNAs from a single contig, currate them now
				if contig != data[2]:
					if HMM:
						NEW = currate_hits(HMM,db,featCount)
						for i in NEW:
							RNA_FINAL[i] = NEW[i]
						HMM = {}

					#Keep going through hits to new contig
					contig = ln[2]
					count = 0

				count += 1
				evalue = float(data[15])
				if evalue <= PARAM["hmmer_min_eval"]:
					model = data[0]
					contig = data[2]
					idTag = contig + '.' + str(count)
					x = Hit(model,contig,idTag)
					x.start = int(data[7])
					x.end = int(data[8])
					x.length = x.end - x.start
					x.score = float(data[14])
					x.eval = evalue
					if data[9] == "+":
						x.strand = "1"
					else:
						x.strand = "-1"
					HMM[count] = x
		if HMM:
			NEW = currate_hits(HMM,db,featCount)
			for i in NEW:
				RNA_FINAL[i] = NEW[i]
	return RNA_FINAL

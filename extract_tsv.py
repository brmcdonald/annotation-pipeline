import sys
import os
import json
import collections as cl
from multiprocessing import Pool



####################################################################################################
usage = """\n

Generates tsv-formatted files from feature JSON files.

<json_feature folder>
<file extension for output files (not including the period)>
<output folder>
<number of cores to use>\n\n"""
argv = sys.argv[1:]

if len(argv) == 4:
	inFolder= argv.pop(0)
	fileExt = argv.pop(0)
	outFolder = argv.pop(0)
	cores = int(argv.pop(0))
else:
	sys.exit(usage)

os.system('mkdir -p {0}'.format(outFolder))
jsons = os.listdir(inFolder)

#Getting DB data from first JSON file and setting up header
ann_header =[]
filler = []
db_len = -1
empty_db = '\t-\t0.0'
with open('{0}/{1}'.format(inFolder,jsons[0]),'r') as f:
	data = json.load(f,object_pairs_hook=cl.OrderedDict)
	for g in data:
		if data[g]['type'] == 'cds' and db_len == -1:
			db_len = len(data[g]['ann'])
			for db in data[g]['ann']:
				ann_header.append(db)
				ann_header.append('{0}_bitscore'.format(db))
				filler.append(empty_db)
			filler = ''.join(filler)
ann_header = '\t'.join(ann_header)
header = 'LocusID\tAnnotation\tStart\tEnd\tPartial\tAbsolute_Start\tStrand\tContig\tType\t{0}\n'.format(ann_header)


def extractTSV(f,header,filler,inFolder,outFolder,fileExt):
	with open('{0}/{1}'.format(inFolder,f),'r') as inF:
		name = f.replace('json',fileExt)
		data = json.load(inF,object_pairs_hook=cl.OrderedDict)
		out = open('{0}/{1}'.format(outFolder,name),'w')
		out.write(header)
		for g in data:
			if data[g]['type'] == 'cds':

				#Compiling data from all database hits
				db_ann = []
				for db in data[g]['ann']:
					if data[g]['ann'][db] == '-':
						db_ann.append('-')
						db_ann.append('0.0')
					else:
						if type(data[g]['ann'][db]) == list:
							score = ','.join(data[g]['score'][db])
							dom = ','.join(data[g]['ann'][db])
						else:
							score = str(data[g]['score'][db])
							dom = data[g]['ann'][db]
						db_ann.append(dom)
						db_ann.append(score)
				db_ann = '\t'.join(db_ann)

				gene_data = (data[g]["id"],data[g]["function"],data[g]["start"],data[g]["end"],data[g]["partial"],
					data[g]["absStart"],data[g]["strand"],data[g]["contig"],data[g]["type"],db_ann)
				gene_data = [str(i) for i in gene_data]
				out.write('\t'.join(gene_data) + '\n')
	
			else:
				gene_data = (data[g]["id"],data[g]["function"],data[g]["start"],data[g]["end"],data[g]["partial"],
					data[g]["absStart"],data[g]["strand"],data[g]["contig"],data[g]["type"],filler)
				gene_data = [str(i) for i in gene_data]
				out.write('\t'.join(gene_data) + '\n')
		out.close()

if cores > 1:
	p=Pool(cores)
	for jFile in jsons:
		p.apply_async(extractTSV, args=(jFile,header,filler,inFolder,outFolder,fileExt))
		pass #for debugging
	p.close()
	p.join()
else:
	for jFile in jsons:
		os.system(jFile)


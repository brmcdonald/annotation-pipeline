import os
import sys
import json
import collections as cl

usage = """\n

Converts a tsv input table to JSON format for the pipeline.

<Input folder>
<tsv input table>
<output file name>\n\n"""
argv = sys.argv[1:]

if len(argv) == 3:
	inFolder = os.path.abspath(argv.pop(0))
	tsvInput = argv.pop(0)
	outFile = argv.pop(0)
else:
	sys.exit(usage)

def tbl2dictlist(f,):
	"""Converts a tab-delimited table into a dictionary: dict[col1] = [other cols]"""
	DATA = cl.OrderedDict()
	with open(f,'r') as info:
		for ln in info:
			ln = ln.replace('\n','')
			if ln:
				ln = ln.split('\t')
				if ln[0]:
					name = ln.pop(0)
					DATA[name] = []
					if len(ln) > 0:
						for i in ln:
							DATA[name].append(i)
	return DATA

DATA = tbl2dictlist(tsvInput)
GENOME = {}

for g in DATA:
	if g and '#' not in g:
		GENOME[DATA[g][0]] = {}
		for idx, i in enumerate(DATA[g]):
			if i == '-':
				DATA[g][idx] = ""

		GENOME[DATA[g][0]]["Name"] = DATA[g][0]
		GENOME[DATA[g][0]]["chromosomeFastas"] = DATA[g][1].split(',')
		GENOME[DATA[g][0]]["chromosomeNames"] = DATA[g][2].split(',')
		GENOME[DATA[g][0]]["plasmidFastas"] = DATA[g][3].split(',')
		GENOME[DATA[g][0]]["plasmidNames"] = DATA[g][4].split(',')
		GENOME[DATA[g][0]]["coordFiles"] = DATA[g][5].split(',')
		GENOME[DATA[g][0]]["genbankFiles"] = DATA[g][6].split(",")
		GENOME[DATA[g][0]]["inputType"] = DATA[g][7]
		GENOME[DATA[g][0]]["orgType"] = DATA[g][8]
		GENOME[DATA[g][0]]["locusTag"] = g
		GENOME[DATA[g][0]]["inputFolder"] = inFolder

out = open(outFile,'w')
json.dump(GENOME,out,indent=4)
out.close()
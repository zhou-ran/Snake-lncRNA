#!/usr/bin/env python
import argparse 
import re
parser = argparse.ArgumentParser(description="Filter the isoforms by the class_code(u,i,x), exon numbers(1) and length(200).\nUse this script after cufflinks/stringtie as a subsitution for next analysis\n\nAuthor: Ran Zhou	Mail: Ranzhou1005@gmail.com",formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument("-i","--input", action="store", dest="gtf",type=str, required=True, 
                    help="the prefix of the gffcompare output files")
parser.add_argument("-o","--out", action="store",type = str, dest="out", default='filter.gtf',\
                    help="the prefix of the output file")
parser.add_argument("-s","--single", action="store_true", dest="single", default=False,\
                     help="plant: yes\tanimal: No \tdefault:plant")
argvs = parser.parse_args()

pool = ['u','i','x']
tmp = []
pattern_name = re.compile(r'transcript_id "(.*?)"')
pattern_cc   = re.compile(r'class_code "(.*?)"')
with open(argvs.gtf+'.tracking') as IN, open(argvs.gtf+'.annotated.gtf') as IN1, open(argvs.out + '.tracking','w') as OUT, open(argvs.out+'.gtf','w') as OUT1:
	OUT.write('\t'.join(['tracking-iso','tracking-gene','gene','transcrip','exon-num','FPKM','TPM','COV','length']) + '\n')
	for i in IN.readlines():
		ALL = i.strip().split('\t')
		if ALL[3] in pool:
			info = re.split(':|\|',ALL[-1])
			#query 0, gene ID 1, transcrip ID 2, exon num 3, FPKM 4, TPM 5, cov 6, length 7
			if int(info[3]) == 1:
				if int(info[7]) >= 200 and float(info[4]) >=2:
					tmp.append(info[2])
					OUT.write('\t'.join(ALL[0:2])+'\t'+'\t'.join(info[1:])+'\n')
				else:pass
			elif int(info[3]) >1:
				if int(info[7]) >= 200 and float(info[4]) >= 1:
					tmp.append(info[2])
					OUT.write('\t'.join(ALL[0:2])+'\t'+'\t'.join(info[1:])+'\n')
				else:pass
			else:
				pass
		else:pass
	for i in IN1.readlines():
		ALL = i.strip().split('\t')
		name = re.findall(pattern_name,i)[0]
		if ALL[2] == 'transcript':
			cc = re.findall(pattern_cc,i)[0]
			if cc in pool:
				if tmp:
					if name in tmp:
						tmp.remove(name)
						OUT1.write(i)
						PASS = True
					else:
						PASS = False
				else:
					break
			else:pass
		else:
			if PASS:
				OUT1.write(i)
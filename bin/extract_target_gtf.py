import argparse 
import re
parser = argparse.ArgumentParser(description="Extract the target gene information from gtf files.\nUse this script after cufflinks/stringtie as a subsitution for next analysis\n\nAuthor: Ran Zhou	Mail: Ranzhou1005@gmail.com",formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument("-t","--target", action="store", dest="gene",type=str, required=True, 
                    help="target genes file, one line for on gene")
parser.add_argument("-g","--gtf", action="store", dest="gtf",type=str, required=True, 
                    help="GTF file")
parser.add_argument("-o","--out", action="store",type = str, dest="out", default='target.gtf',\
                    help="the name of the output file")
argvs = parser.parse_args()
ID = re.compile(r'transcript_id "(.*?)"')

with open(argvs.gene) as IN1, open(argvs.gtf) as IN2, open(argvs.out,'w') as OUT:
	pool = IN1.read().splitlines()
	for i in IN2.readlines():
		if i.startswith('#'):continue
		ALL = i.strip().split('\t')
		if ALL[2] == 'transcript' and ALL[6] != '.':
			iso_ID = re.findall(ID,ALL[-1])[0]
			if pool:
				if iso_ID in pool:
					pool.remove(iso_ID)
					OUT.write('\t'.join(ALL) + '\n')
					cid = True
				else:
					cid = False
			else:break
		else:
			if cid:
				OUT.write('\t'.join(ALL) + '\n')

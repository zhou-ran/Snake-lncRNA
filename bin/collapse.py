import sys
import glob
from itertools import islice
import argparse
import matplotlib
import re
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib_venn import venn3, venn3_circles

parser = argparse.ArgumentParser(description="Merge all potential non-coding transcripts after CPC, CNCI and Pfam.\n\nAuthor: Ran Zhou	Mail: Ranzhou1005@gmail.com",formatter_class= argparse.RawTextHelpFormatter)
parser.add_argument("-i","--input", action="store", dest="file",type=str, required=True, 
                    help="the prefix of the predition output dir")
parser.add_argument("-o","--out", action="store",type = str, dest="out", default='output',\
                    help="the dir of the output file")
parser.add_argument("-f","--fasta", action="store",type = str, dest="fa", required=True,\
                    help="The LncRNA candiate fasta")
parser.add_argument("-g","--gtf", action="store",type = str, dest="gtf", required=True,\
                    help="The LncRNA candiate gtf")

ID = re.compile(r'transcript_id "(.*?)"')
argvs = parser.parse_args()
# cpc, cnci, pfam, fa= sys.argv[1:]
CPC, Pfam, CNCI  = [], [], []

# allfile = glob.glob('data/Lnc_predict/split/*')
allfile = glob.glob(argvs.file+'/*')

def filter(l,line,sline):
	"""
	filter(file, skip line, coding or noncoding)
	"""
	with open(l) as IN:
		non_c, coding = [],[]
		for i in islice(IN,line,None):
			ALL = i.strip().split('\t')
			if ALL[sline] == 'coding':
				coding.append(ALL[0])
			else:
				non_c.append(ALL[0])
	return non_c, coding

def filter_CPC(l,line,sline,orf,info):
	ID = re.compile(r'>(.*?) ')
	SE = re.compile(r'framefinder \((.*?),(.*?)\)')
	with open(l) as IN, open(orf) as IN2:
		non_c, coding = [],[]
		for i in islice(IN,line,None):
			ALL = i.strip().split('\t')
			if ALL[sline] == 'coding':
				coding.append(ALL[0])
			else:
				non_c.append(ALL[0])
		
		for i in IN2.readlines():
			if i.startswith('>'):
				seq_ID = re.findall(ID, i)[0]
				real_ID = seq_ID.split('|')[0]
				# print(real_ID)
				strand = seq_ID[-1]
				tmp = re.findall(SE, i)[0]
				Start, End = int(tmp[0]), int(tmp[1])

				if strand == '+':
					if End == int(LEN_d[real_ID]) and (End - Start)+1 >= 150:
						if seq_ID in non_c:
						# print(seq_ID)
							non_c.remove(seq_ID)
						else:
							pass
					elif (End - Start) >= 600:
						if seq_ID in non_c:
						# print(seq_ID)
							non_c.remove(seq_ID)
						else:pass
				else:
					if Start == 1 and (End - Start) + 1 >= 150:
						if seq_ID in non_c:
							# print(seq_ID)
							non_c.remove(seq_ID)
						else:
							pass
					elif (End - Start) >= 600:
						if seq_ID in non_c:
						# print(seq_ID)
							non_c.remove(seq_ID)
						else:pass
	return non_c, coding




def pfam_filter(l, fa):
	ALLname = []
	coding = []
	with open(l) as IN1, open(fa) as IN2:
		for i in IN2.readlines():
			if i.startswith('>'):
				ALLname.append(i.strip()[1:])
			else:pass
		for i in IN1.readlines():
			if i.startswith('#') or i == '\n':continue
			tmp = i.strip().split()[0]
			name = '.'.join(tmp.split('.')[:-1])
			if name not in coding:
				coding.append('.'.join(tmp.split('.')[:-1]))
			else:pass
		non_c = [x for x in ALLname if x not in coding]
	return non_c
LEN_d = {}
with open(argvs.out + '/final.info') as IN1:
	for i in islice(IN1,1,None):
		ALL = i.strip().split('\t')
		LEN_d[ALL[0]] = ALL[11]

for i in allfile:
	"""
	CPC: ['MSTRG.4160.1|chr01-','1855','noncoding','-1.20735']
	no header
	"""
	CPC  += filter_CPC(i+'/CPC/CPC.LncRNA.txt',0,2,i+'/CPC/ff.fa',LEN_d)[0]
	"""
	CNCI: ['MSTRG.4160.1|chr01-','noncoding','-0.086016','699','798','1855']
	one header
	"""
	CNCI += filter(i+'/CNCI/CNCI.index', 1, 1)[0]
	"""
	except these transcripts emerged in the pfam results
	"""
	Pfam += pfam_filter(i+'/pfam/Pfam_result.txt', i+'/split.fasta')

CPC_s  = set(CPC)
CNCI_s = set(CNCI)
Pfam_s = set(Pfam)


"""
plot a venn diagram
"""
plt.figure(figsize=(4,4))
v = venn3([CPC_s,CNCI_s,Pfam_s],("CPC",'CNCI','Pfam'))
v.get_patch_by_id('100').set_alpha(1.0)
c = venn3_circles([CPC_s,CNCI_s,Pfam_s], linestyle = 'dashed')
plt.title('Coding Potential Prediction')
plt.savefig(argvs.out+'/candiate.pdf')
common_set = CPC_s&CNCI_s&Pfam_s

output = list(common_set)
with open(argvs.fa) as FA, open(argvs.gtf) as GTF, open(argvs.out + '/candiate.fasta','w') as OUT, open(argvs.out + '/candiate.gtf','w') as OUT1, open(argvs.out + '/candiate.list','w') as OUT2:
	OUT2.write('qry_id\tChr'+'\n')
	for i in output:
		OUT2.write('\t'.join(i.split('|'))+'\n')
	"""
	generate the candiate fasta file
	"""
	for i in FA.readlines():
		if i.startswith('>'):
			cid = i.strip()[1:]
			if cid in output:
				OUT.write('>'+cid+'\n')
				PASS = True
			else:
				PASS = False
		else:
			if PASS:
				OUT.write(i)
			else:pass
	"""
	generate the candiate gtf file
	"""
	for i in GTF.readlines():
		ALL = i.strip().split('\t')
		if ALL[2] == 'transcript':
			T_ID = re.findall(ID,ALL[-1])[0]
			T_ID = T_ID +"|"+ ALL[0] + ALL[6]
			if T_ID in output:
				OUT1.write(i)
				PASS = True
			else:
				PASS = False
		else:
			if PASS:
				OUT1.write(i)
			else:
				pass

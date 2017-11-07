import sys
import re
from itertools import islice
iso = re.compile(r'transcript_id "(.*?)";')
file1, file2, fileout = sys.argv[1:]
with open(file1) as IN1, open(file2) as IN2, open(fileout+'LincRNA.gtf','w') as OUT1, open(fileout + 'intronic.gtf','w') as OUT2, open(fileout+'anti-sense.gtf','w') as OUT3:
	dic = {}
	for i in islice(IN1,1,None):
		ALL = i.strip().split('\t')
		dic[ALL[0]] = ALL[4]
	for i in IN2.readlines():
		ALL = i.strip().split('\t')
		tran = re.findall(iso,i)[0]
		if dic[tran] == 'u':
			OUT1.write(i)
		elif dic[tran] == 'i':
			OUT2.write(i)
		elif dic[tran] == 'x':
			OUT3.write(i)
		else:pass


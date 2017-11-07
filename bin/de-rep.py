import sys
import re
filein, fileout = sys.argv[1:]
ID = re.compile(r' transcript_id "(.*?)";')
with open(filein) as IN, open(fileout,'w') as OUT:
	cid = ''
	for i in IN.readlines():
		ALL = i.strip().split('\t')
		if i.startswith('#'):
			OUT.write(i)
		else:
			if ALL[2] == 'transcript':
				tran = re.findall(ID, i)
				if cid == '':
					PASS = True
					cid = tran
					OUT.write(i)
				elif cid == tran:
				# print(i)
					PASS = False
				else:
					cid = tran
					OUT.write(i)
					PASS = True
			else:
				if PASS:
					OUT.write(i)
				else:
					print(cid)
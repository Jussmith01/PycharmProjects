import re

f = open('/home/jujuman/Dropbox/Research/AceticAcidDimerProtonTransfer/aad_10pts_irc.log','r')

ty = []
ls = []
for i in range(0,21):
	ls.append([])

r = re.compile('^\s+(\d+?)((?:\s+?-?\d+?\.\d+)+)\n$')
r2 = re.compile('^\s([CHON]),')
for line in f:
	t = r2.search(line)
	if t:
		ty.append(t.group(1))

	m = r.search(line)
	if m:
		dat = m.group(2).strip().split('  ')
		for i in dat:
			#print (m.group(1))
			ls[int(m.group(1))-1].append(i.strip())	
			#print('DATA: ' + str(dat))
print(ty)

Na = int(len(ty))

o = open('data_irc.dat','w')
o.write('IRCOUT\n')
o.write(str(len(ls))+'\n')

tyl = str(Na) + ','
for i in ty:
	tyl = tyl + str(i) + ','

o.write(tyl+'\n')

count = 0
for i in ls:
	lout = ''
	for j in range(2,3*Na+2):
		lout = lout + str(i[j]) + ','
	lout = lout + i[0] + ',\n'
	print count
	count = count + 1
	o.write(lout)


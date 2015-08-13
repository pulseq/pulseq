#!/usr/bin/env python



print "Testing parser for open MRI format"
print "=================================="

from subprocess import call

def cmp_lines(path_1, path_2):
	"""Compare two files, ignoring line-endings"""
	l1 = l2 = ' '
	with open(path_1, 'r') as f1:
		with open(path_2, 'r') as f2:
			while l1 != '' and l2 != '':
				l1 = f1.readline()
				l2 = f2.readline()
				if l1 != l2:
					return False
	return True

sequences = ['fid','gre','press']
base_dir = '../examples/'
approved_dir = '../examples/approved/'
ok = True

for index, seq in enumerate(sequences):

	with open(base_dir + seq + '.out', 'w') as f:

		call(["./parsemr",base_dir + seq + '.seq'],stdout=f)

	same = cmp_lines(base_dir + seq + '.out',approved_dir + seq + '.out')

	result = "ok" if same else "not ok"
	print "Comparing output {0}: {1}".format(seq,result)

	ok = ok & same

exit(0 if ok else 1)

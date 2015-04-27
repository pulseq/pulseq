#!/usr/bin/env python



print "Testing parser for open MRI format"
print "=================================="

from subprocess import call
import filecmp

sequences = ['fid','gre','press']
base_dir = '../examples/'
approved_dir = '../examples/approved/'
ok = True

for index, seq in enumerate(sequences):

	with open(base_dir + seq + '.out', 'w') as f:

		call(["./parsemr",base_dir + seq + '.seq'],stdout=f)
	
	same = filecmp.cmp(base_dir + seq + '.out',approved_dir + seq + '.out')

	result = "ok" if same else "not ok"
	print "Comparing output {0}: {1}".format(seq,result)
	
	ok = ok & same

exit(0 if ok else 1)


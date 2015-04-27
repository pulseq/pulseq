#!/usr/bin/env python


print "Testing matlab code for open MRI format"
print "======================================="

from subprocess import call
import filecmp

sequences = ['fid','gre','press']
base_dir = '../examples/'
approved_dir = '../examples/approved/'
ok = True

cellStr = '{';
for index, seq in enumerate(sequences):
    cellStr += "'" + base_dir + seq + ".seq',"
cellStr += "}"

funCmd = "parsemr({0}); exit;".format(cellStr)
matlabCmd = 'matlab -nodisplay -nosplash -r "{0}"'.format(funCmd)
call(matlabCmd,shell=True)

for index, seq in enumerate(sequences):
    same = filecmp.cmp(base_dir + seq + '.matlab.out',approved_dir + seq + '.matlab.out')

    result = "ok" if same else "not ok"
    print "Comparing output {0}: {1}".format(seq,result)

    ok = ok & same

exit(0 if ok else 1)
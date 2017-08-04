#!/usr/bin/env python
"""
testmatlab.py
"""

from __future__ import print_function
import subprocess
import os

print("Testing matlab code for open MRI format")
print("=======================================")


def cmp_lines(path_1, path_2):
    """Compare two files, ignoring line-endings"""
    line_1 = line_2 = ' '
    with open(path_1, 'r') as file_1:
        with open(path_2, 'r') as file_2:
            while line_1 != '' and line_2 != '':
                line_1 = file_1.readline()
                line_2 = file_2.readline()
                if line_1 != line_2:
                    return False
    return True


def main():
    """Main"""
    sequences = ['fid', 'gre', 'press']
    binary_sequences = ['gre_binary']
    base_dir = '../examples/'
    approved_dir = '../examples/approved/'
    ok_flag = True

    cell_str = '{'
    for _, seq in enumerate(sequences):
        cell_str += "'" + base_dir + seq + ".seq',"
    for _, seq in enumerate(binary_sequences):
        cell_str += "'" + base_dir + seq + ".bin',"
    cell_str += "}"

    # Cleanup previous results
    for _, seq in enumerate(sequences + binary_sequences):
        try:
            os.remove(base_dir + seq + '.matlab.out')
        except:
            pass

    fun_cmd = "parsemr({0}); exit;".format(cell_str)
    matlab_cmd = 'matlab -nodisplay -nosplash -r ' +\
        '"try; {0}; catch; exit(-1); end"'.format(fun_cmd)
    subprocess.call(matlab_cmd, shell=True)

    for _, seq in enumerate(sequences + binary_sequences):
        same = cmp_lines(base_dir + seq + '.matlab.out',
                         approved_dir + seq + '.matlab.out')

        result = "ok" if same else "not ok"
        print("Comparing output {0}: {1}".format(seq, result))

        ok_flag = ok_flag & same

    exit(0 if ok_flag else 1)


if __name__ == '__main__':
    main()

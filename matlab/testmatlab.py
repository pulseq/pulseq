#!/usr/bin/env python
"""
testmatlab.py
"""

from __future__ import print_function
import subprocess
import os

print("Pulseq Test Suite")
print("=================")


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
    sequences = ['fid', 'gre', 'epi_rs']
    # binary sequences are unsuported now
    #binary_sequences = ['gre_binary']
    binary_sequences = []
    base_dir = '../examples/'
    approved_dir = '../examples/approved/'
    ok_flag = True

    cell_str = '{'
    for seq in sequences:
        cell_str += "'" + base_dir + seq + ".seq',"
    for seq in binary_sequences:
        cell_str += "'" + base_dir + seq + ".bin',"
    cell_str += "}"

    # Cleanup previous results
    for seq in (sequences + binary_sequences):
        try:
            os.remove(base_dir + seq + '.matlab.out')
        except:
            pass

    # run example sequences
    matlab_cmd = 'matlab -nodisplay -nosplash -r ' + \
        '"try; run ' + base_dir + '/run_all; catch; exit(-1); end; exit(0)"' +\
        ' > /dev/null'
    print("Exporting test sequences (this may take a while)...")
    #  print(matlab_cmd)
    subprocess.call(matlab_cmd, shell=True)

    #  fun_cmd = "parsemr({0}); exit;".format(cell_str)
    #  matlab_cmd = 'matlab -nodisplay -nosplash -r ' +\
    #      '"try; {0}; catch; exit(-1); end; exit(0);"'.format(fun_cmd)
    #  print("Executing")
    #  print(matlab_cmd)
    #  subprocess.call(matlab_cmd, shell=True)

    print("Comparing output...")
    for seq in sequences + binary_sequences:
        # compare condensed info
        same = cmp_lines(base_dir + seq + '.matlab.out',
                         approved_dir + seq + '.matlab.out')

        result = "ok" if same else "not ok"
        print("Sequence {0} (condensed): {1}".format(seq, result))

        # compare full .seq files
        same = cmp_lines(base_dir + seq + '.seq',
                         approved_dir + seq + '.seq')

        result = "ok" if same else "not ok"
        print("Sequence {0} (full seq): {1}".format(seq, result))

        ok_flag = ok_flag & same

    exit(0 if ok_flag else 1)


if __name__ == '__main__':
    main()

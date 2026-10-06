"""
window_repeat_content.py

usage: python window_repeat_content.py <windowfile.bed> <repeat.gff> <out.txt>
"""

import sys
from decimal import *

out = open(sys.argv[3], 'w')
out.write('chromosome'+'\t'+'start'+'\t'+'end'+'\t'+'total_repeats'+'\t'+'prop_repeats'+'\n')

windows = []
for line in open(sys.argv[1], 'r'):
        windows.append(line)

repeats = []
for line in open(sys.argv[2], 'r'):
    li=line.strip()
    if not li.startswith("#"):
                chrom = li.split('\t')[0]
                start = li.split('\t')[3]
                end = li.split('\t')[4]
                name = li.split('\t')[8]
                length = abs(int(end) - int(start))
                repeats.append(chrom+'_'+start+'_'+end+'_'+str(length)+'_'+name)

for window in windows:
        total = 0
        size = int(window.split()[2]) - int(window.split()[1])

        matches = [r for r in repeats if r.split('_')[0] == window.split()[0] and int(r.split('_')[1]) > int(window.split()[1]) and int(r.split('_')[2]) < int(window.split()[2])]
        for m in matches:
                length = m.split('_')[3]

                total = total + int(length)

                getcontext().prec = 11
                total_prop = float(total)/float(size)

        if total_prop > 1.0:
                total_prop = 1.0
        print window.split()[0], window.split()[1], window.split()[2], total, total_prop
        out.write(window.split()[0]+'\t'+window.split()[1]+'\t'+window.split()[2]+'\t'+str(total)+'\t'+str(total_prop)+'\n')
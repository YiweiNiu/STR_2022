#!/usr/bin/env python

import os
import re
from glob import glob

#count length
f = open('/home2/niuyw/RefData/Homo_sapiens/hg38.chrom.sizes', 'r')
length = 0
line = f.readline()
while line:
  if re.match('chr\d+\t', line):
    length = length + int(line.split('\t')[1])
  line = f.readline()
f.close()

#count
dir1 = 'snp_het_NyuWa'
dir2 = 'snp_het_1KGP'
fo1 = open('snp_het_NyuWa.tsv', 'wt')
fo2 = open('snp_het_1KGP.tsv', 'wt')


flag = 0
dic1 = {}
for file in glob('%s/*.het' %(dir1)):
  if file == 'snp_het_NyuWa/chrY.het.het':
    continue
  filename = file
  fa = open(filename, 'r')
  if flag == 0:#read name as dic
    line = fa.readline()#title
    line = fa.readline()
    while line:
      dic1[line.split('\t')[0]] = 0
      line = fa.readline()
    flag = 1
    fa.close()
    fa = open(filename, 'r')
  line = fa.readline()
  line = fa.readline()#title
  while line:
    dic1[line.split('\t')[0]] = dic1[line.split('\t')[0]] + int(line.split('\t')[3]) - int(line.split('\t')[1])
    line = fa.readline()
  fa.close()

for a in dic1:
  b = dic1[a]
  fo1.write(a + '\t' + str(b) + '\t' + str(length) + '\t' + str(b/length) + '\n')
fo1.close()


dic2 = {}
for file in glob('%s/*.het' %(dir2)):
  filename = file
  fa = open(filename, 'r')
  if flag == 1:
    line = fa.readline()#title
    line = fa.readline()
    while line:
      dic2[line.split('\t')[0]] = 0
      line = fa.readline()
    flag = 2
    fa.close()
    fa = open(filename, 'r')
  line = fa.readline()
  line = fa.readline()#title
  while line:
    dic2[line.split('\t')[0]] = dic2[line.split('\t')[0]] + int(line.split('\t')[3]) - int(line.split('\t')[1])
    line = fa.readline()
  fa.close()

for a in dic2:
  b = dic2[a]
  fo2.write(a + '\t' + str(b) + '\t' + str(length) + '\t' + str(b/length) + '\n')
fo2.close()



import os
import sys
import numpy as np
import gzip

D_period={}
D_af={}
D_ac={}

with open('hg38_ver13.bed',"r") as f:
   for line in f:
      lin=line.strip().split('\t')
      chrom,start,period = lin[0],lin[1],lin[3]
      D_period[chrom+':'+start] = period

with open('ac',"r") as f:
   for line in f:
      lin=line.strip().split(' ')
      chrom,start,allele,ac = lin
      D_ac[chrom+':'+start+':'+allele] = ac

with open('af',"r") as f:
   for line in f:
      lin=line.strip().split(' ')
      chrom,start,allele,af = lin
      D_af[chrom+':'+start+':'+allele] = af

with open('vep-merge-high.txt',"r") as f, open('lof_allele_new.txt',"w") as fout:
   for line in f:
      lin=line.strip().split('\t')
      pos,allele,tp = lin[0], lin[2], lin[6]

      pos = pos.split('_')[0] +':'+ str(int(pos.split('_')[1])-1)
      key = pos +':'+ allele

      if key in D_af:
         af= float(D_af[key])
         ac= float(D_ac[key])
         for i in [1,0.1,0.01,0.001,0.0001]:
            if af<=i:
               tag = i

         print(line.strip() +'\t'+ tp.split(',')[0] +'\t'+ D_af[key] +'\t'+ str(tag) +'\t'+ D_ac[key] +'\t'+ D_period[pos],file=fout)

import os
import sys
import numpy as np

for i in ['EAS','SAS','EUR','AFR','AMR','nyuwa']:
   with open(i+'.tab',"r") as f, open('expansion-'+i,"w") as fout:
      for line in f:
         lin=line.strip().split('\t')

         if lin[0]!='chrom':
            pos=lin[0]+':'+lin[1]
            acount = lin[5]

            d=[]
            for i in acount.split(','):
               allele=i.split(':')[0]
               count=i.split(':')[1]
               length = len(allele)
               for j in range(0,int(count)):
                  d.append(length)
            
            a=np.array(d)
            percent95 = np.percentile(a,95)
            median = np.median(a) 
            expan_score = float(percent95-median)/float(median)
            
            if expan_score>=2:
               print(pos +'\t'+  str(percent95)+'\t'+ str(median)+'\t'+ str(expan_score),file=fout)

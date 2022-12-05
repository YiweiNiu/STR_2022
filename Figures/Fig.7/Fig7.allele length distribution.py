import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib as mpl
mpl.rcParams.update({'font.size': 15, 'font.family': 'Arial', "savefig.dpi": 300})

pos_draw='chr1:2609035'

with open("total.txt", 'r') as f:
      for line in f:
         lin=line.strip().split("\t")
         if lin[0]!='chrom':
            position = lin[0]+':'+lin[1]
            if position == pos_draw:
               ref, motif, func=lin[-6], lin[13], lin[-3]

for pop in ['nyuwa','EAS','SAS','EUR','AMR','AFR']:
   d=[]
   with open(pop+".tab", 'r') as f:
            for line in f:
               lin=line.strip().split("\t")

               if lin[0]!='chrom':
                  pos = lin[0]+':'+lin[1]
                  if pos == pos_draw:
                     position, acount = lin[0]+':'+lin[1], lin[5].replace('"','')
                     for i in acount.split(','):
                        print(i)
                        for j in range(0, int(i.split(':')[1])):
                           d.append(len(i.split(':')[0]))

   #fig
   fig, ax = plt.subplots(figsize=(4,1.5))
   plt.hist(d,bins=max(d)-min(d),facecolor='#1874CD', log=True)
   plt.axvline(x=len(ref)+0.2,ls=":",c="grey",lw=1) 

   for spine in ['bottom','left']:
       ax.spines[spine].set_linewidth(1.5)
   for spine in ['top','right']:
       ax.spines[spine].set_visible(False)
   plt.tight_layout()

   plt.savefig(pos_draw.replace(':','-')+pop+'.pdf') 

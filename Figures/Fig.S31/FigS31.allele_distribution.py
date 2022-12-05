import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams.update({'font.size': 27, 'font.family': 'Arial', "savefig.dpi": 600})
from matplotlib.pyplot import MultipleLocator
import matplotlib.ticker as mtick
from matplotlib.ticker import MaxNLocator

locus='TBP'

d=[]
with open("extract-eh-process.txt","r") as f:
   for line in f:
         lin = line.strip().split('\t')
         if lin[1]!='LocusId':
            LocusId,path,gt1_filter,gt2_filter = lin[1],int(lin[-6]),int(lin[-5]),int(lin[-4])

            if LocusId==locus:
               if gt1_filter >0:
                  d.append(gt1_filter)
               if gt2_filter >0:
                  d.append(gt2_filter)

               pathfig = path

   print(locus,sorted(d))

   fig, ax = plt.subplots(figsize=(10, 5))
   plt.hist(d,bins=max(d)-min(d)+1,density=False, log=True, color='#1874CD',alpha=0.5)
   plt.xlim(0,max(d)+10)
   plt.axvline(x=pathfig,ls=":",c="red",lw=3) 
   plt.title(locus,fontsize=27)

   plt.xlabel('Repeat number')
   plt.ylabel('Count')

   plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))

   for spine in ['bottom','left']:
         ax.spines[spine].set_linewidth(1.5)
   for spine in ['top','right']:
         ax.spines[spine].set_visible(False)
   plt.tight_layout()

   plt.savefig(locus+'.pdf')


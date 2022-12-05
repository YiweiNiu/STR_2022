import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import math
import matplotlib as mpl

D={}
with open("ld.txt", 'r') as f, open("out", 'w') as fout: 
   for line in f:
      lin=line.strip().split("\t")
      if lin[-1] != "nan":
         snp_pos,str_pos,r2=lin[1],lin[4],lin[-1]

         dist = abs(int(snp_pos)-int(str_pos))
         Bin = str(int(dist)/5000) 

         if Bin not in D:
            D[Bin]=[]
         D[Bin].append(r2)

   for key in D:
      print(key, np.mean(D[key]),file=fout) 

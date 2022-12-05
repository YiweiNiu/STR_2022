
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from math import isnan

df_str = pd.read_csv('STRs/chr21.str_dosage.csv', sep=',')

df_snp = pd.read_csv('SNPs/chr21.raw', sep=' ')
# useless columns
df_snp.drop(df_snp.columns[[0, 2, 3, 4, 5]], axis=1, inplace=True)

# T
df_str = df_str.set_index('site')
df_str = df_str.T

# sort STR table by SNP table
df_str = df_str.reindex(index=df_snp['IID'])
df_snp = df_snp.set_index('IID')

# change names of SNP
x = df_snp.columns.to_list()
x = {i:'chr' + ':'.join(i.split(':')[:2]) for i in x}
df_snp.rename(columns=x, inplace=True)

# corr
tr = 'chr21:14081816'
snp = 'chr21:13847526'
x = df_str[tr]
y = df_snp[snp]

tr = 'chr21:8993159'
snp = 'chr21:13297803'
x = df_str[tr]
y = df_snp[snp]

if isinstance(y, pd.DataFrame):
    cor_lst = []
    for ix,col in y.iteritems():
        cor = x.corr(col)
        cor_lst.append(cor)
    cor = max(cor_lst)
elif isinstance(y, pd.Series):
    cor = x.corr(y)


# cmp with pearsonr
x_lst = x.to_list()
y_lst = y.to_list()
x_n = []
y_n = []
for i,j in enumerate(x_lst):
    if not isnan(j):
        x_n.append(j)
        y_n.append(y_lst[i])

cor, _ = pearsonr(x_n, y_n)





import pandas as pd
import numpy as np


def read_snp_tr_pair(fname=None, chrom=None):
    snp_str_pair = {}
    STRs = {}
    with open(fname, 'rt') as fin:
        for line in fin:
            line = line.strip('\n').split('\t')
            if line[0] != chrom:
                continue
            snp = line[0] + ":" + line[2]
            tr = line[3] + ":" + line[4]
            if snp not in snp_str_pair:
                snp_str_pair[snp] = [tr]
            else:
                snp_str_pair[snp].append(tr)
            STRs[tr] = 0
    return snp_str_pair, STRs

snp_str_pair, STRs = read_snp_tr_pair('gwas_snps.window250k_pSTRs.txt', 'chr21')


df_snp = df_snp = pd.read_csv('SNPs/chr21.raw', sep=' ')
df_snp.drop(df_snp.columns[[0, 2, 3, 4, 5]], axis=1, inplace=True)
# Load STR dosage
df_str = pd.read_csv('STRs/chr21.str_dosage.csv', sep=',')
df_str = df_str.set_index('site')
df_str = df_str.T
# sort STR table by SNP table
df_str = df_str.reindex(index=df_snp['IID'])
df_snp = df_snp.set_index('IID')
# change names of SNP
x = df_snp.columns.to_list()
x = {i:'chr' + ':'.join(i.split(':')[:2]) for i in x}
df_snp.rename(columns=x, inplace=True)


renamer = defaultdict()
for column_name in df_snp.columns[df_snp.columns.duplicated(keep=False)].tolist():
    if column_name not in renamer:
        renamer[column_name] = [column_name+'_0']
    else:
        renamer[column_name].append(column_name +'_'+str(len(renamer[column_name])))
renamer2 = deepcopy(renamer)

df_snp.rename(
    columns=lambda column_name: renamer[column_name].pop(0)
    if column_name in renamer 
    else column_name,
    inplace=True
)

>>> df_snp.columns[df_snp.columns.duplicated()]
Index(['chr21:15205839', 'chr21:21467382', 'chr21:34289445', 'chr21:35056743',
       'chr21:35477121', 'chr21:36105675', 'chr21:36155066', 'chr21:37574227',
       'chr21:38034828', 'chr21:39398611', 'chr21:45984508'],
      dtype='object')


df_cor = pd.concat([df_snp, df_str], axis=1, keys=['df_snp', 'df_str']).corr().loc['df_snp', 'df_str']
dict_cor = df_cor.to_dict("index")




for snp in dict_cor:
    for tr in dict_cor[snp]:


df1 = pd.DataFrame(np.random.rand(10, 4), columns=list('abcd'))
df2 = pd.DataFrame(np.random.rand(10, 3), columns=list('xyz'))










import numpy as np
import gzip

d=[]
with open('trio-vcf/vcf-name.txt', 'r') as f:
    for line in f:
        d.append(line.strip())

for i in d:
    sampleid=i.split('.')[0]
    filename = 'trio-vcf/'+sampleid+'.vcf.gz' 
    f=gzip.open(filename, 'r')  

    with open('stat-gt/'+sampleid+'.txt','w') as fout:
        for line in f:
            if line[0]!='#':
                lin = line.split('\t')
                FILTER=lin[6]

                mend='unkonwn'
                for i in range(9,12):
                    if lin[i].split(':')[3] in ['.','.,.']:
                        mend='lack'

                if mend=='unkonwn':
                    child_a1=int(lin[9].split(':')[3].split(',')[0])
                    child_a2=int(lin[9].split(':')[3].split(',')[1])

                    father_a1=int(lin[10].split(':')[3].split(',')[0])
                    father_a2=int(lin[10].split(':')[3].split(',')[1])

                    mother_a1=int(lin[11].split(':')[3].split(',')[0])
                    mother_a2=int(lin[11].split(':')[3].split(',')[1])

                    if child_a1==father_a1 and child_a2==mother_a1\
                    or child_a1==father_a1 and child_a2==mother_a2 \
                    or child_a1==father_a2 and child_a2==mother_a1 \
                    or child_a1==father_a2 and child_a2==mother_a2 \
                    or child_a1==mother_a1 and child_a2==father_a1 \
                    or child_a1==mother_a1 and child_a2==father_a2 \
                    or child_a1==mother_a2 and child_a2==father_a1 \
                    or child_a1==mother_a2 and child_a2==father_a2 :
                        mend='yes'
                    else:
                        mend='no'                        

                if FILTER=='.':                    
                    print>>fout, lin[0],lin[1],mend
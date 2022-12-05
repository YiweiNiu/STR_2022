#!/usr/bin/env python3

input_file = open("cat_snp_het.tsv", 'rt')
sample_file = open('/home2/niuyw/project/STR/samples/220112_sample_info.GangSTR.txt', 'rt')
output_file = open("population_hete.tsv", 'wt')

pop_dict = {}
line = sample_file.readline()
line = sample_file.readline()
while line:
    line = line.strip('\n')
    li = line.split('\t')
    pop = li[5]
    if pop != 'NA':
        if pop not in pop_dict:
            pop_dict[pop] = set()
        pop_dict[pop].add(li[1])
    line = sample_file.readline()

sum_dict = {}
count_dict = {}
average_dict = {}
for race in pop_dict:
    sum_dict[race] = 0
    count_dict[race] = 0

for line in input_file:
    line = line.strip('\n')
    li = line.split('\t')
    for race in pop_dict:
        if li[0] in pop_dict[race]:
            sum_dict[race] += float(li[3])
            count_dict[race] += 1
            break

output_file.write('population\tsnp_heterozygosity\n')
for race in pop_dict:
    if count_dict[race] != 0:
        output_file.write(race+'\t'+str(sum_dict[race]/count_dict[race])+'\n')

sample_file.close()
input_file.close()
output_file.close()


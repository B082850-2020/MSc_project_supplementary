from subprocess import call
import subprocess
import datetime
import time

import glob
import pandas as pd


# read the liftOver file with hg19 and hg38 coordiantes
df = pd.read_csv('hg19_v2/COVID19_HGI_2021.bed',sep='\t',names=['chr','start','end','hg38_pos'])
df['hg19_pos'] = df['chr'].map(str) + ':' + df['start'].map(str) + '-' + df['end'].map(str)

# read all bed.gz files
msa_mapping=glob.glob('hg19_v2/*bed.gz')
msa_mapping.sort()

# read each file as a dataframe and store into a list
bed_list=[]
for i in range(len(msa_mapping)):
    df1=pd.read_csv(msa_mapping[i],sep='\t',names=['chr','start','end','hg19_pos','score','strand'])
    species=msa_mapping[i].split('_')[2].split('.')[0]
    df1['species']=species
    #df1=df[['chr','start','end','hg19_loc','species']]
    df1=df1.drop_duplicates() #default keep the first occurences
    #print(df1.shape)
    bed_list.append(df1)

# concatenate all dataframe together
all=pd.concat(bed_list,ignore_index=True)

df['hit_freq'] = 0
for i in range(len(df)):
    # spcies mapping
    hit = all[all.hg19_pos==df.hg19_pos[i]].drop_duplicates(subset=['species']).shape[0]
    df['hit_freq'][i] = hit

# read the GWAS txt files as pandas dataframe
gwas_files = glob.glob('GWASdata/*')
gwas_files.sort()
pd_list = []
for i in range(len(gwas_files)):
    df_gwas = pd.read_csv(gwas_files[i],sep='\t')
    pd_list.append(df_gwas)

gwas_name = ['A2_ALL_eur_leave_23andme', 'A2_ALL_eur_leave_ukbb_23andme', 'A2_ALL_leave_23andme', 'A2_ALL_leave_UKBB_23andme', 'B1_ALL_eur_leave_23andme', 'B1_ALL_eur_leave_ukbb_23andme', 'B1_ALL_leave_23andme', 'B1_ALL_leave_UKBB_23andme', 'B2_ALL_eur_leave_23andme', 'B2_ALL_eur_leave_ukbb_23andme', 'B2_ALL_leave_23andme', 'B2_ALL_leave_UKBB_23andme', 'C2_ALL_eur_leave_23andme', 'C2_ALL_eur_leave_ukbb_23andme', 'C2_ALL_leave_23andme', 'C2_ALL_leave_UKBB_23andme']

gwas_list = []
for x in range(len(gwas_name)):
    print(gwas_name[x])
    df2 = pd_list[x][['#CHR','POS','SNP']]
    df2['gwas_name'] = gwas_name[x]
    df2['hg38_pos'] = df2['#CHR'].map(str) + ':' + df2['POS'].map(str)
    df2['hit_freq'] = 0
    for i in range(len(df2)):
        if df2['hg38_pos'][i] != '9:133270497' and df2['hg38_pos'][i] != '9:133270637':
            hit = df[df.hg38_pos == df2.hg38_pos[i]]['hit_freq'].values[0]
            df2['hit_freq'][i] = hit
    df2 = df2[df2['hg38_pos']!='9:133270497']
    df2 = df2[df2['hg38_pos']!='9:133270637']
    df2 = df2[['gwas_name','hg38_pos','hit_freq']]
    gwas_list.append(df2)

all_gwas = pd.concat(gwas_list,ignore_index=True)

all_gwas.to_csv('COVID_variant_count/all_variant_count.csv',sep=",",header=True,index=False)


########################################################################################################
########################################################################################################
########################################################################################################
### variant_counts_v2.py


# load data
print('Loading data...\n\n\n')
f=open('hg19_v2/COVID19_HGI_2021.bed','r')
fc=f.read().split()

# make a dictionary with hg19 coordinates as keys and hg38 coordinates as items
hg_dic={}
for i in range(len(fc)):
    if (i+1)%4==0:
        hg_dic[fc[i-2]]=fc[i]   # {hg19:hg38}

# read all the gz files
## reduce the redundent (v3 update)
zcat='zcat hg19_v2/*bed.gz |sort |uniq |awk \'{print $1,$2,$3,$4}\''   # only first four fields,leave strand and score
genome_var=subprocess.check_output(zcat, shell=True)
genome_var=genome_var.decode('UTF-8').split()

# check the frequency of the hg19 coordinates appear in gz files
var_freq={}
for i in hg_dic:
    c=True
    for j in genome_var:
        if (i in j) & (i!=j):
            c=False
            if i in var_freq:
                var_freq[i]+=1
            else:
                var_freq[i]=1
    if c:       # when the hg19 coordinate does not appear in gz files
        var_freq[i]=0

# nested dictionary
counts = {"A2_ALL_eur_leave_23andme": {} , "A2_ALL_eur_leave_ukbb_23andme": {} , "A2_ALL_leave_23andme": {} , "A2_ALL_leave_UKBB_23andme":{} , "B1_ALL_eur_leave_23andme" : {} , "B1_ALL_eur_leave_ukbb_23andme" : {} , "B1_ALL_leave_23andme" : {} , "B1_ALL_leave_UKBB_23andme" : {} , "B2_ALL_eur_leave_23andme" : {} , "B2_ALL_eur_leave_ukbb_23andme" : {} , "B2_ALL_leave_23andme" : {} , "B2_ALL_leave_UKBB_23andme" : {} , "C2_ALL_eur_leave_23andme" : {} , "C2_ALL_eur_leave_ukbb_23andme" : {} ,"C2_ALL_leave_23andme" : {} , "C2_ALL_leave_UKBB_23andme" : {}}
gwas=list(counts.keys())
gwas.sort()

# read the GWAS txt files as pandas dataframe
gwas_files=glob.glob('GWASdata/*')
gwas_files.sort()
for i in range(len(gwas)):
    df= pd.read_csv(gwas_files[i],sep='\t')
    for j in hg_dic:
        try:    # try if the hg38 coordinate appear in any of the GWAS analysis
            a=df['SNP'].str.find(hg_dic[j]).value_counts()[0]
            counts[gwas[i]][j]=var_freq[j]
        except:
            pass

# write the output file
with open("COVID_variant_output_v3.txt", 'w') as f:
    for key,value in counts.items():
        for variant,freq in value.items():
#        print(key,variant,freq,sep='\t')
            f.write('%s\t%s\t%s\n' % (key, variant,freq))

########################################################################################################
########################################################################################################
########################################################################################################
### variant_counts_v1.py

hg19=[]
hg38=[]
for i in range(len(fc)):
    if i % 2 == 0:
        hg19.append(fc[i])
    else:
        hg38.append(fc[i])
print('Data loaded.\n\n-----------------------\n')

#check hg19 bed.gz files for variant frequencies
print('Checking variant frequencies...\n\n\nThe time now is ' + str(datetime.datetime.now()))
variant_freq=[]
for i in range(len(hg19)):
    zgrep = 'zgrep \'' + hg19[i] + '\' hg19/* | wc -l'
    variant_freq.append(int(subprocess.check_output(zgrep,shell=True))-1)
print('\n\nFinish checking variant frequencies.\n\n-----------------------\n ')

#prepare a dictionary for output
print('Making a dictionary for GWAS and frequency counts\n\n-----------------------\n\n ')
counts = {"A2_ALL_eur_leave_23andme": 0 , "A2_ALL_eur_leave_ukbb_23andme": 0 , "A2_ALL_leave_23andme": 0 , "A2_ALL_leave_UKBB_23andme":0 , "B1_ALL_eur_leave_23andme" : 0 , "B1_ALL_eur_leave_ukbb_23andme" : 0 , "B1_ALL_leave_23andme" : 0 , "B1_ALL_leave_UKBB_23andme" : 0 , "B2_ALL_eur_leave_23andme" : 0 , "B2_ALL_eur_leave_ukbb_23andme" : 0 , "B2_ALL_leave_23andme" : 0 , "B2_ALL_leave_UKBB_23andme" : 0 , "C2_ALL_eur_leave_23andme" : 0 , "C2_ALL_eur_leave_ukbb_23andme" : 0 ,"C2_ALL_leave_23andme" : 0 , "C2_ALL_leave_UKBB_23andme" : 0}
gwas=list(counts.keys())
gwas.sort()

# write frequency counts into the dictionary
print('Writing variant frequencies to corresponding GWAS...\n\n\nThe time now is ' + str(datetime.datetime.now()))
for i in range(len(hg38)):
    for j in range(len(gwas)):
        grep='grep \'' + hg38[i] + '\' GWASdata/* | grep \'' + gwas[j] + '\' | wc -l'
        a=int(subprocess.check_output(grep,shell=True))
        if a!=0:
            counts[gwas[j]]= counts[gwas[j]] + variant_freq[i]
print('\n\nFinish writing variant frequencies.\nFinish time: ' + str(datetime.datetime.now()))

# write the dictionary into a file
print('Writing output...\n\n')
with open("output_new.txt", 'w') as f:
    for key, value in counts.items():
        f.write('%s:%s\n' % (key, value))
print('Output written in output.txt.')

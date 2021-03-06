## module load igmm/apps/python/3.7.3
from subprocess import call
import subprocess
import datetime
import time

import glob
import pandas as pd

# read all bed.gz files
msa_mapping=glob.glob('hg19_v2/*bed.gz')
msa_mapping.sort()

# read each file as a dataframe and store into a list
pd_list=[]
for i in range(len(msa_mapping)):
    df=pd.read_csv(msa_mapping[i],sep='\t',names=['chr','start','end','hg19_loc','score','strand'])
    species=msa_mapping[i].split('_')[2].split('.')[0]
    df['species']=species
    #df1=df[['chr','start','end','hg19_loc','species']]
    df=df.drop_duplicates() #default keep the first occurences
    print(df.shape)
    pd_list.append(df)

# concatenate all dataframe together
all=pd.concat(pd_list,ignore_index=True)

# find the variants with potential duplication
# keep=False so that all the duplications are parsed
dup_map=all[all.duplicated(subset=['hg19_loc','species'],keep=False)]
# make a liftOver input string
dup_map['liftOver_input']=dup_map['chr'].map(str) + ':' + dup_map['start'].map(str) + '-' + dup_map['end'].map(str)
# save the potential duplicates into a file
dup_map.to_csv('liftOver_duplications.tsv',sep="\t",header=True,index=False)

# all the hg19 variants are one-to-one mapped, therefore human variants are not in the potential duplicated subset
# species genome name list
species_names = list(dup_map.drop_duplicates(subset=['species'])['species'])

## count the total number of variants for each species
variant_total={}
for i in range(len(species_names)):
    # all[all['species']=='rheMac3'].shape[0]
    count=all[all['species']==species_names[i]].shape[0]
    variant_total[species_names[i]]=count
print(variant_total)
#{'C57B6J': 1601, 'Rattus': 1630, 'bosTau8': 2491, 'canFam3': 2618, 'felCat8': 2702, 'gorGor3': 5046, 'jacJac1': 1708, 'loxAfr3': 2377, 'micOch1': 1596, 'oryCun2': 1885, 'oviAri3': 2508, 'oviBos': 2471, 'panTro4': 5374, 'ponAbe2': 5216, 'rheMac3': 4758}

# seperate the dataframe by species and write bed files for each species
pd_list2=[]
for i in range(len(species_names)):
    df2=dup_map[dup_map['species']==species_names[i]]
    df2=df2[['chr','start','end','liftOver_input','score','strand']]
    df2.to_csv('liftOver_input/' + species_names[i] + '.bed',sep="\t",header=False,index=False)
    pd_list2.append(df2)
#
    # liftOver_input=list(dup_map.drop_duplicates(subset=['liftOver_input'])['liftOver_input'])
    # output=open('liftOver_input/' + species_names[i] + '.bed','w')
    # for element in liftOver_input:
    #     output.write(element)
    #     output.write('\n')
    # output.close()

# module='module load igmm/apps/R/3.3.0 igmm/apps/BEDTools/2.25.0 igmm/libs/ncurses/6.0 igmm/apps/samtools/1.3 igmm/apps/bcftools/1.3 igmm/apps/vcftools/0.1.13 igmm/libs/ensembl_api/86 igmm/apps/last/847 igmm/compilers/gcc/5.5.0'
# call(module, shell=True)
#
# liftOver='/exports/cmvm/eddie/sbms/groups/young-lab/rob/scripts/hal/halLiftover'
# call(liftOver, shell=True)

# liftOver
for i in range(len(species_names)):
    hal_genome = 'hg19'
    halLiftover = '/exports/cmvm/eddie/sbms/groups/young-lab/rob/scripts/hal/halLiftover /exports/cmvm/eddie/sbms/groups/young-lab/rob/genomes/1509_outgroups.hal ' + species_names[i] + ' liftOver_input/' + species_names[i] + '.bed ' + hal_genome + ' liftOver_output/' + species_names[i] + '_' + hal_genome + '.bed'
    print(str(datetime.datetime.now()))
    print(halLiftover,'\n\n')
    call(halLiftover, shell=True)

# read the liftOver output files
reverse_mapping=glob.glob('liftOver_output/*bed')
reverse_mapping.sort()

# read each file as a dataframe and store all the dataframe in a list
pd_list3=[]
for i in range(len(reverse_mapping)):
    df3=pd.read_csv(reverse_mapping[i],sep='\t',names=['chr','start','end','species_loc','score','strand'])
    df3['species']=species_names[i]
    df3=df3.drop_duplicates()
    pd_list3.append(df3)
# concatenate all the dataframe together
liftOver_all=pd.concat(pd_list3,ignore_index=True)
# make liftOver mapped hg19 locations into a string
liftOver_all['hg19_loc']=liftOver_all['chr'].map(str) + ':' + liftOver_all['start'].map(str) + '-' + liftOver_all['end'].map(str)

# match the original hg19 locations through species locations
# make a new column with original hg19 locations
species_loc=list(liftOver_all['species_loc'])
liftOver_all['hg19_loc_original']=0
for i in range(len(liftOver_all)):
    liftOver_all['hg19_loc_original'][i]=dup_map[dup_map['liftOver_input']==species_loc[i]]['hg19_loc'].to_string().split()[1]

# save the dataframe as tsv file
liftOver_all.to_csv('liftOver_complete.tsv',sep="\t",header=True,index=False)

# seperate the dataframe by species
pd_list4=[]
for i in range(len(species_names)):
    df4=liftOver_all[liftOver_all['species']==species_names[i]]
    df4=df4[['species','species_loc','hg19_loc_original','hg19_loc']]
    pd_list4.append(df4)

# make a nested dictionary {sepcies: {variant: non-human/human_ratio}}
species_dic={'C57B6J': {}, 'Rattus': {}, 'bosTau8': {}, 'canFam3': {}, 'felCat8': {}, 'gorGor3': {}, 'jacJac1': {}, 'loxAfr3': {}, 'micOch1': {}, 'oryCun2': {}, 'oviAri3': {}, 'oviBos': {}, 'panTro4': {}, 'ponAbe2': {}, 'rheMac3': {}}
pd_list5=[]
# for each species
for i in range(len(pd_list4)):
    pd_list4[i]['nh/h_ratio']=0
    hg19_ori_list=list(pd_list4[i].drop_duplicates(subset=['hg19_loc_original'])['hg19_loc_original'])
    # for each orignial hg19 location with dulications
    for j in range(len(hg19_ori_list)):
        df5=pd_list4[i][pd_list4[i]['hg19_loc_original']==hg19_ori_list[j]]
        # count the number of unique liftover reverse mapping
        hg19_count=len(set(list(df5['hg19_loc'])))
        # count the number of unique non-human variants at that hg19 location
        non_human_count=len(set(list(df5['species_loc'])))
        # if 1 means no duplication
        ratio=non_human_count/hg19_count
        # write the dictionary
        species_dic[species_names[i]][hg19_ori_list[j]]=ratio
        # write the ratio into dataframe
        df5['nh/h_ratio']=ratio
        pd_list5.append(df5)
# make dataframe with the ratio column
liftOver_all_ratio=pd.concat(pd_list5,ignore_index=True)
# write the dataframe into a tsv file
liftOver_all_ratio.to_csv('liftOver_complete_ratio.tsv',sep="\t",header=True,index=False)
        #condition=pd_list4[i][pd_list4[i]['hg19_loc_original'].isin([hg19_ori_list[j]])]
        #condition_a=pd_list4[i]['hg19_loc_original'].isin([hg19_ori_list[j]])
        #condition_b=pd_list4[i]['species'].isin([species_names[i]])
        #if condition_a and condition_b:
            #pd_list4[i]['nh/h_ratio']=ratio

## write the dictionary into a file
# {sepcies: {variant: non-human/human_ratio}}
with open("liftOver_ratio.txt", 'w') as f:
    for key,value in species_dic.items():
        for hg19,ratio in value.items():
            f.write('%s\t%s\t%s\n' % (key, hg19, ratio))

## assign score to the nonhuman-human ratio
# ratio >1 == 1, more copies in non-human, loss of duplication in human
# ratio =1 == 0, no change in evolution
# ratio <1 == -1 , less copies in non-human, duplication in human
df5=pd.read_csv('liftOver_ratio.txt',sep='\t',names=['species','hg19_loc_original','nh/h_ratio'])
df5['duplication_score']=0
for i in range(len(df5)):
    if df5['nh/h_ratio'][i]>1:
        df5['duplication_score'][i]=1
    elif df5['nh/h_ratio'][i]==1:
        df5['duplication_score'][i]=0
    else:
        df5['duplication_score'][i]=-1
df5.to_csv('liftOver_ratio_score.tsv',sep="\t",header=True,index=False)

## a dictionary to count the ratio scores
# {species:{ratio_score:count}}
score_dic={'C57B6J': {}, 'Rattus': {}, 'bosTau8': {}, 'canFam3': {}, 'felCat8': {}, 'gorGor3': {}, 'jacJac1': {}, 'loxAfr3': {}, 'micOch1': {}, 'oryCun2': {}, 'oviAri3': {}, 'oviBos': {}, 'panTro4': {}, 'ponAbe2': {}, 'rheMac3': {}}
for i in range(len(species_names)):
    score='grep -v \'duplication_score\' liftOver_ratio_score.tsv | grep ' + species_names[i] + ' | awk \'{print $NF}\' | sort | uniq -c'
    species_score=subprocess.check_output(score, shell=True)
    species_score=species_score.decode('UTF-8').split()
    for j in range(len(species_score)):
        if (j+1)%2==0:
            score_dic[species_names[i]][species_score[j]]=species_score[j-1]

## write the dictionary into a file
with open("liftOver_ratio_score_count.txt", 'w') as f:
    for key,value in score_dic.items():
        for score,count in value.items():
            f.write('%s\t%s\t%s\n' % (key, score, count))

## changed insertion chriteria
## TG/T treated as same locations for both start and end
## vep_input sorted
## print variant differences between COVID19_HGI_2021.bed and txt file
from subprocess import call
import subprocess
import datetime
import time

import glob
import pandas as pd


# load data
print('Loading data...\n\n\n')
f=open('hg19_v2/COVID19_HGI_2021.bed','r')
fc=f.read().split()

# make a dictionary with hg19 coordinates as keys and hg38 coordinates as items
hg_dic={}
for i in range(len(fc)):
    if (i+1)%4==0:
        hg_dic[fc[i-2]]=fc[i]

# read the GWAS txt files as pandas dataframe
gwas_files=glob.glob('GWASdata/*')
gwas_files.sort()
pd_list=[]
for i in range(len(gwas_files)):
    df=pd.read_csv(gwas_files[i],sep='\t')
    pd_list.append(df)

# a function to find diffences between two lists
def Diff(li1, li2):
    li_dif = [i for i in li1 + li2 if i not in li1 or i not in li2]
    return li_dif

# all gwas names, preparing for file names
gwas_name_list=['A2_ALL_eur_leave_23andme', 'A2_ALL_eur_leave_ukbb_23andme', 'A2_ALL_leave_23andme', 'A2_ALL_leave_UKBB_23andme', 'B1_ALL_eur_leave_23andme', 'B1_ALL_eur_leave_ukbb_23andme', 'B1_ALL_leave_23andme', 'B1_ALL_leave_UKBB_23andme', 'B2_ALL_eur_leave_23andme', 'B2_ALL_eur_leave_ukbb_23andme', 'B2_ALL_leave_23andme', 'B2_ALL_leave_UKBB_23andme', 'C2_ALL_eur_leave_23andme', 'C2_ALL_eur_leave_ukbb_23andme', 'C2_ALL_leave_23andme', 'C2_ALL_leave_UKBB_23andme']
#
# loop through gwas files to create vep inputs
for i in range(len(pd_list)):
    gwas_name=gwas_name_list[i]
    # new subset dataframe
    df1=pd_list[i][['#CHR','POS','REF','ALT','SNP','rsid']]
    df1['input']=0  # new input column
    #
    # concatenate input string (chr start end alleles)
    for i in range(len(df1)):
        if len(df1['REF'][i])>1 or df1['ALT'][i]=='-':    # Deletion
            pos_dif=len(df1['REF'][i])-1
            df1['input'][i]=str(df1['#CHR'][i]) + ' ' + str(df1['POS'][i]) + ' ' + str(df1['POS'][i]+pos_dif) + ' ' + df1['REF'][i] + '/' + df1['ALT'][i]
        elif df1['REF'][i]=='-':  # Insertion
            pos_dif=len(df1['ALT'][i])-1
            df1['input'][i]=str(df1['#CHR'][i]) + ' ' + str(df1['POS'][i]) + ' ' + str(df1['POS'][i]-pos_dif) + ' ' + df1['REF'][i] + '/' + df1['ALT'][i]
        else:
            df1['input'][i]=str(df1['#CHR'][i]) + ' ' + str(df1['POS'][i]) + ' ' + str(df1['POS'][i]) + ' ' + df1['REF'][i] + '/' + df1['ALT'][i]
    #
    # put the concatenated input strings of variants that match to COVID19_HGI_2021 hg38 coordinates into a list
    vep_input=[]
    df1['hg38']=df1['#CHR'].map(str) + ':' + df1['POS'].map(str)
    for i in hg_dic:
        for j in range(len(df1)):
            if hg_dic[i]==df1['hg38'][j]:
                vep_input.append(df1['input'][j])
    vep_input.sort()
    #
    # write the list into a txt file
    output=open('vep_input_v2/' + gwas_name + '_vep_input.txt','w')
    for element in vep_input:
        output.write(element)
        output.write('\n')
    output.close()
    print(Diff(list(df1['input']),vep_input))

########################################################################################################
########################################################################################################
########################################################################################################
### vep_v2.py before loop

# # make a dataframe using the GWAS txt file
# df = pd.read_csv('GWASdata/COVID19_HGI_C2_ALL_eur_leave_ukbb_23andme_20210107.txt.gz_1.0E-5.txt',sep='\t')
# df1=df[['#CHR','POS','REF','ALT','SNP','rsid']] # make a subset dataframe
# # make the VEP input strings
# df1['input']=0
# for i in range(len(df1)):
#     if len(df1['REF'][i])>1 or df1['ALT'][i]=='-':    # Deletion
#         pos_dif=len(df1['REF'][i])-1
#         df1['input'][i]=str(df1['#CHR'][i]) + ' ' + str(df1['POS'][i]) + ' ' + str(df1['POS'][i]+pos_dif) + ' ' + df1['REF'][i] + '/' + df1['ALT'][i]
#     elif df1['REF'][i]=='-':  # Insertion
#         pos_dif=len(df1['ALT'][i])-1
#         df1['input'][i]=str(df1['#CHR'][i]) + ' ' + str(df1['POS'][i]) + ' ' + str(df1['POS'][i]-pos_dif) + ' ' + df1['REF'][i] + '/' + df1['ALT'][i]
#     else:
#         df1['input'][i]=str(df1['#CHR'][i]) + ' ' + str(df1['POS'][i]) + ' ' + str(df1['POS'][i]) + ' ' + df1['REF'][i] + '/' + df1['ALT'][i]
#
# # put the concatenated input strings of variants that match to COVID19_HGI_2021 hg38 coordinates into a list
# vep_input=[]
# df1['hg38']=df1['#CHR'].map(str) + ':' + df1['POS'].map(str)
# for i in hg_dic:
#     for j in range(len(df1)):
#         if hg_dic[i]==df1['hg38'][j]:
#             vep_input.append(df1['input'][j])
#
# # write the list into a txt file
# output=open('vep_input.txt','w')
# for element in vep_input:
#      output.write(element)
#      output.write('\n')
# output.close()

########################################################################################################
########################################################################################################
########################################################################################################
### notes

# # Python code t get difference of two lists
# # Not using set()
# def Diff(li1, li2):
#     li_dif = [i for i in li1 + li2 if i not in li1 or i not in li2]
#     return li_dif

# Diff(list(df1['input']),vep_input)
# ['9 133270497 133270498 GA/G', '9 133270637 133270638 AT/A']

# df1[df1['input']=='9 133270497 133270498 GA/G']
#      #CHR        POS REF ALT               SNP          rsid  \
# 413     9  133270497  GA   G  9:133270497:GA:G  rs1291122587
#
#                           input         hg38
# 413  9 133270497 133270498 GA/G  9:133270497
# df1[df1['input']=='9 133270637 133270638 AT/A']
#      #CHR        POS REF ALT               SNP          rsid  \
# 415     9  133270637  AT   A  9:133270637:AT:A  rs1220069967
#
#                           input         hg38
# 415  9 133270637 133270638 AT/A  9:133270637

# make the GWAS name list
# name_list=[]
# for i in range(len(pd_list)):
#     name_list.append('_'.join(gwas_files[i].split('_')[2:8]))
# print(name_list)

# concatenate dataframe strings
# df1['input'] =df1['#CHR'].map(str) + ' ' + df1['POS'].map(str) + ' ' + df1['REF'].map(str) + '/' + df1['ALT'].map(str)

# for i in range(len(df1)):
#     if len(df1['REF'][i])>1 or len(df1['ALT'][i])>1:
#         print(i,df1.loc[[i]])

# example concatenate
# df1['input'][471]=str(df1['#CHR'][471]) + ' ' + str(df1['POS'][471]) + ' ' + str(df1['POS'][471]+1) + ' ' + df1['REF'][471] + '/' + df1['ALT'][471]

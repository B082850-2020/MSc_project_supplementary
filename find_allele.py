from subprocess import call
import subprocess
import glob
import datetime
import time
import os
import sys
import pandas as pd


# read the liftOver file with hg19 and hg38 coordiantes
df = pd.read_csv('hg19/COVID19_HGI_2021.bed',sep='\t',names=['chr','start','end','hg38_pos'])
df['hg19_pos'] = df['chr'].map(str) + ':' + df['start'].map(str) + '-' + df['end'].map(str)

# read the GWAS txt files as pandas dataframe
gwas_files = glob.glob('GWASdata/*')
gwas_files.sort()
pd_list = []
for i in range(len(gwas_files)):
    df_gwas = pd.read_csv(gwas_files[i],sep='\t')
    pd_list.append(df_gwas)

gwas_name=['A2_ALL_eur_leave_23andme', 'A2_ALL_eur_leave_ukbb_23andme', 'A2_ALL_leave_23andme', 'A2_ALL_leave_UKBB_23andme', 'B1_ALL_eur_leave_23andme', 'B1_ALL_eur_leave_ukbb_23andme', 'B1_ALL_leave_23andme', 'B1_ALL_leave_UKBB_23andme', 'B2_ALL_eur_leave_23andme', 'B2_ALL_eur_leave_ukbb_23andme', 'B2_ALL_leave_23andme', 'B2_ALL_leave_UKBB_23andme', 'C2_ALL_eur_leave_23andme', 'C2_ALL_eur_leave_ukbb_23andme', 'C2_ALL_leave_23andme', 'C2_ALL_leave_UKBB_23andme']

all_gwas = pd.concat(pd_list,ignore_index=True)
all_gwas = all_gwas.rename(columns = {'all_inv_var_meta_beta':'meta_beta'})
all_gwas = all_gwas[['#CHR','POS','REF','ALT','SNP','meta_beta']]
# make a hg38 coordinate string for each variant
all_gwas['hg38_pos'] = all_gwas['#CHR'].map(str) + ':' + all_gwas['POS'].map(str)

# find the hg19 coordinate for each variant from the COVID19_HGI_2021.bed file
all_gwas['hg19_pos'] = 0
for i in range(len(all_gwas)):
    if all_gwas['hg38_pos'][i] != '9:133270497' and all_gwas['hg38_pos'][i] != '9:133270637':
        all_gwas['hg19_pos'][i] = df[df['hg38_pos']==all_gwas['hg38_pos'][i]]['hg19_pos'].values[0]



# read the projection file
df1 = pd.read_csv('hg19/' + '0' + '_projections.fa.gz', sep='\t', names=['projection'])
df2 = df1.iloc[::2].reset_index(drop = True)    # reindex and remove the old index
df3 = df1.iloc[1::2].reset_index(drop = True)
df3 = df3.rename(columns = {'projection':'fa_data'})
df1 = pd.concat([df2, df3], axis = 1)

df1['species'] = 0
df1['hg19_pos'] = 0
df1['fa_second_nt'] = 0
for i in range(len(df1)):
    df1['species'][i] = df1['projection'][i].split(':')[0].split('>')[1]
    df1['hg19_pos'][i] = df1['projection'][i].split(':')[1] + ':' + df1['projection'][i].split(':')[2]
    df1['fa_second_nt'][i] = df1['fa_data'][i][1]


hg19_index = list(df1[df1['species']=='hg19'].index)
# len(hg19_index) = 506

df1['REF'] = 0
df1['ALT'] = 0
df1['old_beta'] = 0
for i in range(len(hg19_index)):
    gwas_search = all_gwas[all_gwas['hg19_pos']==df1[df1.index==hg19_index[i]]['hg19_pos'].values[0]].drop_duplicates()
    df1['REF'][hg19_index[i]] = gwas_search['REF'].values[0]
    df1['ALT'][hg19_index[i]] = gwas_search['ALT'].values[0]
    df1['old_beta'][hg19_index[i]] = str(gwas_search['meta_beta'].values[0])

df1 = df1[['species','hg19_pos','fa_data','fa_second_nt','REF','ALT','old_beta']]

c = 0
for i in range(len(df1)):
    if len(df1['fa_data'][i])<=2 and '-' in df1.iloc[i]['fa_data']:
        c += 1
#        print(df1.iloc[i])
        print(df1[df1.index==i])

#         and df1['species'][i]=='hg19'
# len(df1['fa_data'][i])>2  c = 106
# len(df1['fa_data'][i])>3  c = 12  >hg19:chr17:44065900-44065902

## order each variant mapping species by evolutionary distance
for i in range(len(hg19_index)):
    if hg19_index[i] < hg19_index[-1]:
        variant_mapping = df1.iloc[hg19_index[i]:hg19_index[i+1]]
        variant_mapping.species = pd.Categorical(variant_mapping.species,categories=['hg19', 'panTro4','gorGor3','ponAbe2','rheMac3','oryCun2','jacJac1','micOch1','Rattus','C57B6J','oviBos','oviAri3','bosTau8','canFam3','felCat8','loxAfr3'])
        variant_mapping = variant_mapping.sort_values(by = 'species')
        df1.iloc[hg19_index[i]:hg19_index[i+1]] = variant_mapping   #replace the original dataframe with the ordered version
    ## take care of the last variant
    else:
        variant_mapping = df1.iloc[hg19_index[i]:]
        variant_mapping.species = pd.Categorical(variant_mapping.species,categories=['hg19', 'panTro4','gorGor3','ponAbe2','rheMac3','oryCun2','jacJac1','micOch1','Rattus','C57B6J','oviBos','oviAri3','bosTau8','canFam3','felCat8','loxAfr3'])
        variant_mapping = variant_mapping.sort_values(by = 'species')
        df1.iloc[hg19_index[i]:] = variant_mapping   #replace the original dataframe with the ordered version


c0 = 0
c1 = 0
c2 = 0
c3 = 0
df1['new_beta'] = 0
for i in range(len(hg19_index)):
    if len(str(df1['REF'][hg19_index[i]])) == 1 and len(str(df1['ALT'][hg19_index[i]])) == 1 and len(str(df1['fa_data'][hg19_index[i]])) == 2:
        c0+=1   # c0 = 486
        ## ( human variants match to reference alleles )
        if df1['fa_second_nt'][hg19_index[i]] == df1['REF'][hg19_index[i]]:
            #c+=1   # 482
            ## ( chimp variants match to reference alleles )
            if df1['fa_second_nt'][hg19_index[i]+1] == df1['REF'][hg19_index[i]]:
                df1['new_beta'][hg19_index[i]] = df1['old_beta'][hg19_index[i]]
                #print('Hit!')
                c1+=1   # c1 = 248
            ## ( chimp variants match to alternative alleles )
            elif df1['fa_second_nt'][hg19_index[i]+1] == df1['ALT'][hg19_index[i]]:
                df1['new_beta'][hg19_index[i]] = str(-float(df1['old_beta'][hg19_index[i]]))
                #print('Reverse beta score!')
                c2+=1   # c2 = 136
            ## ( chimp variants match to neither )
            else:
                c3+=1   # c3 = 98
                print('Closest species match to neither!')
                print(df1.iloc[hg19_index[i]:hg19_index[i]+2])
                #print(hg19_index[i])
        else:
            print('fa does not match to reference')
            #print(df1[df1.index==hg19_index[i]])

df2=df1[df1.species=='hg19']
df2=df2[df2.new_beta!=0]
df3=df2[['hg19_pos','new_beta']]
df3=df3.reset_index(drop=True)


### updated version
### make sure that it is the chimp allele
### make sure to process the last variant
c0 = 0
c1 = 0
c2 = 0
c3 = 0
c4 = 0
df1['new_beta'] = 0
for i in range(len(hg19_index)):
    if len(str(df1['REF'][hg19_index[i]])) == 1 and len(str(df1['ALT'][hg19_index[i]])) == 1 and len(str(df1['fa_data'][hg19_index[i]])) == 2:
        c0+=1   # c0 = 486
        ## ( human variants match to reference alleles )
        if df1['fa_second_nt'][hg19_index[i]] == df1['REF'][hg19_index[i]]:
            if hg19_index[i] < hg19_index[-1]:
                df1_var = df1.iloc[hg19_index[i]:hg19_index[i+1]]
                if 'panTro4' in list(df1_var.species):
                    chimp_allele = df1_var[df1_var.species=='panTro4']['fa_second_nt'].values[0]
                    #c+=1   # 482
                    ## ( chimp variants match to reference alleles )
                    if chimp_allele == df1['REF'][hg19_index[i]]:
                        df1['new_beta'][hg19_index[i]] = df1['old_beta'][hg19_index[i]]
                        #print('Hit!')
                        c1+=1   # c1 = 248
                    ## ( chimp variants match to alternative alleles )
                    elif chimp_allele == df1['ALT'][hg19_index[i]]:
                        df1['new_beta'][hg19_index[i]] = str(-float(df1['old_beta'][hg19_index[i]]))
                        #print('Reverse beta score!')
                        c2+=1   # c2 = 136
                    ## ( chimp variants match to neither )
                    else:
                        c3+=1   # c3 = 98
                        print('Closest species match to neither!')
                        #print(df1.iloc[hg19_index[i]:hg19_index[i]+2])
                        #print(hg19_index[i])
                else:
                    c4+=1
                    print('chimp is not mapped at ', df1['hg19_pos'][hg19_index[i]])
            ## the last variant
            else:
                df1_var = df1.iloc[hg19_index[i]:]
                if 'panTro4' in list(df1_var.species):
                    chimp_allele = df1_var[df1_var.species=='panTro4']['fa_second_nt'].values[0]
                    #c+=1   # 482
                    ## ( chimp variants match to reference alleles )
                    if chimp_allele == df1['REF'][hg19_index[i]]:
                        df1['new_beta'][hg19_index[i]] = df1['old_beta'][hg19_index[i]]
                        #print('Hit!')
                        c1+=1   # c1 = 248
                    ## ( chimp variants match to alternative alleles )
                    elif chimp_allele == df1['ALT'][hg19_index[i]]:
                        df1['new_beta'][hg19_index[i]] = str(-float(df1['old_beta'][hg19_index[i]]))
                        #print('Reverse beta score!')
                        c2+=1   # c2 = 136
                    ## ( chimp variants match to neither )
                    else:
                        c3+=1   # c3 = 98
                        print('Closest species match to neither!')
                        #print(df1.iloc[hg19_index[i]:hg19_index[i]+2])
                        #print(hg19_index[i])
                else:
                    c4+=1
                    print('chimp is not mapped at ', df1['hg19_pos'][hg19_index[i]])
        else:
            print('human allele does not match to reference')
            #print(df1[df1.index==hg19_index[i]])


########################################################################################################
### make a loop

# read the liftOver file with hg19 and hg38 coordiantes
df = pd.read_csv('hg19/COVID19_HGI_2021.bed',sep='\t',names=['chr','start','end','hg38_pos'])
df['hg19_pos'] = df['chr'].map(str) + ':' + df['start'].map(str) + '-' + df['end'].map(str)

# read the GWAS txt files as pandas dataframe
gwas_files = glob.glob('GWASdata/*')
gwas_files.sort()
pd_list = []
for i in range(len(gwas_files)):
    df_gwas = pd.read_csv(gwas_files[i],sep='\t')
    pd_list.append(df_gwas)

gwas_name=['A2_ALL_eur_leave_23andme', 'A2_ALL_eur_leave_ukbb_23andme', 'A2_ALL_leave_23andme', 'A2_ALL_leave_UKBB_23andme', 'B1_ALL_eur_leave_23andme', 'B1_ALL_eur_leave_ukbb_23andme', 'B1_ALL_leave_23andme', 'B1_ALL_leave_UKBB_23andme', 'B2_ALL_eur_leave_23andme', 'B2_ALL_eur_leave_ukbb_23andme', 'B2_ALL_leave_23andme', 'B2_ALL_leave_UKBB_23andme', 'C2_ALL_eur_leave_23andme', 'C2_ALL_eur_leave_ukbb_23andme', 'C2_ALL_leave_23andme', 'C2_ALL_leave_UKBB_23andme']

all_gwas = pd.concat(pd_list,ignore_index=True)
all_gwas = all_gwas.rename(columns = {'all_inv_var_meta_beta':'meta_beta'})
all_gwas = all_gwas[['#CHR','POS','REF','ALT','SNP','meta_beta']]
# make a hg38 coordinate string for each variant
all_gwas['hg38_pos'] = all_gwas['#CHR'].map(str) + ':' + all_gwas['POS'].map(str)

# find the hg19 coordinate for each variant from the COVID19_HGI_2021.bed file
all_gwas['hg19_pos'] = 0
for i in range(len(all_gwas)):
    if all_gwas['hg38_pos'][i] != '9:133270497' and all_gwas['hg38_pos'][i] != '9:133270637':
        all_gwas['hg19_pos'][i] = df[df['hg38_pos']==all_gwas['hg38_pos'][i]]['hg19_pos'].values[0]

projection_list = []
for x in range(10):
    projection_file = 'hg19/' + str(x) + '_projections.fa.gz'
    print(projection_file)
    # read the projection file
    df1 = pd.read_csv(projection_file, sep='\t', names=['projection'])
    # odd lines have the coordinates
    df2 = df1.iloc[::2].reset_index(drop = True)    # reindex and remove the old index
    # even lines have the fasta data
    df3 = df1.iloc[1::2].reset_index(drop = True)
    df3 = df3.rename(columns = {'projection':'fa_data'}) # rename the column
    # put them back together
    df1 = pd.concat([df2, df3], axis = 1)
    # reformat the dataframe
    df1['species'] = 0
    df1['hg19_pos'] = 0
    df1['fa_second_nt'] = 0     # prepare the second nucleotide to map to the reference allele
    for i in range(len(df1)):
        df1['species'][i] = df1['projection'][i].split(':')[0].split('>')[1]
        df1['hg19_pos'][i] = df1['projection'][i].split(':')[1] + ':' + df1['projection'][i].split(':')[2]
        df1['fa_second_nt'][i] = df1['fa_data'][i][1]
    # indexes for all the human variants
    hg19_index = list(df1[df1['species']=='hg19'].index)
    # len(hg19_index) = 506
    # extract the REF ALT and beta values from all_gwas (GWAS txt files)
    df1['REF'] = 0
    df1['ALT'] = 0
    df1['old_beta'] = 0
    for i in range(len(hg19_index)):
        gwas_search = all_gwas[all_gwas['hg19_pos']==df1[df1.index==hg19_index[i]]['hg19_pos'].values[0]].drop_duplicates()
        df1['REF'][hg19_index[i]] = gwas_search['REF'].values[0]
        df1['ALT'][hg19_index[i]] = gwas_search['ALT'].values[0]
        df1['old_beta'][hg19_index[i]] = str(gwas_search['meta_beta'].values[0])
    # clean df1
    df1 = df1[['species','hg19_pos','fa_data','fa_second_nt','REF','ALT','old_beta']]
    ## order each variant mapping species by evolutionary distance
    for i in range(len(hg19_index)):
        if hg19_index[i] < hg19_index[-1]:
            variant_mapping = df1.iloc[hg19_index[i]:hg19_index[i+1]]
            variant_mapping.species = pd.Categorical(variant_mapping.species,categories=['hg19', 'panTro4','gorGor3','ponAbe2','rheMac3','oryCun2','jacJac1','micOch1','Rattus','C57B6J','oviBos','oviAri3','bosTau8','canFam3','felCat8','loxAfr3'])
            variant_mapping = variant_mapping.sort_values(by = 'species')
            df1.iloc[hg19_index[i]:hg19_index[i+1]] = variant_mapping   #replace the original dataframe with the ordered version
        ## take care of the last variant
        else:
            variant_mapping = df1.iloc[hg19_index[i]:]
            variant_mapping.species = pd.Categorical(variant_mapping.species,categories=['hg19', 'panTro4','gorGor3','ponAbe2','rheMac3','oryCun2','jacJac1','micOch1','Rattus','C57B6J','oviBos','oviAri3','bosTau8','canFam3','felCat8','loxAfr3'])
            variant_mapping = variant_mapping.sort_values(by = 'species')
            df1.iloc[hg19_index[i]:] = variant_mapping   #replace the original dataframe with the ordered version

    ## check human second nucleotide mapping to reference alleles
    ## check if cloest species map to reference or alternative alleles
    c0 = 0
    c1 = 0
    c2 = 0
    c3 = 0
    c4 = 0
    c5 = 0
    df1['new_beta'] = 0
    for i in range(len(hg19_index)):
        if len(str(df1['REF'][hg19_index[i]])) == 1 and len(str(df1['ALT'][hg19_index[i]])) == 1 and len(str(df1['fa_data'][hg19_index[i]])) == 2:
            c0+=1   # c0 = 486
            ## ( human variants match to reference alleles )
            if df1['fa_second_nt'][hg19_index[i]] == df1['REF'][hg19_index[i]]:
                if hg19_index[i] < hg19_index[-1]:
                    df1_var = df1.iloc[hg19_index[i]:hg19_index[i+1]]
                    if 'panTro4' in list(df1_var.species):
                        chimp_allele = df1_var[df1_var.species=='panTro4']['fa_second_nt'].values[0]
                        #c+=1   # 482
                        ## ( chimp variants match to reference alleles )
                        if chimp_allele == df1['REF'][hg19_index[i]]:
                            df1['new_beta'][hg19_index[i]] = df1['old_beta'][hg19_index[i]]
                            #print('Hit!')
                            c1+=1   # c1 = 248
                        ## ( chimp variants match to alternative alleles )
                        elif chimp_allele == df1['ALT'][hg19_index[i]]:
                            df1['new_beta'][hg19_index[i]] = str(-float(df1['old_beta'][hg19_index[i]]))
                            #print('Reverse beta score!')
                            c2+=1   # c2 = 136
                        ## ( chimp variants match to neither )
                        else:
                            c3+=1   # c3 = 98
                            print('Closest species match to neither!')
                            #print(df1.iloc[hg19_index[i]:hg19_index[i]+2])
                            #print(hg19_index[i])
                    else:
                        c4+=1
                        print('chimp is not mapped at ' + df1['hg19_pos'][hg19_index[i]])

                ## the last variant
                else:
                    df1_var = df1.iloc[hg19_index[i]:]
                    if 'panTro4' in list(df1_var.species):
                        chimp_allele = df1_var[df1_var.species=='panTro4']['fa_second_nt'].values[0]
                        #c+=1   # 482
                        ## ( chimp variants match to reference alleles )
                        if chimp_allele == df1['REF'][hg19_index[i]]:
                            df1['new_beta'][hg19_index[i]] = df1['old_beta'][hg19_index[i]]
                            #print('Hit!')
                            c1+=1   # c1 = 248
                        ## ( chimp variants match to alternative alleles )
                        elif chimp_allele == df1['ALT'][hg19_index[i]]:
                            df1['new_beta'][hg19_index[i]] = str(-float(df1['old_beta'][hg19_index[i]]))
                            #print('Reverse beta score!')
                            c2+=1   # c2 = 136
                        ## ( chimp variants match to neither )
                        else:
                            c3+=1   # c3 = 98
                            print('Closest species match to neither!')
                            #print(df1.iloc[hg19_index[i]:hg19_index[i]+2])
                            #print(hg19_index[i])
                    else:
                        c4+=1
                        print('chimp is not mapped at ' + df1['hg19_pos'][hg19_index[i]])
            else:
                c5+=1
                print('human allele does not match to reference')
    df2 = df1[df1.species=='hg19']
    df2 = df2[df2.new_beta!=0]
    df3 = df2[['hg19_pos','new_beta']]
    df3 = df3.reset_index(drop=True) #reset index and remove original index
    projection_list.append(df3)

    print('human allele does not match to reference: ' + str(c5))
    print('no liftover from human to chimp: ' + str(c4))
    print('human allele match to reference alleles: ' + str(c0 - c5))
    print('chimp allele match to reference alleles: ' + str(c1) + '\t Rate: ' + str(c1/c0*100))
    print('chimp allele match to alternative alleles: ' + str(c2) + '\t Rate: ' + str(c2/c0*100))
    print('chimp allele match to neither: '+ str(c3) + '\t Rate: ' + str(c3/c0*100) + '\n')


all_beta = pd.concat(projection_list,ignore_index=True)
output_dir = 'COVID_projection/projection_v2/'
all_beta.to_csv(output_dir + 'all_projection_beta.csv',sep=",",header=True,index=False)

gwas_beta_list = []
for x in range(len(gwas_name)):
    print(gwas_name[x])
    df1 = pd_list[x]
    df1['gwas_name'] = gwas_name[x]
    df1['hg38_pos'] = df1['#CHR'].map(str) + ':' + df1['POS'].map(str)
    df1['hg19_pos'] = 0
    df1['new_beta'] = 0
    for i in range(len(df1)):
        if df1['hg38_pos'][i] != '9:133270497' and df1['hg38_pos'][i] != '9:133270637':
            df1['hg19_pos'][i] = df[df['hg38_pos']==df1['hg38_pos'][i]]['hg19_pos'].values[0]
            if df1['hg19_pos'][i] in list(all_beta['hg19_pos']):
                df1['new_beta'][i] = all_beta[all_beta['hg19_pos']==df1['hg19_pos'][i]]['new_beta'].values[0]
    df1=df1[['gwas_name','hg38_pos','hg19_pos','new_beta']]
    gwas_beta_list.append(df1)
    df1.to_csv(output_dir + 'output_v1/' + gwas_name[x] + '_projection_beta.csv',sep=",",header=True,index=False)

all_gwas_beta = pd.concat(gwas_beta_list,ignore_index=True)
all_gwas_beta.to_csv(output_dir + 'output_v1/' + 'all_gwas_projection_beta.csv',sep=",",header=True,index=False)

## without all the 0 values
gwas_beta_list = []
for x in range(len(gwas_name)):
    df1 = pd.read_csv(output_dir + 'output_v1/' + gwas_name[x] + '_projection_beta.csv',sep=',')
    df1 = df1[df1.new_beta!=0]
    gwas_beta_list.append(df1)
    df1.to_csv(output_dir + 'output_v2/' + gwas_name[x] + '_projection_beta.csv',sep=",",header=True,index=False)

all_gwas_beta = pd.concat(gwas_beta_list,ignore_index=True)
all_gwas_beta.to_csv(output_dir + 'output_v2/' + 'all_gwas_projection_beta.csv',sep=",",header=True,index=False)


########################################################################################################
########################################################################################################
### For checking the matching of reference alleles


c = 0
for i in range(len(hg19_index)):
    if df1['fa_first_nt'][hg19_index[i]] not in df1['REF'][hg19_index[i]]:
        #c += 1     #C = 414 /412
        if df1['fa_first_nt'][hg19_index[i]] == df1['ALT'][hg19_index[i]]:
            c += 1      # C = 120 /118
            print(df1[df1.index==hg19_index[i]])

for i in range(10):
    print(df1[df1.index==hg19_index[i]+1])

c = 0
for i in range(len(df1)):
    if len(str(df1['REF'][i])) > 1 or len(str(df1['ALT'][i])) > 1:
        c+=1
        print(df1[df1.index==i])
        # c = 12

c=0
#df1['new_beta'] = 0
for i in range(len(hg19_index)):
    if df1['fa_first_nt'][hg19_index[i]] == df1['REF'][hg19_index[i]]:
        c+=1
        if df1['fa_first_nt'][hg19_index[i]+1] == df1['REF'][hg19_index[i]]:
            df1['new_beta'][hg19_index[i]] = df1['old_beta'][hg19_index[i]]
            print('Hit!')
        elif df1['fa_first_nt'][hg19_index[i]+1] == df1['ALT'][hg19_index[i]]:
            df1['new_beta'][hg19_index[i]] = str(-float(df1['old_beta'][hg19_index[i]]))
            print('Reverse beta score!')
        else:
            print('No hit for closest species!')
            print(hg19_index[i])

c=0
#df1['new_beta'] = 0
for i in range(len(hg19_index)):
    if df1['fa_second_nt'][hg19_index[i]] != df1['REF'][hg19_index[i]]:
        c+=1
        print(df1[df1.index==hg19_index[i]])



########################################################################################################
########################################################################################################
########################################################################################################
### for cat (felCat8)

# read the liftOver file with hg19 and hg38 coordiantes
df = pd.read_csv('hg19/COVID19_HGI_2021.bed',sep='\t',names=['chr','start','end','hg38_pos'])
df['hg19_pos'] = df['chr'].map(str) + ':' + df['start'].map(str) + '-' + df['end'].map(str)

# read the GWAS txt files as pandas dataframe
gwas_files = glob.glob('GWASdata/*')
gwas_files.sort()
pd_list = []
for i in range(len(gwas_files)):
    df_gwas = pd.read_csv(gwas_files[i],sep='\t')
    pd_list.append(df_gwas)

gwas_name=['A2_ALL_eur_leave_23andme', 'A2_ALL_eur_leave_ukbb_23andme', 'A2_ALL_leave_23andme', 'A2_ALL_leave_UKBB_23andme', 'B1_ALL_eur_leave_23andme', 'B1_ALL_eur_leave_ukbb_23andme', 'B1_ALL_leave_23andme', 'B1_ALL_leave_UKBB_23andme', 'B2_ALL_eur_leave_23andme', 'B2_ALL_eur_leave_ukbb_23andme', 'B2_ALL_leave_23andme', 'B2_ALL_leave_UKBB_23andme', 'C2_ALL_eur_leave_23andme', 'C2_ALL_eur_leave_ukbb_23andme', 'C2_ALL_leave_23andme', 'C2_ALL_leave_UKBB_23andme']

all_gwas = pd.concat(pd_list,ignore_index=True)
all_gwas = all_gwas.rename(columns = {'all_inv_var_meta_beta':'meta_beta'})
all_gwas = all_gwas[['#CHR','POS','REF','ALT','SNP','meta_beta']]
# make a hg38 coordinate string for each variant
all_gwas['hg38_pos'] = all_gwas['#CHR'].map(str) + ':' + all_gwas['POS'].map(str)

# find the hg19 coordinate for each variant from the COVID19_HGI_2021.bed file
all_gwas['hg19_pos'] = 0
for i in range(len(all_gwas)):
    if all_gwas['hg38_pos'][i] != '9:133270497' and all_gwas['hg38_pos'][i] != '9:133270637':
        all_gwas['hg19_pos'][i] = df[df['hg38_pos']==all_gwas['hg38_pos'][i]]['hg19_pos'].values[0]

projection_list = []
for x in range(10):
    projection_file = 'hg19/' + str(x) + '_projections.fa.gz'
    print(projection_file)
    # read the projection file
    df1 = pd.read_csv(projection_file, sep='\t', names=['projection'])
    # odd lines have the coordinates
    df2 = df1.iloc[::2].reset_index(drop = True)    # reindex and remove the old index
    # even lines have the fasta data
    df3 = df1.iloc[1::2].reset_index(drop = True)
    df3 = df3.rename(columns = {'projection':'fa_data'}) # rename the column
    # put them back together
    df1 = pd.concat([df2, df3], axis = 1)
    # reformat the dataframe
    df1['species'] = 0
    df1['hg19_pos'] = 0
    df1['fa_second_nt'] = 0     # prepare the second nucleotide to map to the reference allele
    for i in range(len(df1)):
        df1['species'][i] = df1['projection'][i].split(':')[0].split('>')[1]
        df1['hg19_pos'][i] = df1['projection'][i].split(':')[1] + ':' + df1['projection'][i].split(':')[2]
        df1['fa_second_nt'][i] = df1['fa_data'][i][1]
    # indexes for all the human variants
    hg19_index = list(df1[df1['species']=='hg19'].index)
    # len(hg19_index) = 506
    # extract the REF ALT and beta values from all_gwas (GWAS txt files)
    df1['REF'] = 0
    df1['ALT'] = 0
    df1['old_beta'] = 0
    for i in range(len(hg19_index)):
        gwas_search = all_gwas[all_gwas['hg19_pos']==df1[df1.index==hg19_index[i]]['hg19_pos'].values[0]].drop_duplicates()
        df1['REF'][hg19_index[i]] = gwas_search['REF'].values[0]
        df1['ALT'][hg19_index[i]] = gwas_search['ALT'].values[0]
        df1['old_beta'][hg19_index[i]] = str(gwas_search['meta_beta'].values[0])
    # clean df1
    df1 = df1[['species','hg19_pos','fa_data','fa_second_nt','REF','ALT','old_beta']]
    ## order each variant mapping species by evolutionary distance
    for i in range(len(hg19_index)):
        if hg19_index[i] < hg19_index[-1]:
            variant_mapping = df1.iloc[hg19_index[i]:hg19_index[i+1]]
            variant_mapping.species = pd.Categorical(variant_mapping.species,categories=['hg19', 'panTro4','gorGor3','ponAbe2','rheMac3','oryCun2','jacJac1','micOch1','Rattus','C57B6J','oviBos','oviAri3','bosTau8','canFam3','felCat8','loxAfr3'])
            variant_mapping = variant_mapping.sort_values(by = 'species')
            df1.iloc[hg19_index[i]:hg19_index[i+1]] = variant_mapping   #replace the original dataframe with the ordered version
        ## take care of the last variant
        else:
            variant_mapping = df1.iloc[hg19_index[i]:]
            variant_mapping.species = pd.Categorical(variant_mapping.species,categories=['hg19', 'panTro4','gorGor3','ponAbe2','rheMac3','oryCun2','jacJac1','micOch1','Rattus','C57B6J','oviBos','oviAri3','bosTau8','canFam3','felCat8','loxAfr3'])
            variant_mapping = variant_mapping.sort_values(by = 'species')
            df1.iloc[hg19_index[i]:] = variant_mapping   #replace the original dataframe with the ordered version

    ## check human second nucleotide mapping to reference alleles
    ## check if cloest species map to reference or alternative alleles
    c0 = 0
    c1 = 0
    c2 = 0
    c3 = 0
    c4 = 0
    c5 = 0
    df1['new_beta'] = 0
    for i in range(len(hg19_index)):
        if len(str(df1['REF'][hg19_index[i]])) == 1 and len(str(df1['ALT'][hg19_index[i]])) == 1 and len(str(df1['fa_data'][hg19_index[i]])) == 2:
            c0+=1   # c0 = 486
            ## ( human variants match to reference alleles )
            if df1['fa_second_nt'][hg19_index[i]] == df1['REF'][hg19_index[i]]:
                if hg19_index[i] < hg19_index[-1]:
                    df1_var = df1.iloc[hg19_index[i]:hg19_index[i+1]]
                    if 'felCat8' in list(df1_var.species):
                        cat_allele = df1_var[df1_var.species=='felCat8']['fa_second_nt'].values[0]
                        #c+=1   # 482
                        ## ( cat variants match to reference alleles )
                        if cat_allele == df1['REF'][hg19_index[i]]:
                            df1['new_beta'][hg19_index[i]] = df1['old_beta'][hg19_index[i]]
                            #print('Hit!')
                            c1+=1   # c1 = 248
                        ## ( cat variants match to alternative alleles )
                        elif cat_allele == df1['ALT'][hg19_index[i]]:
                            df1['new_beta'][hg19_index[i]] = str(-float(df1['old_beta'][hg19_index[i]]))
                            #print('Reverse beta score!')
                            c2+=1   # c2 = 136
                        ## ( cat variants match to neither )
                        else:
                            c3+=1   # c3 = 98
                            #print('cat allele match to neither')
                            #print(df1.iloc[hg19_index[i]:hg19_index[i]+2])
                            #print(hg19_index[i])
                    else:
                        c4+=1
                        #print('cat is not mapped at ' + df1['hg19_pos'][hg19_index[i]])

                ## the last variant
                else:
                    df1_var = df1.iloc[hg19_index[i]:]
                    if 'felCat8' in list(df1_var.species):
                        cat_allele = df1_var[df1_var.species=='felCat8']['fa_second_nt'].values[0]
                        #c+=1   # 482
                        ## ( cat variants match to reference alleles )
                        if cat_allele == df1['REF'][hg19_index[i]]:
                            df1['new_beta'][hg19_index[i]] = df1['old_beta'][hg19_index[i]]
                            #print('Hit!')
                            c1+=1   # c1 = 248
                        ## ( cat variants match to alternative alleles )
                        elif cat_allele == df1['ALT'][hg19_index[i]]:
                            df1['new_beta'][hg19_index[i]] = str(-float(df1['old_beta'][hg19_index[i]]))
                            #print('Reverse beta score!')
                            c2+=1   # c2 = 136
                        ## ( cat variants match to neither )
                        else:
                            c3+=1   # c3 = 98
                            #print('cat allele match to neither')
                            #print(df1.iloc[hg19_index[i]:hg19_index[i]+2])
                            #print(hg19_index[i])
                    else:
                        c4+=1
                        #print('cat is not mapped at ' + df1['hg19_pos'][hg19_index[i]])
            else:
                c5+=1
                print('human allele does not match to reference')

    df2 = df1[df1.species=='hg19']
    df2 = df2[df2.new_beta!=0]
    df3 = df2[['hg19_pos','new_beta']]
    df3 = df3.reset_index(drop=True) #reset index and remove original index
    projection_list.append(df3)

    print('human allele does not match to reference: ' + str(c5))
    print('no liftover from human to cat: ' + str(c4))
    print('human allele match to reference alleles: ' + str(c0 - c5))
    print('cat allele match to reference alleles: ' + str(c1) + '\t Rate: ' + str(c1/c0*100))
    print('cat allele match to alternative alleles: ' + str(c2) + '\t Rate: ' + str(c2/c0*100))
    print('cat allele match to neither: '+ str(c3) + '\t Rate: ' + str(c3/c0*100) + '\n')


all_beta = pd.concat(projection_list,ignore_index=True)
output_dir = 'COVID_projection/cat_projection_v1/'
all_beta.to_csv(output_dir + 'all_projection_beta.csv',sep=",",header=True,index=False)

gwas_beta_list = []
for x in range(len(gwas_name)):
    print(gwas_name[x])
    df1 = pd_list[x]
    df1['gwas_name'] = gwas_name[x]
    df1['hg38_pos'] = df1['#CHR'].map(str) + ':' + df1['POS'].map(str)
    df1['hg19_pos'] = 0
    df1['new_beta'] = 0
    for i in range(len(df1)):
        if df1['hg38_pos'][i] != '9:133270497' and df1['hg38_pos'][i] != '9:133270637':
            df1['hg19_pos'][i] = df[df['hg38_pos']==df1['hg38_pos'][i]]['hg19_pos'].values[0]
            if df1['hg19_pos'][i] in list(all_beta['hg19_pos']):
                df1['new_beta'][i] = all_beta[all_beta['hg19_pos']==df1['hg19_pos'][i]]['new_beta'].values[0]
    df1=df1[['gwas_name','hg38_pos','hg19_pos','new_beta']]
    #df1.to_csv(output_dir + 'output_v1/' + gwas_name[x] + '_projection_beta.csv',sep=",",header=True,index=False)
    df1 = df1[df1.new_beta!=0]  ## without all the 0 values
    gwas_beta_list.append(df1)
    df1.to_csv(output_dir + 'output_v1/' + gwas_name[x] + '_projection_beta.csv',sep=",",header=True,index=False)

all_gwas_beta = pd.concat(gwas_beta_list,ignore_index=True)
all_gwas_beta.to_csv(output_dir + 'output_v1/' + 'all_gwas_projection_beta.csv',sep=",",header=True,index=False)


########################################################################################################
########################################################################################################
########################################################################################################
### for all other species

# read the liftOver file with hg19 and hg38 coordiantes
df = pd.read_csv('hg19/COVID19_HGI_2021.bed',sep='\t',names=['chr','start','end','hg38_pos'])
df['hg19_pos'] = df['chr'].map(str) + ':' + df['start'].map(str) + '-' + df['end'].map(str)

# read the GWAS txt files as pandas dataframe
gwas_files = glob.glob('GWASdata/*')
gwas_files.sort()
pd_list = []
for i in range(len(gwas_files)):
    df_gwas = pd.read_csv(gwas_files[i],sep='\t')
    pd_list.append(df_gwas)

gwas_name=['A2_ALL_eur_leave_23andme', 'A2_ALL_eur_leave_ukbb_23andme', 'A2_ALL_leave_23andme', 'A2_ALL_leave_UKBB_23andme', 'B1_ALL_eur_leave_23andme', 'B1_ALL_eur_leave_ukbb_23andme', 'B1_ALL_leave_23andme', 'B1_ALL_leave_UKBB_23andme', 'B2_ALL_eur_leave_23andme', 'B2_ALL_eur_leave_ukbb_23andme', 'B2_ALL_leave_23andme', 'B2_ALL_leave_UKBB_23andme', 'C2_ALL_eur_leave_23andme', 'C2_ALL_eur_leave_ukbb_23andme', 'C2_ALL_leave_23andme', 'C2_ALL_leave_UKBB_23andme']

all_gwas = pd.concat(pd_list,ignore_index=True)
all_gwas = all_gwas.rename(columns = {'all_inv_var_meta_beta':'meta_beta'})
all_gwas = all_gwas[['#CHR','POS','REF','ALT','SNP','meta_beta']]
# make a hg38 coordinate string for each variant
all_gwas['hg38_pos'] = all_gwas['#CHR'].map(str) + ':' + all_gwas['POS'].map(str)

# find the hg19 coordinate for each variant from the COVID19_HGI_2021.bed file
all_gwas['hg19_pos'] = 0
for i in range(len(all_gwas)):
    if all_gwas['hg38_pos'][i] != '9:133270497' and all_gwas['hg38_pos'][i] != '9:133270637':
        all_gwas['hg19_pos'][i] = df[df['hg38_pos']==all_gwas['hg38_pos'][i]]['hg19_pos'].values[0]

species_an = ['panTro4','gorGor3','ponAbe2','rheMac3','oryCun2','jacJac1','micOch1','Rattus','C57B6J','oviBos','oviAri3','bosTau8','felCat8','canFam3','loxAfr3']
species_cn = ['chimpanzee','gorilla', 'orangutan', 'macaque', 'rabbit', 'egyptian_jerboa', 'prairie_vole', 'rat', 'mouse', 'musk_ox', 'sheep', 'cow','cat', 'dog', 'elephant']
species_names = {}
for i in range(len(species_an)):
    species_names[species_an[i]] = species_cn[i]

for ancestor in species_an:
    mkdir = 'mkdir -p COVID_projection/' + species_names[ancestor] + '_projection_v1/output_v1'
    silent_call = call(mkdir, shell=True)
    #rm = 'rm -r COVID_projection/'+ ancestor + '_projection_v1'
    #silent_call = call(rm, shell=True)


for ancestor in species_an:
    print(ancestor + ':' + species_names[ancestor] + '\n')
    output_dir = 'COVID_projection/' + species_names[ancestor] + '_projection_v1/'
    f = open('COVID_projection/' + species_names[ancestor] + '_projection_v1/' + species_names[ancestor] + '_find_allele_log.txt', 'w')
    projection_list = []
    for x in range(10):
        projection_file = 'hg19/' + str(x) + '_projections.fa.gz'
        print(projection_file)
        # read the projection file
        df1 = pd.read_csv(projection_file, sep='\t', names=['projection'])
        # odd lines have the coordinates
        df2 = df1.iloc[::2].reset_index(drop = True)    # reindex and remove the old index
        # even lines have the fasta data
        df3 = df1.iloc[1::2].reset_index(drop = True)
        df3 = df3.rename(columns = {'projection':'fa_data'}) # rename the column
        # put them back together
        df1 = pd.concat([df2, df3], axis = 1)
        # reformat the dataframe
        df1['species'] = 0
        df1['hg19_pos'] = 0
        df1['fa_second_nt'] = 0     # prepare the second nucleotide to map to the reference allele
        for i in range(len(df1)):
            df1['species'][i] = df1['projection'][i].split(':')[0].split('>')[1]
            df1['hg19_pos'][i] = df1['projection'][i].split(':')[1] + ':' + df1['projection'][i].split(':')[2]
            df1['fa_second_nt'][i] = df1['fa_data'][i][1]
        # indexes for all the human variants
        hg19_index = list(df1[df1['species']=='hg19'].index)
        # len(hg19_index) = 506
        # extract the REF ALT and beta values from all_gwas (GWAS txt files)
        df1['REF'] = 0
        df1['ALT'] = 0
        df1['old_beta'] = 0
        for i in range(len(hg19_index)):
            gwas_search = all_gwas[all_gwas['hg19_pos']==df1[df1.index==hg19_index[i]]['hg19_pos'].values[0]].drop_duplicates()
            df1['REF'][hg19_index[i]] = gwas_search['REF'].values[0]
            df1['ALT'][hg19_index[i]] = gwas_search['ALT'].values[0]
            df1['old_beta'][hg19_index[i]] = str(gwas_search['meta_beta'].values[0])
        # clean df1
        df1 = df1[['species','hg19_pos','fa_data','fa_second_nt','REF','ALT','old_beta']]
        ## order each variant mapping species by evolutionary distance
        for i in range(len(hg19_index)):
            if hg19_index[i] < hg19_index[-1]:
                variant_mapping = df1.iloc[hg19_index[i]:hg19_index[i+1]]
                variant_mapping.species = pd.Categorical(variant_mapping.species,categories=['hg19', 'panTro4','gorGor3','ponAbe2','rheMac3','oryCun2','jacJac1','micOch1','Rattus','C57B6J','oviBos','oviAri3','bosTau8','canFam3','felCat8','loxAfr3'])
                variant_mapping = variant_mapping.sort_values(by = 'species')
                df1.iloc[hg19_index[i]:hg19_index[i+1]] = variant_mapping   #replace the original dataframe with the ordered version
            ## take care of the last variant
            else:
                variant_mapping = df1.iloc[hg19_index[i]:]
                variant_mapping.species = pd.Categorical(variant_mapping.species,categories=['hg19', 'panTro4','gorGor3','ponAbe2','rheMac3','oryCun2','jacJac1','micOch1','Rattus','C57B6J','oviBos','oviAri3','bosTau8','canFam3','felCat8','loxAfr3'])
                variant_mapping = variant_mapping.sort_values(by = 'species')
                df1.iloc[hg19_index[i]:] = variant_mapping   #replace the original dataframe with the ordered version

        ## check human second nucleotide mapping to reference alleles
        ## check if cloest species map to reference or alternative alleles
        c0 = 0
        c1 = 0
        c2 = 0
        c3 = 0
        c4 = 0
        c5 = 0
        df1['new_beta'] = 0
        for i in range(len(hg19_index)):
            if len(str(df1['REF'][hg19_index[i]])) == 1 and len(str(df1['ALT'][hg19_index[i]])) == 1 and len(str(df1['fa_data'][hg19_index[i]])) == 2:
                c0+=1   # c0 = 486
                ## ( human variants match to reference alleles )
                if df1['fa_second_nt'][hg19_index[i]] == df1['REF'][hg19_index[i]]:
                    if hg19_index[i] < hg19_index[-1]:
                        df1_var = df1.iloc[hg19_index[i]:hg19_index[i+1]]
                        if ancestor in list(df1_var.species):
                            ancestor_allele = df1_var[df1_var.species==ancestor]['fa_second_nt'].values[0]
                            #c+=1   # 482
                            ## ( cat variants match to reference alleles )
                            if ancestor_allele == df1['REF'][hg19_index[i]]:
                                df1['new_beta'][hg19_index[i]] = df1['old_beta'][hg19_index[i]]
                                #print('Hit!')
                                c1+=1   # c1 = 248
                            ## ( cat variants match to alternative alleles )
                            elif ancestor_allele == df1['ALT'][hg19_index[i]]:
                                df1['new_beta'][hg19_index[i]] = str(-float(df1['old_beta'][hg19_index[i]]))
                                #print('Reverse beta score!')
                                c2+=1   # c2 = 136
                            ## ( cat variants match to neither )
                            else:
                                c3+=1   # c3 = 98
                                #print('cat allele match to neither')
                                #print(df1.iloc[hg19_index[i]:hg19_index[i]+2])
                                #print(hg19_index[i])
                        else:
                            c4+=1
                            #print('cat is not mapped at ' + df1['hg19_pos'][hg19_index[i]])

                    ## the last variant
                    else:
                        df1_var = df1.iloc[hg19_index[i]:]
                        if ancestor in list(df1_var.species):
                            ancestor_allele = df1_var[df1_var.species==ancestor]['fa_second_nt'].values[0]
                            #c+=1   # 482
                            ## ( cat variants match to reference alleles )
                            if ancestor_allele == df1['REF'][hg19_index[i]]:
                                df1['new_beta'][hg19_index[i]] = df1['old_beta'][hg19_index[i]]
                                #print('Hit!')
                                c1+=1   # c1 = 248
                            ## ( cat variants match to alternative alleles )
                            elif ancestor_allele == df1['ALT'][hg19_index[i]]:
                                df1['new_beta'][hg19_index[i]] = str(-float(df1['old_beta'][hg19_index[i]]))
                                #print('Reverse beta score!')
                                c2+=1   # c2 = 136
                            ## ( cat variants match to neither )
                            else:
                                c3+=1   # c3 = 98
                                #print('cat allele match to neither')
                                #print(df1.iloc[hg19_index[i]:hg19_index[i]+2])
                                #print(hg19_index[i])
                        else:
                            c4+=1
                            #print('cat is not mapped at ' + df1['hg19_pos'][hg19_index[i]])
                else:
                    c5+=1
                    print('human allele does not match to reference')

        df2 = df1[df1.species=='hg19']
        df2 = df2[df2.new_beta!=0]
        df3 = df2[['hg19_pos','new_beta']]
        df3 = df3.reset_index(drop=True) #reset index and remove original index
        projection_list.append(df3)

        f.write('human allele does not match to reference: ' + str(c5) + '\n')
        f.write('no liftover from human to '+ ancestor + ': ' + str(c4) + '\n')
        f.write('human allele match to reference alleles: ' + str(c0 - c5) + '\n')
        f.write(ancestor + ' allele match to reference alleles: ' + str(c1) + '\t Rate: ' + str(c1/c0*100) + '\n')
        f.write(ancestor + ' allele match to alternative alleles: ' + str(c2) + '\t Rate: ' + str(c2/c0*100) + '\n')
        f.write(ancestor + ' allele match to neither: '+ str(c3) + '\t Rate: ' + str(c3/c0*100) + '\n\n')
    f.close()

    output_dir = 'COVID_projection/' + species_names[ancestor] + '_projection_v1/'
    all_beta = pd.concat(projection_list,ignore_index=True)
    peint(all_beta)
    all_beta.to_csv(output_dir + 'all_projection_beta.csv',sep=",",header=True,index=False)

    gwas_beta_list = []
    for x in range(len(gwas_name)):
        print(gwas_name[x])
        df1 = pd_list[x]
        df1['gwas_name'] = gwas_name[x]
        df1['hg38_pos'] = df1['#CHR'].map(str) + ':' + df1['POS'].map(str)
        df1['hg19_pos'] = 0
        df1['new_beta'] = 0
        for i in range(len(df1)):
            if df1['hg38_pos'][i] != '9:133270497' and df1['hg38_pos'][i] != '9:133270637':
                df1['hg19_pos'][i] = df[df['hg38_pos']==df1['hg38_pos'][i]]['hg19_pos'].values[0]
                if df1['hg19_pos'][i] in list(all_beta['hg19_pos']):
                    df1['new_beta'][i] = all_beta[all_beta['hg19_pos']==df1['hg19_pos'][i]]['new_beta'].values[0]
        df1=df1[['gwas_name','hg38_pos','hg19_pos','new_beta']]
        #df1.to_csv(output_dir + 'output_v1/' + gwas_name[x] + '_projection_beta.csv',sep=",",header=True,index=False)
        df1 = df1[df1.new_beta!=0]  ## without all the 0 values
        gwas_beta_list.append(df1)
        df1.to_csv(output_dir + 'output_v1/' + gwas_name[x] + '_projection_beta.csv',sep=",",header=True,index=False)

    all_gwas_beta = pd.concat(gwas_beta_list,ignore_index=True)
    all_gwas_beta.to_csv(output_dir + 'output_v1/' + 'all_gwas_projection_beta.csv',sep=",",header=True,index=False)


########################################################################################################

species_an = ['panTro4','gorGor3','ponAbe2','rheMac3','oryCun2','jacJac1','micOch1','Rattus','C57B6J','oviBos','oviAri3','bosTau8','felCat8','canFam3','loxAfr3']
species_cn = ['chimpanzee','gorilla', 'orangutan', 'macaque', 'rabbit', 'egyptian_jerboa', 'prairie_vole', 'rat', 'mouse', 'musk_ox', 'sheep', 'cow','cat', 'dog', 'elephant']
species_names = {}
for i in range(len(species_an)):
    species_names[species_an[i]] = species_cn[i]

f = open('COVID_projection/find_allele_ratio.tsv', 'w')
f.write('species' + '\t' + 'ref' + '\t' + 'alt' + '\t' + 'neither' + '\t' + 'ref/alt' + '\n')
for ancestor in species_an:
    output_dir = 'COVID_projection/' + species_names[ancestor] + '_projection_v1/'
    df1_el = pd.read_csv(output_dir + species_names[ancestor] + '_find_allele_log.txt', sep='\t', names=['freq_info','rate_info'])
    # reformat the dataframe
    df1_el['name'] = 0
    df1_el['freq'] = 0
    df1_el['rate'] = 0
    df1_el['total_freq'] = 0
    for i in range(len(df1_el)):
        df1_el['name'][i] = df1_el.freq_info[i].split(':')[0]
        df1_el['freq'][i] = df1_el.freq_info[i].split(':')[1]
        if str(df1_el.rate_info[i]) != 'nan':
            df1_el['rate'][i] = df1_el.rate_info[i].split(':')[1]
    # take subset of the dataframe
    df2_el = df1_el[['name','freq','rate']]
    # count the sum mapping
    total_ref = sum(list(df2_el.iloc[3::6].freq))
    total_alt = sum(list(df2_el.iloc[4::6].freq))
    total_neither = sum(list(df2_el.iloc[5::6].freq))
    print(species_names[ancestor] + '\t' + str(total_ref) + '\t' + str(total_alt) + '\t' + str(total_neither) + '\t' + str(total_ref/total_alt) + '\n')
    f.write(species_names[ancestor] + '\t' + str(total_ref) + '\t' + str(total_alt) + '\t' + str(total_neither) + '\t' + str(total_ref/total_alt) + '\n')

f.close()





print(ancestor + ':' + species_names[ancestor])
projection_list = []
for x in range(10):
    projection_file = 'hg19/' + str(x) + '_projections.fa.gz'
    print(projection_file)
    # read the projection file
    df1 = pd.read_csv(projection_file, sep='\t', names=['projection'])
    # odd lines have the coordinates
    df2 = df1.iloc[::2].reset_index(drop = True)    # reindex and remove the old index
    # even lines have the fasta data
    df3 = df1.iloc[1::2].reset_index(drop = True)
    df3 = df3.rename(columns = {'projection':'fa_data'}) # rename the column
    # put them back together
    df1 = pd.concat([df2, df3], axis = 1)
    # reformat the dataframe
    df1['species'] = 0
    df1['hg19_pos'] = 0
    df1['fa_second_nt'] = 0     # prepare the second nucleotide to map to the reference allele
    for i in range(len(df1)):
        df1['species'][i] = df1['projection'][i].split(':')[0].split('>')[1]
        df1['hg19_pos'][i] = df1['projection'][i].split(':')[1] + ':' + df1['projection'][i].split(':')[2]
        df1['fa_second_nt'][i] = df1['fa_data'][i][1]
    # indexes for all the human variants
    hg19_index = list(df1[df1['species']=='hg19'].index)
    # len(hg19_index) = 506
    # extract the REF ALT and beta values from all_gwas (GWAS txt files)
    df1['REF'] = 0
    df1['ALT'] = 0
    df1['old_beta'] = 0
    for i in range(len(hg19_index)):
        gwas_search = all_gwas[all_gwas['hg19_pos']==df1[df1.index==hg19_index[i]]['hg19_pos'].values[0]].drop_duplicates()
        df1['REF'][hg19_index[i]] = gwas_search['REF'].values[0]
        df1['ALT'][hg19_index[i]] = gwas_search['ALT'].values[0]
        df1['old_beta'][hg19_index[i]] = str(gwas_search['meta_beta'].values[0])
    # clean df1
    df1 = df1[['species','hg19_pos','fa_data','fa_second_nt','REF','ALT','old_beta']]
    ## order each variant mapping species by evolutionary distance
    for i in range(len(hg19_index)):
        if hg19_index[i] < hg19_index[-1]:
            variant_mapping = df1.iloc[hg19_index[i]:hg19_index[i+1]]
            variant_mapping.species = pd.Categorical(variant_mapping.species,categories=['hg19', 'panTro4','gorGor3','ponAbe2','rheMac3','oryCun2','jacJac1','micOch1','Rattus','C57B6J','oviBos','oviAri3','bosTau8','canFam3','felCat8','loxAfr3'])
            variant_mapping = variant_mapping.sort_values(by = 'species')
            df1.iloc[hg19_index[i]:hg19_index[i+1]] = variant_mapping   #replace the original dataframe with the ordered version
        ## take care of the last variant
        else:
            variant_mapping = df1.iloc[hg19_index[i]:]
            variant_mapping.species = pd.Categorical(variant_mapping.species,categories=['hg19', 'panTro4','gorGor3','ponAbe2','rheMac3','oryCun2','jacJac1','micOch1','Rattus','C57B6J','oviBos','oviAri3','bosTau8','canFam3','felCat8','loxAfr3'])
            variant_mapping = variant_mapping.sort_values(by = 'species')
            df1.iloc[hg19_index[i]:] = variant_mapping   #replace the original dataframe with the ordered version

    ## check human second nucleotide mapping to reference alleles
    ## check if cloest species map to reference or alternative alleles
    c0 = 0
    c1 = 0
    c2 = 0
    c3 = 0
    c4 = 0
    c5 = 0
    df1['new_beta'] = 0
    for i in range(len(hg19_index)):
        if len(str(df1['REF'][hg19_index[i]])) == 1 and len(str(df1['ALT'][hg19_index[i]])) == 1 and len(str(df1['fa_data'][hg19_index[i]])) == 2:
            c0+=1   # c0 = 486
            ## ( human variants match to reference alleles )
            if df1['fa_second_nt'][hg19_index[i]] == df1['REF'][hg19_index[i]]:
                if hg19_index[i] < hg19_index[-1]:
                    df1_var = df1.iloc[hg19_index[i]:hg19_index[i+1]]
                    if ancestor in list(df1_var.species):
                        ancestor_allele = df1_var[df1_var.species==ancestor]['fa_second_nt'].values[0]
                        #c+=1   # 482
                        ## ( cat variants match to reference alleles )
                        if ancestor_allele == df1['REF'][hg19_index[i]]:
                            df1['new_beta'][hg19_index[i]] = df1['old_beta'][hg19_index[i]]
                            #print('Hit!')
                            c1+=1   # c1 = 248
                        ## ( cat variants match to alternative alleles )
                        elif ancestor_allele == df1['ALT'][hg19_index[i]]:
                            df1['new_beta'][hg19_index[i]] = str(-float(df1['old_beta'][hg19_index[i]]))
                            #print('Reverse beta score!')
                            c2+=1   # c2 = 136
                        ## ( cat variants match to neither )
                        else:
                            c3+=1   # c3 = 98
                            #print('cat allele match to neither')
                            #print(df1.iloc[hg19_index[i]:hg19_index[i]+2])
                            #print(hg19_index[i])
                    else:
                        c4+=1
                        #print('cat is not mapped at ' + df1['hg19_pos'][hg19_index[i]])

                ## the last variant
                else:
                    df1_var = df1.iloc[hg19_index[i]:]
                    if ancestor in list(df1_var.species):
                        ancestor_allele = df1_var[df1_var.species==ancestor]['fa_second_nt'].values[0]
                        #c+=1   # 482
                        ## ( cat variants match to reference alleles )
                        if ancestor_allele == df1['REF'][hg19_index[i]]:
                            df1['new_beta'][hg19_index[i]] = df1['old_beta'][hg19_index[i]]
                            #print('Hit!')
                            c1+=1   # c1 = 248
                        ## ( cat variants match to alternative alleles )
                        elif ancestor_allele == df1['ALT'][hg19_index[i]]:
                            df1['new_beta'][hg19_index[i]] = str(-float(df1['old_beta'][hg19_index[i]]))
                            #print('Reverse beta score!')
                            c2+=1   # c2 = 136
                        ## ( cat variants match to neither )
                        else:
                            c3+=1   # c3 = 98
                            #print('cat allele match to neither')
                            #print(df1.iloc[hg19_index[i]:hg19_index[i]+2])
                            #print(hg19_index[i])
                    else:
                        c4+=1
                        #print('cat is not mapped at ' + df1['hg19_pos'][hg19_index[i]])
            else:
                c5+=1
                print('human allele does not match to reference')

    df2 = df1[df1.species=='hg19']
    df2 = df2[df2.new_beta!=0]
    df3 = df2[['hg19_pos','new_beta']]
    df3 = df3.reset_index(drop=True) #reset index and remove original index
    projection_list.append(df3)

#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
from pandas import read_csv
import urllib.request
import time


# In[2]:


def get_assemblies(phages_list, path):
    '''
    This function download genomes from flirting list of Bacteriophages to concrete path
    '''
    #import assembly_summary_refseq file to dataframe
    assembly_sum = pd.read_csv('/home/amanda/assembly_summary_refseq.txt', sep='\t', skiprows=1) 
    #names of columns
    assembly_sum.columns = [
        'assembly_accession',
        'bioproject','biosample',
        'wgs_master','refseq_category',
        'taxid','species_taxid','organism_name',
        'infraspecific_name','isolate','version_status',
        'assembly_level','release_type','genome_rep',
        'seq_rel_date','asm_name','submitter','gbrs_paired_asm',
        'paired_asm_comp','ftp_path','excluded_from_refseq',
        'relation_to_type_material','asm_not_live_date'
    ]

    for assembly in phages_list:
        # searching specific genomes from list
        my_df = assembly_sum[(assembly_sum['assembly_accession'] == assembly)]
        #Process the newly created file and download genomes from NCBI website
        ftp = my_df['ftp_path'].tolist() #making ftp list --> path to download genomes
        asm = my_df['asm_name'].to_list() #making asm list --> asm necessary as part of suffix

        for elem in ftp:
            for i in asm:
                file_in = assembly + '.fna.gz' #gzip format
                fullfilename = os.path.join(path, file_in) #directory and file name
                file_suffix=elem+'/'+assembly+'_'+i+'_genomic.fna.gz'
                try:   
                    if os.path.isfile(fullfilename): #if genome is in directory, skip it and continue the rest of them
                        print(file_in, " already exist")
                        continue
                    else:
                        response = urllib.request.urlretrieve(file_suffix, fullfilename) #download genomes
                        print("Download:", file_in)
                        time.sleep(1)

                except:        
                    print("Skipping", file_in, " - doesn't exist.") #If there is an error or the ftp server doesn't have the genome, skip it


# In[ ]:





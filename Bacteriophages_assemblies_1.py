#!/usr/bin/env python
# coding: utf-8

# In[158]:


## IMPORT ##

import pandas as pd
from pandas import read_csv
from Bio import Entrez
import os
import urllib
import requests


# In[2]:


seq = pd.read_csv('sequences.csv') 
seq


# In[3]:


#new dataframe with individual columns
new_seq = seq[['Assembly', 'Species', 'Molecule_type', 'Family', 'Host', 'GenBank_Title']].copy()


# In[4]:


new_seq


# In[5]:


#sorting values by Host and Family
new_seq.sort_values(by=['Host', 'Family'])


# In[6]:


#searching Family Siphoviridae which Host is Lactococcus lactis
siphoviridae = new_seq.loc [(new_seq['Family'] == 'Siphoviridae') & (new_seq['Host'] == 'Lactococcus lactis')] 
siphoviridae


# In[7]:


asb = siphoviridae["Assembly"].tolist()


# In[8]:


asb


# In[9]:


with open('asb.txt', 'w') as f:
    for item in asb:
        f.write("%s\n" % item)


# In[10]:


#from https://dmnfarrell.github.io/bioinformatics/assemblies-genbank-python
def get_assembly_summary(id):
    """Get esummary for an entrez id"""
    from Bio import Entrez
    esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
    esummary_record = Entrez.read(esummary_handle)
    return esummary_record

def get_assemblies(term, download=True, path='assemblies'):
    """Download genbank assemblies for a given search term.
    Args:
        term: search term, usually organism name
        download: whether to download the results
        path: folder to save to
    """
    #provide your own mail here
    Entrez.email = "amandakoryl@gmail.com"
    handle = Entrez.esearch(db="assembly", term=term, retmax='100')
    record = Entrez.read(handle)
    ids = record['IdList']
    print (f'found {len(ids)} ids')
    links = []
    for id in ids:
        #get summary
        summary = get_assembly_summary(id)
        #get ftp link
        url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
        if url == '':
            continue
        label = os.path.basename(url)
        link = os.path.join(url,label+'_genomic.fna.gz')
        print (link)
        links.append(link)
        if download == True:
            #download link
            urllib.request.urlretrieve(link, f'{label}.fna.gz')
    return links


# In[11]:


get_ipython().run_cell_magic('bash', '', '\n\nwget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt\n\n')


# In[150]:


get_ipython().run_cell_magic('bash', '', '\ngrep GCF_000910615.1 assembly_summary_refseq.txt \\\n    | awk \'BEGIN{FS="\\t"}{if($12=="Complete Genome"){print $20}}\' \\\n    | awk \'BEGIN{OFS=FS="/"}{print $0,$NF"_genomic.fna.gz"}\' \\\n    > urls.txt\n\n#Now you can go through your urls file\n\nIFS=$\'\\n\'; for NEXT in $(cat urls.txt); do wget "$NEXT" iconv -f ISO-8859-1 -t UTF-8 ; done')


# In[144]:


import codecs

doc = codecs.open('assembly_summary_refseq.txt','rU','UTF-8') #open for reading with "universal" type set

df = pd.read_csv(doc, sep='\t')


# In[145]:


df


# In[ ]:


# assembly_accession	


# In[185]:


def get_assemblies(asb, download=True):
    from Bio import Entrez
    Entrez.email = 'amandakoryl@gmail.com'
    for id in asb:
        handle = Entrez.esearch(db="assembly", term=id)
        record = Entrez.read(handle)
    for i in record['IdList']:
    # Get Assembly Summary
        esummary_handle = Entrez.esummary(db="assembly", id=i, report="full")
        esummary_record = Entrez.read(esummary_handle)
    links = []
    for id in asb:

        handle1 = Entrez.efetch(db="nucleotide", id=id, retmode='text')
        record = Entrez.read(handle1)
        url = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
        if url == '':
            print("ups")
            break
        label = os.path.basename(url)
        print(label)
        link = os.path.join(url,label+'_genomic.fna.gz')
        print(link)
        links.append(link)
        if download == True:
        #download link
            urllib.request.urlretrieve(link, f'{label}.fna.gz')
        return links




    
    


# In[186]:


get_assemblies(asb, download=True)


# In[273]:



import re
first = []
with open('assembly_summary_refseq.txt', "r", encoding = "utf-8") as f:
    lines = f.readlines()[1:]
    for line in lines:
        line = str(line)
        re.sub(r'(\t)+', '\t', line)
        line.replace("#", "")
        first.append(line)

        


# In[274]:


print(first[:100])


# In[239]:



assembly_summary_file="ftp://ftp.ncbi.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
#Download the file using wget sysyem call
subprocess.call("wget "+assembly_summary_file, shell=True)
#Reformat the file to pandas-friendly format
subprocess.call("sed -i '1d' assembly_summary_refseq.txt",shell=True)
subprocess.call("sed -i 's/^# //' assembly_summary_refseq.txt", shell=True)
#Read the file as a dataframe - using read_table
#Use read_table if the column separator is tab

assembly_sum = pd.read_csv('assembly_summary_refseq.txt', dtype='unicode', header=None, delimiter='\t+', engine='python')
assembly_sum.columns = ["assembly_accession", "bioproject", "biosample", "wgs_master", "refseq_category", "taxid", "species_taxid", "organism_name", "infraspecific_name", "isolate", "version_status", "assembly_level", "release_type", "genome_rep", "seq_rel_date", "asm_name", "submitter", "gbrs_paired_asm", "paired_asm_comp", "ftp_path", "excluded_from_refseq", "relation_to_type_material"]


# In[240]:


assembly_sum


# In[ ]:


def get_this(asb, outfile='refseq_genome.txt'):
for id in asb:
        my_df=assembly_sum[(assembly_sum['assembly_accession'] == id) &
                           ((assembly_sum['refseq_category'] == 'reference genome') |
                            (assembly_sum['refseq_category'] == 'representative genome')
                           )]
        my_df=my_df[['ftp_path','assembly_accession','asm_name']]
        #Process the newly created file and download genomes from NCBI website
        my_df.to_csv(outfile,mode='w',index=False,header=None)
        process_url_file(outfile)
    return


# In[238]:


get_this(asb, outfile='refseq_genome.txt')


# In[ ]:





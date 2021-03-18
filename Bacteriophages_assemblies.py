#!/usr/bin/env python
# coding: utf-8

# In[1]:


## IMPORT ##

import pandas as pd
from pandas import read_csv
from Bio import Entrez
import os
import urllib


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


# In[ ]:


get_ipython().run_cell_magic('bash', '', '\n\nwget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt\n\n\ngrep -f asb.txt assembly_summary_refseq.txt \\\n    | awk \'BEGIN{FS="\\t"}{if($12=="Complete Genome"){print $20}}\' \\\n    | awk \'BEGIN{OFS=FS="/"}{print $0,$NF"_genomic.fna.gz"}\' \\\n    > urls.txt\n\n#Now you can go through your urls file\nIFS=$\'\\n\'; for NEXT in $(cat urls.txt); do wget "$NEXT"; done')


# In[ ]:





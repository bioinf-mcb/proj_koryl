
### IMPORT ###

import os
from Bio import Entrez
import subprocess
import pandas as pd
import re
import sys
import requests
from bs4 import BeautifulSoup

####



def folder():
    '''
    This function create new folder in directory
    '''
    parent_dir = input("Write parent dir: ")
    folder_name = input("Write new folder name: ")
    path = os.path.join(parent_dir, folder_name)
    os.mkdir(path)
    print("Directory '% s' created" % path)
    
def directory():
    os.chdir(input("Write path to change directory: "))
    print('Current directory is: ', os.getcwd())

def acc_list(acc_path):
    with open(acc_path, 'r', encoding='UTF-8') as f:
        lines = f.readlines()
        lines = [line.rstrip() for line in lines]
    for i in range(0, len(lines), 7000):
        yield lines[i:i + 7000]
        
def get_ncbi_genomes():

    Entrez.email = input('Write your e-mail: ')
    list_of_lists = acc_list(input("Write path to the accession file: "))
    filename = input("Write path with new filename (in FASTA format): ")

    for elem in list_of_lists:    
        net_handle = Entrez.efetch(
            db="nucleotide", id=elem, rettype="fasta", retmode="text"
        )
        out_handle = open(filename, "a")
        out_handle.write(net_handle.read())
        out_handle.close()
        net_handle.close()
        print("Saved")

def pyfasta():
    n =input("Write on how many files do you want to split the multifasta file: ")
    fasta_file = input("Write multifasta filename: ")
    subprocess.call('pyfasta split -n{} {}'.format(n, fasta_file), shell=True)
    
def bacphlip():
    
    fasta_file = input("Write genomes fasta filename: ")
    subprocess.call('nohup bacphlip -i {} -f --multi_fasta &'.format(fasta_file), shell=True)

def metadata():
    
    metadata_df = pd.read_csv(input("Write *.csv filename with metadata: "))
    bacphlip_df = pd.read_csv(input("Write BACPHLIP filename: "), sep = '\t',names = ['Accession', 'Virulent', 'Temperate'], header= None)
    delete_row = bacphlip_df[bacphlip_df["Virulent"]=='Virulent'].index
    bacphlip_df = bacphlip_df.drop(delete_row)
    df = (metadata_df.merge(bacphlip_df, on='Accession', how='inner'))
    df.to_csv('sequences.csv')
    
def fastaParser(infile):
    seqs = []
    headers = []
    with open(input("Write multifasta file: ")) as f:
        sequence = ""
        header = None
        for line in f:
            if line.startswith('>'):
                headers.append(line[1:-1])
                if header:
                    seqs.append([sequence])
                sequence = ""
                header = line[1:]
            else:
                sequence += line.rstrip()
        seqs.append([sequence])
    return headers, seqs

   
def count():
    
    fasta = sys.argv[1]
    headers, seqs = fastaParser(fasta)
    flat_seqs = [item for sublist in seqs for item in sublist]
    genomes = list(zip(headers, flat_seqs))
    metadata_df = pd.read_csv(input("Write *.csv filename with metadata: "))
    ncount = [seq.count('N') for header, seq in genomes]
    accession = re.findall(r'[A-Z]+\w?\d{5,8}\.\d', str(headers))
    df = pd.DataFrame({'Accession':accession, 'Ncount': ncount})
    new_df = df.merge(metadata_df,  how='right', on=['Accession'])
    new_df.to_csv(input('Write filename with new data: '))
    

def complete():
    
    fasta = sys.argv[1]
    headers, seqs = fastaParser(fasta)
    genomes = headers
    metadata_df = pd.read_csv(input("Write *.csv filename with metadata: "))
    complete = [c for c in genomes if "complete" in c]
    accession_complete = re.findall(r'[A-Z]+\w?\d{5,8}\.\d', str(complete))
    df_complete = pd.DataFrame({'Accession':accession_complete, 'Completness': 'complete'})
    partial = [p for p in genomes if 'partial' in p]
    accession_partial = re.findall(r'[A-Z]+\w?\d{5,8}\.\d', str(partial))
    df_partial = pd.DataFrame({'Accession':accession_partial, 'Completness': 'partial'})
    strings = ['complete', 'partial']
    unknown = [s for s in genomes if not any(x in s for x in strings)]
    accession_unknown = re.findall(r'[A-Z]+\w?\d{5,8}\.\d', str(unknown))
    df_unknown = pd.DataFrame({'Accession':accession_unknown, 'Completness': 'unknown'})
    frames = [df_complete, df_partial, df_unknown]
    result = pd.concat(frames)
    new_df = result.merge(metadata_df,  how='right', on=['Accession'])
    new_df.to_csv(input('Write filename with new data: '))


def length():
    
    fasta = sys.argv[1]
    headers, seqs = fastaParser(fasta)
    flat_seqs = [item for sublist in seqs for item in sublist]
    genomes = list(zip(headers, flat_seqs))
    length = [len(seq) for header, seq in genomes]
    metadata_df = pd.read_csv(input("Write *.csv filename with metadata: "))
    accession = re.findall(r'[A-Z]+\w?\d{5,8}\.\d', str(headers))
    df = pd.DataFrame({'Accession':accession, 'Length': length})
    new_df = df.merge(metadata_df,  how='right', on=['Accession'])
    new_df.to_csv(input('Write filename with new data: '))
        
                  
    
def phanotate():
    import subprocess
    filename = input("Write FASTA genomes filename: ")
    phanotate_directory = input("Write phanotate directory: ")
    phanotate_filename = input("Write Phanotate filename: ")
    subprocess.call('nohup {} {} >> {} &'.format(phanotate_directory, filename, phanotate_filename), shell=True)

def removed():
    
    url_list = []
    accession = []
    metadata_df = pd.read_csv(input("Write *.csv filename with metadata: "))
    with open(input("Write accession filename: "), 'r') as f:
        for lines in f:
            lines = f.read().splitlines()
            for line in lines:
                url = "https://www.ncbi.nlm.nih.gov/nuccore/" + line
                url_list.append(url)
                accession.append(line)

    the_word = 'Record removed'

    for url in url_list:
        r = requests.get(url, allow_redirects=False)
        soup = BeautifulSoup(r.content.lower(), 'lxml')
        words = soup.find_all(text=lambda text: text and the_word.lower() in text)
        count = len(words)
#         print('\nUrl: {}\ncontains {} of words: {}'.format(url, count, the_word))
        for elem in accession:
            if count > 0:
                df = pd.DataFrame({'Accession':elem, 'Removed': 'yes'}, index=[0])
            else:
                df1 = pd.DataFrame({'Accession':elem, 'Removed': 'no'}, index=[0])
    frames = [df, df1]
    result = pd.concat(frames)
    new_df = result.merge(metadata_df,  how='right', on=['Accession'])
    new_df.to_csv(input('Write filename with new data: '))


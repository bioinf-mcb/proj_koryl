import os
import pandas as pd
import requests
from bs4 import BeautifulSoup
import subprocess
from pandas import ExcelFile
from Bio import Entrez
from Bio import SeqIO
from Bio import GenBank
import re
from datetime import datetime
import argparse
import shutil
import numpy as np



def check(accessions, metadata):
    with open(accessions, 'r', encoding='UTF-8') as f:
        lines = f.readlines()
        lines = [line.rstrip() for line in lines]
    m = pd.read_csv(metadata)
    acc =  m['Accession'].tolist()
    set_accessions = (set(lines) == set(acc))
    if set_accessions == True:
        print('************************')
        print('The number of Accessions in the metadata and the list of accession file from NCBI Virus matches')
    else:
        print('_______________________')
        print('This is list of missed accessions in *.acc accession file: ',list(set(lines) - set(acc)))
        print('This is list of missed accessions in *.csv metadata file: ',list(set(acc) - set(lines)))

def acc_list(acc_path):
    with open(acc_path, 'r', encoding='UTF-8') as f:
        lines = f.readlines()
        lines = [line.rstrip() for line in lines]
    for i in range(0, len(lines), 7000):
        yield lines[i:i + 7000]
        
        
def fastaParser(infile):
    seqs = []
    headers = []
    with open(infile) as f:
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

        
def get_ncbi_genomes(args):
    
    accessions = args.input
    metadata = args.metadata
    directory = args.directory
    filename = args.output
    check(accessions, metadata)
    Entrez.email = args.email
    ictv = args.ictv
    genbank = args.genbank
    

    list_of_lists = acc_list(args.input)       
    folder = os.path.join(directory + 'database_'+ datetime.now().strftime('%Y-%m-%d'))
    os.makedirs(folder)
    old_directory = os.getcwd()
    new_directory = folder
    for file in os.listdir(old_directory):
        if file.endswith(".acc"):
            source1 = os.path.join(old_directory, accessions)
            destination1 = os.path.join(new_directory, accessions)
            shutil.copyfile(source1, destination1)
        elif file.endswith(".csv"):
            source2 = os.path.join(old_directory, metadata)
            destination2 = os.path.join(new_directory, metadata)
            shutil.copyfile(source2, destination2)
        elif file.endswith(".xlsx"):
            source3 = os.path.join(old_directory, ictv)
            destination3 = os.path.join(new_directory, ictv)
            shutil.copyfile(source3, destination3)
        else: 
            continue
    print('************************\nMoved files to new directory\n************************')

    for elem in list_of_lists:
        net_handle = Entrez.efetch(
            db="nucleotide", id=elem, rettype="fasta", retmode="text"
        )
        out_handle = open(os.path.join(new_directory, filename), "a+")
        out_handle.write(net_handle.read())
        out_handle.close()
        net_handle.close()
        print("Part of genomes saved")
    print('************************\nAll genomes saved\n************************')

    print('************************\nStarting processing \n************************')

    metadata_df = pd.read_csv(os.path.join(new_directory, metadata))

    metadata_accession = metadata_df['Accession'].tolist()


       
    headers, seqs = fastaParser(os.path.join(new_directory, filename))
    flat_seqs = [item for sublist in seqs for item in sublist]
    genomes = list(zip(headers, flat_seqs))
    ncount = [seq.count('N') for header, seq in genomes]
    accession = re.findall(r'[A-Z]+\w?\d{5,8}\.\d', str(headers))
    complete = [item[0] for item in genomes if 'complete' in item[0]]
    accession_complete = re.findall(r'[A-Z]+\w?\d{5,8}\.\d', str(complete))
    partial = [item[0] for item in genomes if 'partial' in item[0]]
    accession_partial = re.findall(r'[A-Z]+\w?\d{5,8}\.\d', str(partial))
    length = [len(seq) for header, seq in genomes]
    short_dict = {}
    for elem in length:
        if elem > 3000:
            short_dict = {k:'no' for k in accession}
        else:
            short_dict = {k:'yes' for k in accession}
            
            
    
    ####DICTIONARIES
    partial_d = {k:'partial' for k in accession_partial}
    complete_d = {k:'complete' for k in accession_complete}
    length_d = dict(zip(accession, length))
    ncount_d = dict(zip(accession, ncount))
#   d = {k:None for k in metadata_accession}
    updated_completness = complete_d.copy()
    updated_completness.update(partial_d)
    
    print('************************\nChecking all removed genomes\n************************')
    
    url_list = []
    
    for elem in metadata_accession:
        url = "https://www.ncbi.nlm.nih.gov/nuccore/" + elem
        url_list.append(url)
    
    the_word = 'Record removed'
    count_list = []
    
    for url in url_list:
        r = requests.get(url, allow_redirects=False)
        soup = BeautifulSoup(r.content.lower(), 'lxml')
        words = soup.find_all(text=lambda text: text and the_word.lower() in text)
        count = len(words) 
        count_list.append(count)
        print('\nUrl: {}\ncontains {} of words: {}'.format(url, count, the_word))
    #removed_df =  pd.DataFrame(list(zip(metadata_accession, count_list)), columns =['Accession', 'Removed'])
    #####REMOVED DICTIONARY
    removed_d = (dict(zip(metadata_accession, count_list)))
    
    print('************************\nAdding information to metadata\n************************')
    
    metadata_df['Nuc_length'] =  metadata_df['Accession'].map(length_d)
    metadata_df['Ncount'] =  metadata_df['Accession'].map(ncount_d)
    metadata_df['Removed'] =  metadata_df['Accession'].map(removed_d)
    metadata_df['Completness'] =  metadata_df['Accession'].map(updated_completness)
    metadata_df['Short_sequence'] = metadata_df['Accession'].map(short_dict)
    
    
    partial_df = metadata_df.loc[metadata_df['Completness'] == 'partial']
    partial_list = partial_df['Accession'].tolist()
    
    print('************************\nDownload genbank file and change information about completness\n************************')
    

    list_of_lists = partial_list

    

    for elem in list_of_lists:
        if len(elem) > 0:
            net_handle = Entrez.efetch(
                db="nucleotide", id=elem, rettype="gb", retmode="text"
            )
            out_handle = open(os.path.join(new_directory, genbank), "a+")
            out_handle.write(net_handle.read())
            out_handle.close()
            net_handle.close()
            print(elem, "Saved")


            with open(os.path.join(new_directory, genbank)) as handle:
                for record in GenBank.parse(handle):
                    if 'full length' in record.comment:
                        print(record.accession, "OK")
                        for col in metadata_df:
                            metadata_df['Completness'] = metadata_df['Completness'].replace('partial', 'complete')
        else:
            continue

                    

    metadata_df['accession_underscore'] = [x.split('.')[0] for x in metadata_df['Accession']]
    metadata_df['accession_underscore'] = metadata_df['accession_underscore'].astype(str)
    taxa = pd.ExcelFile(destination3)
    df = taxa.parse("VMRb36")
    new_df = df.loc[df['Host Source'] == 'bacteria']
    ictv_df = new_df[['Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Genome composition', 'Host Source', 
                 'Virus GENBANK accession', 'Virus REFSEQ accession']].copy()
    ictv_df = ictv_df.rename(columns = {'Species': 'Species_ictv', 'Family': 'Family_ictv', 'Genus': 'Genus_ictv', 'Realm' : 'Realm_ictv', 
                                        'Kingdom' : 'Kingdom_ictv', 'Phylum' : 'Phylum_ictv', 'Class' : 'Class_ictv', 'Order' : 'Order_ictv', 
                                        'Genome composition' : 'Genome_composition_ictv', 'Host Source' : 'Host_source_ictv', 
                                        'Virus REFSEQ accession':'Virus_REFSEQ_accession', 'Virus GENBANK accession':'Virus_GENBANK_accession'})
    df1 = metadata_df.merge(ictv_df, left_on='accession_underscore', right_on = 'Virus_REFSEQ_accession', how='left')
    metadata_df_all = df1.merge(ictv_df, left_on='accession_underscore', right_on = 'Virus_GENBANK_accession', how='left')
    
    metadata_df_all.Kingdom_ictv_x = metadata_df_all.Kingdom_ictv_x.combine_first(metadata_df_all.Kingdom_ictv_y)
    metadata_df_all.Phylum_ictv_x = metadata_df_all.Phylum_ictv_x.combine_first(metadata_df_all.Phylum_ictv_y)
    metadata_df_all.Class_ictv_x = metadata_df_all.Class_ictv_x.combine_first(metadata_df_all.Class_ictv_y)
    metadata_df_all.Order_ictv_x = metadata_df_all.Order_ictv_x.combine_first(metadata_df_all.Order_ictv_y)
    metadata_df_all.Family_ictv_x = metadata_df_all.Family_ictv_x.combine_first(metadata_df_all.Family_ictv_y)
    metadata_df_all['Genome_composition_ictv_x'] = metadata_df_all['Genome_composition_ictv_x'].combine_first(metadata_df_all['Genome_composition_ictv_y'])
    metadata_df_all['Host_source_ictv_x'] = metadata_df_all['Host_source_ictv_x'].combine_first(metadata_df_all['Host_source_ictv_y'])
    metadata_df_all['Virus_GENBANK_accession_x'] = metadata_df_all['Virus_GENBANK_accession_x'].combine_first(metadata_df_all['Virus_GENBANK_accession_y'])
    metadata_df_all['Virus_REFSEQ_accession_x'] = metadata_df_all['Virus_REFSEQ_accession_x'].combine_first(metadata_df_all['Virus_REFSEQ_accession_y'])
    metadata_df_all.Realm_ictv_x = metadata_df_all.Realm_ictv_x.combine_first(metadata_df_all.Realm_ictv_y)
    metadata_df_all.Genus_ictv_x = metadata_df_all.Genus_ictv_x.combine_first(metadata_df_all.Genus_ictv_y)
    metadata_df_all.Species_ictv_x = metadata_df_all.Species_ictv_x.combine_first(metadata_df_all.Species_ictv_y)
    
    metadata_df_all.drop([ 
             'Realm_ictv_y', 'Kingdom_ictv_y', 'Phylum_ictv_y', 'Class_ictv_y', 'Order_ictv_y', 'Genus_ictv_y', 
             'Genome_composition_ictv_y', 'Host_source_ictv_y', 'Virus_GENBANK_accession_y', 'Virus_REFSEQ_accession_y', 
             'Family_ictv_y', 'Species_ictv_y'], axis=1, inplace=True)
    metadata_df_all.columns = metadata_df_all.columns.str.rstrip('_x')
    
    df4 = metadata_df_all.groupby(['Accession'])
    dictio = {}
    for x,y in df4:
    
        if len(y['Genus_ictv'].unique()) > 1:
            y['Genus_ictv'] = np.where((len(y['Genus_ictv'].unique()) > 1), 
                                       y['Genus_ictv'].unique()[0]+'/'+y['Genus_ictv'].unique()[1], 
                                       y['Genus_ictv'])
            accession = list(y['Accession'])
            genus = list(y['Genus_ictv'])
            dictio = dict(zip(accession, genus))
        for key in dictio.keys():
            metadata_df_all["Genus_ictv"] = np.where((metadata_df_all["Accession"] == key), dictio.values(), metadata_df_all["Genus_ictv"])
    metadata_df_all = metadata_df_all.drop_duplicates(subset='Accession', keep = 'first')
    
    long_sequences = []  # Setup an empty list
    for record in SeqIO.parse(os.path.join(new_directory, filename), "fasta"):
        if len(record.seq) > 3000:
        # Add this record to our list
            long_sequences.append(record)
    SeqIO.write(long_sequences,os.path.join(new_directory, filename), "fasta")
    metadata_df_all = metadata_df_all.drop(metadata_df_all.loc[metadata_df_all['Length']<3000].index, inplace=False)
    metadata_df_all = metadata_df_all.reset_index()

    metadata_df_all.to_csv(destination2)
    

    



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Download phage database with metadata')
    parser.add_argument('-d', '--directory', help='Create path')
    parser.add_argument('-i', '--input', type=str, help='Accession file from NCBI Virus')
    parser.add_argument('-e', '--email', type=str, help='E-mail address to tell NCBI who you are')
    parser.add_argument('-o', '--output',  type=str, help='Directory for all genomes file')
    parser.add_argument('-c','--metadata', type=str, help='Metadata filename from NCBI Virus')
    parser.add_argument('-t', '--ictv',  type=str, help='ICTV filename')
    parser.add_argument('-g', '--genbank',  type=str, help='GenBank file to get information about completness')

    args = parser.parse_args()
    get_ncbi_genomes(args)




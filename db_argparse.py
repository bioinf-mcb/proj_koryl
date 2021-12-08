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


        
def get_ncbi_genomes(args):
    
    accessions = args.input
    metadata = args.metadata
    directory = args.directory
    filename = args.output
    check(accessions, metadata)
    Entrez.email = args.email
    ictv = 'ICTV.xlsx'

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
    metadata_df = pd.read_csv(destination2)
    taxa = pd.ExcelFile(destination3)
    df = taxa.parse("VMRb36")
    new_df = df.loc[df['Host Source'] == 'bacteria']
    metadata_df['Accession_underscore'] = [x.split('.')[0] for x in metadata_df['Accession']]
    metadata_df['Accession_underscore'] = metadata_df['Accession_underscore'].astype(str)
    new = new_df[['Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Genome composition', 
                  'Host Source', 'Virus GENBANK accession', 'Virus REFSEQ accession']].copy()
    df1 = metadata_df.merge(new, left_on='Accession_underscore', right_on = 'Virus REFSEQ accession', how='left')
    df2 = df1.merge(new, left_on='Accession_underscore', right_on = 'Virus GENBANK accession', how='left')
    df2.Kingdom_x = df2.Kingdom_x.combine_first(df2.Kingdom_y)
    df2.Phylum_x = df2.Phylum_x.combine_first(df2.Phylum_y)
    df2.Class_x = df2.Class_x.combine_first(df2.Class_y)
    df2.Order_x = df2.Order_x.combine_first(df2.Order_y)
    df2.Family_x = df2.Family_x.combine_first(df2.Family_y)
    df2['Genome composition_x'] = df2['Genome composition_x'].combine_first(df2['Genome composition_y'])
    df2['Host Source_x'] = df2['Host Source_x'].combine_first(df2['Host Source_y'])
    df2['Virus GENBANK accession_x'] = df2['Virus GENBANK accession_x'].combine_first(df2['Virus GENBANK accession_y'])
    df2['Virus REFSEQ accession_x'] = df2['Virus REFSEQ accession_x'].combine_first(df2['Virus REFSEQ accession_y'])
    df2.Realm_x = df2.Realm_x.combine_first(df2.Realm_y)
    df2.Genus_x = df2.Genus_x.combine_first(df2.Genus_y)
    df2.Species_x = df2.Species_x.combine_first(df2.Species_y)
    df2.drop([ 
              'Realm_y', 'Kingdom_y', 'Phylum_y', 'Class_y', 'Order_y', 'Genus_y', 
              'Genome composition_y', 'Host Source_y', 'Virus GENBANK accession_y', 
              'Virus REFSEQ accession_y', 
              'Family_y', 'Species_y'], axis=1, inplace=True)
    df2 = df2.rename(columns = {'Species': 'Species_ncbi', 'Family': 'Family_ncbi', 'Genus': 'Genus_ncbi', 'Realm_x' : 'Realm_ictv', 
                                'Kingdom_x' : 'Kingdom_ictv', 'Phylum_x' : 'Phylum_ictv', 'Class_x' : 'Class_ictv', 'Order_x' : 'Order_ictv', 
                                'Family_x' : 'Family_ictv', 'Genus_x' : 'Genus', 'Species_x' : 'Species', 
                                'Genome composition_x' : 'Genome composition_ictv', 'Host Source_x' : 'Host source_ictv', 
                                'Virus GENBANK accession_x' : 'Virus GENBANK accession', 'Virus REFSEQ accession_x' : 'Virus REFSEQ accession'})
    url_list = []
    accession = []
   
    with open(accessions, 'r') as f:
        for lines in f:
            lines = f.read().splitlines()
            for line in lines:
                url = "https://www.ncbi.nlm.nih.gov/nuccore/" + line
                url_list.append(url)
                accession.append(line)
    the_word = 'Record removed'
    count_list = []

    for url in url_list:
        r = requests.get(url, allow_redirects=False)
        soup = BeautifulSoup(r.content.lower(), 'lxml')
        words = soup.find_all(text=lambda text: text and the_word.lower() in text)
        count = len(words)
        count_list.append(count)
        print('\nUrl: {}\ncontains {} of words: {}'.format(url, count, the_word))
    for acc in accession:
        removed_df =  pd.DataFrame(list(zip(accession, count_list)), columns =['Accession', 'Removed'])
    
    records = SeqIO.parse(os.path.join(new_directory, filename), 'fasta')
    seq = [record for record in records]
    ncount = [seq[0].seq.count('N')]
    accession = []
    accession_complete = []
    accession_partial = []
    length = []
    data_short = {}
    data_long = {}
    for seq in records:
        accession = re.findall(r'[A-Z]+\w?\d{5,8}\.\d', str(seq.id))
        nuc = str(seq.seq)
        length = [len(nuc) for elem in nuc]
        complete = [c for c in seq.id if 'complete' in c]
        accession_complete = re.findall(r'[A-Z]+\w?\d{5,8}\.\d', str(complete))
        partial = [p for p in seq.id if not 'complete' in p]
        accession_partial = re.findall(r'[A-Z]+\w?\d{5,8}\.\d', str(partial))
      
        if len(nuc) < 3000:
            data_short = {'Accession': accession, 'Short_sequence': 'yes'}
        else:
            data_long = {'Accession': accession, 'Short_sequence': 'no'}
    dict_short = list(zip(data_short, data_long))
    df_processing_data_short = pd.DataFrame.from_dict(dict_short, orient='index')
    df_processing_data_short.transpose()
    data_partial = {'Accession': accession_partial, 'Completness': 'partial'}
    df_processing_data_partial = pd.DataFrame.from_dict(data_partial, orient='index')
    df_processing_data_partial.transpose() 
    data_complete =  {'Accession': accession_complete, 'Completness': 'complete'}
    df_processing_data_complete = pd.DataFrame.from_dict(data_complete, orient='index')
    df_processing_data_complete.transpose()
    data_length = {'Accession': accession, 'Length': length}
    df_processing_data_length = pd.DataFrame.from_dict(data_length, orient='index') 
    df_processing_data_length.transpose()
    data_ncount = {'Accession': accession, 'N_count':  ncount}
    df_processing_data_ncount = pd.DataFrame.from_dict(data_ncount, orient='index')
    df_processing_data_ncount.transpose()
    frames = [df_processing_data_partial, df_processing_data_complete, df_processing_data_length, df_processing_data_ncount, df2, removed_df, df_processing_data_short]
    all_df = pd.concat(frames)
    new_df = metadata_df.merge(all_df,  how='left', on=['Accession'])
    new_df.to_csv(destination2)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Download phage database with metadata')
    parser.add_argument('-d', '--directory', help='Create path')
    parser.add_argument('-i', '--input', type=str, help='Accession file from NCBI Virus')
    parser.add_argument('-e', '--email', type=str, help='E-mail address to tell NCBI who you are')
    parser.add_argument('-o', '--output',  type=str, help='Directory for all genomes file')
    parser.add_argument('-c','--metadata', type=str, help='Metadata filename from NCBI Virus')
    parser.add_argument('-f', '--ictv',  type=str, help='ICTV filename')

    args = parser.parse_args()
    get_ncbi_genomes(args)



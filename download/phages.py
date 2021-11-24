import os
import time
import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
import sys
import numpy as np
import requests
from bs4 import BeautifulSoup
import subprocess
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import GenBank
import re
import itertools


def directory():

    os.chdir(input("Write path to change directory: "))
    print('Current directory is: ', os.getcwd())


def folder():
    """
    This function create new folder in directory
    """
    parent_dir = input('Write parent dir: ')
    folder_name = input("Write new folder name: ")
    path = os.path.join(parent_dir, folder_name)
    os.mkdir(path)
    print("Directory '% s' created" % path)


def acc_list(acc_path):
    with open(acc_path, 'r', encoding='UTF-8') as f:
        lines = f.readlines()
        lines = [line.rstrip() for line in lines]
    for i in range(0, len(lines), 7000):
        yield lines[i:i + 7000]


def get_ncbi_genomes():


    Entrez.email = input('Write your e-mail: ')
    list_of_lists = acc_list(input("Write path to the accession file: "))
    filename = input("Write new filename (in FASTA format): ")

    for elem in list_of_lists:
        net_handle = Entrez.efetch(
            db="nucleotide", id=elem, rettype="fasta", retmode="text"
        )
        out_handle = open(filename, "a")
        out_handle.write(net_handle.read())
        out_handle.close()
        net_handle.close()
        print("Saved")


def bachplip_metadata():

    metadata_df = pd.read_csv(input("Write *.csv filename with metadata: "))
    bacphlip_df = pd.read_csv(input("Write BACPHLIP filename: "), sep='\t',
                              names=['Accession', 'Virulent', 'Temperate'],
                              header=None)
    delete_row = bacphlip_df[bacphlip_df["Virulent"]=='Virulent'].index
    bacphlip_df = bacphlip_df.drop(delete_row)
    df = (metadata_df.merge(bacphlip_df, on='Accession', how='inner'))
    df.to_csv(input('Write *.csv filename with new metadadata: '))


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


def count_n():

    fasta = sys.argv[1]
    headers, seqs = fastaParser(fasta)
    flat_seqs = [item for sublist in seqs for item in sublist]
    genomes = list(zip(headers, flat_seqs))
    metadata_df = pd.read_csv(input("Write *.csv filename with metadata: "))
    ncount = [seq.count('N') for header, seq in genomes]
    accession = re.findall(r'[A-Z]+\w?\d{5,8}\.\d', str(headers))
    df = pd.DataFrame({'Accession': accession, 'Ncount': ncount})
    new_df = df.merge(metadata_df, how='right', on=['Accession'])
    new_df.to_csv(input('Write *.csv filename with new metadadata: '))


def complete():

    fasta = sys.argv[1]
    headers, seqs = fastaParser(fasta)
    genomes = headers
    metadata_df = pd.read_csv(input("Write *.csv filename with metadata: "))
    complete = [c for c in genomes if "complete" in c]
    accession_complete = re.findall(r'[A-Z]+\w?\d{5,8}\.\d', str(complete))
    df_complete = pd.DataFrame({'Accession': accession_complete, 'Completness': 'complete'})
    partial = [p for p in genomes if not 'complete' in p]
    accession_partial = re.findall(r'[A-Z]+\w?\d{5,8}\.\d', str(partial))
    df_partial = pd.DataFrame({'Accession': accession_partial, 'Completness': 'partial'})
    frames = [df_complete, df_partial]
    result = pd.concat(frames)
    new_df = result.merge(metadata_df,  how='right', on=['Accession'])
    new_df.to_csv(input('Write *.csv filename with new metadadata: '))


def length():

    fasta = sys.argv[1]
    headers, seqs = fastaParser(fasta)
    flat_seqs = [item for sublist in seqs for item in sublist]
    genomes = list(zip(headers, flat_seqs))
    length = [len(seq) for header, seq in genomes]
    metadata_df = pd.read_csv(input("Write *.csv filename with metadata: "))
    accession = re.findall(r'[A-Z]+\w?\d{5,8}\.\d', str(headers))
    #     print(accession)
    df = pd.DataFrame({'Accession': accession, 'Length': length})
    new_df = df.merge(metadata_df, how='right', on=['Accession'])
    new_df.to_csv(input('Write *.csv filename with new metadadata: '))


def removed_check():

    url_list = []
    accession = []
    #
    with open(input("Write accession filename: "), 'r') as f:
        for lines in f:
            lines = f.read().splitlines()
            for line in lines:
                url = "https://www.ncbi.nlm.nih.gov/nuccore/" + line
                url_list.append(url)
                accession.append(line)

    the_word = 'Record removed'
    with open('remove.txt', 'a+') as fout:

        for url in url_list:
            r = requests.get(url, allow_redirects=False)
            soup = BeautifulSoup(r.content.lower(), 'lxml')
            words = soup.find_all(text=lambda text: text and the_word.lower() in text)
            count = len(words)
            print('\nUrl: {}\ncontains {} of words: {}'.format(url, count, the_word))
            fout.write("{} {} {}\n".format(url, '\t', count))


def remove():

    result = pd.read_csv('remove.txt', sep='\t', names=['Url', 'Remove'])
    result['Accession'] = [x.split('/')[4] for x in result['Url']]
    # result['Accession'] = result['Accession'].astype(str)
    # result['Remove'] = result['Remove'].astype(str)
    result.drop('Url', axis=1, inplace=True)
    result['Accession'] = result['Accession'].str.strip()
    result.to_csv('remove.txt')


def remove_to_csv():

    result = pd.read_csv('remove.txt')
    metadata_df = pd.read_csv(input('Write *.csv filename with metadadata: '))
    # result['Remove'].dtypes
    new_df = metadata_df.merge(result, how='left', on=['Accession'])
    # new_df = result.merge(metadata_df, on='Accession', how='outer')
    new_df.to_csv(input('Write *.csv filename with new metadadata: '))


def pyfasta():

    n = input("Write on how many files do you want to split the multifasta file: ")
    fasta_file = input("Write multifasta filename: ")
    subprocess.call('pyfasta split -n{} {}'.format(n, fasta_file), shell=True)


def partial():

    p = pd.read_csv(input('Write *.csv filename with metadata: '))
    lis = p.loc[p['Completness'] == 'partial']
    partial_list = lis['Accession'].tolist()
    return partial_list


def partial_complete():

    url_list = []
    partial_list = partial()
    for elem in partial_list:
        url = "https://www.ncbi.nlm.nih.gov/nuccore/" + elem
        url_list.append(url)
    the_word = 'full length'
    for url in url_list:
        r = requests.get(url, allow_redirects=False)
        soup = BeautifulSoup(r.content.lower(), 'lxml')
        words = soup.find_all(text=lambda text: text and the_word.lower() in text)
        count = len(words)
        print('\nUrl: {}\ncontains {} of words: {}'.format(url, count, the_word))


def download_genbank():

    Entrez.email = input('Write your e-mail: ')
    list_of_lists = partial()
    filename = input("Write new filename (in *.gb format): ")

    for elem in list_of_lists:
        net_handle = Entrez.efetch(
            db="nucleotide", id=elem, rettype="gb", retmode="text"
        )
        out_handle = open(filename, "a")
        out_handle.write(net_handle.read())
        out_handle.close()
        net_handle.close()
        print("Saved")


def partial_to_complete():

    df = pd.read_csv(input("Write *.csv filename with metadata: "))
    with open(input('Write *.gb filename: ')) as handle:
        for record in GenBank.parse(handle):
            if 'full length' in record.comment:
                print(record.accession, "OK")
                for col in df:
                    df['Completness'] = df['Completness'].replace('partial', 'complete')
    df.to_csv(input('Write *.csv filename with new metadadata: '))


def bacphlip():

    fasta_file = input("Write genomes fasta filename: ")
    subprocess.call('nohup bacphlip -i {} -f --multi_fasta &'.format(fasta_file), shell=True)


def phanotate():

    filename = input("Write FASTA genomes filename: ")
    phanotate_directory = input("Write phanotate directory: ")
    phanotate_filename = input("Write Phanotate filename: ")
    subprocess.call('nohup {} {} >> {} &'.format(phanotate_directory, filename, phanotate_filename), shell=True)


def remove_fasta():

    records = SeqIO.parse(open(input('Write *.fasta genomes filename: ')), 'fasta')
    with open(input('Write *.fasta filename: '), 'w') as f:
        for seq in records:
            nuc = str(seq.seq)
            if not len(nuc) < 3000:
                SeqIO.write([seq], f, "fasta")


def new_phanotate():

    record = SeqIO.parse(input('Write *.fasta genomes filename: '), "fasta")
    phanotate_list = []
    genomes_list = []

    with open(input('Write *.txt Phanotate filename: ')) as f:
        lines = f.read().split('\n')[:-1]
        for c, line in enumerate(lines):
            if line.startswith('#id'):
                phanotate_list.append(re.findall(r'[A-Z]+\w?\d{5,8}\.\d', line))
                flat_list = list(set(itertools.chain(*phanotate_list)))
    for elem in record:
        genomes_list.append(elem.id)
    not_same = [i for i in genomes_list if not i in flat_list]

    not_analyzed = []  # Setup an empty list
    for records in SeqIO.parse(input('Write *.fasta genomes filename'), "fasta"):
        if records.id in not_same:
            # Add this record to our list
            not_analyzed.append(records)

    print("Found %i sequences" % len(not_analyzed))

    SeqIO.write(not_analyzed, input('Write new *.fasta filename: '), "fasta")


def new_bacphlip():

    genomes_list = []
    df = pd.read_csv(input('Write *.bacphlip filename: '), sep='\t')
    bacphlip_list = df['Unnamed: 0'].to_list()
    record = SeqIO.parse(input('Write *.fasta genomes filename'), "fasta")
    for elem in record:
        genomes_list.append(elem.id)
    not_same =  [i for i in genomes_list if not i in bacphlip_list]

    not_analyzed = []  # Setup an empty list
    for records in record:
        if records.id in not_same:
            # Add this record to our list
            not_analyzed.append(records)

    print("Found %i sequences" % len(not_analyzed))

    SeqIO.write(not_analyzed, input('Write new *.fasta filename: '), "fasta")


def n_remove():

    records = SeqIO.parse(open(input('Write *.fasta genomes filename')), 'fasta')
    with open(input('Write new *.fasta genomes filename: '), 'w') as f:
        for seq in records:
            nuc = str(seq.seq)
            if nuc.count('N') < 100:
                print('N count: ', nuc.count('N'), 'ID: ', seq.id)
                SeqIO.write([seq], f, "fasta")


def ictv_to_df():

    ictv = pd.ExcelFile(input('Write ICTV *.xlsc filename: '))
    df = ictv.parse("VMRb36")
    new_df = df.loc[df['Host Source'] == 'bacteria']
    new_df.to_csv('ICTV.csv')


def ictv_to_csv():

    pp = pd.read_csv(input("Write *.csv filename with metadata: "))
    pp['Accession_noversion'] = [x.split('.')[0] for x in pp['Accession']]
    pp['Accession_noversion'] = pp['Accession_noversion'].astype(str)
    i = pd.read_csv('ICTV.csv')
    new = i[['Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Genome composition',
             'Host Source', 'Virus GENBANK accession', 'Virus REFSEQ accession']].copy()
    df1 = pp.merge(new, left_on=['Accession_noversion'], right_on=['Virus GENBANK accession'],
                   how='left', suffixes=('', '_y', '_x'))
    df1.drop(df1.filter(regex='_y$').columns.tolist(), axis=1, inplace=True)
    df2 = df1.merge(new, left_on=['Accession_noversion'], right_on=['Virus REFSEQ accession'],
                    how='right', suffixes=('', '_y'))
    df2.to_csv(input('Write *.csv filename with new metadadata: '))



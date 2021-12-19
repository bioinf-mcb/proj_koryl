
######################### IMPORTS
import os
import pandas as pd
import requests
from bs4 import BeautifulSoup
from pandas import ExcelFile
from Bio import Entrez
from Bio import SeqIO
from Bio import GenBank
import re
from datetime import datetime
import argparse
import shutil
import numpy as np
#####################################


def check(accessions, metadata):
    '''
    
    

    Parameters
    ----------
    accessions : *.acc filename
        Accession file from NCBI Virus
    metadata : *.csv
        Metadata file from NCBI Virus

    -------
    Help function to check if the number of accessions 
    of a file with extension *.acc is the 
    same as the number of rows in the metadata file

    '''
    with open(accessions, 'r', encoding='UTF-8') as f: #open file with accessions ids
        lines = f.readlines() #return a list with each line
        lines = [line.rstrip() for line in lines] #remove any unexpected characters (like space on the end of line)
    m = pd.read_csv(metadata) #open metadata dataframe
    acc =  m['Accession'].tolist() #Get all accessions ids from "Accession" column to list
    set_accessions = (set(lines) == set(acc)) #check if the number of accessions from *.acc file is the same as number of accessions in metadata
    if set_accessions == True: #if lists are same
        print('************************')
        print('The number of Accessions in the metadata and the list of accession file from NCBI Virus matches')
    else: #if lists are not same
        print('_______________________')
        print('This is list of missed accessions in *.acc accession file: ',list(set(lines) - set(acc))) #print what accession is additional in the metadata 
        print('This is list of missed accessions in *.csv metadata file: ',list(set(acc) - set(lines))) #print what accession is additional in the *.acc accession file

def acc_list(acc):
    '''
    Parameters
    ----------
    acc : 
    ------
     A helper to split the accession list as Entrez 
     can download to 10,000 records.
    '''
    with open(acc, 'r', encoding='UTF-8') as f: #open accession file
        lines = f.readlines()
        lines = [line.rstrip() for line in lines]
    for i in range(0, len(lines), 7000):
        yield lines[i:i + 7000] #split every 7000 elements from list
        
        
def fastaParser(infile):
    '''
    ----------
        --output
    Auxiliary function. It takes the file of 
    fasta genomes and reads headlines 
    and sequences    

    Returns
    -------
    '''
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
    '''
    Main function. It takes genomes and adds 
    information to metadata such as sequence 
    length, completeness, and information 
    from ICTV.
    
    '''
    ####### Name variables for ease of use ######
    accessions = args.input
    metadata = args.metadata
    directory = args.directory
    filename = args.output
    check(accessions, metadata)
    Entrez.email = args.email
    ictv = args.ictv
    genbank = args.genbank
    
    ##############################################
    
    list_of_lists = acc_list(args.input) #Load split accession items
    folder = os.path.join(directory + 'database_'+ datetime.now().strftime('%Y-%m-%d')) #create PATH to new folder in path with all databases
    os.makedirs(folder) #create this folder
    old_directory = os.getcwd() #check your current directory
    new_directory = folder 
    for file in os.listdir(old_directory): #copy files to new directory
        if file.endswith(".acc"):
            source_acc = os.path.join(old_directory, accessions)
            destination_acc = os.path.join(new_directory, accessions)
            shutil.copyfile(source_acc, destination_acc)
        elif file.endswith(".csv"):
            source_csv = os.path.join(old_directory, metadata)
            destination_csv = os.path.join(new_directory, metadata)
            shutil.copyfile(source_csv, destination_csv)
        elif file.endswith(".xlsx"):
            source_ictv = os.path.join(old_directory, ictv)
            destination_ictv = os.path.join(new_directory, ictv)
            shutil.copyfile(source_ictv, destination_ictv)
        else: 
            continue
    print('************************\nMoved files to new directory\n************************')

    for elem in list_of_lists: #for every 7000 genomes
        net_handle = Entrez.efetch(
            db="nucleotide", id=elem, rettype="fasta", retmode="text"
        )
        out_handle = open(os.path.join(new_directory, filename), "a+") #open new file
        out_handle.write(net_handle.read())
        out_handle.close()
        net_handle.close()
        print("Part of genomes saved")
    print('All genomes saved\n************************')

    print('Starting processing \n************************')

    metadata_df = pd.read_csv(os.path.join(new_directory, metadata)) #read dataframe with metadata

    metadata_accession = metadata_df['Accession'].tolist() #save all rows from metadata column "Accession" to list


       
    headers, seqs = fastaParser(os.path.join(new_directory, filename)) #read headers and sequences from multigenome fasta file
    flat_seqs = [item for sublist in seqs for item in sublist]
    genomes = list(zip(headers, flat_seqs))
    ncount = [seq.count('N') for header, seq in genomes] #count how many "N"'s in every sequence
    accession = re.findall(r'[A-Z]+\w?\d{5,8}\.\d', str(headers)) #find accessions ids in headers 
    complete = [item[0] for item in genomes if 'complete' in item[0]] #make list with headers if "complete" is in element
    accession_complete = re.findall(r'[A-Z]+\w?\d{5,8}\.\d', str(complete)) #find accessions ids in headers 
    partial = [item[0] for item in genomes if 'partial' in item[0]] #create list with headers if "partial" is in element
    accession_partial = re.findall(r'[A-Z]+\w?\d{5,8}\.\d', str(partial)) #find accessions ids in headers 
    length = [len(seq) for header, seq in genomes] #make a list with length of sequences
            
    #### DICTIONARIES -> They'll be add to metadata dataframe
    partial_d = {k:'partial' for k in accession_partial}
    complete_d = {k:'complete' for k in accession_complete}
    length_d = dict(zip(accession, length))
    ncount_d = dict(zip(accession, ncount))
    updated_completness = complete_d.copy()
    updated_completness.update(partial_d)
    
    print('************************\nChecking all removed genomes')
    
    url_list = [] #empty list to contain urls
    
    for elem in metadata_accession: #for every accession record
        url = "https://www.ncbi.nlm.nih.gov/nuccore/" + elem #create path to website
        url_list.append(url) #add url path to list
    
    the_word = 'Record removed' 
    count_list = [] #empty list to contain information if record was removed from the NCBI site
    
    for url in url_list:
        r = requests.get(url, allow_redirects=False)
        soup = BeautifulSoup(r.content.lower(), 'lxml')
        words = soup.find_all(text=lambda text: text and the_word.lower() in text) #find string from the_word value
        count = len(words) #how many strings from the_word was on the NCBI site
        count_list.append(count) #add count value to list
        print('\nUrl: {}\ncontains {} of words: {}'.format(url, count, the_word))
    #removed_df =  pd.DataFrame(list(zip(metadata_accession, count_list)), columns =['Accession', 'Removed'])
    #####REMOVED DICTIONARY
    removed_d = (dict(zip(metadata_accession, count_list))) #dictionary with accessions and information if record was removed 
    
    print('************************\nAdding information to metadata')
    ############ ADD ALL DICTIONARIES TO METADATA
    metadata_df['Nuc_length'] =  metadata_df['Accession'].map(length_d)
    metadata_df['Ncount'] =  metadata_df['Accession'].map(ncount_d)
    metadata_df['Removed'] =  metadata_df['Accession'].map(removed_d)
    metadata_df['Completness'] =  metadata_df['Accession'].map(updated_completness)

    
    
    partial_df = metadata_df.loc[metadata_df['Completness'] == 'partial'] #check if there is any partial genome in data
    partial_list = partial_df['Accession'].tolist()
    
    print('************************\nDownload genbank file and change information about completness')
    

    list_of_lists = partial_list

    
    ##### IF THERE IS ANY "PARTIAL" GENOME DOWNLOAD GENBANK FILE AND CHECK IF IT'S GOT "full length" in comment section.
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

    print("Adding ICTV to metadata")                

    metadata_df['accession_underscore'] = [x.split('.')[0] for x in metadata_df['Accession']] #make a list of accession without version
    metadata_df['accession_underscore'] = metadata_df['accession_underscore'].astype(str)
    taxa = pd.ExcelFile(destination_ictv) #read Excel file from the ICTV site
    df = taxa.parse("VMRb36") 
    new_df = df.loc[df['Host Source'] == 'bacteria'] #create dataframe only with data where Host_source is bacteria
    ictv_df = new_df[['Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Genome composition', 'Host Source', 
                 'Virus GENBANK accession', 'Virus REFSEQ accession']].copy() #this columns will be add to metadata
    ictv_df = ictv_df.rename(columns = {'Species': 'Species_ictv', 'Family': 'Family_ictv', 'Genus': 'Genus_ictv', 'Realm' : 'Realm_ictv', 
                                        'Kingdom' : 'Kingdom_ictv', 'Phylum' : 'Phylum_ictv', 'Class' : 'Class_ictv', 'Order' : 'Order_ictv', 
                                        'Genome composition' : 'Genome_composition_ictv', 'Host Source' : 'Host_source_ictv', 
                                        'Virus REFSEQ accession':'Virus_REFSEQ_accession', 'Virus GENBANK accession':'Virus_GENBANK_accession'}) #rename columns to be sure later if it's ictv data
    #### merge dataframes to add ictv dataframe
    df1 = metadata_df.merge(ictv_df, left_on='accession_underscore', right_on = 'Virus_REFSEQ_accession', how='left')
    metadata_df_all = df1.merge(ictv_df, left_on='accession_underscore', right_on = 'Virus_GENBANK_accession', how='left')
    ##############################################################################################################
    ############## extract information from columns to have all ICTV data in one
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
    metadata_df_all = metadata_df_all.rename(columns = {'Species': 'Species_ncbi', 'Family': 'Family_ncbi', 'Genus': 'Genus_ncbi'})
    
    metadata_df_all.drop([ 
             'Realm_ictv_y', 'Kingdom_ictv_y', 'Phylum_ictv_y', 'Class_ictv_y', 'Order_ictv_y', 'Genus_ictv_y', 
             'Genome_composition_ictv_y', 'Host_source_ictv_y', 'Virus_GENBANK_accession_y', 'Virus_REFSEQ_accession_y', 
             'Family_ictv_y', 'Species_ictv_y'], axis=1, inplace=True) #drop columns with "_y" suffixes 
    metadata_df_all.columns = metadata_df_all.columns.str.rstrip('_x') #drop suffix "_x" 
    ####### IN ICTV THERE ARE SOME DUPLICATES
    df4 = metadata_df_all.groupby(['Accession']) #group metadata dataframe 
    dictio = {} #empty dictionary to contain information about duplicates
    for x,y in df4:
        if len(y['Genus_ictv'].unique()) > 1: #some duplicates have different Genus_ictv values
            y['Genus_ictv'] = np.where((len(y['Genus_ictv'].unique()) > 1), 
                                       y['Genus_ictv'].unique()[0]+'/'+y['Genus_ictv'].unique()[1], 
                                       y['Genus_ictv']) #change this duplicates with adding both information into one row
            accession = list(y['Accession'])
            genus = list(y['Genus_ictv'])
            dictio = dict(zip(accession, genus))
        for key in dictio.keys():
            metadata_df_all["Genus_ictv"] = np.where((metadata_df_all["Accession"] == key), dictio.values(), metadata_df_all["Genus_ictv"]) #change duplicated Genus_ictv value in metadata dataframe
    metadata_df_all = metadata_df_all.drop_duplicates(subset='Accession', keep = 'first') #drop duplicates and save first values
    ############### DELETE SHORT SEQUENCES FROM MULTIGENOME FASTA FILE AND FROM METADATA DATAFRAME
    long_sequences = []  # Setup an empty list
    for record in SeqIO.parse(os.path.join(new_directory, filename), "fasta"):
        if len(record.seq) > 3000:
        # Add this record to our list
            long_sequences.append(record)
    SeqIO.write(long_sequences,os.path.join(new_directory, filename), "fasta")
    metadata_df_all = metadata_df_all.drop(metadata_df_all.loc[metadata_df_all['Length']<3000].index, inplace=False)
    metadata_df_all = metadata_df_all.reset_index()

    metadata_df_all.to_csv(destination_csv) #save metadata
    

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Command-line program for retrieving bacteriophage genomes from NCBI Virus and metadata. Additionally, adding information from ICTV and information about completeness, sequence length and checking that all genomes are still available on the NCBI Virus website')
    parser.add_argument('-d', '--directory', help='Enter the path to the database folder')
    parser.add_argument('-i', '--input', type=str, help='Enter the file name in *.acc format. Download from NCBI Virus website')
    parser.add_argument('-e', '--email', type=str, help='NCBI requires your e-mail address every time to know who you are')
    parser.add_argument('-o', '--output',  type=str, help='Enter the file name in *.fasta format to retrieve bacteriophage genomes')
    parser.add_argument('-c','--metadata', type=str, help='Enter the  metadata file name in *.csv format. Download from NCBI Virus website')
    parser.add_argument('-t', '--ictv',  type=str, help='Enter the name of the ICTV file in the *.xlsx format. Download this file from the ICTV website.')
    parser.add_argument('-g', '--genbank',  type=str, help='Enter the file name in * .gb format. Genbank records of incomplete genomes will be kept here to check if they are full sequence length')

    args = parser.parse_args()
    get_ncbi_genomes(args)




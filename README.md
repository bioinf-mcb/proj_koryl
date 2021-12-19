## ABOUT PHAGE DATABASE PIPELINE

Command-line program for retrieving bacteriophage genomes from NCBI Virus and metadata. Additionally, adding information from ICTV and information about completeness, sequence length and checking that all genomes are still available on the NCBI Virus website.
The program was written in Python 3.8.3.

## NECESSARY INSTALLATIONS

#### PANDAS:

```
conda install pandas
```
#### BS4:
```
conda install -c conda-forge bs4
```
or
```
conda install -c conda-forge/label/cf202003 bs4
```
#### BIOPYTHON:
```
conda install -c anaconda biopython
```
#### ARGPARSE:
```
conda install -c conda-forge argparse
```
#### REGEX:
Download one of the following:
```
conda install -c conda-forge regex
```
```
conda install -c conda-forge/label/gcc7 regex
```
```
conda install -c conda-forge/label/cf201901 regex
```
```
conda install -c conda-forge/label/cf202003 regex
```
#### SHUTIL:
```
pip install shutil
```
#### NUMPY:
```
pip install numpy
```
#### REQUESTS:
```
conda install -c anaconda requests
```

## BEFORE RUNNING PROGRAM

Before you run the program, go to the website: 
https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Bacteriophage,%20all%20taxids&Completeness_s=complete
On the NCBI website click "Download" at the top. Then select in the column "Accession" -> "Nucleotide". Press "Next" twice. In the third step, select "Accession with version". And then download file with all accessions ids.

Additionally, press "Download" again at the top of the page. Now select the "Current table view result" column and press "CSV format". Press "Next" twice. In the third step, press "Select all" and "Accession with version". Now you can download the metadata files.

### Now you should upload these files to your virtual machine.

## HOW TO RUN PROGRAM

Download this repository.
If you don't have git downloaded on your virtual machine, please do the following:

```
sudo apt install git-all
```
Change the current working directory to the location where you want the cloned directory.

Then clone repository:
```
git clone https://github.com/bioinf-mcb/proj_koryl.git
```
Change the path to the proj_koryl folder

```
cd proj_koryl
```

Change the permission to access python file
```
chmod u+x manage_db.py
```

## RUN PROGRAM

To check the arguments you need to run the program, you can type

```
python manage_db.py -h
```

```
'-d', '--directory', help='Enter the path to the database folder'
'-i', '--input', help='Enter the file name in *.acc format. Download from NCBI Virus website'
'-e', '--email', help='NCBI requires your e-mail address every time to know who you are'
'-o', '--output', help='Enter the file name in *.fasta format to retrieve bacteriophage genomes'
'-c','--metadata', help='Enter the  metadata file name in *.csv format. Download from NCBI Virus website'
'-t', '--ictv', help='Enter the name of the ICTV file in the *.xlsx format. Download this file from the ICTV website.'
'-g', '--genbank', help='Enter the file name in * .gb format. Genbank records of incomplete genomes will be kept here to check if they are full sequence length'
```

To run the program, for example, you should enter:
```
python manage_db.py -d /home/MCB/your_login_name/path_to_database/ -i sequences.csv -e your_email@domain.com -c sequences.csv -t ICTV.xlsx -g partial.gb -o genomes.fasta
```









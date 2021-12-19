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









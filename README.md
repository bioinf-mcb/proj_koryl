## ABOUT PHAGE DATABASE PIPELINE

The Phage Database pipeline is used to retrieve bacteriophage (phage) genomes. 
You need to download the GenBank accession or RefSeq assembly .csv file from the NCBI Virus website 
(https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Bacteriophage,%20all%20taxids). 
With this list you can download bacteriophage genomes from the NCBI FTP server. 
Then, using external command-line programs, the genomes can be annotated.

## NECESSARY INSTALLATIONS

#### PANDAS:

pip install pandas
conda install pandas

#### SAMTOOLS:
https://github.com/samtools/samtools

```
conda install -c conda-forge -c bioconda samtools=1.12
```

or
```
conda create -n samtools112 -c conda-forge -c bioconda samtools=1.12
```

#### PYFASTA

```
conda install -c bioconda pyfasta
````
or
```
conda install -c bioconda/label/cf201901 pyfasta
```


#### PHANOTATE:
https://github.com/deprekate/PHANOTATE
```
pip3 install phanotate
```
or
```
git clone https://github.com/deprekate/PHANOTATE.git
cd PHANOTATE
python3 setup.py install
```

## PIPELINE EXECUTION

Pipeline is written in Python 3.8 and the individual functions need to be imported.
```
import download_genomes
from download_genomes import create_folder, assembly, download_summary, orf, phage_list, phanotate, samtools, inventory
```

#### Create folders in database:
```
download_genomes.create_folder.folder(directory,folder's name)
```

#### Download assembly_summary_refseq/assembly_summary genbank

```
download_genomes.download_summary.download_refseq_summary(directory/filename)
download_genomes.download_summary.download_genbank_summary(directory/filename)
```

#### Create list of accesion/assembly ids

```
download_genomes.phage_list.create_refseq_list('NCBI Virus *.csv file directory, directory/filename)
download_genomes.phage_list.create_genbank_list('NCBI Virus *.csv file directory, directory/filename)
```
#### Download genomes from FTP NCBI server

```
download_genomes.assembly.download_refseq/download_genbank(folder path where genomes will be downloaded, assembly_summary_refseq/assembly_summary genbank file directory, file path with list of accesions/assemblies)
```

#### SAMTOOLS

```
download_genomes.samtools.to_one_file(downloaded genomes directory)
```
### PYFASTA
Split multiple fasta file into more smaller files to run Phanotate on a few cores
```
pyfasta split -n8 genomes.fasta
```
### PHANOTATE

In terminal create new file:
run_phanotate.sh

```
#!/bin/bash
 
input="${1}"   
INPATH= write  inpath directory
OUTPATH= write outpath directory

/home/*/PHANOTATE/phanotate.py "${INPATH}/genomes.${input}.fasta" >> "${OUTPATH}/Phanotate.${input}.txt" 
echo "${INPUT} done."
```
##### Run PHANOTATE
```
echo 1 2 3 4 | xargs -n1 -d ' ' -P4 -I{} ./tmp.sh {}
```
###### echo 
recalls the number of the new fast file created with pyfasta
###### -P
recalls the number of cores the phanotate should run on

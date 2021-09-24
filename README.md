## ABOUT PHAGE DATABASE PIPELINE

The Phage Database pipeline is used to retrieve bacteriophage (phage) genomes. 
It includes Entrez Python library, Bachplip tool and Phanotate tool.
Firstly, it is necessary to download file with accessions of interest (*.acc file) and file with metadata (*.csv)

## NECESSARY INSTALLATIONS

#### PANDAS:
```
pip install pandas
```
or 
```
conda install pandas
```

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
from download_genomes import folder, directory, acc_list, get_ncbi_genomes, pyfasta, bacphlip, metadata
```

#### Create folders in database:
```
download_genomes.folder()
```
This definition will ask you to write parent dir and name of a new folder

#### Change current directory:
```
download_genomes.directory()
```
This definition will ask you to write path to change directory

#### Download assembly_summary_refseq/assembly_summary genbank

```
download_genomes.download_summary.download_refseq_summary(directory/filename)
download_genomes.download_summary.download_genbank_summary(directory/filename)
```

#### Download genomes:
```
download_genomes.get_ncbi_genomes()
```
This definition will ask you to write you e-mail adress (it is necessary to tell NCBI who you are), filename with all accessions and genomes filename. 

### PYFASTA
Split multiple fasta file into more smaller files to run Phanotate on a few cores
In terminal:
```
pyfasta split -n8 genomes.fasta
```
or in Python:
```
download_genomes.pyfasta()
```
The definition will ask you to write on how many files do you want to split the multifasta file and genomes filename.


#### BACPHLIP
Bacphlip is a tool to predict if phage is virulent or temperate
Write in terminal: 
```
bacplip -i multifasta_file --multi_fasta
```
or in Python:
```
download_genomes.bacphlip()
```
The definition will ask you to write multifasta filename
###### Add bacphlip virulent prediction to metadata
```
download_genomes.metadata()
```

#### PHANOTATE

In terminal create new file:
phanotate.sh

```
#!/bin/bash
 
input="${1}"   
INPATH= write  inpath directory
OUTPATH= write outpath directory

{write path to the phanotate.py file}/PHANOTATE/phanotate.py "${INPATH}/genomes.${input}.fasta" >> "${OUTPATH}/Phanotate.${input}.txt" 
echo "${INPUT} done."
```
##### Run PHANOTATE
```
chmod u+x phanotate.sh
nohup echo -n 1 2 3 4 | xargs -n1 -d ' ' -P4 -I{} ./phanotate.sh {} &
```
###### echo 
recalls the number of the new fast file created with pyfasta
###### -P
recalls the number of cores the phanotate should run on

import os
import subprocess
def to_one_file(path):
  ''' this function changes directory and with samtools create one fasta file with all downloaded genomes and compressed *.bgzf file, and *.fai file'''
    os.chdir(path)
    print("Current Working Directory: " , os.getcwd())
    subprocess.Popen('gunzip GCF*.fna.gz; cat GCF*.fna>> genomy.fasta; rm GCF*.fna; bgzip -c genomy.fasta >> genomy.fasta.bgzf; samtools faidx genomy.fasta', shell=True, stderr=subprocess.PIPE)

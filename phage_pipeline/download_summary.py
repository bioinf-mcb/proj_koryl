import urllib.request
def download_summary():
  '''function to download assembly_summary_refseq from NCBI FTP server'''
  urllib.request.urlretrieve("https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt", "assembly_summary_refseq.txt")

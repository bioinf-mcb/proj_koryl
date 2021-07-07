import subprocess, sys

def phanotate(cmd):
    ''' function to annotate genomes where open reading frames are beneficial paths, while gaps and overlaps are penalized paths'''
    p = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
    out = p.stderr.read(1)
    sys.stdout.write(out)
    sys.stdout.flush()

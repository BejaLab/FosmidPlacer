
import ncbi_genome_download as ngd
import subprocess as sp
import re
from os.path import join, isfile
from os import listdir
from glob import glob

fna_re = re.compile('.+_genomic[.]fna')
faa_re = re.compile('.+[.]faa')

accession = snakemake.input["accession"]

dirname_list = glob(f'data/*/*/{accession}')
assert len(dirname_list) == 1, 'Expected one directory, got: ' + len(dirname_list)
dirname = dirname_list.pop()

gz_files = glob(join(dirname, '*.gz'))
if gz_files: sp.run([ 'gunzip' ] + gz_files)

files = listdir(dirname)
fna_list = [ f for f in files if isfile(join(dirname, f)) and fna_re.search(f)]
assert len(fna_list) == 1, 'Expected one *_genomic.fna file, got: ' + len(fna_list)
fna = fna_list.pop()

faa = [ f for f in files if isfile(join(dirname, f)) and faa_re.search(f)]

if not faa:
    sp.call([ 'gmsn.pl', '--format', 'GFF', '--faa', '--prok', fna[0] ], cwd = dirname)

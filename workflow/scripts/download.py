
from ncbi_genome_download import download
from subprocess import run
from glob import glob

accession_dict = { 'genbank': [], 'refseq': [] }
with open(input[0]) as fp:
    for accession in fp:
        accession = accession.strip()
            prefix, skip = accession.split('_', 1)
            section = 'refseq' if prefix == 'GCF' else 'genbank'
            if not glob(f'data/{section}/*/{accession}/*.fna*'):
                accession_dict[section].append(accession)

for section, accessions in accession_dict.items():
    if accessions:
        download(section = section, file_formats = 'fasta,protein-fasta', assembly_accessions = accessions, output = 'data', parallel = threads)
    gzip_files = glob(f'data/{section}/*/*/*.gz')
    if gzip_files:
        run([ 'gunzip' ] + gzip_files)

with open(output[0], 'w') as fp:
        fp.write('OK')

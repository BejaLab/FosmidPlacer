
threads=$1

parallel -j"$threads" \
	cd "data/*/*/{}/" \;
       	[ -s "*.faa" ] '||' gmsn.pl --format GFF --prok --faa --output *_genomic.fna \
	:::: data/gtdb_{genbank,refseq}.txt
echo OK > data/predict.ok

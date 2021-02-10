
threads=$1

mkdir -p data/orthologs

parallel usearch -db {1} -ublast {2} -evalue 1e-10 -accel .4 -threads 1 -blast6out {2//}/{1/.}.ublast ::: data/samples/*.faa ::: data/{genbank,refseq}/*/*/*.faa
parallel --colsep _ sed -E 's/\\s.+//' 'data/*/*/{2}_{3}/*.ublast' \| sort -u \| xargs samtools faidx 'data/*/*/{2}_{3}/*.faa' \> data/orthologs/{1}_{2}_{3}.faa :::: data/gtdb_list.txt

echo OK > data/usearch.ok

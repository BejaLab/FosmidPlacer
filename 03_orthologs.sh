mkdir -p orthologs
parallel sed -E 's/\\s.+//' blast/*-{/.}.ublast \| sort -u \| xargs samtools faidx {} \> orthologs/{/} ::: gtdbtk/*.faa
(cd orthologs; ln -fs ../contigs/*.faa ./)
proteinortho -p=ublast -conn=1 -alpha=0.8 -project=ublast orthologs/*.faa

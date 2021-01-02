mkdir -p blast
parallel usearch -makeudb_ublast {} -output {}.udb ::: contigs/*.faa
parallel [ -f blast/{/.1}-{/.2}.ublast ] '||' usearch -db {1}.udb -ublast {2} -evalue 1e-10 -accel .4 -threads 1 -blast6out blast/{1/.}-{2/.}.ublast ::: contigs/*.faa ::: gtdbtk/*.faa

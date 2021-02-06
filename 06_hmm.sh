find gtdbtk contigs           -name '*.faa' | parallel hmmsearch -o /dev/null --tblout {2.}+{1/.}.tblout -A {2.}+{1/.}.stk {1} {2} ::: ../hmm/*.hmm :::: -
find gtdbtk contigs -size +1c -name '*.stk' | parallel esl-reformat -o {.}.fasta fasta {}
find gtdbtk contigs -size +1c -name '*.stk' | parallel --colsep [+] usearch -db ../hmm/{2.}.faa.udb -ublast {1}+{2.}.fasta -id 0.5 -evalue 1e-10 -accel .4 -threads 1 -blast6out {1}+{2.}.ublast

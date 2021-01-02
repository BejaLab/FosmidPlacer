
source taxa.cfg

parallel mkdir -p phylogeny/{/.} \; xargs samtools faidx {} \< gtdbtk.txt \> phylogeny/{/.}/gene.trim_w_gaps ::: "$GTDBTK_DATA_PATH/$domain"_msa_marker_genes_all/*.faa 

parallel mafft --reorder --localpair --maxiterate 1000 {} \> {.}.mafft ::: phylogeny/*/gene.faa
parallel trimal -in {} -out {.}.trim -automated1                       ::: phylogeny/*/gene.mafft
parallel trimal -in {} -out {.}.trim -noallgaps                        ::: phylogeny/*/gene.trim_w_gaps
parallel [ -s {}.treefile ] '||' iqtree -bb 1000 -nt 1 -s {}           ::: phylogeny/*/gene.trim

run_treeshrink.py -i phylogeny -t gene.trim.treefile -a gene.trim -O treeshrink --force
cat phylogeny/*/treeshrink.treefile > treeshrink.trees
astral-constrained -t 2 -i treeshrink.trees -o astral-constrained.tree -j clade.reduced.tree
astral             -t 2 -i treeshrink.trees -o astral.tree

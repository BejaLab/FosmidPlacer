# FosmidPlacer

Phylogenetic placement for environmental DNA fragments.

## Input

As input (in folder `input`) the pipeline takes a set of fasta files with protein sequences from the environmental clones. Also required is a full installation of the [GTDB database](https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/auxillary_files/gtdbtk_r95_data.tar.gz). On top of the default GTDB installation, we need the marker alignments (for [Archaea](https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/genomic_files_all/ar122_msa_marker_genes_all_r95.tar.gz) and for [Bacteria](https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/genomic_files_all/bac120_msa_marker_genes_all_r95.tar.gz)) that should be extracted in the root directory of GTDB.

Configuration file in `config/config.json` should provide the following required metadata:

* `taxon` named taxon in GTDB classification to which the fosmids are known to belong, the lower the rank the better
* `gtdb_path` path to the GTDB database installation
* `gtdb_prefix` either `bac120` (Bacteria) or `ar122` (Archaea), or some other prefix in the future releases

## Output

The main output from the pipeline is the file `data/results.jtree`, an extensively-annotated rooted tree in jtree format of [`treeio`](https://bioconductor.org/packages/release/bioc/html/treeio.html). For each node the following data is provided:

* `branch.length.erable` branch length calculated by ERaBLE (same as the branch length)
* `split` the bipartition corresponding to the node
* `all.tips.new` a flag indicating whether all of the tips from the node are new with respect to the gtdb reference tree
* `gtdb.split` correponding bipartition in GTDB (if there was such a biparitition)

For internal nodes only:

* `branch.length.astral` for internal nodes, length of the branch in coalescence units as calculated by ASTRAL
* `QC` the total number of quartets around the branch
* `EN` effective number of genes
* for other ASTRAL annotations see [ASTRAL's tutorial](https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md#newick-annotations)

For tips only:

* `genes.phylo` number of genes used for phylogeny
* full set of metadata provided by GTDB

## Dependencies

The workflow's dependencies are not hooked to conda and should be installed separately. Use the rules `dependencies-cli` and `dependencies-r` to check.

## Known issues

The rule `ncbi_genome_download` utilizes the python package `ncbi_genome_download` to download data from NCBI Assembly for the reference assemblies. In case of discrepancy between GTDB and NCBI (e.g. RefSeq deleted), this step will get stuck and require manual resolution. In such cases the missing files (`*_genomic.fna.gz`) should be downloaded e.g. from the GenBank section manually and put in the corresponding folder under `data/refseq/{accession}/`.

There is currently no mechanism to resolve the rare situations of the fosmids containing genes from GTDB's reference marker list which introduces an incomplete redundancy in the species inference step.

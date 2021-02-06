
packages <- c("dplyr", "tidyr", "treeio", "phytools", "ape", "tools", "seqinr", "phangorn", "bioformatr")
invisible(lapply(packages, function(pkg) suppressPackageStartupMessages(library(pkg, character.only = T))))

options(stringsAsFactors = F)

metadata <- Sys.getenv("GTDBTK_DATA_PATH") %>%
	file.path("metadata/genome_metadata.tsv") %>%
	read.table(sep = "\t", header = T, quote = "", na.strings = c("na", ""), comment.char = "") %>%
	mutate(prefix = ifelse(grepl("^GCF", accession), "RS_", "GB_"), label = paste0(prefix, accession))

get.outgroup <- function(tree) {
	left_right <- child(tree, rootnode(tree))
	left_tips  <- offspring(tree, left_right[1], tiponly = T, self_include = T)
	right_tips <- offspring(tree, left_right[2], tiponly = T, self_include = T)
	tree$tip.label[if (length(left_tips) < length(right_tips)) left_tips else right_tips]
}
get.descendants <- function(tree, node, labels.ignore = c()) {
	mapply(offspring, .node = node, MoreArgs = list(.data = tree, tiponly = T, self_include = T), SIMPLIFY = F) %>%
		lapply(function(x) tree$tip.label[x]) %>%
		lapply(function(x) x[! x %in% labels.ignore]) %>%
		lapply(sort) %>%
		sapply(paste, collapse = ";") %>%
		`[<-`(. == "", NA)
}
get.splits <- function(tree, nodes, labels.ignore = c()) {
	as.splits(tree) %>%
		as.matrix %>%
		data.frame(node = 1:nrow(.)) %>%
		gather(label, split, -node) %>%
		filter(! label %in% labels.ignore) %>%
		arrange(label) %>%
		group_by(node, split) %>%
		summarize(label = paste(label, collapse = ";"), .groups = "drop_last") %>%
		summarize(left  = pmin(label[split == 0], label[split == 1]), right = pmax(label[split == 0], label[split == 1]), .groups = "drop") %>%
		`[`(match(nodes, .$node),) %>%
		mutate(split = ifelse(is.na(left), NA, paste(left, right, sep = "|"))) %>%
		pull
}
all.tips.new <- function(tree, node, new.tips) {
	mapply(offspring, .node = node, MoreArgs = list(.data = tree, tiponly = T, self_include = T), SIMPLIFY = F) %>%
		lapply(function(x) tree$tip.label[x] %in% new.tips) %>%
		sapply(all)
}

txt.fnames <- Sys.glob("../hmm/*.txt")
txt <- lapply(txt.fnames, read.table, col.name = c("Seq.Name", "Group")) %>%
	setNames(basename(txt.fnames) %>% file_path_sans_ext) %>%
	bind_rows(.id = "hmm") %>%
	mutate(Group = paste("hmm", hmm, Group, sep = "."))

ublast.fnames <- Sys.glob(c("gtdbtk/*+*.ublast", "contigs/*+*.ublast")) %>%
	file.info %>%
	filter(size > 0) %>%
	rownames
hmm.hits <- lapply(ublast.fnames, read.table, sep = "\t", col.names = c("qseqid", "Seq.Name", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")) %>%
	setNames(ublast.fnames) %>%
	bind_rows(.id = "fname") %>%
	distinct(qseqid, .keep_all = T) %>%
	extract(fname, into = c("label", "hmm"), regex = ".+/(.+)[+](.+).ublast") %>%
	filter(pident > 50) %>%
	left_join(txt, by = "Seq.Name") %>%
	group_by(label, Group) %>%
	summarize(n = n()) %>%
	spread(Group, n)

clade <- read.jtree("clade.jtree")

gene.files <- Sys.glob("phylogeny/*/treeshrink.trim")
gene.names <- file_path_sans_ext(gene.files)
gene.lens <- gene.files %>%
	lapply(read.fasta) %>%
	lapply(`[`, 1) %>%
	lapply(unlist) %>%
	lapply(length) %>%
	setNames(gene.names)
gene.dists <- paste0(gene.names, ".treefile") %>%
	lapply(read.tree) %>%
	do.call(c,.) %>%
	lapply(cophenetic) %>%
	lapply(as.table) %>%
	setNames(gene.names)

contigs <- Sys.glob("contigs/*.faa")
contig.names <- basename(contigs) %>% file_path_sans_ext
genes.all <- lapply(contigs, read.fasta, as.string = T) %>%
	sapply(length) %>%
	setNames(contig.names)
genes.phylo <- lapply(gene.dists, colnames) %>%
	unlist %>% table %>% c

# Write distance matrices for ERABLE
lapply(gene.names, function(gene.name) c(
	"",
	paste(ncol(gene.dists[[gene.name]]), gene.lens[[gene.name]]),
	capture.output(write.table(gene.dists[[gene.name]], sep = "\t", col.names = F, quote = F))
)) %>% {c(length(.), unlist(.))} %>%
	cat(file = "phylogeny.trees.dists", sep = "\n")

astral <- read.astral("astral.tree")
astral@data <- bind_rows(astral@data, data.frame(node = 1:Ntip(astral)))
astral.backbone <- unroot(astral@phylo)
astral.backbone$edge.length <- NULL
write.tree(astral.backbone, "astral.backbone.tree")

# Run ERABLE
system("erable -i phylogeny.trees.dists -t astral.backbone.tree -o erable.tree")

erable.tree <- read.tree("erable.tree.length.nwk") %>% midpoint
erable.tree$edge.length[erable.tree$edge.length < 0] <- 0

clade.data <- as_tibble(clade) %>%
	filter(node > Ntip(clade)) %>%
	mutate(gtdb.root = parent == rootnode(clade@phylo)) %>%
	mutate(gtdb.descendants = get.descendants(clade@phylo, node)) %>%
	mutate(gtdb.split       = get.splits(clade@phylo, node)) %>%
	select(gtdb.descendants, gtdb.split, gtdb.root, gtdb.support = support, gtdb.taxon = taxon, gtdb.branch.length = branch.length)

astral.data <- as_tibble(astral) %>%
	filter(!is.na(branch.length), !is.nan(branch.length)) %>%
	rename(branch.length.astral = branch.length) %>%
	mutate(gtdb.split = get.splits(astral@phylo, node, contig.names)) %>%
	mutate(split      = get.splits(astral@phylo, node)) %>%
	left_join(clade.data, by = "gtdb.split") %>%
	group_by(gtdb.split) %>%
	mutate(gtdb.skip = any(gtdb.split == split) & gtdb.split != split) %>%
	mutate(gtdb.support = ifelse(gtdb.skip, NA, gtdb.support), gtdb.root = ifelse(gtdb.skip, NA, gtdb.root)) %>%
	mutate(EN = as.integer(EN), QC = as.integer(QC)) %>%
	select(-node, -parent, -label, -gtdb.skip, -gtdb.taxon)
gtdb.taxa <- as_tibble(erable.tree) %>%
	mutate(gtdb.descendants = get.descendants(erable.tree, node, contig.names), descendants = get.descendants(erable.tree, node)) %>%
	left_join(clade.data, by = "gtdb.descendants") %>%
	group_by(gtdb.descendants) %>%
	mutate(gtdb.skip = any(gtdb.descendants == descendants) & gtdb.descendants != descendants) %>%
	filter(!gtdb.skip, !is.na(gtdb.taxon)) %>%
	select(descendants, gtdb.descendants, gtdb.taxon)

as_tibble(erable.tree) %>%
	mutate(split        = get.splits(erable.tree, node)) %>%
	mutate(descendants  = get.descendants(erable.tree, node)) %>%
	mutate(all.tips.new = all.tips.new(erable.tree, node, contig.names)) %>%
	mutate(branch.length.erable = branch.length) %>%
	left_join(astral.data, by = "split") %>%
	left_join(gtdb.taxa,   by = "descendants") %>%
	left_join(hmm.hits,    by = "label") %>%
	left_join(metadata,    by = "label") %>%
	mutate(genes.all = genes.all[label], genes.phylo = genes.phylo[label]) %>%
	mutate_if(is.character, list(~na_if(.,""))) %>%
	distinct(node, .keep_all = T) %>%
	as.treedata %>%
	write.jtree("erable.jtree")


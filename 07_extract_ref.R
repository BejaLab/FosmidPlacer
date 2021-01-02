
packages <- c("dplyr", "tidyr", "treeio", "phytools", "ape", "phangorn", "tools", "seqinr", "bioformatr")
invisible(lapply(packages, function(pkg) suppressPackageStartupMessages(library(pkg, character.only = T))))

options(stringsAsFactors = F)

metadata <- Sys.getenv("GTDBTK_DATA_PATH") %>%
        file.path("metadata/genome_metadata.tsv") %>%
        read.table(sep = "\t", header = T, quote = "", fill = T, na.strings = c("na", "")) %>%
	mutate(prefix = ifelse(grepl("^GCF", accession), "RS_", "GB_"), label = paste0(prefix, accession)) %>%
	select(label, ncbi_organism_name)

find.outgroup <- function(tree) {
	left_right <- child(tree, rootnode(tree))
	left_tips  <- offspring(tree, left_right[1], tiponly = T, self_include = T)
	right_tips <- offspring(tree, left_right[2], tiponly = T, self_include = T)
	tree$tip.label[if (length(left_tips) < length(right_tips)) left_tips else right_tips]
}
get.splits <- function(tree, labels.ignore = c()) {
	as.splits(tree) %>%
		as.matrix %>%
		data.frame(node = 1:nrow(.)) %>%
		gather(label, split, -node) %>%
		filter(! label %in% labels.ignore) %>%
		arrange(label) %>%
		group_by(node, split) %>%
		summarize(label = paste(label, collapse = ";"), .groups = "drop_last") %>%
		summarize(left  = pmin(label[split == 0], label[split == 1]), right = pmax(label[split == 0], label[split == 1]), .groups = "drop") %>%
		mutate(split = paste(left, right, sep = "|"), .groups = "drop") %>%
		select(node, split)
}
all.tips.new <- function(node, tree, new.tips) {
	mapply(offspring, .node = node, MoreArgs = list(.data = tree, tiponly = T, self_include = T), SIMPLIFY = F) %>%
		lapply(function(x) tree$tip.label[x] %in% new.tips) %>%
		sapply(all)
}


clade <- read.jtree("clade.jtree")

tblout.fnames <- Sys.glob(c("gtdbtk/*.tblout", "contigs/*.tblout"))
tblout <- lapply(tblout.fnames, read.hmmer.tblout) %>%
	setNames(tblout.fnames) %>%
	bind_rows(.id = "fname") %>%
	filter(full.sequence.E.value < 1e-10) %>%
	mutate(label = basename(fname)) %>%
	mutate(label = gsub("[+].+", "", label)) %>%
	select(label, query.name) %>%
	group_by(label, query.name) %>%
	summarize(hits = n()) %>%
	spread(query.name, hits)

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

new.labels <- astral@phylo$tip.label %>% `[`(! . %in% clade@phylo$tip.label)

clade.data <- as_tibble(clade) %>%
	filter(node > Ntip(clade)) %>%
	mutate(gtdb.root = parent == rootnode(clade@phylo)) %>%
	left_join(get.splits(clade@phylo), by = "node") %>%
	select(split.old = split, gtdb.root, gtdb.support = support, gtdb.taxon = taxon, gtdb.branch.length = branch.length)

astral.data <- as_tibble(astral) %>%
	filter(!is.na(branch.length), !is.nan(branch.length)) %>%
	rename(branch.length.astral = branch.length) %>%
	left_join(get.splits(astral@phylo, new.labels), by = "node") %>%
	rename(split.old = split) %>%
	left_join(get.splits(astral@phylo), by = "node") %>%
	left_join(clade.data, by = c("split.old")) %>%
	group_by(split.old) %>%
	mutate(gtdb.skip = any(split == split.old) & split != split.old) %>%
	mutate(gtdb.support = ifelse(gtdb.skip, NA, gtdb.support), gtdb.taxon = ifelse(gtdb.skip, NA, gtdb.taxon), gtdb.root = ifelse(gtdb.skip, NA, gtdb.root)) %>%
	select(-node, -parent, -label, -gtdb.skip)

as_tibble(erable.tree) %>%
	left_join(get.splits(erable.tree), by = "node") %>%
	mutate(branch.length.erable = branch.length) %>%
	left_join(astral.data, by = "split") %>%
	left_join(tblout, by = "label") %>%
	left_join(metadata, by = "label") %>%
	mutate(all.tips.new = all.tips.new(node, erable.tree, new.labels)) %>%
	mutate_if(is.character, list(~na_if(.,""))) %>%
	distinct(node, .keep_all = T) %>%
	as.treedata %>%
	write.jtree("erable.jtree")

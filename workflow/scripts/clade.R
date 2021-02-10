#!/usr/bin/Rscript

taxon  <- snakemake@params$taxon
gtdb_tree <- snakemake@input$tree
marker_genes <- snakemake@input$marker_genes

suppressPackageStartupMessages(library(treeio))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ape))

ref.tree <- read.tree(gtdb_tree)
nodes <- as_tibble(ref.tree) %>%
	filter(node > Ntip(ref.tree), node != rootnode(ref.tree)) %>%
	separate(label, into = c("support", "taxon"), sep = ":", fill = "right", convert = T) %>%
	select(node, support, taxon)
clade.node <- filter(nodes, taxon %in% !!taxon) %>% pull(node)
if (length(clade.node) == 0) {
	write("Taxon not found on the GTDB reference tree", stderr())
	q(status = 1)
}
if (length(clade.node) > 1) {
	clade.node <- getMRCA(ref.tree, clade.node)
}
clade <- as_tibble(ref.tree) %>%
	left_join(nodes, by = "node") %>%
	mutate(label = gsub("^(RS|GB)_", "", label)) %>%
	as.treedata %>%
	tree_subset(node = clade.node, levels_back = 0)

cat(clade@phylo$tip.label, sep = "\n", file = "data/gtdb_list.txt")

write.jtree(clade, "data/clade.jtree")
clade.reduced <- clade@phylo
edges <- as_tibble(clade) %>%
	filter(node > Ntip(clade), node != rootnode(clade)) %>%
	left_join(data.frame(parent = clade@phylo$edge[,1], node = clade@phylo$edge[,2], edge.id = 1:nrow(clade@phylo$edge)), by = c("node","parent")) %>%
	separate(label, into = c("support","taxon"), sep = ":", fill = "right", convert = T) %>%
	arrange(edge.id) %>%
	mutate(keep = (!is.na(taxon) & support > 89 | support > 96)) %>%
	select(edge.id, keep)

clade.reduced$edge.length[edges$edge.id] <- edges$keep
clade.reduced <- di2multi(clade.reduced)
clade.reduced$edge.length <- NULL
clade.reduced$node.label  <- NULL
write.tree(clade.reduced, file = "data/clade.reduced.tree")

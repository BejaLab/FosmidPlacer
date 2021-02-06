#!/usr/bin/Rscript

dir <- commandArgs(T)
if (length(dir) > 0 & dir.exists(dir[1])) pwd(dir[1])

source("taxa.cfg")

if (! exists("taxon") ) {
	write("taxon not specified in taxa.cfg", stderr())
	q(status = 1)
}

suppressPackageStartupMessages(library(treeio))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ape))

ref.tree <- read.tree("gtdbtk.tree")
nodes <- as_tibble(ref.tree) %>%
	filter(node > Ntip(ref.tree), node != rootnode(ref.tree)) %>%
	separate(label, into = c("support", "taxon"), sep = ":", fill = "right", convert = T) %>%
	select(node, support, taxon)
clade.node <- filter(nodes, taxon %in% !!taxon) %>%
	pull(node) %>%
	{ifelse(length(.) > 1, getMRCA(ref.tree, .), .)}

clade <- as_tibble(ref.tree) %>%
	left_join(nodes, by = "node") %>%
	as.treedata %>%
	tree_subset(node = clade.node, levels_back = 0)

cat(clade@phylo$tip.label, file = "gtdbtk.txt", sep = "\n")
write.jtree(clade, "clade.jtree")

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
write.tree(clade.reduced, file = "clade.reduced.tree")

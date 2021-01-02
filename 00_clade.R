suppressPackageStartupMessages(library(treeio))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ape))

source("taxa.cfg")

if (! exists("taxon") | ! exists("domain") ) {
	write("taxon and domain must be specified in taxa.cfg", stderr())
	q(status = 1)
}

ref.tree.fname <- paste0(domain, "_r95.tree")

ref.tree <- Sys.getenv("GTDBTK_DATA_PATH") %>%
	file.path(ref.tree.fname) %>%
	read.tree
nodes <- as_tibble(ref.tree) %>%
	filter(node > Ntip(ref.tree), node != rootnode(ref.tree)) %>%
	separate(label, into = c("support", "taxon"), sep = ":", fill = "right", convert = T) %>%
	select(node, support, taxon)
node <- filter(nodes, taxon == !!taxon) %>% pull(node)

clade <- as_tibble(ref.tree) %>%
	left_join(nodes, by = "node") %>%
	as.treedata %>%
	tree_subset(node = node, levels_back = 0)

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

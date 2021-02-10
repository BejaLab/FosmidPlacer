library(dplyr)
library(tidyr)
library(seqinr)

contigs <- Sys.glob("input/*.faa")

fnames <- Sys.glob("data/gtdb/*.faa_selected") %>%
	`[`(file.size(.) > 0) %>%
	c(contigs)

faa <- lapply(fnames, read.fasta, as.string = T, seqtype = "AA") %>%
	setNames(basename(fnames)) %>%
	lapply(unlist) %>%
	lapply(data.frame) %>%
	bind_rows(.id = "faa") %>%
	mutate(rownames(.)) %>%
	setNames(c("faa","sequence","gene"))

if (!dir.exists("data/phylogeny")) dir.create("data/phylogeny")

orthogroups <- read.table("data/ublast.proteinortho.tsv", comment.char = "", sep = "\t", header = T, na.strings = "*") %>%
	mutate(og = sprintf("og%03d", 1:n())) %>%
	filter(.[basename(contigs)] %>% is.na %>% `!` %>% rowSums %>% `>`(0)) %>%
	filter(X..Species == Genes, Genes > 4) %>%
	gather(faa, gene, -X..Species, -Genes, -Alg..Conn., -og, na.rm = T) %>%
	left_join(faa, by = c("faa", "gene")) %>%
	mutate(fasta = sprintf(">%s %s\n%s", sub("\\.faa[a-zA-Z_0-9]*$", "", faa), gene, sequence)) %>%
	group_by(og) %>%
	group_map(function(og, meta) {
		og.name <- paste0("orthogroup_", meta$og)
		og.dir <- file.path("data/phylogeny", og.name)
		if (!dir.exists(og.dir)) dir.create(og.dir, recursive = T)
		pull(og, fasta) %>% cat(file = file.path(og.dir, "gene.faa"), sep = "\n")
		return(og.name)
	}) %>% unlist

fileConn <- file("data/filter_orthogroups.txt")
writeLines(orthogroups, fileConn)
close(fileConn)

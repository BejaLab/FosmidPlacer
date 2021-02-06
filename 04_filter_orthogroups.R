library(dplyr)
library(tidyr)
library(seqinr)
library(tools)

contigs <- Sys.glob("contigs/*.faa") %>% basename

fnames <- Sys.glob("orthologs/*.faa") %>% `[`(file.size(.) > 0)
faa <- lapply(fnames, read.fasta, as.string = T, seqtype = "AA") %>%
	setNames(basename(fnames)) %>%
	lapply(unlist) %>%
	lapply(data.frame) %>%
	bind_rows(.id = "faa") %>%
	mutate(rownames(.)) %>%
	setNames(c("faa","sequence","gene"))
if (!dir.exists("phylogeny")) dir.create("phylogeny")
orthogroups <- read.table("ublast.proteinortho.tsv", comment.char="", sep = "\t", header = T, na.strings = "*") %>%
	mutate(og = sprintf("og%03d", 1:n())) %>%
	filter(.[contigs] %>% is.na %>% `!` %>% rowSums %>% `>`(0)) %>%
	filter(X..Species == Genes, Genes > 4) %>%
	gather(faa, gene, -X..Species, -Genes, -Alg..Conn., -og, na.rm = T) %>%
	left_join(faa, by = c("faa", "gene")) %>%
	mutate(fasta = sprintf(">%s %s\n%s", file_path_sans_ext(faa), gene, sequence)) %>%
	group_by(og) %>%
	group_map(function(og, meta) {
		og.dir <- file.path("phylogeny", meta$og)
		if (!dir.exists(og.dir)) dir.create(og.dir, recursive = T)
		pull(og, fasta) %>% cat(file = file.path(og.dir, "gene.faa"), sep = "\n")
	})


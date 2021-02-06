
mkdir -p gtdbtk

predict_genes() {
	local base=$1
	mkdir -p "gtdbtk/$base"
	cd "gtdbtk/$base"
	gmsn.pl --format GFF --prok --output /dev/stdout "../$base.fna" | awk -F\\t '{sub(/ .+/,"",$1)}{sub(/gene_id=/,sprintf("ID=%s_gmsn_",$1),$9)}1' OFS=\\t |
		tee "$base.gff" | gffread -H -y "../$base.faa" -g "../$base.fna" -
}
export -f predict_genes

# NB: make paths more flexible
parallel \
	ln -fs ../../../gtdbtk/{}.fna gtdbtk/ \; \
	[ -e ../../gtdbtk/{}.faa ] '&&' ln -fs ../../../gtdbtk/{}.faa gtdbtk/ \
	:::: gtdbtk.txt
parallel [ -s gtdbtk/{}.faa ] '||' predict_genes {} :::: gtdbtk.txt


OKBLUE = '\033[94m'
ENDC = '\033[0m'
FAIL = '\033[91m'

all.packages <- installed.packages()
dep.packages <- snakemake@params$packages

installed <- dep.packages[dep.packages %in% rownames(all.packages)]
not.installed <- dep.packages[! dep.packages %in% rownames(all.packages)]

package.versions <- data.frame(
	package = installed,
	version = all.packages[installed, "Version"]
)

if (length(not.installed) > 0) {
	write(FAIL, stderr())
	write("The following R packages are required:", stderr())
	write(paste(" - ", not.installed), stderr())
	write(ENDC, stderr())
	q(status = 1)
}

fileConn <- file("data/dependencies-r.txt")
writeLines(with(package.versions, paste0(package, ": v. ", version)), fileConn)
close(fileConn)

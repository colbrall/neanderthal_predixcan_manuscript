# source code from DECIPHER package https://rdrr.io/bioc/DECIPHER/
# hopefully to avoid having to install/import entire package
# outputs R dendrogram object in Newick tree format 

WriteDendrogram <- function(x,
	file="",
	quoteLabels=TRUE,
	convertBlanks=!quoteLabels,
	internalLabels=TRUE,
	digits=10,
	append=FALSE) {

	# error checking
	if (class(x) != "dendrogram")
		stop("x is not a dendrogram.")
	if (!is.logical(quoteLabels))
		stop("quoteLabels must be a logical.")
	if (!is.logical(convertBlanks))
		stop("convertBlanks must be a logical.")
	if (!is.numeric(digits))
		stop("digits must be a numeric.")
	if (floor(digits)!=digits)
		stop("digits must be a whole number.")
	if (digits < 1)
		stop("digits must be at least 1.")
	if (!is.logical(append))
		stop("append must be a logical.")
	if (!is.logical(internalLabels))
		stop("internalLabels must be a logical.")

	if (is.character(file)) {
		if (file == "") {
			file <- stdout()
		} else if (substring(file, 1L, 1L) == "|") {
			file <- pipe(substring(file, 2L), "w")
			on.exit(close(file))
		} else {
			file <- file(file, "w")
			on.exit(close(file))
		}
	}

	getLab <- function(LAB) {
		if (is.null(LAB))
			return("")
		lab <- gsub("'", "''", LAB, fixed=TRUE)
		if (convertBlanks)
			lab <- gsub(" ", "_", lab, fixed=TRUE)
		if (quoteLabels)
			lab <- paste('"', lab, '"', sep="")
		return(lab)
	}

	.dendrogram2newick <- function(x, height=attr(x, "height"), root=TRUE) {
		if (is.leaf(x)) {
			cat(getLab(attr(x, "label")),
				":",
				round(height - attr(x, "height"),
					digits=digits),
				sep="",
				file=file,
				append=TRUE)
		} else {
			cat("(",
				file=file,
				append=TRUE)
			for (i in seq_along(x)) {
				.dendrogram2newick(x[[i]],
					attr(x, "height"),
					root=FALSE)
				if (i < length(x))
					cat(",",
						file=file,
						append=TRUE)
			}
			if (root) {
				cat(");\n",
					file=file,
					append=TRUE)
			} else {
				cat(")",
					ifelse(internalLabels,
						getLab(attr(x, "edgetext")),
						""),
					":",
					round(height - attr(x, "height"),
						digits=digits),
					sep="",
					file=file,
					append=TRUE)
			}
		}
	}

	if (!append) # overwrite the file
		cat("", file=file)
	.dendrogram2newick(x)
	invisible(NULL)
}

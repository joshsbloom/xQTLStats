#' Query variant effects for a set of strains
#'
#' Given a variant lookup object (from a qs2 file) and a vector of strain names,
#' returns a data.table of all variants that are polymorphic in that strain subset.
#' The result left-joins variant annotations with coding effects, so coding columns
#' (GENEID, CONSEQUENCE, etc.) are filled with NA for non-coding variants. Variants
#' affecting overlapping ORFs will have multiple rows (one per gene).
#'
#' Per-strain dosage columns (0/1/2/NA) are included, along with subset allele
#' frequency and alt allele count.
#'
#' @param lookup List from \code{qs2::qs_read('variant_lookup.qs')}, containing
#'   \code{Z}, \code{missing_mask}, \code{strain_ids}, \code{variants},
#'   \code{variant_ids}, and \code{coding}.
#' @param strain.names Character vector of strain names to query.
#' @return A data.table with one row per variant (or per variant x gene for
#'   coding variants in overlapping ORFs). Columns from \code{variants} and
#'   \code{coding} are merged, with coding columns set to NA for non-coding
#'   variants. Additional columns: \code{subset_AC}, \code{subset_AF}, and
#'   one column per strain with integer dosage (0/1/2/NA).
#' @export
query_variants <- function(lookup, strain.names) {

    si <- match(strain.names, lookup$strain_ids)
    missing <- strain.names[is.na(si)]
    si <- si[!is.na(si)]
    if (length(si) == 0) stop("No matching strains found in lookup$strain_ids")
    if (length(missing) > 0) {
        message("  Note: ", length(missing), " strain(s) not found in lookup: ",
                paste(head(missing, 5), collapse = ', '),
                if (length(missing) > 5) ', ...' else '')
    }
    matched.names <- lookup$strain_ids[si]

    # Subset Z to these strains
    Z.sub <- lookup$Z[si, , drop = FALSE]

    # Restore NAs from missing mask
    if (!is.null(lookup$missing_mask)) {
        mask.sub <- lookup$missing_mask[si, , drop = FALSE]
    } else {
        mask.sub <- NULL
    }

    # Identify polymorphic variants in this subset
    acnt <- as.integer(Matrix::colSums(Z.sub))
    if (!is.null(mask.sub)) {
        non_missing <- as.integer(nrow(Z.sub) - Matrix::colSums(mask.sub))
    } else {
        non_missing <- rep(nrow(Z.sub), ncol(Z.sub))
    }
    max_dosage <- 2L * non_missing
    if (length(si) == 1L) {
        # Single strain: keep all non-REF variants (the other parent is the reference)
        varies <- acnt > 0L & non_missing > 0L
    } else {
        # Multiple strains: keep variants polymorphic within the subset
        varies <- acnt > 0L & acnt < max_dosage & non_missing > 0L
    }

    # Subset to polymorphic columns
    poly.idx <- lookup$variant_ids[varies]
    Z.poly <- Z.sub[, varies, drop = FALSE]
    acnt <- acnt[varies]
    non_missing <- non_missing[varies]

    # Densify and restore NAs for polymorphic variants only
    Z.dense <- as.matrix(Z.poly)
    if (!is.null(mask.sub)) {
        mask.poly <- mask.sub[, varies, drop = FALSE]
        Z.dense[as.matrix(mask.poly) == 1] <- NA_integer_
    }

    # Build per-strain genotype columns (transposed: variants as rows)
    gt.dt <- data.table::as.data.table(t(Z.dense))
    data.table::setnames(gt.dt, matched.names)
    data.table::set(gt.dt, j = 'idx', value = poly.idx)

    # Subset allele stats
    data.table::set(gt.dt, j = 'subset_AC', value = as.integer(acnt))
    data.table::set(gt.dt, j = 'subset_AF', value = acnt / (2L * non_missing))

    # Join variant annotations — avoid data.table [ keyed lookup (dispatch issues in package namespace).
    var_rows <- match(poly.idx, lookup$variants$idx)
    vars_dt <- lookup$variants[var_rows, ]
    gt_cols <- setdiff(names(gt.dt), names(vars_dt))
    result <- cbind(vars_dt, gt.dt[, gt_cols, with = FALSE])

    # Left-join coding annotations by idx
    coding_match <- match(result$idx, lookup$coding$idx)
    coding_cols <- setdiff(names(lookup$coding), 'idx')
    for (cc in coding_cols) {
        data.table::set(result, j = cc, value = lookup$coding[[cc]][coding_match])
    }

    # Sort by genomic position
    chr.order <- paste0('chr', as.roman(1:16))
    data.table::set(result, j = 'chr_f', value = factor(result$chr, levels = chr.order))
    data.table::setorder(result, chr_f, pos)
    data.table::set(result, j = 'chr_f', value = NULL)

    data.table::setkey(result, 'idx')
    result
}

#' Extract variants within QTL peak intervals for parent strains
#'
#' Given a GRanges object from \code{\link{findQTLPeaks}} and a variant lookup
#' object, queries polymorphic variants between the parent strains within each
#' QTL confidence interval. Returns a named list with one element per QTL peak.
#'
#' @param peaks GRanges object from \code{\link{findQTLPeaks}}.
#' @param lookup List from \code{qs2::qs_read('variant_lookup.qs')}, the same
#'   input used by \code{\link{query_variants}}.
#' @param strain.names Character vector of strain names (parent strains) to
#'   query in the lookup. If NULL (default), uses the unique values from
#'   the \code{high.parent} and \code{low.parent} metadata columns of peaks.
#' @return A named list with one element per QTL peak. Each element is a list
#'   containing:
#'   \describe{
#'     \item{peak}{A one-row GRanges for the interval.}
#'     \item{variants}{A data.table of polymorphic variants within the interval
#'       (output of \code{\link{query_variants}} filtered to the region).
#'       NULL if no variants found.}
#'   }
#'   List names are formatted as \code{chr:start-end}.
#' @export
extractQTLVariants <- function(peaks, lookup, strain.names = NULL) {

    if (!requireNamespace("GenomicRanges", quietly = TRUE))
        stop("GenomicRanges is required. Install with: BiocManager::install('GenomicRanges')")

    if (length(peaks) == 0) return(list())

    # Determine strains to query
    if (is.null(strain.names)) {
        strain.names <- unique(c(
            GenomicRanges::mcols(peaks)$high.parent,
            GenomicRanges::mcols(peaks)$low.parent
        ))
    }

    # Resolve strain row indices in Z (once)
    si <- match(strain.names, lookup$strain_ids)
    missing <- strain.names[is.na(si)]
    si <- si[!is.na(si)]
    if (length(si) == 0) stop("No matching strains found in lookup$strain_ids")
    if (length(missing) > 0) {
        message("  Note: ", length(missing), " strain(s) not found in lookup: ",
                paste(head(missing, 5), collapse = ', '),
                if (length(missing) > 5) ', ...' else '')
    }
    matched.names <- lookup$strain_ids[si]

    # Build a GRanges of all variants (point ranges) for interval overlap.
    # variant_ids and variants$idx are aligned: row i in variants = column i in Z.
    var_gr <- GenomicRanges::GRanges(
        seqnames = lookup$variants$chr,
        ranges = IRanges::IRanges(start = lookup$variants$pos, width = 1L)
    )

    # Ensure peaks is a proper GRanges, not a list
    if (is.list(peaks) && !methods::is(peaks, "GRanges")) {
        stop("peaks must be a GRanges object, not a list. ",
             "If you have a list of GRanges, call extractQTLVariants on each element.")
    }

    # Single findOverlaps call: interval tree, O(n log n + k) total
    hits <- GenomicRanges::findOverlaps(peaks, var_gr)
    peak_hits <- S4Vectors::queryHits(hits)
    var_hits  <- S4Vectors::subjectHits(hits)

    message("  findOverlaps: ", length(hits), " variant hits across ", length(peaks), " peaks")

    # Group variant column indices by peak
    hits_by_peak <- split(var_hits, peak_hits)

    result <- vector("list", length(peaks))
    peak_names <- paste0(
        as.character(GenomicRanges::seqnames(peaks)), ":",
        GenomicRanges::start(peaks), "-", GenomicRanges::end(peaks)
    )

    for (i in seq_along(peaks)) {
        pk <- peaks[i]
        col_idx <- hits_by_peak[[as.character(i)]]

        if (is.null(col_idx) || length(col_idx) == 0) {
            message("  Peak ", i, ": no overlap hits")
            result[[i]] <- list(peak = pk, variants = NULL)
            next
        }

        # Extract small submatrix: matched strains x interval variants
        Z.sub <- lookup$Z[si, col_idx, drop = FALSE]

        # Filter to non-REF (or polymorphic) variants in this subset
        acnt <- as.integer(Matrix::colSums(Z.sub))
        if (!is.null(lookup$missing_mask)) {
            mask.sub <- lookup$missing_mask[si, col_idx, drop = FALSE]
            non_missing <- as.integer(nrow(Z.sub) - Matrix::colSums(mask.sub))
        } else {
            mask.sub <- NULL
            non_missing <- rep(nrow(Z.sub), length(col_idx))
        }

        if (length(si) == 1L) {
            varies <- acnt > 0L & non_missing > 0L
        } else {
            max_dosage <- 2L * non_missing
            varies <- acnt > 0L & acnt < max_dosage & non_missing > 0L
        }

        message("  Peak ", i, ": ", length(col_idx), " overlap variants, ",
                sum(varies), " pass filter, length(si)=", length(si))

        if (!any(varies)) {
            result[[i]] <- list(peak = pk, variants = NULL)
            next
        }

        # Densify only the passing variants
        Z.poly <- Z.sub[, varies, drop = FALSE]
        Z.dense <- as.matrix(Z.poly)
        if (!is.null(mask.sub)) {
            mask.poly <- mask.sub[, varies, drop = FALSE]
            Z.dense[as.matrix(mask.poly) == 1] <- NA_integer_
        }

        poly_col_idx <- col_idx[varies]
        poly_var_ids <- lookup$variant_ids[poly_col_idx]

        # Build genotype columns
        gt.dt <- data.table::as.data.table(t(Z.dense))
        data.table::setnames(gt.dt, matched.names)
        data.table::set(gt.dt, j = 'idx', value = poly_var_ids)
        data.table::set(gt.dt, j = 'subset_AC', value = as.integer(acnt[varies]))
        nm <- non_missing[varies]
        data.table::set(gt.dt, j = 'subset_AF', value = acnt[varies] / (2L * nm))

        # Join with annotations — avoid data.table [ keyed lookup (dispatch issues in package namespace).
        var_rows <- match(poly_var_ids, lookup$variants$idx)
        vars_dt <- lookup$variants[var_rows, ]
        # gt.dt already has idx; drop it before cbind to avoid duplication
        gt_cols <- setdiff(names(gt.dt), names(vars_dt))
        interval_vars <- cbind(vars_dt, gt.dt[, gt_cols, with = FALSE])
        # Left-join coding annotations by idx
        coding_match <- match(interval_vars$idx, lookup$coding$idx)
        coding_cols <- setdiff(names(lookup$coding), 'idx')
        for (cc in coding_cols) {
            data.table::set(interval_vars, j = cc, value = lookup$coding[[cc]][coding_match])
        }
        message("  Peak ", i, ": final rows=", nrow(interval_vars))
        data.table::setorder(interval_vars, pos)

        result[[i]] <- list(
            peak = pk,
            variants = if (nrow(interval_vars) > 0) interval_vars else NULL
        )
    }

    names(result) <- peak_names
    result
}

#' Flatten nested QTL variant results into a single data.table
#'
#' Combines the nested list output of multiple \code{\link{extractQTLVariants}}
#' calls into a single data.table. Per-strain genotype columns (which may differ
#' across intervals) are collapsed into a single \code{genotypes} string column
#' (e.g. \code{"ACD=2"} or \code{"ACD=2;BYa=0"}), so all rows share the same
#' schema. Peak metadata (parent strains, LOD, position) is added per row.
#'
#' @param qtl_vars Named list of \code{\link{extractQTLVariants}} results,
#'   typically one per cross/comparison.
#' @param lookup Variant lookup object (same as passed to \code{\link{extractQTLVariants}}).
#'   Used to identify which columns are strain genotypes vs annotations.
#' @return A data.table with variant annotations, genotype strings, and peak
#'   metadata columns: comparison, peak, high.parent, low.parent, peak.pos, peak.LOD.
#'   Returns an empty data.table if no variants are found.
#' @export
flattenQTLVariants <- function(qtl_vars, lookup) {

    if (!requireNamespace("GenomicRanges", quietly = TRUE))
        stop("GenomicRanges is required. Install with: BiocManager::install('GenomicRanges')")

    strain_ids <- lookup$strain_ids

    rows <- list()
    for (comp in names(qtl_vars)) {
        for (interval_name in names(qtl_vars[[comp]])) {
            entry <- qtl_vars[[comp]][[interval_name]]
            if (is.null(entry$variants)) next
            v <- data.table::copy(entry$variants)

            pk <- entry$peak

            # Identify strain genotype columns by matching against lookup strain IDs
            strain_cols <- intersect(names(v), strain_ids)

            # Collapse per-strain dosages into a single string column
            if (length(strain_cols) > 0) {
                gt_strings <- apply(
                    as.matrix(v[, strain_cols, with = FALSE]), 1,
                    function(row) paste(strain_cols, row, sep = "=", collapse = ";")
                )
                data.table::set(v, j = "genotypes", value = gt_strings)
                for (sc in strain_cols) data.table::set(v, j = sc, value = NULL)
            }

            data.table::set(v, j = "comparison", value = comp)
            data.table::set(v, j = "peak", value = interval_name)
            data.table::set(v, j = "high.parent",
                            value = GenomicRanges::mcols(pk)$high.parent)
            data.table::set(v, j = "low.parent",
                            value = GenomicRanges::mcols(pk)$low.parent)
            data.table::set(v, j = "peak.pos",
                            value = GenomicRanges::mcols(pk)$peak.pos)
            data.table::set(v, j = "peak.LOD",
                            value = GenomicRanges::mcols(pk)$peak.LOD)
            rows[[length(rows) + 1]] <- v
        }
    }

    if (length(rows) == 0) return(data.table::data.table())
    data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
}

#' Flatten nested QTL CNV results into a single data.table
#'
#' Combines the nested list output of multiple \code{\link{extractQTLCNVs}}
#' calls into a single data.table with peak metadata columns added per row.
#'
#' @param qtl_cnvs Named list of \code{\link{extractQTLCNVs}} results,
#'   typically one per cross/comparison.
#' @return A data.table with CNV call columns (chr, start, end, Strain,
#'   CNV_type, CNV_size, normalized_RD, e-values, q0) plus peak metadata:
#'   comparison, peak, high.parent, low.parent, peak.pos, peak.LOD.
#'   Returns an empty data.table if no CNVs are found.
#' @export
flattenQTLCNVs <- function(qtl_cnvs) {

    if (!requireNamespace("GenomicRanges", quietly = TRUE))
        stop("GenomicRanges is required.")

    rows <- list()
    for (comp in names(qtl_cnvs)) {
        for (interval_name in names(qtl_cnvs[[comp]])) {
            entry <- qtl_cnvs[[comp]][[interval_name]]
            if (is.null(entry$cnvs)) next
            v <- data.table::as.data.table(entry$cnvs)

            pk <- entry$peak
            data.table::set(v, j = "comparison", value = comp)
            data.table::set(v, j = "peak", value = interval_name)
            data.table::set(v, j = "high.parent",
                            value = GenomicRanges::mcols(pk)$high.parent)
            data.table::set(v, j = "low.parent",
                            value = GenomicRanges::mcols(pk)$low.parent)
            data.table::set(v, j = "peak.pos",
                            value = GenomicRanges::mcols(pk)$peak.pos)
            data.table::set(v, j = "peak.LOD",
                            value = GenomicRanges::mcols(pk)$peak.LOD)
            rows[[length(rows) + 1]] <- v
        }
    }

    if (length(rows) == 0) return(data.table::data.table())
    data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
}

#' Load the FungAMR antifungal resistance database
#'
#' Reads the first sheet (Supp_table_1_FungAMR) from the FungAMR Excel file
#' (Bédard et al. 2025) and returns it as a data.table with cleaned column
#' names. The database contains 35,792 curated antifungal resistance mutations
#' across fungal species, with S. cerevisiae systematic names mapped where
#' available.
#'
#' @param path Path to the FungAMR xlsx file.
#'   Default: \code{"/data1/yeast/reference/fungAMR/41564_2025_2084_MOESM3_ESM.xlsx"}
#' @return A data.table with columns including: species, gene_or_protein_name,
#'   drug, mutation, confidence_score, sc_systematic_name, and others.
#' @export
loadFungAMR <- function(path = "/data1/yeast/reference/fungAMR/41564_2025_2084_MOESM3_ESM.xlsx") {

    if (!requireNamespace("readxl", quietly = TRUE))
        stop("readxl is required. Install with: install.packages('readxl')")

    df <- readxl::read_excel(path, sheet = 1)

    # Convert to data.table
    dt <- data.table::as.data.table(df)

    # Clean column names: lowercase, replace spaces/special chars with underscores
    clean_names <- tolower(names(dt))
    clean_names <- gsub("[^a-z0-9_]", "_", clean_names)
    clean_names <- gsub("_+", "_", clean_names)
    clean_names <- gsub("^_|_$", "", clean_names)
    # First column is an unnamed row index from the xlsx
    if (clean_names[1] == "1") clean_names[1] <- "row_id"
    data.table::setnames(dt, clean_names)

    dt
}

#' Load yeast gene coordinates from SGD GFF file
#'
#' Parses an SGD GFF file and returns a GRanges of gene features with
#' systematic names. Cached after first load within a session.
#'
#' @param gff_path Path to the SGD GFF file (plain or gzipped).
#'   Default: \code{"/data1/yeast/reference/yeast_ref/saccharomyces_cerevisiae.gff.gz"}
#' @return GRanges with one range per gene and a \code{systematic_name} metadata column.
#' @export
loadSGDGenes <- function(gff_path = "/data1/yeast/reference/yeast_ref/saccharomyces_cerevisiae.gff.gz") {

    if (!requireNamespace("GenomicRanges", quietly = TRUE) ||
        !requireNamespace("IRanges", quietly = TRUE))
        stop("GenomicRanges and IRanges are required.")

    # Read GFF, skip comment lines
    lines <- readLines(gff_path)
    lines <- lines[!grepl("^#", lines) & nchar(lines) > 0]

    # Parse tab-delimited fields
    fields <- strsplit(lines, "\t")
    types <- vapply(fields, function(x) x[3], character(1))

    # Keep only gene features
    gene_idx <- which(types == "gene")
    if (length(gene_idx) == 0) stop("No gene features found in GFF file")

    chr   <- vapply(fields[gene_idx], function(x) x[1], character(1))
    start <- as.integer(vapply(fields[gene_idx], function(x) x[4], character(1)))
    end   <- as.integer(vapply(fields[gene_idx], function(x) x[5], character(1)))
    strand <- vapply(fields[gene_idx], function(x) x[7], character(1))
    attrs <- vapply(fields[gene_idx], function(x) x[9], character(1))

    # Extract systematic name from ID= attribute
    sys_name <- sub(".*ID=([^;]+).*", "\\1", attrs)

    GenomicRanges::GRanges(
        seqnames = chr,
        ranges = IRanges::IRanges(start = start, end = end),
        strand = strand,
        systematic_name = sys_name
    )
}

#' Annotate QTL peaks with overlapping genes from a gene set
#'
#' Given a GRanges from \code{\link{findQTLPeaks}} and a character vector of
#' yeast systematic names (e.g. from FungAMR), identifies which genes from
#' the set overlap each QTL interval. Returns the peaks GRanges with an
#' added \code{gene_hits} metadata column containing comma-separated
#' systematic names of overlapping genes (or NA if none).
#'
#' @param peaks GRanges object from \code{\link{findQTLPeaks}}.
#' @param gene_set Character vector of yeast systematic names to test
#'   (e.g. \code{c("YER145C", "YLR404W", "YGR060W")}).
#' @param gene_gr GRanges of all yeast genes with a \code{systematic_name}
#'   metadata column. If NULL (default), loaded via \code{\link{loadSGDGenes}}.
#' @param gff_path Path to SGD GFF file, passed to \code{\link{loadSGDGenes}}
#'   if \code{gene_gr} is NULL.
#' @return The input \code{peaks} GRanges with an added \code{gene_hits}
#'   metadata column (character). Intervals with no overlapping genes from
#'   the set have NA.
#' @export
annotateQTLWithGenes <- function(peaks, gene_set, gene_gr = NULL,
                                  gff_path = "/data1/yeast/reference/yeast_ref/saccharomyces_cerevisiae.gff.gz") {

    if (!requireNamespace("GenomicRanges", quietly = TRUE))
        stop("GenomicRanges is required.")

    if (length(peaks) == 0) {
        GenomicRanges::mcols(peaks)$gene_hits <- character(0)
        return(peaks)
    }

    # Load gene coordinates if not provided
    if (is.null(gene_gr)) {
        gene_gr <- loadSGDGenes(gff_path)
    }

    # Subset to genes in the target set
    in_set <- GenomicRanges::mcols(gene_gr)$systematic_name %in% gene_set
    gene_gr_sub <- gene_gr[in_set]

    if (length(gene_gr_sub) == 0) {
        GenomicRanges::mcols(peaks)$gene_hits <- rep(NA_character_, length(peaks))
        return(peaks)
    }

    # Find overlaps: which peaks overlap which genes from the set
    hits <- GenomicRanges::findOverlaps(peaks, gene_gr_sub, ignore.strand = TRUE)
    peak_hits <- S4Vectors::queryHits(hits)
    gene_names <- GenomicRanges::mcols(gene_gr_sub)$systematic_name[S4Vectors::subjectHits(hits)]

    # Collapse gene names per peak
    gene_hits <- rep(NA_character_, length(peaks))
    if (length(hits) > 0) {
        by_peak <- split(gene_names, peak_hits)
        for (pk_idx in names(by_peak)) {
            gene_hits[as.integer(pk_idx)] <- paste(unique(by_peak[[pk_idx]]), collapse = ",")
        }
    }

    GenomicRanges::mcols(peaks)$gene_hits <- gene_hits
    peaks
}

#' Test enrichment of a gene set within QTL intervals by circular permutation
#'
#' Tests whether QTL intervals overlap more genes from a target gene set than
#' expected by chance. The null distribution is generated by circularly
#' permuting QTL intervals around each chromosome, preserving interval sizes
#' and per-chromosome structure. For each permutation a random offset is
#' drawn per chromosome and all intervals on that chromosome are shifted
#' (wrapping at chromosome boundaries). The test statistic is the number of
#' unique genes from the set overlapping any interval.
#'
#' @param peaks GRanges object from \code{\link{findQTLPeaks}}.
#' @param gene_set Character vector of yeast systematic names to test.
#' @param gene_gr GRanges of all yeast genes with a \code{systematic_name}
#'   metadata column. If NULL (default), loaded via \code{\link{loadSGDGenes}}.
#' @param n_perm Number of circular permutations (default 10000).
#' @param chr_sizes Named integer vector of chromosome lengths. If NULL,
#'   uses sacCer3 sizes (chrI–chrXVI).
#' @param gff_path Path to SGD GFF file, passed to \code{\link{loadSGDGenes}}
#'   if \code{gene_gr} is NULL.
#' @return A list with:
#'   \describe{
#'     \item{observed}{Number of unique gene set genes overlapping the peaks.}
#'     \item{observed_genes}{Character vector of the overlapping gene names.}
#'     \item{null_distribution}{Integer vector of length \code{n_perm} with
#'       the count of overlapping genes for each permutation.}
#'     \item{p_value}{Proportion of permutations with count >= observed
#'       (one-sided, greater).}
#'     \item{n_perm}{Number of permutations performed.}
#'     \item{n_gene_set}{Number of genes from gene_set found in gene_gr.}
#'   }
#' @export
testQTLGeneEnrichment <- function(peaks, gene_set, gene_gr = NULL,
                                   n_perm = 10000, chr_sizes = NULL,
                                   gff_path = "/data1/yeast/reference/yeast_ref/saccharomyces_cerevisiae.gff.gz") {

    if (!requireNamespace("GenomicRanges", quietly = TRUE) ||
        !requireNamespace("IRanges", quietly = TRUE))
        stop("GenomicRanges and IRanges are required.")

    # Load gene coordinates if not provided
    if (is.null(gene_gr)) {
        gene_gr <- loadSGDGenes(gff_path)
    }

    # Default sacCer3 chromosome sizes (nuclear only)
    if (is.null(chr_sizes)) {
        chr_sizes <- c(
            chrI = 230218L,   chrII = 813184L,   chrIII = 316620L,
            chrIV = 1531933L, chrV = 576874L,     chrVI = 270161L,
            chrVII = 1090940L, chrVIII = 562643L, chrIX = 439888L,
            chrX = 745751L,   chrXI = 666816L,    chrXII = 1078177L,
            chrXIII = 924431L, chrXIV = 784333L,  chrXV = 1091291L,
            chrXVI = 948066L
        )
    }

    # Subset gene_gr to the target set
    in_set <- GenomicRanges::mcols(gene_gr)$systematic_name %in% gene_set
    gene_gr_sub <- gene_gr[in_set]
    n_gene_set <- length(gene_gr_sub)

    if (length(peaks) == 0 || n_gene_set == 0) {
        return(list(
            observed = 0L,
            observed_genes = character(0),
            null_distribution = integer(n_perm),
            p_value = 1.0,
            n_perm = n_perm,
            n_gene_set = n_gene_set
        ))
    }

    # Observed overlaps
    obs_hits <- GenomicRanges::findOverlaps(peaks, gene_gr_sub, ignore.strand = TRUE)
    observed_genes <- unique(GenomicRanges::mcols(gene_gr_sub)$systematic_name[S4Vectors::subjectHits(obs_hits)])
    observed <- length(observed_genes)

    # Pre-extract peak info for vectorized permutation
    peak_chr <- as.character(GenomicRanges::seqnames(peaks))
    peak_start <- GenomicRanges::start(peaks)
    peak_width <- GenomicRanges::width(peaks)
    n_peaks <- length(peaks)

    # Chromosomes that have peaks and their sizes
    peak_chroms <- unique(peak_chr)
    chr_size_vec <- chr_sizes[peak_chr]  # per-peak chromosome size

    # ── Vectorized circular permutation ──────────────────────────────
    # Build ALL permuted intervals at once (n_perm * n_peaks rows),
    # then do a single findOverlaps call against the gene set.
    # This avoids 10,000 separate GRanges constructions + overlap calls.

    # For each permutation, draw one offset per chromosome, then map to per-peak offsets
    # Matrix: n_perm rows x length(peak_chroms) columns
    offset_mat <- matrix(0L, nrow = n_perm, ncol = length(peak_chroms))
    colnames(offset_mat) <- peak_chroms
    for (j in seq_along(peak_chroms)) {
        offset_mat[, j] <- sample.int(chr_sizes[peak_chroms[j]], n_perm, replace = TRUE)
    }

    # Map per-chromosome offsets to per-peak offsets: n_perm x n_peaks
    chr_idx <- match(peak_chr, peak_chroms)
    # Each row p: offset for peak i is offset_mat[p, chr_idx[i]]
    # Vectorize: expand to full vectors of length n_perm * n_peaks
    perm_idx <- rep(seq_len(n_perm), each = n_peaks)   # permutation ID per row
    peak_idx <- rep(seq_len(n_peaks), times = n_perm)   # peak index per row

    offsets_expanded <- offset_mat[cbind(perm_idx, chr_idx[peak_idx])]
    starts_expanded <- as.integer((peak_start[peak_idx] + offsets_expanded - 1L) %% chr_size_vec[peak_idx] + 1L)
    ends_expanded <- starts_expanded + peak_width[peak_idx] - 1L
    # Clamp to chromosome end for intervals that wrap
    chsize_expanded <- chr_size_vec[peak_idx]
    ends_expanded <- pmin(ends_expanded, chsize_expanded)
    chr_expanded <- peak_chr[peak_idx]

    # Single GRanges for all permuted intervals
    all_perm_gr <- GenomicRanges::GRanges(
        seqnames = chr_expanded,
        ranges = IRanges::IRanges(start = starts_expanded, end = ends_expanded)
    )

    # Single findOverlaps call
    all_hits <- GenomicRanges::findOverlaps(all_perm_gr, gene_gr_sub, ignore.strand = TRUE)
    q_hits <- S4Vectors::queryHits(all_hits)
    s_hits <- S4Vectors::subjectHits(all_hits)

    # Map query hits back to permutation index and count unique genes per permutation
    hit_perm <- perm_idx[q_hits]
    # Unique gene per permutation: paste perm + subject to deduplicate
    perm_gene <- paste(hit_perm, s_hits, sep = "_")
    perm_gene_unique <- !duplicated(perm_gene)
    null_dist <- tabulate(hit_perm[perm_gene_unique], nbins = n_perm)

    p_value <- (sum(null_dist >= observed) + 1) / (n_perm + 1)

    list(
        observed = observed,
        observed_genes = observed_genes,
        null_distribution = null_dist,
        p_value = p_value,
        n_perm = n_perm,
        n_gene_set = n_gene_set
    )
}

#' Extract CNV calls within QTL peak intervals for parent strains
#'
#' Given a GRanges from \code{\link{findQTLPeaks}} and a CNV GRanges from
#' \code{\link{loadCNVnator}}, finds CNV calls that overlap each QTL interval
#' and belong to the parental strains indicated in the peaks metadata.
#' Returns a named list with one element per QTL peak.
#'
#' @param peaks GRanges object from \code{\link{findQTLPeaks}}.
#' @param cnv_gr GRanges from \code{\link{loadCNVnator}} with a \code{Strain}
#'   metadata column.
#' @param strain.names Character vector of strain names to query. If NULL
#'   (default), uses the unique values from \code{high.parent} and
#'   \code{low.parent} metadata columns of peaks.
#' @return A named list with one element per QTL peak. Each element is a list
#'   containing:
#'   \describe{
#'     \item{peak}{A one-row GRanges for the interval.}
#'     \item{cnvs}{A data.frame of CNV calls within the interval for the
#'       parent strains. NULL if no CNVs found.}
#'   }
#'   List names are formatted as \code{chr:start-end}.
#' @export
extractQTLCNVs <- function(peaks, cnv_gr, strain.names = NULL) {

    if (!requireNamespace("GenomicRanges", quietly = TRUE))
        stop("GenomicRanges is required.")

    if (length(peaks) == 0) return(list())

    # Determine strains to query
    if (is.null(strain.names)) {
        strain.names <- unique(c(
            GenomicRanges::mcols(peaks)$high.parent,
            GenomicRanges::mcols(peaks)$low.parent
        ))
    }

    # Subset cnv_gr to parental strains
    cnv_strains <- GenomicRanges::mcols(cnv_gr)$Strain
    in_parents <- cnv_strains %in% strain.names
    if (!any(in_parents)) {
        missing <- setdiff(strain.names, unique(cnv_strains))
        if (length(missing) > 0) {
            message("  Note: ", length(missing), " strain(s) not found in cnv_gr: ",
                    paste(head(missing, 5), collapse = ', '))
        }
    }
    cnv_sub <- cnv_gr[in_parents]

    if (length(cnv_sub) == 0) {
        result <- vector("list", length(peaks))
        peak_names <- paste0(
            as.character(GenomicRanges::seqnames(peaks)), ":",
            GenomicRanges::start(peaks), "-", GenomicRanges::end(peaks)
        )
        for (i in seq_along(peaks)) {
            result[[i]] <- list(peak = peaks[i], cnvs = NULL)
        }
        names(result) <- peak_names
        return(result)
    }

    # Single findOverlaps: peaks vs parental CNVs
    hits <- GenomicRanges::findOverlaps(peaks, cnv_sub, ignore.strand = TRUE)
    peak_hits <- S4Vectors::queryHits(hits)
    cnv_hits <- S4Vectors::subjectHits(hits)

    # Group CNV indices by peak
    hits_by_peak <- split(cnv_hits, peak_hits)

    result <- vector("list", length(peaks))
    peak_names <- paste0(
        as.character(GenomicRanges::seqnames(peaks)), ":",
        GenomicRanges::start(peaks), "-", GenomicRanges::end(peaks)
    )

    for (i in seq_along(peaks)) {
        pk <- peaks[i]
        cnv_idx <- hits_by_peak[[as.character(i)]]

        if (is.null(cnv_idx) || length(cnv_idx) == 0) {
            result[[i]] <- list(peak = pk, cnvs = NULL)
            next
        }

        # Extract matching CNVs as a data.frame
        matched_cnvs <- cnv_sub[cnv_idx]
        cnv_df <- data.frame(
            chr = as.character(GenomicRanges::seqnames(matched_cnvs)),
            start = GenomicRanges::start(matched_cnvs),
            end = GenomicRanges::end(matched_cnvs),
            as.data.frame(GenomicRanges::mcols(matched_cnvs)),
            stringsAsFactors = FALSE
        )

        result[[i]] <- list(peak = pk, cnvs = cnv_df)
    }

    names(result) <- peak_names
    result
}

#' Load CNVnator results and convert to a GRanges object
#'
#' Reads a CNVnator results TSV file (one row per CNV call per strain),
#' parses the coordinates column, converts chromosome names from
#' \code{chromosome1}–\code{chromosome16} to \code{chrI}–\code{chrXVI},
#' and returns a GRanges with all original columns as metadata.
#' Optionally saves the result as a qs2 file.
#'
#' @param path Path to the CNVnator TSV file.
#' @param save_qs2 If not NULL, path to save the GRanges as a qs2 file.
#' @return GRanges with metadata columns: Strain, CNV_type, CNV_size,
#'   normalized_RD, e_val1–e_val4, q0.
#' @export
loadCNVnator <- function(path, save_qs2 = NULL) {

    if (!requireNamespace("GenomicRanges", quietly = TRUE) ||
        !requireNamespace("IRanges", quietly = TRUE))
        stop("GenomicRanges and IRanges are required.")

    dt <- data.table::fread(path, sep = "\t", header = TRUE)

    # Parse coordinates column: "chromosome1:3001-27000"
    coords <- dt$coordinates
    chr_raw <- sub(":.*", "", coords)
    range_str <- sub(".*:", "", coords)
    start_pos <- as.integer(sub("-.*", "", range_str))
    end_pos <- as.integer(sub(".*-", "", range_str))

    # Convert chromosome names: chromosome1 → chrI, chromosome16 → chrXVI
    chr_num <- as.integer(sub("chromosome", "", chr_raw))
    chr_roman <- paste0("chr", as.roman(chr_num))

    gr <- GenomicRanges::GRanges(
        seqnames = chr_roman,
        ranges = IRanges::IRanges(start = start_pos, end = end_pos),
        Strain = dt$Strain,
        CNV_type = dt$CNV_type,
        CNV_size = dt$CNV_size,
        normalized_RD = dt$normalized_RD,
        e_val1 = dt[[6]],
        e_val2 = dt[[7]],
        e_val3 = dt[[8]],
        e_val4 = dt[[9]],
        q0 = dt$q0
    )

    if (!is.null(save_qs2)) {
        if (!requireNamespace("qs2", quietly = TRUE))
            stop("qs2 is required to save. Install with: install.packages('qs2')")
        qs2::qs_save(gr, save_qs2)
        message("Saved GRanges to: ", save_qs2)
    }

    gr
}

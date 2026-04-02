library(xQTLStats)
library(qs2)
library(ggpubr)
library(parallel)
library(pbmcapply)

# ══════════════════════════════════════════════════════════════════════
# Fluconazole xQTL Analysis
# ══════════════════════════════════════════════════════════════════════
# This script processes GATK ASEReadCounter output from bulk-segregant
# yeast experiments and performs QTL mapping for fluconazole resistance.
#
# Workflow:
#   1. Load genotype reference and phase allele counts per sample
#   2. Filter low-coverage samples
#   3. Smooth allele frequencies per chromosome (loess)
#   4. Contrast drug vs control to identify QTL
#   5. Find significant peaks and extract candidate variants
#   6. Extract parental CNVs within QTL intervals
#   7. Annotate peaks with known resistance genes (FungAMR)
#   8. Test enrichment of resistance genes in QTL intervals
# ══════════════════════════════════════════════════════════════════════

# ── Configuration ──────────────────────────────────────────────────
# Sparse genotype matrix (dgCMatrix) for ~1011 yeast strains x ~1.2M SNPs.
# Used to phase allele counts relative to parental genotypes.
Z.qs.file = '/media/hoffman2/PUBLIC_SHARED/yeast/reference/3000/VCFs_SNPs/Z_1011_filtered.qs'
# Hoffman2 path:
# Z.qs.file = '/u/project/kruglyak/PUBLIC_SHARED/yeast/reference/3000/VCFs_SNPs/Z_1011_filtered.qs'

# Pre-annotated variant lookup (coding effects, phyloP, etc.)
# Used downstream for extracting candidate variants within QTL peaks.
vlookup = '/data1/yeast/reference/3000/VCFs_SNPs/variant_lookup.qs'
# Hoffman2 path:
# vlookup = '/u/project/kruglyak/PUBLIC_SHARED/yeast/reference/3000/VCFs_SNPs/variant_lookup.qs'

# Directory containing per-sample GATK ASEReadCounter .txt files
count_dir = '/data2/xQTL/FluconazoleV1_041825/count_files/'

# Sample key file: maps sample IDs to parent strains and conditions
key.file = '/home/jbloom/Dropbox/code/ComplexTraitTools/xQTLStats/examples/FluconazoleV1_041825/flatkey.csv'

# Number of cores for parallel steps (adjust to your machine)
ncores = 16

# ── Load reference data ───────────────────────────────────────────
# Genetic maps for physical→genetic position interpolation
gmaps = readRDS(system.file('reference', 'yeast_gmaps.RDS', package = 'xQTLStats'))
gmap  = gmaps[['A']]
uchr  = paste0('chr', as.character(as.roman(1:16)))

Z.obj   = qs2::qs_read(Z.qs.file)
samples = readr::read_csv(key.file, show_col_types = FALSE)

# ── 1. Pre-extract per-cross Z subsets ────────────────────────────
# The full Z matrix is too large to safely fork via mclapply (copy-on-write
# doesn't help if child processes touch different rows). Extract a minimal
# Z.obj for each unique biparental cross (just the 2 parent rows) before
# forking, then free the full matrix.
unique.crosses = unique(paste(samples$parent_a, samples$parent_alpha, sep = '|'))
cross.Z = lapply(unique.crosses, function(cr) {
    p.names = strsplit(cr, '\\|')[[1]]
    list(
        Z           = Z.obj$Z[p.names, ],
        variant_ids = Z.obj$variant_ids
    )
})
names(cross.Z) = unique.crosses

rm(Z.obj); gc()

# ── 2. Build phased count tables (parallel across samples) ────────
# For each sample, reads the GATK ASEReadCounter output, identifies
# segregating sites between the two parents, and phases allele counts
# so that 'ref' and 'alt' correspond to parent_a and parent_alpha.
countdfs = pbmclapply(1:nrow(samples), function(n) {
    s       = samples$'Sample ID'[n]
    sname   = paste0(s, '.txt')
    p.names = c(samples$parent_a[n], samples$parent_alpha[n])
    cr.key  = paste(p.names, collapse = '|')
    makeCountTablesGATK2(count_dir, sname, p.names, cross.Z[[cr.key]], gmap)
}, mc.cores = ncores)
names(countdfs) = samples$'Sample ID'

# ── 3. Filter low-coverage samples ───────────────────────────────
# Mean per-site read depth < 1 indicates a failed or very low-coverage
# library. These are excluded from downstream analysis.
tdepth.v = sapply(countdfs, function(x) sum(x$ref + x$alt) / nrow(x))
keep     = names(countdfs)[tdepth.v > 1]
cat("Samples passing depth filter:", length(keep), "/", length(countdfs), "\n")

# ── 4. Smooth allele frequency differences (parallel) ────────────
# Per sample: bin allele counts, fit loess per chromosome with
# AICc-optimized span, and compute smoothed AF + standard error.
afds = pbmclapply(keep, function(snn) {
    calcAFD(countdfs[[snn]], experiment.name = snn,
            sample.size = 1e4, sel.strength = 0.9)
}, mc.cores = ncores)
names(afds) = keep

# ── 5. Generate per-sample AF plots (parallel) ───────────────────
plots = pbmclapply(keep, function(snn) {
    plotIndividualExperiment(afds[[snn]], snn)
}, mc.cores = ncores)
names(plots) = keep

# ── 6. Inspect individual samples ────────────────────────────────
# Check which samples passed the depth filter
samples$oldLabel_parent_a[tdepth.v > 1]

# Visualize a single sample's allele frequency across the genome
plots[[7]]

# ── 7a. Single-sample fitness QTL (vs expected AF = 0.5) ─────────
# Test whether a single sample's allele frequency deviates from the
# null expectation of 0.5 (no selection). Useful for detecting
# viability/fitness QTL without a paired control.
resultsFitness = calcContrastStats(
    results = list(afds[[s.L]]),
    L = paste0('_', s.L),
    R = 0.5
)
plotSummarySignedShade(
    resultsFitness, use.LOD = TRUE,
    p1.name = samples$parent_a[samples$"Sample ID" == s.L],
    p2.name = samples$parent_alpha[samples$"Sample ID" == s.L]
)

# ── 7b. Two-sample contrast: drug vs control ─────────────────────
# Compare a drug-treated sample against a matched control from the
# same cross. The contrast removes shared fitness QTL and isolates
# drug-specific resistance loci.
s.L = keep[82]  # drug-treated sample
s.R = keep[87]  # matched control

results = calcContrastStats(
    results = list(afds[[s.L]], afds[[s.R]]),
    L = paste0('_', s.L),
    R = paste0('_', s.R)
)

# ── 8. Visualize contrast results ────────────────────────────────
# Four-panel figure:
#   Top row:    individual AF profiles for drug (L) and control (R)
#   Bottom row: AF difference (contrast) and signed significance
library(patchwork)

# Helper: look up parent names for a sample
p1_of <- function(sid) samples$parent_a[samples$"Sample ID" == sid]
p2_of <- function(sid) samples$parent_alpha[samples$"Sample ID" == sid]

pA = plotIndividualExperiment(afds[[s.L]], s.L,
        p1.name = p1_of(s.L), p2.name = p2_of(s.L))
pB = plotIndividualExperiment(afds[[s.R]], s.R,
        p1.name = p1_of(s.R), p2.name = p2_of(s.R))
pC = plotContrast(results, s.L, s.R,
        p1.name = p1_of(s.L), p2.name = p2_of(s.L))
pD = plotSummarySignedShade(results,
        p1.name = p1_of(s.L), p2.name = p2_of(s.L)) +
    ggtitle(paste(s.L, '-', s.R))

(pA | pB) / (pC | pD)
# ggsave('zWGsquare.png', width = 20, height = 20, units = "cm")

# ── 9. Find QTL peaks ────────────────────────────────────────────
# Identify LOD peaks above threshold and define confidence intervals
# using a LOD drop method. Returns a GRanges with peak position,
# LOD score, and which parent allele is enriched.
qpeaks = findQTLPeaks(results, LOD.threshold = 5, LOD.drop = 3,
    p1.name = p1_of(s.L), p2.name = p2_of(s.L))

# ── 10. Repeat for a second cross ────────────────────────────────
s.L = keep[46]
s.R = keep[49]

results2 = calcContrastStats(
    results = list(afds[[s.L]], afds[[s.R]]),
    L = paste0('_', s.L),
    R = paste0('_', s.R)
)

qpeaks2 = findQTLPeaks(results2, LOD.threshold = 5, LOD.drop = 3,
    p1.name = p1_of(s.L), p2.name = p2_of(s.L))

# ── 11. Extract candidate variants within QTL intervals ──────────
# Collect peaks from multiple crosses into a named list, then query
# the variant lookup for segregating variants within each interval.
# Parent strain names on the GRanges determine which strains are queried;
# strains not in the lookup (e.g. S288C = reference) are skipped.
peaks = list()
peaks[['x1']] = qpeaks
peaks[['x2']] = qpeaks2

# Load the pre-annotated variant lookup
lookup = qs2::qs_read(vlookup)

# Extract variants per cross, per peak interval
qtl_vars = list()
for (i in names(peaks)) {
    qtl_vars[[i]] = extractQTLVariants(peaks[[i]], lookup)
}

# ── 12. Flatten to a single table ────────────────────────────────
# Combine all crosses/intervals into one data.table. Strain genotype
# columns (which differ across crosses) are collapsed into a single
# 'genotypes' string column (e.g. "ACD=2"). Peak metadata (parent
# strains, LOD, position) are added per row.
qtl_vars.df = flattenQTLVariants(qtl_vars, lookup)

# ── 13. Extract parental CNVs within QTL intervals ───────────────
# Load CNVnator structural variant calls for ~3000 yeast strains.
# Chromosome names are converted from chromosome1→chrI format.
# For each QTL interval, extract CNV calls (deletions/duplications)
# belonging to the parental strains of that cross. This identifies
# structural variants that may underlie the QTL signal.
cnv.file = "/media/hoffman2/PUBLIC_SHARED/yeast/reference/3000/full3039Matrix.CNVnatorResults.tsv"
# Hoffman2 path:
# cnv.file = "/u/project/kruglyak/PUBLIC_SHARED/yeast/reference/3000/full3039Matrix.CNVnatorResults.tsv"
cnv_gr = loadCNVnator(cnv.file)

qtl_cnv = list()
for (i in names(peaks)) {
    qtl_cnv[[i]] = extractQTLCNVs(peaks[[i]], cnv_gr)
}

# ── 14. Flatten CNV results to a single table ────────────────────
# Same structure as flattenQTLVariants: one row per CNV call, with
# comparison, peak interval, parent strain info, and LOD added.
qtl_cnv.df = flattenQTLCNVs(qtl_cnv)

# ── 15. Load FungAMR antifungal resistance database ──────────────
# Bédard et al. 2025 — curated database of ~36k fungal antifungal
# resistance mutations with S. cerevisiae systematic names mapped
# where available. The sc_systematic_name column can contain
# comma-separated names when a mutation maps to multiple genes.
fungamr = loadFungAMR()
# Hoffman2 path:
# fungamr = loadFungAMR("/u/project/kruglyak/PUBLIC_SHARED/yeast/reference/fungAMR/41564_2025_2084_MOESM3_ESM.xlsx")

# ── 16. Build azole resistance gene set ──────────────────────────
# Extract all S. cerevisiae systematic names associated with any azole
# drug in FungAMR. The sc_systematic_name field can be comma-separated
# (e.g. "YER145C,YLR404W"), so split and deduplicate.
azole_genes = na.omit(unique(fungamr$sc_systematic_name[grepl('zole|ZOLE', fungamr$drug)]))
azole_genes = unique(trimws(unlist(strsplit(azole_genes, ','))))

# ── 17. Annotate QTL peaks with azole gene overlaps ─────────────
# For each QTL interval, check whether any known azole resistance genes
# fall within the interval. Adds a 'gene_hits' column to the GRanges
# with comma-separated systematic names (NA if none).
# On Hoffman2, pass gff_path explicitly:
# qpeaks.annotated = annotateQTLWithGenes(qpeaks, azole_genes,
#     gff_path = '/u/project/kruglyak/PUBLIC_SHARED/yeast/reference/yeast_ref/saccharomyces_cerevisiae.gff.gz')
qpeaks.annotated = annotateQTLWithGenes(qpeaks, azole_genes)

# This can also be applied across crosses by annotating each peak set:
# qpeaks2.annotated = annotateQTLWithGenes(qpeaks2, azole_genes)

# ── 18. Test enrichment of azole genes in QTL intervals ──────────
# Circular genome permutation test: are QTL intervals enriched for
# azole resistance genes beyond what's expected by chance?
# Strategy: circularly shift all intervals per chromosome by a random
# offset (preserving interval sizes and number), count gene set overlaps,
# repeat n_perm times to build a null distribution.
# Returns observed count, overlapping gene names, null distribution,
# and a one-sided p-value.
enrich_result = testQTLGeneEnrichment(qpeaks.annotated, azole_genes, n_perm = 10000)

cat("Observed azole genes in QTL intervals:", enrich_result$observed,
    "out of", enrich_result$n_gene_set, "in gene set\n")
cat("Genes:", paste(enrich_result$observed_genes, collapse = ", "), "\n")
cat("Permutation p-value:", enrich_result$p_value, "\n")

# Visualize the null distribution
# hist(enrich_result$null_distribution, main = "Azole gene enrichment",
#      xlab = "Number of azole genes overlapping permuted intervals")
# abline(v = enrich_result$observed, col = "red", lwd = 2)

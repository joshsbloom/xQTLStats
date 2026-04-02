## xQTLStats: R package for simulation and analysis of bulk-segregant QTL experiments in yeast

xQTLStats supports the full workflow for extreme QTL (xQTL) mapping in yeast biparental crosses: simulating haploid cross progeny with realistic meiotic recombination, phasing bulk-segregant sequencing data, performing QTL mapping via allele frequency difference (AFD) statistics with loess smoothing, and annotating QTL intervals with candidate variants, structural variants, and known resistance gene sets.

### Installation

```r
devtools::install_github("joshsbloom/xQTLStats", ref = "main")
```

### Core Modules

| File | Layer | Key Functions |
|---|---|---|
| `R/xQTL_sim.R` | Simulation | `getCrossVCF`, `simHaploidSegsFN`, `simXQTLExperiment`, `makeCountTablesGATK`, `makeCountTablesGATK2` |
| `R/xQTL_stats.R` | Statistical analysis | `calcAFD`, `calcMetaAFD`, `calcContrastStats`, `findQTLPeaks`, `plotSummary`, `plotSummarySignedShade`, `testQTLGeneEnrichment` |
| `R/xQTL_analysis.R` | AlphaSimR integration | `createFounderPop` |
| `R/query_variants.R` | Variant annotation | `query_variants`, `extractQTLVariants`, `extractQTLCNVs`, `annotateQTLWithGenes`, `loadFungAMR`, `loadCNVnator` |

### Analysis Workflow

1. **Align + count**: Nextflow pipeline (`bwa mem` + `sambamba` + `GATK ASEReadCounter`)
2. **Phase**: Assign allele counts to parental haplotypes using a sparse genotype matrix (`makeCountTablesGATK2`)
3. **Smooth**: Loess-smoothed allele frequency differences per chromosome (`calcAFD`)
4. **Contrast**: Compare drug vs control (or vs expected AF) to compute LOD scores (`calcContrastStats`)
5. **Peak calling**: Identify QTL peaks with LOD drop confidence intervals (`findQTLPeaks`)
6. **Variant extraction**: Query segregating SNPs/indels within QTL intervals for parental strains (`extractQTLVariants`)
7. **CNV extraction**: Query parental structural variants within QTL intervals (`extractQTLCNVs`)
8. **Gene set annotation**: Annotate peaks with overlapping genes from curated sets like FungAMR (`annotateQTLWithGenes`)
9. **Enrichment testing**: Circular genome permutation test for gene set enrichment in QTL intervals (`testQTLGeneEnrichment`)

### Examples

- [FluconazoleV1_041825](examples/FluconazoleV1_041825) — Full worked example: Nextflow alignment pipeline, R analysis script covering steps 1–9 above for fluconazole resistance mapping in biparental yeast crosses
- [do_simulations.R](examples/do_simulations.R) — Simulating xQTL experiments with known genetic architectures

### Key Dependencies

- **data.table**, **Matrix**: core data structures for genotype matrices and variant tables
- **ggplot2**, **ggpubr**, **patchwork**: visualization
- **GenomicRanges**, **IRanges** (Suggests): interval operations for peak calling, variant extraction, and enrichment testing
- **qs2** (Suggests): fast serialization for large genotype objects
- **Meiosis** (Suggests): meiotic crossover simulation
- **readxl** (Suggests): loading FungAMR database

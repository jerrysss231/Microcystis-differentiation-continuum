# Population-genetic statistics for Microcystis with PopGenome
#
# This workflow follows the PopGenome analysis used by Stanojković et al.
# (Nature Communications, 2024; doi:10.1038/s41467-024-46459-6), adapted
# to the C1-C5 Microcystis lineage framework and the project directory layout.

suppressPackageStartupMessages({
  library(PopGenome)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

PROJECT_ROOT <- normalizePath(
  Sys.getenv("MICROCYSTIS_PROJECT_ROOT", unset = "."),
  mustWork = FALSE
)

data_path <- function(...) file.path(PROJECT_ROOT, "data", ...)
result_path <- function(...) file.path(PROJECT_ROOT, "results", ...)

vcf_file <- data_path("snp", "core_alignment_final.fasta.vcf")
lineage_file <- result_path(
  "population_structure",
  "Master_Clade_Assignments.csv"
)
out_dir <- result_path("population_genetics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(vcf_file, lineage_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing input file(s):\n", paste(missing_files, collapse = "\n"))
}

clean_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("_genomic_fna$", "", x)
  x <- sub("\\.fna$", "", x)
  x <- sub("\\.fa(sta)?$", "", x)
  x
}

read_vcf_header <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(close(con), add = TRUE)

  contig_lines <- character(0)
  sample_names <- NULL

  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "##contig=<")) {
      contig_lines <- c(contig_lines, line)
    }
    if (startsWith(line, "#CHROM")) {
      fields <- strsplit(line, "\\t", fixed = FALSE)[[1]]
      sample_names <- fields[-seq_len(9)]
      break
    }
  }

  if (length(contig_lines) != 1) {
    stop(
      "Expected one VCF contig in the header, found ",
      length(contig_lines), "."
    )
  }
  if (is.null(sample_names) || length(sample_names) == 0) {
    stop("No VCF sample names found in the #CHROM header line.")
  }

  contig_id <- sub(".*ID=([^,>]+).*", "\\1", contig_lines)
  length_text <- sub(".*length=([0-9]+).*", "\\1", contig_lines)
  if (identical(length_text, contig_lines)) {
    stop("VCF contig header does not contain a sequence length.")
  }

  list(
    contig = contig_id,
    length = as.numeric(length_text),
    samples = sample_names
  )
}

make_windows <- function(sequence_length, width, jump) {
  start <- seq(1, sequence_length, by = jump)
  stop <- start + width
  keep <- stop < sequence_length
  tibble(
    start = start[keep],
    stop = stop[keep],
    mid = start[keep] + width / 2
  )
}

add_windows <- function(values, windows, label) {
  values <- as.data.frame(values, check.names = FALSE)
  if (nrow(values) != nrow(windows)) {
    stop(
      label, " returned ", nrow(values), " windows, but ",
      nrow(windows), " window coordinates were generated."
    )
  }
  bind_cols(windows, values)
}

pair_names <- function(x) {
  x <- gsub("/", "_", x, fixed = TRUE)
  x <- gsub("Clade_", "C", x, fixed = TRUE)
  x
}

vcf_header <- read_vcf_header(vcf_file)
lineages <- read.csv(lineage_file, stringsAsFactors = FALSE)
if (!all(c("Strain_ID", "Clade") %in% names(lineages))) {
  stop("Lineage file must contain Strain_ID and Clade columns.")
}

lineages <- lineages %>%
  transmute(
    Strain_ID = clean_id(Strain_ID),
    Lineage = sub("^Clade_", "C", as.character(Clade))
  ) %>%
  distinct(Strain_ID, .keep_all = TRUE)

vcf_samples <- tibble(
  VCF_sample = vcf_header$samples,
  Strain_ID = clean_id(vcf_header$samples)
)

sample_map <- vcf_samples %>%
  left_join(lineages, by = "Strain_ID")

if (any(is.na(sample_map$Lineage))) {
  missing_ids <- sample_map$VCF_sample[is.na(sample_map$Lineage)]
  stop(
    "No lineage assignment for VCF sample(s): ",
    paste(missing_ids, collapse = ", ")
  )
}

population_list <- split(sample_map$VCF_sample, sample_map$Lineage)
population_list <- population_list[
  order(as.numeric(sub("^C", "", names(population_list))))
]

cat("Reading VCF with PopGenome...\n")
genome <- readVCF(
  vcf_file,
  numcols = 50000,
  tid = vcf_header$contig,
  frompos = 1,
  topos = vcf_header$length,
  include.unknown = TRUE
)
genome <- set.populations(genome, population_list, diploid = FALSE)

# Genome-wide diversity and differentiation: 50-kb windows, 12.5-kb step.
window_50 <- 50000
jump_50 <- 12500
windows_50 <- make_windows(vcf_header$length, window_50, jump_50)

cat("Calculating nucleotide diversity, FST and Dxy...\n")
genome_50 <- sliding.window.transform(
  genome,
  width = window_50,
  jump = jump_50,
  type = 2
)
fst_stats <- F_ST.stats(genome_50, mode = "nucleotide")

pi_50 <- fst_stats@nuc.diversity.within / window_50
fst_50 <- t(fst_stats@nuc.F_ST.pairwise)
dxy_50 <- t(fst_stats@nuc.diversity.between / window_50)

colnames(pi_50) <- sub("^pop[._]?", "C", colnames(pi_50), ignore.case = TRUE)
colnames(fst_50) <- pair_names(colnames(fst_50))
colnames(dxy_50) <- pair_names(colnames(dxy_50))

pi_windows <- add_windows(pi_50, windows_50, "Nucleotide diversity")
fst_windows <- add_windows(fst_50, windows_50, "FST")
dxy_windows <- add_windows(dxy_50, windows_50, "Dxy")

write.csv(pi_windows, file.path(out_dir, "nucleotide_diversity_50kb.csv"), row.names = FALSE)
write.csv(fst_windows, file.path(out_dir, "pairwise_FST_50kb.csv"), row.names = FALSE)
write.csv(dxy_windows, file.path(out_dir, "pairwise_Dxy_50kb.csv"), row.names = FALSE)

pi_long <- pi_windows %>%
  pivot_longer(-c(start, stop, mid), names_to = "Lineage", values_to = "pi") %>%
  filter(is.finite(pi))

pi_summary <- pi_long %>%
  group_by(Lineage) %>%
  summarise(
    n_windows = n(),
    mean_pi = mean(pi),
    median_pi = median(pi),
    .groups = "drop"
  )
write.csv(pi_summary, file.path(out_dir, "nucleotide_diversity_summary.csv"), row.names = FALSE)

kw_pi <- kruskal.test(pi ~ Lineage, data = pi_long)
writeLines(
  capture.output(kw_pi),
  file.path(out_dir, "nucleotide_diversity_kruskal_wallis.txt")
)

pairwise_summary <- tibble(
  Pair = intersect(colnames(fst_50), colnames(dxy_50)),
  mean_FST = colMeans(fst_50[, intersect(colnames(fst_50), colnames(dxy_50)), drop = FALSE], na.rm = TRUE),
  mean_Dxy = colMeans(dxy_50[, intersect(colnames(fst_50), colnames(dxy_50)), drop = FALSE], na.rm = TRUE)
)
write.csv(
  pairwise_summary,
  file.path(out_dir, "pairwise_FST_Dxy_summary.csv"),
  row.names = FALSE
)

# Neutrality statistics: 10-kb windows, 2.5-kb step.
window_10 <- 10000
jump_10 <- 2500
windows_10 <- make_windows(vcf_header$length, window_10, jump_10)

cat("Calculating Tajima's D and Fu and Li's F...\n")
genome_10 <- sliding.window.transform(
  genome,
  width = window_10,
  jump = jump_10,
  type = 2
)
genome_10 <- neutrality.stats(genome_10)

tajima_10 <- genome_10@Tajima.D
fu_li_f_10 <- genome_10@Fu.Li.F
colnames(tajima_10) <- names(population_list)
colnames(fu_li_f_10) <- names(population_list)

tajima_windows <- add_windows(tajima_10, windows_10, "Tajima's D")
fu_li_f_windows <- add_windows(fu_li_f_10, windows_10, "Fu and Li's F")

write.csv(tajima_windows, file.path(out_dir, "tajima_D_10kb.csv"), row.names = FALSE)
write.csv(fu_li_f_windows, file.path(out_dir, "fu_li_F_10kb.csv"), row.names = FALSE)

neutrality_tests <- bind_rows(lapply(names(population_list), function(lineage) {
  tajima_values <- tajima_windows[[lineage]]
  fu_values <- fu_li_f_windows[[lineage]]
  tajima_values <- tajima_values[is.finite(tajima_values)]
  fu_values <- fu_values[is.finite(fu_values)]

  tajima_test <- t.test(tajima_values, mu = 0)
  fu_test <- t.test(fu_values, mu = 0)

  tibble(
    Lineage = lineage,
    Statistic = c("Tajima_D", "Fu_Li_F"),
    n_windows = c(length(tajima_values), length(fu_values)),
    mean = c(mean(tajima_values), mean(fu_values)),
    t = c(unname(tajima_test$statistic), unname(fu_test$statistic)),
    df = c(unname(tajima_test$parameter), unname(fu_test$parameter)),
    p_value = c(tajima_test$p.value, fu_test$p.value)
  )
}))
write.csv(
  neutrality_tests,
  file.path(out_dir, "neutrality_one_sample_t_tests.csv"),
  row.names = FALSE
)

write.csv(sample_map, file.path(out_dir, "sample_lineage_map.csv"), row.names = FALSE)
cat("PopGenome analysis complete. Results saved to: ", out_dir, "\n", sep = "")

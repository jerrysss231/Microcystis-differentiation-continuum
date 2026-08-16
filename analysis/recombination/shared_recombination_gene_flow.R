PROJECT_ROOT <- normalizePath(
  Sys.getenv("MICROCYSTIS_PROJECT_ROOT", unset = "."),
  mustWork = FALSE
)

data_path <- function(...) file.path(PROJECT_ROOT, "data", ...)
result_path <- function(...) file.path(PROJECT_ROOT, "results", ...)

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(readxl)
  library(IRanges)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(tibble)
  library(scales)
})

gff_path <- data_path("recombination", "gubbins.recombination_predictions.gff")
clade_path <- result_path("population_structure", "Master_Clade_Assignments.csv")
stats_path <- data_path("recombination", "Gubbins_Recombination_Stats_140.xlsx")

out_dir <- result_path("shared_recombination")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

strict_within_requires_other_same_lineage <- TRUE
base_family <- "Times New Roman"

clean_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("_genomic_fna$", "", x)
  x <- gsub("\\.fna$", "", x)
  x <- gsub("\\.fa$", "", x)
  x
}

extract_attr <- function(attr, key) {
  m <- str_match(attr, paste0("(^|;)", key, "=\"?([^\";]+)\"?"))
  m[, 3]
}

parse_taxa <- function(attr) {
  taxa_raw <- extract_attr(attr, "taxa")
  if (length(taxa_raw) == 0 || is.na(taxa_raw) || taxa_raw == "") return(character(0))
  taxa <- unlist(strsplit(taxa_raw, "[,[:space:]]+"))
  taxa <- clean_id(taxa[taxa != ""])
  unique(taxa)
}

order_lineages <- function(x) {
  num <- suppressWarnings(as.numeric(gsub("[^0-9]", "", x)))
  x[order(num, x)]
}

calc_unique_bases_from_df <- function(starts, ends) {
  if (length(starts) == 0) return(0)
  ir <- IRanges::IRanges(start = as.integer(starts), end = as.integer(ends))
  sum(BiocGenerics::width(IRanges::reduce(ir)))
}

clade_map <- read.csv(clade_path, stringsAsFactors = FALSE)
clade_candidates <- c("Clade", "fastBAPS_Optimised", "fastBAPS_Classic")
clade_col <- clade_candidates[clade_candidates %in% colnames(clade_map)][1]
if (length(clade_col) == 0) stop("No lineage-assignment column found in clade map.")
clade_map$Clade <- clade_map[[clade_col]]

required_clade_cols <- c("Strain_ID", "Clade")
missing_clade_cols <- setdiff(required_clade_cols, colnames(clade_map))
if (length(missing_clade_cols) > 0) stop("Missing lineage-map columns: ", paste(missing_clade_cols, collapse = ", "))

clade_map <- clade_map %>%
  mutate(Strain_clean = clean_id(Strain_ID)) %>%
  filter(!is.na(Strain_clean), Strain_clean != "", !is.na(Clade), Clade != "") %>%
  distinct(Strain_clean, Clade)

gubbins_stats <- read_excel(stats_path)
names(gubbins_stats) <- make.names(names(gubbins_stats))
required_stats_cols <- c("Node", "Genome.Length")
missing_stats_cols <- setdiff(required_stats_cols, colnames(gubbins_stats))
if (length(missing_stats_cols) > 0) stop("Missing recombination-statistics columns: ", paste(missing_stats_cols, collapse = ", "))

gubbins_stats <- gubbins_stats %>%
  mutate(Node_clean = clean_id(Node), Genome.Length = as.numeric(Genome.Length)) %>%
  filter(!is.na(Node_clean), Node_clean != "", !is.na(Genome.Length), Genome.Length > 0)

genome_len_lookup <- setNames(gubbins_stats$Genome.Length, gubbins_stats$Node_clean)
valid_strains <- intersect(clade_map$Strain_clean, names(genome_len_lookup))
if (length(valid_strains) == 0) stop("No strains are shared between lineage and genome-length inputs.")

clade_map2 <- clade_map %>% filter(Strain_clean %in% valid_strains)
strain2clade <- setNames(clade_map2$Clade, clade_map2$Strain_clean)
all_strains <- names(strain2clade)
all_lineages <- order_lineages(unique(clade_map2$Clade))

cat("Matched strains:", length(all_strains), "\n")
cat("Lineages:", length(all_lineages), "\n")

gff <- suppressWarnings(read_tsv(gff_path, comment = "#", col_names = FALSE, show_col_types = FALSE))
if (ncol(gff) < 9) stop("Gubbins GFF must contain at least nine columns.")
gff <- gff[, 1:9]
colnames(gff) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

gff2 <- gff %>%
  mutate(event_id = row_number(), start = as.numeric(start), end = as.numeric(end), taxa = map(attributes, parse_taxa)) %>%
  mutate(taxa = map(taxa, ~ intersect(.x, all_strains)), n_taxa = map_int(taxa, length)) %>%
  filter(!is.na(start), !is.na(end), end >= start, n_taxa >= 1)

within_collector <- vector("list", nrow(gff2))
outside_collector <- vector("list", nrow(gff2))
pair_collector <- vector("list", nrow(gff2))

for (i in seq_len(nrow(gff2))) {
  ev_start <- gff2$start[i]
  ev_end <- gff2$end[i]
  taxa_i <- gff2$taxa[[i]]
  if (length(taxa_i) == 0) next

  w_strains <- character(0)
  o_strains <- character(0)
  p_strains <- character(0)
  p_lineages <- character(0)

  for (s in taxa_i) {
    Ls <- strain2clade[[s]]
    other_taxa <- setdiff(taxa_i, s)
    if (length(other_taxa) == 0) next
    other_lineages <- strain2clade[other_taxa]
    same_lineage_taxa <- other_taxa[other_lineages == Ls]
    diff_lineages <- unique(other_lineages[other_lineages != Ls])

    within_flag <- if (strict_within_requires_other_same_lineage) length(same_lineage_taxa) > 0 else any(other_lineages == Ls, na.rm = TRUE)
    if (within_flag) {
      w_strains <- c(w_strains, s)
      p_strains <- c(p_strains, s)
      p_lineages <- c(p_lineages, Ls)
    }
    if (length(diff_lineages) > 0) {
      o_strains <- c(o_strains, s)
      for (Lj in diff_lineages) {
        p_strains <- c(p_strains, s)
        p_lineages <- c(p_lineages, Lj)
      }
    }
  }

  if (length(w_strains) > 0) within_collector[[i]] <- data.frame(strain = w_strains, start = ev_start, end = ev_end)
  if (length(o_strains) > 0) outside_collector[[i]] <- data.frame(strain = o_strains, start = ev_start, end = ev_end)
  if (length(p_strains) > 0) pair_collector[[i]] <- data.frame(strain = p_strains, target_lineage = p_lineages, start = ev_start, end = ev_end)
}

within_df <- bind_rows(within_collector)
outside_df <- bind_rows(outside_collector)
pair_df <- bind_rows(pair_collector)

within_bp_tbl <- within_df %>% group_by(strain) %>% summarise(Within_bp = calc_unique_bases_from_df(start, end), .groups = "drop")
outside_bp_tbl <- outside_df %>% group_by(strain) %>% summarise(Outside_bp = calc_unique_bases_from_df(start, end), .groups = "drop")

results <- data.frame(
  Strain = all_strains,
  Clade = strain2clade[all_strains],
  Genome_Length = genome_len_lookup[all_strains],
  stringsAsFactors = FALSE
) %>%
  left_join(within_bp_tbl, by = c("Strain" = "strain")) %>%
  left_join(outside_bp_tbl, by = c("Strain" = "strain")) %>%
  mutate(
    Within_bp = ifelse(is.na(Within_bp), 0, Within_bp),
    Outside_bp = ifelse(is.na(Outside_bp), 0, Outside_bp),
    Within = Within_bp / Genome_Length * 100,
    Outside = Outside_bp / Genome_Length * 100
  )

lineage_summary <- results %>%
  group_by(Clade) %>%
  summarise(
    Count = n(),
    Average_within = mean(Within, na.rm = TRUE),
    Average_outside = mean(Outside, na.rm = TRUE),
    Within_min = min(Within, na.rm = TRUE),
    Within_max = max(Within, na.rm = TRUE),
    Outside_min = min(Outside, na.rm = TRUE),
    Outside_max = max(Outside, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(results, file.path(out_dir, "Shared_Recombination_per_strain.csv"), row.names = FALSE)
write.csv(lineage_summary, file.path(out_dir, "Shared_Recombination_lineage_summary.csv"), row.names = FALSE)

pair_bp_tbl <- pair_df %>%
  group_by(strain, target_lineage) %>%
  summarise(Shared_bp = calc_unique_bases_from_df(start, end), .groups = "drop")

strain_lineage_pct <- expand.grid(Strain = all_strains, Target_Lineage = all_lineages, stringsAsFactors = FALSE) %>%
  left_join(pair_bp_tbl, by = c("Strain" = "strain", "Target_Lineage" = "target_lineage")) %>%
  mutate(
    Shared_bp = ifelse(is.na(Shared_bp), 0, Shared_bp),
    Genome_Length = genome_len_lookup[Strain],
    Shared_pct = Shared_bp / Genome_Length * 100,
    Source_Lineage = strain2clade[Strain]
  )

pair_mat <- matrix(0, nrow = length(all_lineages), ncol = length(all_lineages), dimnames = list(all_lineages, all_lineages))
for (Li in all_lineages) {
  for (Lj in all_lineages) {
    if (Li == Lj) {
      pair_mat[Li, Lj] <- strain_lineage_pct %>%
        filter(Source_Lineage == Li, Target_Lineage == Li) %>%
        summarise(v = mean(Shared_pct, na.rm = TRUE)) %>% pull(v)
    } else {
      v1 <- strain_lineage_pct %>% filter(Source_Lineage == Li, Target_Lineage == Lj) %>% summarise(v = mean(Shared_pct, na.rm = TRUE)) %>% pull(v)
      v2 <- strain_lineage_pct %>% filter(Source_Lineage == Lj, Target_Lineage == Li) %>% summarise(v = mean(Shared_pct, na.rm = TRUE)) %>% pull(v)
      pair_mat[Li, Lj] <- mean(c(v1, v2), na.rm = TRUE)
    }
  }
}

write.csv(as.data.frame(pair_mat) %>% rownames_to_column("Lineage"), file.path(out_dir, "Shared_Recombination_pairwise_matrix.csv"), row.names = FALSE)
resist_mat <- 100 - pair_mat
diag(resist_mat) <- 0
write.csv(as.data.frame(resist_mat) %>% rownames_to_column("Lineage"), file.path(out_dir, "GeneFlow_resistance_matrix.csv"), row.names = FALSE)

pair_long <- as.data.frame(pair_mat) %>%
  rownames_to_column("Lineage1") %>%
  pivot_longer(-Lineage1, names_to = "Lineage2", values_to = "SharedFraction")
p2 <- ggplot(pair_long, aes(x = Lineage2, y = Lineage1, size = SharedFraction, color = SharedFraction)) +
  geom_point(alpha = 0.85) +
  scale_size_continuous(range = c(1, 12)) +
  labs(x = NULL, y = NULL, title = "Gene flow among lineages") +
  theme_bw(base_size = 14) + theme(panel.grid = element_blank())
ggsave(file.path(out_dir, "shared_recombination_bubble_plot.png"), p2, width = 8, height = 7, dpi = 300)

plot_df <- results %>%
  select(Strain, Clade, Within, Outside) %>%
  pivot_longer(cols = c(Within, Outside), names_to = "Type", values_to = "GenomeFraction") %>%
  mutate(Type = factor(Type, levels = c("Outside", "Within"), labels = c("Between lineages", "Within lineages")))

wilcox_test <- wilcox.test(results$Outside, results$Within)
p_label <- if (wilcox_test$p.value < 0.001) "p < 0.001" else paste0("p = ", formatC(wilcox_test$p.value, format = "f", digits = 3))

p_box <- ggplot(plot_df, aes(x = Clade, y = GenomeFraction, fill = Clade)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ Type, scales = "free") +
  labs(x = NULL, y = "Genome fraction subjected to recombination (%)", caption = paste0("Wilcoxon rank-sum test: ", p_label)) +
  theme_classic(base_size = 12, base_family = base_family) + theme(legend.position = "none")
ggsave(file.path(out_dir, "within_between_recombination_boxplot.png"), p_box, width = 10, height = 5, dpi = 300)

cat("Shared-recombination analysis complete. Results saved to: ", out_dir, "\n")

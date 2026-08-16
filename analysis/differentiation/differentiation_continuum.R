PROJECT_ROOT <- normalizePath(
  Sys.getenv("MICROCYSTIS_PROJECT_ROOT", unset = "."),
  mustWork = FALSE
)

data_path <- function(...) file.path(PROJECT_ROOT, "data", ...)
result_path <- function(...) file.path(PROJECT_ROOT, "results", ...)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(forcats)
})

infile <- data_path("continuum", "pairwise_continuum_master_table.csv")
outdir <- result_path("continuum")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "pairwise_metric_matrices"), showWarnings = FALSE, recursive = TRUE)

df <- read.csv(infile, stringsAsFactors = FALSE, check.names = FALSE)
cat("Input master table:\n")
print(df)

required_cols <- c(
  "Pair", "mean_FST", "mean_Dxy", "mean_ANI_distance",
  "mean_genecontent_distance", "mean_gene_flow"
)
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

parse_pair <- function(x) {
  parts <- strsplit(x, "_")[[1]]
  if (length(parts) != 2) stop("Bad Pair format: ", x)
  parts
}

get_all_clades <- function(pair_vec) {
  clades <- unique(unlist(strsplit(pair_vec, "_")))
  clades[order(as.numeric(gsub("[^0-9]", "", clades)), clades)]
}

build_pair_matrix <- function(df, value_col, pair_col = "Pair", clades = NULL, diag_value = 0) {
  if (is.null(clades)) clades <- get_all_clades(df[[pair_col]])
  mat <- matrix(
    NA_real_, nrow = length(clades), ncol = length(clades),
    dimnames = list(clades, clades)
  )
  diag(mat) <- diag_value
  for (i in seq_len(nrow(df))) {
    pp <- parse_pair(df[[pair_col]][i])
    mat[pp[1], pp[2]] <- df[[value_col]][i]
    mat[pp[2], pp[1]] <- df[[value_col]][i]
  }
  mat
}

upper_vec <- function(mat) mat[upper.tri(mat)]
permute_matrix_labels <- function(mat, perm) mat[perm, perm]

matrix_perm_correlation <- function(mat1, mat2, method = "spearman", nperm = 9999, seed = 123) {
  stopifnot(all(dim(mat1) == dim(mat2)))
  stopifnot(all(rownames(mat1) == rownames(mat2)))
  stopifnot(all(colnames(mat1) == colnames(mat2)))
  set.seed(seed)
  obs <- suppressWarnings(cor(
    upper_vec(mat1), upper_vec(mat2), method = method,
    use = "pairwise.complete.obs"
  ))
  labs <- rownames(mat2)
  null_vals <- numeric(nperm)
  for (i in seq_len(nperm)) {
    perm <- sample(labs, length(labs), replace = FALSE)
    mat2_perm <- permute_matrix_labels(mat2, perm)
    null_vals[i] <- suppressWarnings(cor(
      upper_vec(mat1), upper_vec(mat2_perm), method = method,
      use = "pairwise.complete.obs"
    ))
  }
  p_two <- (sum(abs(null_vals) >= abs(obs)) + 1) / (nperm + 1)
  list(observed = obs, p_value = p_two, null = null_vals)
}

safe_scale <- function(x) {
  if (sd(x, na.rm = TRUE) == 0 || all(is.na(x))) rep(0, length(x)) else as.numeric(scale(x))
}

safe_rescale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (is.na(rng[1]) || is.na(rng[2]) || diff(rng) == 0) rep(0.5, length(x)) else (x - rng[1]) / diff(rng)
}

theme_paper <- function(base_family = "", base_size = 12) {
  theme_bw(base_family = base_family, base_size = base_size) +
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.border = element_rect(linewidth = 0.6, color = "black"),
      axis.line = element_blank(),
      axis.ticks = element_line(linewidth = 0.4, color = "black"),
      axis.ticks.length = unit(0.15, "cm"),
      axis.text = element_text(color = "black"), axis.title = element_text(face = "bold"),
      legend.position = "right", legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      plot.margin = margin(8, 10, 8, 8)
    )
}

plot_null_hist <- function(null_vals, observed, title, outfile) {
  p <- ggplot(data.frame(x = null_vals), aes(x = x)) +
    geom_histogram(bins = 40, fill = "grey75", color = "white") +
    geom_vline(xintercept = observed, color = "firebrick", linewidth = 1) +
    labs(x = "Permuted correlation", y = "Count", title = title,
         subtitle = paste0("Observed rho = ", round(observed, 3))) +
    theme_paper(base_family = "", base_size = 12)
  ggsave(outfile, p, width = 6.5, height = 4.8, dpi = 300)
}

all_clades <- get_all_clades(df$Pair)
FST_mat <- build_pair_matrix(df, "mean_FST", clades = all_clades, diag_value = 0)
Dxy_mat <- build_pair_matrix(df, "mean_Dxy", clades = all_clades, diag_value = 0)
ANI_mat <- build_pair_matrix(df, "mean_ANI_distance", clades = all_clades, diag_value = 0)
GCD_mat <- build_pair_matrix(df, "mean_genecontent_distance", clades = all_clades, diag_value = 0)
GF_mat <- build_pair_matrix(df, "mean_gene_flow", clades = all_clades, diag_value = 0)
RESIST_mat <- max(GF_mat, na.rm = TRUE) - GF_mat
diag(RESIST_mat) <- 0

matrix_outputs <- list(
  FST_matrix = FST_mat,
  Dxy_matrix = Dxy_mat,
  ANI_distance_matrix = ANI_mat,
  GeneContent_distance_matrix = GCD_mat,
  GeneFlow_matrix = GF_mat,
  GeneFlow_resistance_matrix = RESIST_mat
)
for (nm in names(matrix_outputs)) {
  write.csv(
    as.data.frame(matrix_outputs[[nm]]) %>% rownames_to_column("Clade"),
    file.path(outdir, "pairwise_metric_matrices", paste0(nm, ".csv")),
    row.names = FALSE
  )
}

res_FST_RESIST <- matrix_perm_correlation(FST_mat, RESIST_mat, nperm = 9999, seed = 1)
res_Dxy_RESIST <- matrix_perm_correlation(Dxy_mat, RESIST_mat, nperm = 9999, seed = 2)
res_ANI_RESIST <- matrix_perm_correlation(ANI_mat, RESIST_mat, nperm = 9999, seed = 3)
perm_results <- data.frame(
  Comparison = c("FST_vs_GeneFlowResistance", "Dxy_vs_GeneFlowResistance", "ANI_vs_GeneFlowResistance"),
  rho = c(res_FST_RESIST$observed, res_Dxy_RESIST$observed, res_ANI_RESIST$observed),
  p_value = c(res_FST_RESIST$p_value, res_Dxy_RESIST$p_value, res_ANI_RESIST$p_value)
)
write.csv(perm_results, file.path(outdir, "permutation_matrix_correlation_results.csv"), row.names = FALSE)

plot_null_hist(res_FST_RESIST$null, res_FST_RESIST$observed,
               "Permutation test: FST vs gene-flow resistance",
               file.path(outdir, "permutation_null_distribution_FST_vs_RESIST.png"))
plot_null_hist(res_Dxy_RESIST$null, res_Dxy_RESIST$observed,
               "Permutation test: Dxy vs gene-flow resistance",
               file.path(outdir, "permutation_null_distribution_Dxy_vs_RESIST.png"))
plot_null_hist(res_ANI_RESIST$null, res_ANI_RESIST$observed,
               "Permutation test: ANI distance vs gene-flow resistance",
               file.path(outdir, "permutation_null_distribution_ANI_vs_RESIST.png"))

cont_df <- df %>%
  select(Pair, mean_FST, mean_Dxy, mean_ANI_distance, mean_genecontent_distance, mean_gene_flow) %>%
  mutate(
    mean_gene_flow_rev = -mean_gene_flow,
    z_FST = safe_scale(mean_FST),
    z_Dxy = safe_scale(mean_Dxy),
    z_ANI_distance = safe_scale(mean_ANI_distance),
    z_genecontent = safe_scale(mean_genecontent_distance),
    z_gene_flow_rev = safe_scale(mean_gene_flow_rev)
  ) %>%
  rowwise() %>%
  mutate(Continuum_score_eq = mean(c(z_FST, z_Dxy, z_ANI_distance, z_genecontent, z_gene_flow_rev), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Continuum_score_eq_01 = safe_rescale01(Continuum_score_eq))

pca_input <- cont_df %>%
  select(Pair, mean_FST, mean_Dxy, mean_ANI_distance, mean_genecontent_distance, mean_gene_flow_rev)
X <- as.matrix(pca_input[, -1])
rownames(X) <- pca_input$Pair
pca <- prcomp(X, center = TRUE, scale. = TRUE)
scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
scores$Pair <- rownames(scores)
loadings <- as.data.frame(pca$rotation[, 1:2, drop = FALSE])
loadings$Variable <- rownames(loadings)
var_expl <- summary(pca)$importance[2, 1:2] * 100
write.csv(scores, file.path(outdir, "pairwise_continuum_PCA_scores.csv"), row.names = FALSE)
write.csv(loadings, file.path(outdir, "pairwise_continuum_PCA_loadings.csv"), row.names = FALSE)

scores2 <- scores %>% left_join(cont_df %>% select(Pair, Continuum_score_eq), by = "Pair")
pc1_cor <- suppressWarnings(cor(scores2$PC1, scores2$Continuum_score_eq, method = "spearman", use = "pairwise.complete.obs"))
if (!is.na(pc1_cor) && pc1_cor < 0) {
  scores2$PC1_continuum <- -scores2$PC1
  loadings$PC1 <- -loadings$PC1
} else {
  scores2$PC1_continuum <- scores2$PC1
}
scores2 <- scores2 %>%
  mutate(
    Continuum_score_pca = safe_scale(PC1_continuum),
    Continuum_score_pca_01 = safe_rescale01(Continuum_score_pca)
  )
cont_df <- cont_df %>%
  left_join(scores2 %>% select(Pair, PC1, PC2, Continuum_score_pca, Continuum_score_pca_01), by = "Pair") %>%
  mutate(
    Continuum_score_consensus = rowMeans(cbind(Continuum_score_eq, Continuum_score_pca), na.rm = TRUE),
    Continuum_score_consensus_01 = safe_rescale01(Continuum_score_consensus)
  ) %>%
  arrange(desc(Continuum_score_consensus)) %>%
  mutate(Rank = row_number())

n_pairs <- nrow(cont_df)
cont_df <- cont_df %>%
  mutate(
    Continuum_group = case_when(
      Rank <= ceiling(n_pairs / 3) ~ "High divergence",
      Rank <= ceiling(2 * n_pairs / 3) ~ "Intermediate divergence",
      TRUE ~ "Low divergence"
    ),
    Continuum_group = factor(Continuum_group,
                             levels = c("Low divergence", "Intermediate divergence", "High divergence"))
  )
write.csv(cont_df, file.path(outdir, "pairwise_continuum_with_scores.csv"), row.names = FALSE)
ranked_df <- cont_df %>%
  select(Rank, Pair, Continuum_group, mean_FST, mean_Dxy, mean_ANI_distance,
         mean_genecontent_distance, mean_gene_flow, Continuum_score_eq,
         Continuum_score_pca, Continuum_score_consensus, Continuum_score_eq_01,
         Continuum_score_pca_01, Continuum_score_consensus_01) %>%
  arrange(Rank)
write.csv(ranked_df, file.path(outdir, "pairwise_continuum_ranked.csv"), row.names = FALSE)

cont_df$Pair <- gsub("_", "\u2013", cont_df$Pair)
arrow_scale <- 1.8
load_plot <- loadings %>% mutate(xend = PC1 * arrow_scale, yend = PC2 * arrow_scale)
plot_scores <- cont_df %>% select(Pair, PC1, PC2, Continuum_score_consensus_01, Continuum_group)
p_pca <- ggplot() +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey70") +
  geom_point(data = plot_scores,
             aes(x = PC1, y = PC2, fill = Continuum_score_consensus_01),
             shape = 21, color = "black", size = 4.2, stroke = 0.45) +
  ggrepel::geom_text_repel(data = plot_scores, aes(x = PC1, y = PC2, label = Pair), max.overlaps = Inf) +
  geom_segment(data = load_plot, aes(x = 0, y = 0, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.18, "cm")), linewidth = 0.65, color = "grey30") +
  ggrepel::geom_text_repel(data = load_plot, aes(x = xend, y = yend, label = Variable),
                           segment.color = NA, max.overlaps = Inf) +
  scale_fill_gradient2(low = "#3B7DDD", mid = "white", high = "#C53A32",
                       midpoint = median(plot_scores$Continuum_score_consensus_01, na.rm = TRUE),
                       limits = c(0, 1), name = "Continuum score") +
  labs(x = paste0("PC1 (", round(var_expl[1], 1), "%)"),
       y = paste0("PC2 (", round(var_expl[2], 1), "%)"),
       title = "Pairwise differentiation continuum") +
  coord_equal() + theme_paper()
ggsave(file.path(outdir, "pairwise_continuum_PCA_publication.png"), p_pca, width = 8.4, height = 6.8, dpi = 400)
ggsave(file.path(outdir, "pairwise_continuum_PCA_publication.pdf"), p_pca, width = 8.4, height = 6.8)

bar_df <- cont_df %>% arrange(Continuum_score_consensus) %>% mutate(Pair = factor(Pair, levels = Pair))
group_colors <- c("High divergence" = "#B2182B", "Intermediate divergence" = "#878787", "Low divergence" = "#2166AC")
p_bar <- ggplot(bar_df, aes(x = Pair, y = Continuum_score_consensus_01, fill = Continuum_group)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.25) +
  coord_flip(expand = FALSE) + scale_fill_manual(values = group_colors) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(x = NULL, y = "Composite continuum score") + theme_paper(base_size = 9)
ggsave(file.path(outdir, "Fig_panel_A_barplot.png"), p_bar, width = 89, height = 80, units = "mm", dpi = 600)

heat_df <- cont_df %>%
  select(Pair, z_FST, z_Dxy, z_ANI_distance, z_genecontent, z_gene_flow_rev) %>%
  pivot_longer(cols = -Pair, names_to = "Metric", values_to = "Zscore")
heat_df$Metric <- factor(
  heat_df$Metric,
  levels = c("z_FST", "z_Dxy", "z_ANI_distance", "z_genecontent", "z_gene_flow_rev"),
  labels = c("FST", "Dxy", "ANI distance", "Gene-content\ndistance", "Reversed\ngene flow")
)
p_heat <- ggplot(heat_df, aes(x = Metric, y = Pair, fill = Zscore)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_distiller(type = "div", palette = "RdBu", direction = -1, name = "Z-score") +
  labs(x = NULL, y = NULL) + theme_paper(base_size = 9) +
  theme(panel.border = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(outdir, "Fig_panel_B_heatmap.png"), p_heat, width = 89, height = 80, units = "mm", dpi = 600)

cat("Continuum analysis complete. Results saved to: ", outdir, "\n")

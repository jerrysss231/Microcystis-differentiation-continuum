PROJECT_ROOT <- normalizePath(
  Sys.getenv("MICROCYSTIS_PROJECT_ROOT", unset = "."),
  mustWork = FALSE
)

data_path <- function(...) file.path(PROJECT_ROOT, "data", ...)
result_path <- function(...) file.path(PROJECT_ROOT, "results", ...)

req_pkgs <- c("gradientForest", "extendedForest", "terra", "dplyr", "readxl", "ggplot2")
missing_pkgs <- req_pkgs[!sapply(req_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "))
invisible(lapply(req_pkgs, library, character.only = TRUE))

metadata_path <- data_path("metadata", "Microcystis_final_metadata.xlsx")
gea_dir <- result_path("gea")
out_dir <- gea_dir
future_tif <- data_path("climate", "wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp245_2041-2060.tif")

predictor_vars <- c("BIO1", "BIO2", "BIO8", "BIO12", "BIO15", "BIO19")
gf_ntree <- 500
gf_corr_threshold <- 0.5

extract_accession_key <- function(x) {
  x <- as.character(x)
  m1_pos <- regexpr("GC[AF]_[0-9]+\\.[0-9]+", x)
  if (m1_pos[1] > 0) return(sub("^GC[AF]_", "", regmatches(x, m1_pos)))
  m2_pos <- regexpr("GC[AF]_[0-9]+_[0-9]+", x)
  if (m2_pos[1] > 0) {
    key <- sub("^GC[AF]_", "", regmatches(x, m2_pos))
    return(sub("_([0-9]+)$", ".\\1", key))
  }
  NA_character_
}

read_lfmm_matrix <- function(file_path) {
  x <- as.matrix(utils::read.table(file_path, header = FALSE, stringsAsFactors = FALSE))
  storage.mode(x) <- "numeric"
  x
}

need_files <- c(
  file.path(gea_dir, "Environment_matrix_used.csv"),
  file.path(gea_dir, "Matched_samples_used.csv"),
  file.path(gea_dir, "LFMM2_candidate_SNPs.csv"),
  file.path(gea_dir, "LFMM2_Full_SNP_environment_results.csv"),
  file.path(gea_dir, "snp_all_for_lea.lfmm_imputed.lfmm")
)
missing_files <- need_files[!file.exists(need_files)]
if (length(missing_files) > 0) stop("Missing input files:\n", paste(missing_files, collapse = "\n"))

env_df <- read.csv(file.path(gea_dir, "Environment_matrix_used.csv"), stringsAsFactors = FALSE, check.names = FALSE)
if (!("strain_key" %in% colnames(env_df))) stop("Environment_matrix_used.csv is missing strain_key.")
rownames(env_df) <- trimws(env_df$strain_key)
missing_env_cols <- setdiff(predictor_vars, colnames(env_df))
if (length(missing_env_cols) > 0) stop("Missing environmental predictors: ", paste(missing_env_cols, collapse = ", "))
env_gea <- env_df[, predictor_vars, drop = FALSE]

matched_df <- read.csv(file.path(gea_dir, "Matched_samples_used.csv"), stringsAsFactors = FALSE)
matched_samples <- trimws(matched_df$strain_key)
if (!all(matched_samples %in% rownames(env_gea))) stop("Matched sample IDs are absent from the environmental matrix.")
env_gea <- env_gea[match(matched_samples, rownames(env_gea)), , drop = FALSE]
rownames(env_gea) <- matched_samples

lfmm_sig_df <- read.csv(file.path(gea_dir, "LFMM2_candidate_SNPs.csv"), stringsAsFactors = FALSE)
if (!("SNP_ID" %in% colnames(lfmm_sig_df))) stop("LFMM2_candidate_SNPs.csv is missing SNP_ID.")
lfmm_sig_snp <- lfmm_sig_df$SNP_ID

Y_lfmm <- read_lfmm_matrix(file.path(gea_dir, "snp_all_for_lea.lfmm_imputed.lfmm"))
if (nrow(Y_lfmm) != length(matched_samples)) stop("SNP matrix row count does not match the sample list.")
rownames(Y_lfmm) <- matched_samples

lfmm_full <- read.csv(file.path(gea_dir, "LFMM2_Full_SNP_environment_results.csv"), stringsAsFactors = FALSE)
all_snp_ids <- unique(lfmm_full$SNP_ID)
if (length(all_snp_ids) != ncol(Y_lfmm)) stop("SNP matrix columns do not match LFMM2 result IDs.")
colnames(Y_lfmm) <- all_snp_ids

if (!all(lfmm_sig_snp %in% colnames(Y_lfmm))) stop("Some LFMM2 candidate SNPs are missing from the imputed matrix.")
Y_gf <- Y_lfmm[, lfmm_sig_snp, drop = FALSE]
if (!identical(rownames(Y_gf), rownames(env_gea))) stop("Gradient Forest response and predictor sample orders differ.")

gf_input <- cbind(env_gea[, predictor_vars, drop = FALSE], as.data.frame(Y_gf, check.names = FALSE))
write.csv(cbind(strain_key = rownames(gf_input), gf_input), file.path(out_dir, "GF_input_LFMM2_candidate_SNPs.csv"), row.names = FALSE)

gf_maxLevel <- max(1, floor(log2(0.368 * nrow(gf_input) / 2)))
gf_fit <- gradientForest(
  data = gf_input,
  predictor.vars = predictor_vars,
  response.vars = colnames(Y_gf),
  ntree = gf_ntree,
  maxLevel = gf_maxLevel,
  corr.threshold = gf_corr_threshold,
  trace = TRUE
)
saveRDS(gf_fit, file.path(out_dir, "GF_model_LFMM2_candidate_SNPs.rds"))

gf_imp_raw <- gf_fit$overall.imp
if (is.matrix(gf_imp_raw) || is.data.frame(gf_imp_raw)) {
  gf_imp <- data.frame(Variable = rownames(gf_imp_raw), Overall_Importance = as.numeric(gf_imp_raw[, 1]))
} else {
  gf_imp <- data.frame(Variable = names(gf_imp_raw), Overall_Importance = as.numeric(gf_imp_raw))
}
gf_imp <- gf_imp %>% arrange(desc(Overall_Importance))
write.csv(gf_imp, file.path(out_dir, "GF_variable_importance.csv"), row.names = FALSE)

p_imp <- ggplot(gf_imp, aes(x = reorder(Variable, Overall_Importance), y = Overall_Importance)) +
  geom_col() + coord_flip() + theme_bw(base_size = 12) +
  labs(x = "Environmental variable", y = "Overall importance")
ggsave(file.path(out_dir, "GF_variable_importance.png"), p_imp, width = 6, height = 5)

gf_current_pred <- as.data.frame(predict(gf_fit, newdata = env_gea[, predictor_vars, drop = FALSE]))
rownames(gf_current_pred) <- rownames(env_gea)
write.csv(cbind(strain_key = rownames(gf_current_pred), gf_current_pred), file.path(out_dir, "GF_transformed_current_environment.csv"), row.names = FALSE)

if (!file.exists(future_tif)) stop("Future climate raster not found: ", future_tif)
future_rast <- terra::rast(future_tif)
normalize_accession_key <- function(x) sub("^0+", "", trimws(as.character(x)))
metadata_full <- read_excel(metadata_path) %>%
  mutate(
    strain_key_raw = trimws(sapply(`Assembly Accession`, extract_accession_key)),
    strain_key_norm = normalize_accession_key(strain_key_raw),
    latitude = suppressWarnings(as.numeric(latitude)),
    longitude = suppressWarnings(as.numeric(longitude))
  ) %>%
  filter(!is.na(strain_key_raw), !is.na(latitude), !is.na(longitude))

idx <- match(normalize_accession_key(rownames(env_gea)), metadata_full$strain_key_norm)
if (any(is.na(idx))) stop("Some Gradient Forest samples could not be matched to metadata coordinates.")
metadata_coords <- metadata_full[idx, , drop = FALSE]
rownames(metadata_coords) <- rownames(env_gea)
pts <- terra::vect(metadata_coords[, c("longitude", "latitude")], geom = c("longitude", "latitude"), crs = "EPSG:4326")
future_vals <- as.data.frame(terra::extract(future_rast, pts))

future_env <- data.frame(
  strain_key = rownames(metadata_coords),
  BIO1 = future_vals[, 2], BIO2 = future_vals[, 3], BIO8 = future_vals[, 9],
  BIO12 = future_vals[, 13], BIO15 = future_vals[, 16], BIO19 = future_vals[, 20]
)
rownames(future_env) <- future_env$strain_key
future_env_use <- future_env[rownames(env_gea), predictor_vars, drop = FALSE]
write.csv(future_env, file.path(out_dir, "future_env_ssp245_2041_2060.csv"), row.names = FALSE)

gf_future_pred <- as.data.frame(predict(gf_fit, newdata = future_env_use))
rownames(gf_future_pred) <- rownames(future_env_use)
write.csv(cbind(strain_key = rownames(gf_future_pred), gf_future_pred), file.path(out_dir, "GF_transformed_future_environment.csv"), row.names = FALSE)

gf_offset <- sqrt(rowSums((as.matrix(gf_future_pred) - as.matrix(gf_current_pred))^2))
gf_offset_df <- data.frame(strain_key = rownames(env_gea), Genomic_Offset = gf_offset)
write.csv(gf_offset_df, file.path(out_dir, "GF_future_genomic_offset.csv"), row.names = FALSE)

p_offset <- ggplot(gf_offset_df, aes(x = Genomic_Offset)) +
  geom_histogram(bins = 30) + theme_bw(base_size = 12) +
  labs(x = "Genomic offset", y = "Count")
ggsave(file.path(out_dir, "GF_future_genomic_offset_histogram.png"), p_offset, width = 6, height = 5)

gf_both <- rbind(cbind(gf_current_pred, Scenario = "Current"), cbind(gf_future_pred, Scenario = "Future"))
gf_both_num <- gf_both[, setdiff(colnames(gf_both), "Scenario"), drop = FALSE]
gf_both_pca <- prcomp(gf_both_num, center = TRUE, scale. = TRUE)
gf_both_scores <- as.data.frame(gf_both_pca$x[, 1:2, drop = FALSE])
gf_both_scores$Scenario <- gf_both$Scenario
gf_both_scores$strain_key <- c(rownames(gf_current_pred), rownames(gf_future_pred))
write.csv(gf_both_scores, file.path(out_dir, "GF_current_future_PCA_scores.csv"), row.names = FALSE)

p_both <- ggplot(gf_both_scores, aes(x = PC1, y = PC2, color = Scenario)) +
  geom_point(size = 2, alpha = 0.8) + theme_bw(base_size = 12) +
  labs(x = "PC1", y = "PC2")
ggsave(file.path(out_dir, "GF_current_vs_future_PCA.png"), p_both, width = 6, height = 5)

cat("Gradient Forest analysis complete. Results saved to: ", out_dir, "\n")

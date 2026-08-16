PROJECT_ROOT <- normalizePath(
  Sys.getenv("MICROCYSTIS_PROJECT_ROOT", unset = "."),
  mustWork = FALSE
)

data_path <- function(...) file.path(PROJECT_ROOT, "data", ...)
result_path <- function(...) file.path(PROJECT_ROOT, "results", ...)

packages <- c("vegan", "geosphere", "readxl", "dplyr", "ggplot2", "gridExtra", "tidyr")
missing_pkgs <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "))
invisible(lapply(packages, library, character.only = TRUE))

extract_accession_key <- function(x) {
  x <- as.character(x)
  m1 <- regmatches(x, regexpr("GC[AF]_[0-9]+\\.[0-9]+", x))
  if (length(m1) > 0 && !is.na(m1) && nchar(m1) > 0) return(sub("^GC[AF]_", "", m1))
  m2 <- regmatches(x, regexpr("GC[AF]_[0-9]+_[0-9]+", x))
  if (length(m2) > 0 && !is.na(m2) && nchar(m2) > 0) {
    key <- sub("^GC[AF]_", "", m2)
    return(sub("_([0-9]+)$", ".\\1", key))
  }
  NA_character_
}

run_module_pca <- function(df, vars, prefix, n_pc = 2) {
  subdf <- df %>% select(all_of(vars)) %>% mutate(across(everything(), as.numeric))
  pca_obj <- prcomp(subdf, center = TRUE, scale. = TRUE)
  pcs <- as.data.frame(pca_obj$x[, 1:n_pc, drop = FALSE])
  colnames(pcs) <- paste0(prefix, "_PC", seq_len(n_pc))
  list(
    pca = pca_obj,
    pcs = pcs,
    var_explained = summary(pca_obj)$importance[2, ],
    cum_explained = summary(pca_obj)$importance[3, ]
  )
}

get_positive_pcoa_axes <- function(dist_obj) {
  pcoa <- cmdscale(dist_obj, k = attr(dist_obj, "Size") - 1, eig = TRUE, add = TRUE)
  eig_used <- pcoa$eig[seq_len(ncol(pcoa$points))]
  as.data.frame(pcoa$points[, eig_used > 0, drop = FALSE])
}

check_same_size <- function(...) {
  sizes <- sapply(list(...), function(x) attr(x, "Size"))
  if (length(unique(sizes)) != 1) stop("Distance matrices do not have the same size.")
}

check_no_na_dist <- function(x, name) {
  if (any(is.na(as.vector(x)))) stop(paste0("NA detected in distance object: ", name))
}

metadata_file <- data_path("metadata", "Microcystis_final_metadata.xlsx")
ani_file <- data_path("ani", "fastani_all_vs_all.txt.matrix")
phylip_file <- data_path("snp", "gubbins.filtered_polymorphic_sites.phylip")
panaroo_file <- data_path("pangenome", "gene_presence_absence.Rtab")
outdir <- result_path("landscape_genomics")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

bio_vars <- paste0("BIO", 1:19)
srad_vars <- paste0("wc2.1_2.5m_srad_", sprintf("%02d", 1:12))
uvb_vars <- c(
  "56459_UVB1_Annual_Mean_UV-B",
  "56460_UVB2_UV-B_Seasonality",
  "56461_UVB3_Mean_UV-B_of_Highest_Month",
  "56462_UVB4_Mean_UV-B_of_Lowest_Month",
  "56463_UVB5_Sum_of_UV-B_Radiation_of_Highest_Quarter",
  "56464_UVB6_Sum_of_UV-B_Radiation_of_Lowest_Quarter"
)
all_env_vars <- c(bio_vars, srad_vars, uvb_vars)

metadata <- read_excel(metadata_file) %>%
  mutate(
    strain_id = `Assembly Accession`,
    strain_key = sapply(strain_id, extract_accession_key)
  ) %>%
  filter(!is.na(latitude), !is.na(longitude)) %>%
  filter(if_all(all_of(all_env_vars), ~ !is.na(.))) %>%
  filter(!is.na(strain_key))

ani_lines <- readLines(ani_file)
n_ani <- as.integer(trimws(ani_lines[1]))
ani_names_raw <- sapply(ani_lines[2:(n_ani + 1)], function(x) basename(trimws(strsplit(x, "\t")[[1]][1])))
ani_keys <- sapply(ani_names_raw, extract_accession_key)
ani_mat <- matrix(NA, nrow = n_ani, ncol = n_ani, dimnames = list(ani_names_raw, ani_names_raw))
for (i in seq_len(n_ani)) {
  parts <- strsplit(trimws(ani_lines[i + 1]), "\t")[[1]]
  vals <- suppressWarnings(as.numeric(parts[-1]))
  if (length(vals) > 0) ani_mat[i, seq_along(vals)] <- vals
}
for (i in seq_len(n_ani)) {
  for (j in seq_len(n_ani)) {
    if (is.na(ani_mat[i, j]) && !is.na(ani_mat[j, i])) ani_mat[i, j] <- ani_mat[j, i]
  }
}
diag(ani_mat) <- 100
ani_lookup <- data.frame(ani_name = ani_names_raw, ani_key = ani_keys, stringsAsFactors = FALSE) %>%
  filter(!is.na(ani_key)) %>% distinct(ani_key, .keep_all = TRUE)

phylip_lines <- readLines(phylip_file)
header <- strsplit(trimws(phylip_lines[1]), "\\s+")[[1]]
n_snp_samples <- as.integer(header[1])
snp_names_raw <- character(n_snp_samples)
snp_keys <- character(n_snp_samples)
snp_seqs <- character(n_snp_samples)
for (i in seq_len(n_snp_samples)) {
  parts <- strsplit(trimws(phylip_lines[i + 1]), "\\s+")[[1]]
  snp_names_raw[i] <- parts[1]
  snp_keys[i] <- extract_accession_key(parts[1])
  snp_seqs[i] <- paste(parts[-1], collapse = "")
}
snp_lookup <- data.frame(snp_name = snp_names_raw, snp_key = snp_keys, stringsAsFactors = FALSE) %>%
  filter(!is.na(snp_key)) %>% distinct(snp_key, .keep_all = TRUE)

panaroo_raw <- read.delim(panaroo_file, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
gene_ids <- panaroo_raw[[1]]
panaroo_mat <- panaroo_raw[, -1, drop = FALSE]
panaroo_colnames_raw <- colnames(panaroo_mat)
panaroo_keys <- sapply(panaroo_colnames_raw, extract_accession_key)
panaroo_lookup <- data.frame(panaroo_name = panaroo_colnames_raw, panaroo_key = panaroo_keys, stringsAsFactors = FALSE) %>%
  filter(!is.na(panaroo_key)) %>% distinct(panaroo_key, .keep_all = TRUE)
panaroo_mat2 <- as.data.frame(panaroo_mat, check.names = FALSE)
for (j in seq_len(ncol(panaroo_mat2))) panaroo_mat2[[j]] <- as.numeric(panaroo_mat2[[j]])
rownames(panaroo_mat2) <- gene_ids

common_keys <- Reduce(intersect, list(metadata$strain_key, ani_lookup$ani_key, snp_lookup$snp_key, panaroo_lookup$panaroo_key))
if (length(common_keys) < 10) stop("Too few common genomes after alignment.")

meta_sub <- metadata %>%
  filter(strain_key %in% common_keys) %>%
  distinct(strain_key, .keep_all = TRUE) %>%
  arrange(match(strain_key, common_keys))
ani_sub_names <- ani_lookup %>% filter(ani_key %in% common_keys) %>% arrange(match(ani_key, common_keys)) %>% pull(ani_name)
snp_sub_names <- snp_lookup %>% filter(snp_key %in% common_keys) %>% arrange(match(snp_key, common_keys)) %>% pull(snp_name)
snp_sub_keys <- snp_lookup %>% filter(snp_key %in% common_keys) %>% arrange(match(snp_key, common_keys)) %>% pull(snp_key)
panaroo_sub_names <- panaroo_lookup %>% filter(panaroo_key %in% common_keys) %>% arrange(match(panaroo_key, common_keys)) %>% pull(panaroo_name)

ani_sub <- ani_mat[ani_sub_names, ani_sub_names]
rownames(ani_sub) <- colnames(ani_sub) <- common_keys

snp_sub_seqs <- snp_seqs[match(snp_sub_names, snp_names_raw)]
snp_pdist_mat <- matrix(0, nrow = length(snp_sub_names), ncol = length(snp_sub_names))
rownames(snp_pdist_mat) <- colnames(snp_pdist_mat) <- snp_sub_keys
for (i in 1:(length(snp_sub_names) - 1)) {
  si <- strsplit(snp_sub_seqs[i], "")[[1]]
  for (j in (i + 1):length(snp_sub_names)) {
    sj <- strsplit(snp_sub_seqs[j], "")[[1]]
    valid <- si != "-" & sj != "-" & si != "N" & sj != "N"
    d <- if (sum(valid) == 0) NA else sum(si[valid] != sj[valid]) / sum(valid)
    snp_pdist_mat[i, j] <- d
    snp_pdist_mat[j, i] <- d
  }
}
diag(snp_pdist_mat) <- 0
if (any(is.na(snp_pdist_mat))) stop("NA detected in SNP distance matrix.")

panaroo_sub <- panaroo_mat2[, panaroo_sub_names, drop = FALSE]
colnames(panaroo_sub) <- common_keys
gene_content_mat <- t(as.matrix(panaroo_sub))
gene_content_mat[gene_content_mat > 0] <- 1
gene_content_mat[gene_content_mat <= 0] <- 0

ani_dist <- as.dist(1 - ani_sub / 100)
snp_dist <- as.dist(snp_pdist_mat)
gene_dist <- vegdist(gene_content_mat, method = "jaccard", binary = TRUE)
coords <- as.matrix(meta_sub[, c("longitude", "latitude")])
geo_dist <- as.dist(distm(coords, fun = distGeo) / 1000)

bio_pca_res <- run_module_pca(meta_sub, bio_vars, prefix = "BIO", n_pc = 2)
srad_pca_res <- run_module_pca(meta_sub, srad_vars, prefix = "SRAD", n_pc = 2)
uvb_pca_res <- run_module_pca(meta_sub, uvb_vars, prefix = "UVB", n_pc = 2)
env_pcs <- bind_cols(bio_pca_res$pcs, srad_pca_res$pcs, uvb_pca_res$pcs)
env_dist <- dist(env_pcs, method = "euclidean")

write.csv(round(bio_pca_res$pca$rotation[, 1:2], 4), file.path(outdir, "BIO_PCA_loadings.csv"))
write.csv(round(srad_pca_res$pca$rotation[, 1:2], 4), file.path(outdir, "SRAD_PCA_loadings.csv"))
write.csv(round(uvb_pca_res$pca$rotation[, 1:2], 4), file.path(outdir, "UVB_PCA_loadings.csv"))

check_same_size(ani_dist, snp_dist, gene_dist, geo_dist, env_dist)
check_no_na_dist(ani_dist, "ani_dist")
check_no_na_dist(snp_dist, "snp_dist")
check_no_na_dist(gene_dist, "gene_dist")
check_no_na_dist(geo_dist, "geo_dist")
check_no_na_dist(env_dist, "env_dist")

mantel_rows <- list(
  c("ANI", "Geography", mantel(ani_dist, geo_dist, method = "pearson", permutations = 9999)$statistic, mantel(ani_dist, geo_dist, method = "pearson", permutations = 9999)$signif),
  c("ANI", "Environment", mantel(ani_dist, env_dist, method = "pearson", permutations = 9999)$statistic, mantel(ani_dist, env_dist, method = "pearson", permutations = 9999)$signif),
  c("SNP", "Geography", mantel(snp_dist, geo_dist, method = "pearson", permutations = 9999)$statistic, mantel(snp_dist, geo_dist, method = "pearson", permutations = 9999)$signif),
  c("SNP", "Environment", mantel(snp_dist, env_dist, method = "pearson", permutations = 9999)$statistic, mantel(snp_dist, env_dist, method = "pearson", permutations = 9999)$signif),
  c("GeneContent", "Geography", mantel(gene_dist, geo_dist, method = "pearson", permutations = 9999)$statistic, mantel(gene_dist, geo_dist, method = "pearson", permutations = 9999)$signif),
  c("GeneContent", "Environment", mantel(gene_dist, env_dist, method = "pearson", permutations = 9999)$statistic, mantel(gene_dist, env_dist, method = "pearson", permutations = 9999)$signif)
)
mantel_results <- as.data.frame(do.call(rbind, mantel_rows), stringsAsFactors = FALSE)
names(mantel_results) <- c("Genetic_layer", "Comparison", "r", "p_value")
mantel_results$r <- as.numeric(mantel_results$r)
mantel_results$p_value <- as.numeric(mantel_results$p_value)
write.csv(mantel_results, file.path(outdir, "mantel_results_threeway.csv"), row.names = FALSE)

pm_ani_env_geo <- mantel.partial(ani_dist, env_dist, geo_dist, method = "pearson", permutations = 9999)
pm_snp_env_geo <- mantel.partial(snp_dist, env_dist, geo_dist, method = "pearson", permutations = 9999)
pm_gene_env_geo <- mantel.partial(gene_dist, env_dist, geo_dist, method = "pearson", permutations = 9999)
pm_ani_geo_env <- mantel.partial(ani_dist, geo_dist, env_dist, method = "pearson", permutations = 9999)
pm_snp_geo_env <- mantel.partial(snp_dist, geo_dist, env_dist, method = "pearson", permutations = 9999)
pm_gene_geo_env <- mantel.partial(gene_dist, geo_dist, env_dist, method = "pearson", permutations = 9999)

partial_mantel_results <- data.frame(
  Genetic_layer = c("ANI", "SNP", "GeneContent", "ANI", "SNP", "GeneContent"),
  Test = c(rep("Environment | Geography", 3), rep("Geography | Environment", 3)),
  r = c(pm_ani_env_geo$statistic, pm_snp_env_geo$statistic, pm_gene_env_geo$statistic,
        pm_ani_geo_env$statistic, pm_snp_geo_env$statistic, pm_gene_geo_env$statistic),
  p_value = c(pm_ani_env_geo$signif, pm_snp_env_geo$signif, pm_gene_env_geo$signif,
              pm_ani_geo_env$signif, pm_snp_geo_env$signif, pm_gene_geo_env$signif)
)
write.csv(partial_mantel_results, file.path(outdir, "partial_mantel_results_threeway.csv"), row.names = FALSE)

env_df <- bind_cols(meta_sub %>% select(strain_key, longitude, latitude), env_pcs) %>%
  mutate(
    Lon = as.numeric(scale(longitude)),
    Lat = as.numeric(scale(latitude)),
    Lon2 = Lon^2,
    Lat2 = Lat^2,
    LonLat = Lon * Lat
  )
module_vars <- env_df %>% select(BIO_PC1, BIO_PC2, SRAD_PC1, SRAD_PC2, UVB_PC1, UVB_PC2)
space_vars <- env_df %>% select(Lon, Lat, Lon2, Lat2, LonLat)

run_distance_analysis <- function(dist_obj, dist_name, env_df, module_vars, space_vars, outdir) {
  formula_adonis <- as.formula(paste0("dist_obj ~ ", paste(colnames(module_vars), collapse = " + "), " + ", paste(colnames(space_vars), collapse = " + ")))
  adonis_res <- adonis2(formula_adonis, data = env_df, permutations = 9999, by = "margin")
  formula_cap <- as.formula(paste0("dist_obj ~ ", paste(colnames(module_vars), collapse = " + "), " + Condition(", paste(colnames(space_vars), collapse = " + "), ")"))
  cap_res <- capscale(formula_cap, data = env_df, add = TRUE)
  cap_overall <- anova.cca(cap_res, permutations = 9999)
  cap_axis <- anova.cca(cap_res, by = "axis", permutations = 9999)
  axes <- get_positive_pcoa_axes(dist_obj)
  varpart_res <- varpart(axes, module_vars, space_vars)
  capture.output(adonis_res, file = file.path(outdir, paste0("adonis2_", dist_name, ".txt")))
  capture.output(cap_overall, file = file.path(outdir, paste0("dbrda_", dist_name, "_overall.txt")))
  capture.output(cap_axis, file = file.path(outdir, paste0("dbrda_", dist_name, "_axis.txt")))
  capture.output(varpart_res, file = file.path(outdir, paste0("varpart_", dist_name, ".txt")))
  invisible(list(adonis = adonis_res, cap = cap_res, varpart = varpart_res))
}

run_distance_analysis(ani_dist, "ani_threeway", env_df, module_vars, space_vars, outdir)
run_distance_analysis(snp_dist, "snp_threeway", env_df, module_vars, space_vars, outdir)
run_distance_analysis(gene_dist, "genecontent_threeway", env_df, module_vars, space_vars, outdir)

bio_load <- round(bio_pca_res$pca$rotation[, 1:2], 6)
top_pc1 <- rownames(bio_load)[order(abs(bio_load[, "PC1"]), decreasing = TRUE)][1:5]
top_pc2 <- rownames(bio_load)[order(abs(bio_load[, "PC2"]), decreasing = TRUE)][1:5]
candidate_bio <- unique(c(top_pc1, top_pc2))
bio_env_df <- cbind(env_df, meta_sub[, candidate_bio, drop = FALSE])
for (v in candidate_bio) bio_env_df[[v]] <- as.numeric(bio_env_df[[v]])

run_single_bio_adonis <- function(dist_obj, dist_name, vars, data_df) {
  do.call(rbind, lapply(vars, function(v) {
    fml <- as.formula(paste0("dist_obj ~ `", v, "` + Lon + Lat + Lon2 + Lat2 + LonLat"))
    fit <- adonis2(fml, data = data_df, permutations = 9999, by = "margin")
    data.frame(
      Genetic_layer = dist_name,
      BIO_variable = v,
      R2 = fit[v, "R2"],
      F = fit[v, "F"],
      p_value = fit[v, "Pr(>F)"]
    )
  }))
}

bio_single_results <- rbind(
  run_single_bio_adonis(ani_dist, "ANI", candidate_bio, bio_env_df),
  run_single_bio_adonis(snp_dist, "SNP", candidate_bio, bio_env_df),
  run_single_bio_adonis(gene_dist, "GeneContent", candidate_bio, bio_env_df)
)
write.csv(bio_single_results, file.path(outdir, "BIO_candidate_single_variable_adonis.csv"), row.names = FALSE)

cat("Landscape-genomic analysis complete. Results saved to: ", outdir, "\n")

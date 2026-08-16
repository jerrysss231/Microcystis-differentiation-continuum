PROJECT_ROOT <- normalizePath(
  Sys.getenv("MICROCYSTIS_PROJECT_ROOT", unset = "."),
  mustWork = FALSE
)

data_path <- function(...) file.path(PROJECT_ROOT, "data", ...)
result_path <- function(...) file.path(PROJECT_ROOT, "results", ...)

packages <- c("ape", "phytools", "phangorn", "vegan", "ggplot2", "TreeDist")
missing_pkgs <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "))
invisible(lapply(packages, library, character.only = TRUE))

path_gtdb <- data_path("trees", "gtdb_rooted_tree.treefile")
path_snp <- data_path("trees", "gubbins.final_tree_pruned.tre")
path_core <- data_path("trees", "SpeciesTree_final.treefile")
out_dir <- result_path("tree_concordance")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

clean_tree <- function(tree) {
  tree$tip.label <- sapply(tree$tip.label, extract_accession_key)
  tree <- drop.tip(tree, tree$tip.label[is.na(tree$tip.label)])
  keep.tip(tree, unique(tree$tip.label))
}

tree_gtdb <- clean_tree(read.tree(path_gtdb))
tree_snp <- clean_tree(read.tree(path_snp))
tree_core <- clean_tree(read.tree(path_core))

common_tips <- Reduce(intersect, list(tree_gtdb$tip.label, tree_snp$tip.label, tree_core$tip.label))
if (length(common_tips) < 5) stop("Too few shared tips after accession matching.")
common_sorted <- sort(common_tips)

tree_gtdb_sub <- keep.tip(tree_gtdb, common_sorted)
tree_snp_sub <- keep.tip(tree_snp, common_sorted)
tree_core_sub <- keep.tip(tree_core, common_sorted)

u_core <- unroot(tree_core_sub)
u_snp <- unroot(tree_snp_sub)
u_gtdb <- unroot(tree_gtdb_sub)
n_tip <- length(common_sorted)
max_rf <- 2 * (n_tip - 3)

calc_rf_stats <- function(t1, t2, comparison_name) {
  rf_raw <- RF.dist(t1, t2)
  rf_norm <- rf_raw / max_rf
  data.frame(
    Comparison = comparison_name,
    N_tips = length(t1$tip.label),
    RF_distance = rf_raw,
    Max_RF_distance = max_rf,
    RF_normalized = rf_norm,
    RF_percent = rf_norm * 100,
    RF_similarity_percent = (1 - rf_norm) * 100
  )
}

rf_df <- rbind(
  calc_rf_stats(u_core, u_snp, "Core_vs_SNP"),
  calc_rf_stats(u_core, u_gtdb, "Core_vs_GTDB"),
  calc_rf_stats(u_snp, u_gtdb, "SNP_vs_GTDB")
)
write.csv(rf_df, file.path(out_dir, "Tree_RF_distance_summary.csv"), row.names = FALSE)

calc_treedist_stats <- function(t1, t2, comparison_name) {
  mci_raw <- MutualClusteringInfo(t1, t2, normalize = FALSE)
  mci_norm <- MutualClusteringInfo(t1, t2, normalize = TRUE)
  cid_raw <- ClusteringInfoDistance(t1, t2, normalize = FALSE)
  cid_norm <- ClusteringInfoDistance(t1, t2, normalize = TRUE)
  spi_raw <- SharedPhylogeneticInfo(t1, t2, normalize = FALSE)
  spi_norm <- SharedPhylogeneticInfo(t1, t2, normalize = TRUE)
  data.frame(
    Comparison = comparison_name,
    N_tips = length(t1$tip.label),
    MCI_raw = as.numeric(mci_raw),
    MCI_normalized = as.numeric(mci_norm),
    MCI_percent = as.numeric(mci_norm) * 100,
    CID_raw = as.numeric(cid_raw),
    CID_normalized = as.numeric(cid_norm),
    CID_percent = as.numeric(cid_norm) * 100,
    SPI_raw = as.numeric(spi_raw),
    SPI_normalized = as.numeric(spi_norm),
    SPI_percent = as.numeric(spi_norm) * 100
  )
}

treedist_df <- rbind(
  calc_treedist_stats(u_core, u_snp, "Core_vs_SNP"),
  calc_treedist_stats(u_core, u_gtdb, "Core_vs_GTDB"),
  calc_treedist_stats(u_snp, u_gtdb, "SNP_vs_GTDB")
)
write.csv(treedist_df, file.path(out_dir, "Tree_TreeDist_similarity_summary.csv"), row.names = FALSE)

dist_core <- cophenetic.phylo(tree_core_sub)[common_sorted, common_sorted]
dist_snp <- cophenetic.phylo(tree_snp_sub)[common_sorted, common_sorted]
dist_gtdb <- cophenetic.phylo(tree_gtdb_sub)[common_sorted, common_sorted]

mantel_core_snp <- mantel(as.dist(dist_core), as.dist(dist_snp), method = "pearson", permutations = 9999)
mantel_core_gtdb <- mantel(as.dist(dist_core), as.dist(dist_gtdb), method = "pearson", permutations = 9999)
mantel_snp_gtdb <- mantel(as.dist(dist_snp), as.dist(dist_gtdb), method = "pearson", permutations = 9999)

mantel_df <- data.frame(
  Comparison = c("Core_vs_SNP", "Core_vs_GTDB", "SNP_vs_GTDB"),
  Mantel_r = c(mantel_core_snp$statistic, mantel_core_gtdb$statistic, mantel_snp_gtdb$statistic),
  P_value = c(mantel_core_snp$signif, mantel_core_gtdb$signif, mantel_snp_gtdb$signif)
)
write.csv(mantel_df, file.path(out_dir, "Tree_Mantel_test_summary.csv"), row.names = FALSE)

plot_tanglegram <- function(t1, t2, name1, name2, filename) {
  co_trees <- cophylo(midpoint.root(t1), midpoint.root(t2), rotate = TRUE)
  png(file.path(out_dir, filename), width = 2400, height = 1800, res = 300)
  par(mar = c(1, 1, 3, 1))
  plot(co_trees, link.type = "curved", link.lwd = 1.2,
       link.col = make.transparent("steelblue", 0.5), fsize = 0.6)
  title(main = paste(name1, "vs", name2), cex.main = 1.5)
  dev.off()
}

plot_tanglegram(tree_core_sub, tree_snp_sub, "Core gene ML tree", "Whole-genome SNP tree", "Tanglegram_Core_vs_SNP.png")
plot_tanglegram(tree_core_sub, tree_gtdb_sub, "Core gene ML tree", "GTDB-Tk reference tree", "Tanglegram_Core_vs_GTDB.png")
plot_tanglegram(tree_snp_sub, tree_gtdb_sub, "Whole-genome SNP tree", "GTDB-Tk reference tree", "Tanglegram_SNP_vs_GTDB.png")

plot_tree_distance_scatter <- function(d1, d2, name1, name2, mantel_obj, filename) {
  idx <- upper.tri(d1)
  df <- data.frame(x = d1[idx], y = d2[idx])
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(alpha = 0.15, size = 0.8) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
    theme_bw(base_size = 12) +
    labs(
      title = paste(name1, "vs", name2),
      x = paste(name1, "cophenetic distance"),
      y = paste(name2, "cophenetic distance"),
      subtitle = paste0("Mantel r = ", round(as.numeric(mantel_obj$statistic), 3),
                        ", p = ", signif(mantel_obj$signif, 3))
    )
  ggsave(file.path(out_dir, filename), p, width = 6.5, height = 5)
}

plot_tree_distance_scatter(dist_core, dist_snp, "Core", "SNP", mantel_core_snp, "TreeDistanceScatter_Core_vs_SNP.png")
plot_tree_distance_scatter(dist_core, dist_gtdb, "Core", "GTDB", mantel_core_gtdb, "TreeDistanceScatter_Core_vs_GTDB.png")
plot_tree_distance_scatter(dist_snp, dist_gtdb, "SNP", "GTDB", mantel_snp_gtdb, "TreeDistanceScatter_SNP_vs_GTDB.png")

summary_df <- merge(rf_df, treedist_df, by = c("Comparison", "N_tips"), all = TRUE)
summary_df <- merge(summary_df, mantel_df, by = "Comparison", all = TRUE)
write.csv(summary_df, file.path(out_dir, "Tree_Comparison_All_Summary.csv"), row.names = FALSE)

cat("Tree concordance analysis complete. Shared tips: ", n_tip, ". Results saved to: ", out_dir, "\n", sep = "")

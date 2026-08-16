PROJECT_ROOT <- normalizePath(
  Sys.getenv("MICROCYSTIS_PROJECT_ROOT", unset = "."),
  mustWork = FALSE
)

data_path <- function(...) file.path(PROJECT_ROOT, "data", ...)
result_path <- function(...) file.path(PROJECT_ROOT, "results", ...)

# ============================================================
# Test associations between lineage assignment and geography
# ============================================================
req_pkgs <- c("readxl", "dplyr", "ggplot2", "maps")
missing_pkgs <- req_pkgs[!sapply(req_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "))
library(readxl)
library(dplyr)
library(ggplot2)
library(maps)

meta_file <- data_path("metadata", "Microcystis_final_metadata_with_clade.xlsx")
out_dir <- result_path("lineage_geography")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Reading metadata...\n")
df <- read_excel(meta_file)

clade_candidates <- c("Clade", "fastBAPS_Optimised", "fastBAPS_Classic")
clade_col <- clade_candidates[clade_candidates %in% colnames(df)][1]
if (length(clade_col) == 0) stop("metadata must contain a lineage-assignment column.")
need_cols <- c("latitude", "longitude", "Continent", "Country_Code", "Assembly Accession")
miss_cols <- setdiff(need_cols, colnames(df))
if (length(miss_cols) > 0) {
  stop("Metadata is missing required columns: ", paste(miss_cols, collapse = ", "))
}

df <- df %>%
  mutate(
    latitude = as.numeric(latitude),
    longitude = as.numeric(longitude),
    abs_latitude = abs(latitude),
    Clade = as.factor(.data[[clade_col]]),
    Continent = as.factor(Continent),
    Country_Code = as.factor(Country_Code)
  ) %>%
  filter(!is.na(Clade), !is.na(latitude), !is.na(longitude))

cat("Samples retained:", nrow(df), "\n")
cat("Number of lineages:", length(unique(df$Clade)), "\n")

write.csv(
  df,
  file.path(out_dir, "metadata_used_for_clade_geography.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("Plotting global sampling map...\n")
world_df <- map_data("world")

p_map <- ggplot() +
  geom_polygon(
    data = world_df,
    aes(x = long, y = lat, group = group),
    fill = "grey90",
    color = "white",
    linewidth = 0.2
  ) +
  geom_point(
    data = df,
    aes(x = longitude, y = latitude, color = Clade),
    size = 2.2,
    alpha = 0.85
  ) +
  coord_quickmap() +
  theme_bw(base_size = 12) +
  labs(
    title = "Geographic distribution of fastBAPS clades",
    x = NULL,
    y = NULL,
    color = "Clade"
  )

ggsave(
  file.path(out_dir, "fastBAPS_clade_world_map.png"),
  p_map,
  width = 11,
  height = 6,
  dpi = 300
)

cat("Testing lineage-by-continent association...\n")
tab_cont <- table(df$Clade, df$Continent)
write.csv(
  as.data.frame.matrix(tab_cont),
  file.path(out_dir, "Clade_by_Continent_table.csv"),
  fileEncoding = "UTF-8"
)
print(tab_cont)

chisq_cont <- suppressWarnings(chisq.test(tab_cont))
expected_min <- min(chisq_cont$expected)
cat("Minimum expected count in continent table:", expected_min, "\n")

if (any(chisq_cont$expected < 5)) {
  cat("Using Fisher exact test with simulated p-value because expected counts are small...\n")
  fisher_cont <- fisher.test(tab_cont, simulate.p.value = TRUE, B = 9999)

  sink(file.path(out_dir, "Clade_Continent_association_test.txt"))
  cat("Fisher's Exact Test (simulated p-value)\n\n")
  print(fisher_cont)
  sink()
} else {
  sink(file.path(out_dir, "Clade_Continent_association_test.txt"))
  cat("Chi-squared test\n\n")
  print(chisq_cont)
  sink()
}

cat("Testing geographic-coordinate differences among lineages...\n")
kw_lat <- kruskal.test(latitude ~ Clade, data = df)
kw_abslat <- kruskal.test(abs_latitude ~ Clade, data = df)
kw_lon <- kruskal.test(longitude ~ Clade, data = df)

sink(file.path(out_dir, "Clade_geographic_variable_tests.txt"))
cat("Kruskal-Wallis test: latitude ~ clade\n")
print(kw_lat)
cat("\nKruskal-Wallis test: absolute latitude ~ clade\n")
print(kw_abslat)
cat("\nKruskal-Wallis test: longitude ~ clade\n")
print(kw_lon)
sink()

cat("Plotting geographic summaries...\n")
p_lat <- ggplot(df, aes(x = Clade, y = latitude, fill = Clade)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
  theme_bw(base_size = 12) +
  labs(
    title = paste0("Latitude by clade (Kruskal-Wallis p = ", signif(kw_lat$p.value, 3), ")"),
    x = "Clade",
    y = "Latitude"
  ) +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "Latitude_by_clade.png"), p_lat, width = 7, height = 5, dpi = 300)

p_abslat <- ggplot(df, aes(x = Clade, y = abs_latitude, fill = Clade)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
  theme_bw(base_size = 12) +
  labs(
    title = paste0("Absolute latitude by clade (Kruskal-Wallis p = ", signif(kw_abslat$p.value, 3), ")"),
    x = "Clade",
    y = "Absolute latitude"
  ) +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "Absolute_latitude_by_clade.png"), p_abslat, width = 7, height = 5, dpi = 300)

p_lon <- ggplot(df, aes(x = Clade, y = longitude, fill = Clade)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
  theme_bw(base_size = 12) +
  labs(
    title = paste0("Longitude by clade (Kruskal-Wallis p = ", signif(kw_lon$p.value, 3), ")"),
    x = "Clade",
    y = "Longitude"
  ) +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "Longitude_by_clade.png"), p_lon, width = 7, height = 5, dpi = 300)

summary_geo <- df %>%
  group_by(Clade) %>%
  summarise(
    N = n(),
    mean_lat = mean(latitude, na.rm = TRUE),
    median_lat = median(latitude, na.rm = TRUE),
    mean_abs_lat = mean(abs_latitude, na.rm = TRUE),
    median_abs_lat = median(abs_latitude, na.rm = TRUE),
    mean_lon = mean(longitude, na.rm = TRUE),
    median_lon = median(longitude, na.rm = TRUE),
    n_continent = n_distinct(Continent),
    n_country = n_distinct(Country_Code),
    .groups = "drop"
  )

write.csv(
  summary_geo,
  file.path(out_dir, "Clade_geographic_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("\nAnalysis complete. Main output files:\n")
cat("1. fastBAPS_clade_world_map.png\n")
cat("2. Clade_by_Continent_table.csv\n")
cat("3. Clade_Continent_association_test.txt\n")
cat("4. Clade_geographic_variable_tests.txt\n")
cat("5. Latitude_by_clade.png / Absolute_latitude_by_clade.png / Longitude_by_clade.png\n")
cat("6. Clade_geographic_summary.csv\n")

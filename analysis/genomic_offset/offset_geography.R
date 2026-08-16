PROJECT_ROOT <- normalizePath(
  Sys.getenv("MICROCYSTIS_PROJECT_ROOT", unset = "."),
  mustWork = FALSE
)

data_path <- function(...) file.path(PROJECT_ROOT, "data", ...)
result_path <- function(...) file.path(PROJECT_ROOT, "results", ...)

pkgs <- c("readxl", "dplyr", "ggplot2", "stringr", "readr")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "))
invisible(lapply(pkgs, library, character.only = TRUE))

has_maps <- require("maps", character.only = TRUE, quietly = TRUE)
if (!has_maps) message("Package 'maps' not installed. World map plot will be skipped.")

metadata_file <- data_path("metadata", "Microcystis_final_metadata.xlsx")
offset_file <- result_path("gea", "GF_future_genomic_offset.csv")
outdir <- result_path("gea", "offset_geography")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

extract_accession_key <- function(x) {
  x <- as.character(x)
  m1 <- str_extract(x, "GC[AF]_[0-9]+\\.[0-9]+")
  if (!is.na(m1)) {
    key <- sub("^GC[AF]_", "", m1)
    parts <- strsplit(key, "\\.")[[1]]
    return(paste0(as.character(as.integer(parts[1])), ".", parts[2]))
  }
  m2 <- str_extract(x, "GC[AF]_[0-9]+_[0-9]+")
  if (!is.na(m2)) {
    key <- sub("^GC[AF]_", "", m2)
    key <- sub("_([0-9]+)$", ".\\1", key)
    parts <- strsplit(key, "\\.")[[1]]
    return(paste0(as.character(as.integer(parts[1])), ".", parts[2]))
  }
  NA_character_
}

meta <- read_excel(metadata_file)
if (!"strain_key" %in% colnames(meta)) {
  if (!"Assembly Accession" %in% colnames(meta)) stop("Metadata must contain strain_key or Assembly Accession.")
  meta <- meta %>% mutate(strain_key = sapply(`Assembly Accession`, extract_accession_key))
}
meta <- meta %>% mutate(latitude = as.numeric(latitude), longitude = as.numeric(longitude))

offset_df <- read_csv(offset_file, show_col_types = FALSE)
if (!all(c("strain_key", "Genomic_Offset") %in% colnames(offset_df))) {
  stop("Offset file must contain strain_key and Genomic_Offset.")
}
offset_df <- offset_df %>% mutate(strain_key = as.character(strain_key), Genomic_Offset = as.numeric(Genomic_Offset))

merged <- meta %>%
  inner_join(offset_df, by = "strain_key") %>%
  filter(!is.na(latitude), !is.na(longitude), !is.na(Genomic_Offset)) %>%
  mutate(abs_latitude = abs(latitude))

write.csv(merged, file.path(outdir, "merged_metadata_offset.csv"), row.names = FALSE, fileEncoding = "UTF-8")
capture.output(summary(merged$Genomic_Offset), file = file.path(outdir, "offset_summary.txt"))

cor_lat <- cor.test(merged$latitude, merged$Genomic_Offset, method = "spearman", exact = FALSE)
cor_abslat <- cor.test(merged$abs_latitude, merged$Genomic_Offset, method = "spearman", exact = FALSE)

sink(file.path(outdir, "offset_latitude_correlation.txt"))
cat("Spearman correlation: offset vs latitude\n")
print(cor_lat)
cat("\nSpearman correlation: offset vs absolute latitude\n")
print(cor_abslat)
sink()

p_lat <- ggplot(merged, aes(x = latitude, y = Genomic_Offset)) +
  geom_point(alpha = 0.75, size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  theme_bw(base_size = 12) +
  labs(x = "Latitude", y = "Genomic offset",
       title = paste0("Offset vs latitude\nSpearman rho = ", round(cor_lat$estimate, 3), ", p = ", signif(cor_lat$p.value, 3)))
ggsave(file.path(outdir, "offset_vs_latitude.png"), p_lat, width = 7, height = 5, dpi = 300)

p_abslat <- ggplot(merged, aes(x = abs_latitude, y = Genomic_Offset)) +
  geom_point(alpha = 0.75, size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  theme_bw(base_size = 12) +
  labs(x = "Absolute latitude", y = "Genomic offset",
       title = paste0("Offset vs absolute latitude\nSpearman rho = ", round(cor_abslat$estimate, 3), ", p = ", signif(cor_abslat$p.value, 3)))
ggsave(file.path(outdir, "offset_vs_abs_latitude.png"), p_abslat, width = 7, height = 5, dpi = 300)

if ("Continent" %in% colnames(merged)) {
  cont_df <- merged %>% filter(!is.na(Continent), Continent != "")
  if (nrow(cont_df) > 0) {
    cont_order <- cont_df %>% group_by(Continent) %>% summarise(med = median(Genomic_Offset), .groups = "drop") %>% arrange(med) %>% pull(Continent)
    cont_df$Continent <- factor(cont_df$Continent, levels = cont_order)
    p_cont <- ggplot(cont_df, aes(x = Continent, y = Genomic_Offset)) +
      geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15, alpha = 0.65, size = 1.8) +
      theme_bw(base_size = 12) + labs(x = "Continent", y = "Genomic offset") +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
    ggsave(file.path(outdir, "offset_by_continent.png"), p_cont, width = 8, height = 5, dpi = 300)
    cont_summary <- cont_df %>% group_by(Continent) %>% summarise(n = n(), mean = mean(Genomic_Offset), median = median(Genomic_Offset), sd = sd(Genomic_Offset), .groups = "drop")
    write.csv(cont_summary, file.path(outdir, "offset_by_continent_summary.csv"), row.names = FALSE)
  }
}

threshold <- quantile(merged$Genomic_Offset, probs = 0.90, na.rm = TRUE)
high_offset <- merged %>% filter(Genomic_Offset >= threshold) %>% arrange(desc(Genomic_Offset))
keep_cols <- c("strain_key", "Assembly Accession", "Organism Scientific Name", "Organism Qualifier", "geo_loc_name", "latitude", "longitude", "Country_Code", "Continent", "Genomic_Offset")
keep_cols <- keep_cols[keep_cols %in% colnames(high_offset)]
write.csv(high_offset[, keep_cols, drop = FALSE], file.path(outdir, "high_offset_samples.csv"), row.names = FALSE)

if (has_maps) {
  world_df <- map_data("world")
  p_map <- ggplot() +
    geom_polygon(data = world_df, aes(x = long, y = lat, group = group), fill = "grey90", color = "white", linewidth = 0.2) +
    geom_point(data = merged, aes(x = longitude, y = latitude, color = Genomic_Offset), size = 2.2, alpha = 0.85) +
    scale_color_gradient(low = "skyblue", high = "red") + coord_quickmap() + theme_bw(base_size = 12) +
    labs(x = NULL, y = NULL, color = "Genomic\noffset", title = "Geographic distribution of genomic offset")
  ggsave(file.path(outdir, "offset_world_map.png"), p_map, width = 11, height = 6, dpi = 300)
}

compact_cols <- c("strain_key", "Organism Scientific Name", "geo_loc_name", "latitude", "longitude", "Country_Code", "Continent", "Genomic_Offset")
compact_cols <- compact_cols[compact_cols %in% colnames(merged)]
write.csv(merged[, compact_cols, drop = FALSE], file.path(outdir, "metadata_offset_compact.csv"), row.names = FALSE)

cat("Offset geography analysis complete. Results saved to: ", outdir, "\n")

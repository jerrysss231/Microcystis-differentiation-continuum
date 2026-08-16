# Microcystis population structure: fastBAPS, snapclust and ANI
# Inputs are resolved relative to MICROCYSTIS_PROJECT_ROOT.

suppressPackageStartupMessages({
  library(fastbaps)
  library(adegenet)
  library(ape)
  library(pheatmap)
})

PROJECT_ROOT <- normalizePath(
  Sys.getenv("MICROCYSTIS_PROJECT_ROOT", unset = "."),
  mustWork = FALSE
)

data_path <- function(...) file.path(PROJECT_ROOT, "data", ...)
result_path <- function(...) file.path(PROJECT_ROOT, "results", ...)

fasta_file <- data_path("snp", "gubbins.filtered_polymorphic_sites.fasta")
ani_matrix_file <- data_path("ani", "fastani_all_vs_all.txt.matrix")

out_dir <- result_path("population_structure")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

output_csv <- file.path(out_dir, "Master_Clade_Assignments.csv")
output_aic_plot <- file.path(out_dir, "Snapclust_AIC_Curve.pdf")
output_ani_heatmap <- file.path(out_dir, "FastANI_Heatmap.pdf")

required_files <- c(fasta_file, ani_matrix_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing input file(s):\n", paste(missing_files, collapse = "\n"))
}

cat("Running fastBAPS...\n")
sparse_data <- import_fasta_sparse_nt(fasta_file)
baps_fit <- fast_baps(sparse_data)
part_baps <- best_baps_partition(sparse_data, baps_fit)

cat("Running snapclust...\n")
dna_obj <- read.dna(fasta_file, format = "fasta")
gen_data <- DNAbin2genind(dna_obj)

aic_res <- snapclust.choose.k(max = 10, x = gen_data, IC = AIC)

pdf(output_aic_plot, width = 6, height = 5)
plot(
  aic_res,
  type = "b",
  pch = 19,
  xlab = "Number of clusters (K)",
  ylab = "AIC"
)
dev.off()

best_k <- as.numeric(names(which.min(aic_res)))
snap_fit <- snapclust(gen_data, k = best_k)
part_snap <- snap_fit$group

strain_names <- names(part_baps)
master_df <- data.frame(
  Strain_ID = strain_names,
  Clade = paste0("Clade_", part_baps[strain_names]),
  fastBAPS_Optimised = paste0("Clade_", part_baps[strain_names]),
  Snapclust = paste0("Clade_", part_snap[strain_names]),
  stringsAsFactors = FALSE
)

write.csv(master_df, output_csv, row.names = FALSE)

cat("Rendering ANI heatmap...\n")
tryCatch({
  ani_data <- read.table(
    ani_matrix_file,
    header = FALSE,
    row.names = 1,
    fill = TRUE
  )
  ani_matrix <- as.matrix(ani_data)
  storage.mode(ani_matrix) <- "numeric"

  pdf(output_ani_heatmap, width = 10, height = 10)
  pheatmap(
    ani_matrix,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    display_numbers = FALSE,
    border_color = NA,
    main = "fastANI genome similarity (%)"
  )
  dev.off()
}, error = function(e) {
  warning("ANI heatmap was skipped: ", conditionMessage(e))
})

cat(
  "Population-structure analysis completed. Preferred snapclust K = ",
  best_k, "\n", sep = ""
)

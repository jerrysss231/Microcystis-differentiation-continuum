PROJECT_ROOT <- normalizePath(
  Sys.getenv("MICROCYSTIS_PROJECT_ROOT", unset = "."),
  mustWork = FALSE
)

data_path <- function(...) file.path(PROJECT_ROOT, "data", ...)
result_path <- function(...) file.path(PROJECT_ROOT, "results", ...)

required_packages <- c("dplyr", "ggplot2", "FSA", "readxl", "tidyr", "ggsignif")
missing_pkgs <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0) {
  stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "))
}
invisible(lapply(required_packages, library, character.only = TRUE))

font_choice <- Sys.getenv("MICROCYSTIS_FONT", unset = "Times New Roman")

stats_path <- data_path("recombination", "Gubbins_Recombination_Stats_140.xlsx")
clade_path <- result_path("population_structure", "Master_Clade_Assignments.csv")
output_dir <- result_path("recombination_parameters")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stats <- readxl::read_excel(stats_path)
clade_map <- read.csv(clade_path)
clade_candidates <- c("Clade", "fastBAPS_Optimised", "fastBAPS_Classic")
clade_col <- clade_candidates[clade_candidates %in% colnames(clade_map)][1]
if (length(clade_col) == 0) stop("No lineage-assignment column found in clade map.")
clade_map$Clade <- clade_map[[clade_col]]

stats$Node_clean <- gsub("_genomic_fna$", "", stats$Node)
clade_map$Strain_clean <- gsub("_genomic_fna$", "", clade_map$Strain_ID)

data <- stats %>%
  left_join(
    clade_map %>% select(Strain_clean, Clade),
    by = c("Node_clean" = "Strain_clean")
  ) %>%
  filter(!is.na(Clade)) %>%
  mutate(Clade_Label = gsub("Clade_", "", Clade))

data$Clade_Label <- factor(
  data$Clade_Label,
  levels = sort(as.numeric(unique(data$Clade_Label)))
)

param_summary <- data %>%
  group_by(Clade_Label) %>%
  summarise(
    Count = n(),
    Rho_Theta_Mean = round(mean(`rho/theta`, na.rm = TRUE), 3),
    Rho_Theta_Min = round(min(`rho/theta`, na.rm = TRUE), 3),
    Rho_Theta_Max = round(max(`rho/theta`, na.rm = TRUE), 3),
    Rho_Theta_Range = paste0(Rho_Theta_Min, " - ", Rho_Theta_Max),
    r_m_Mean = round(mean(`r/m`, na.rm = TRUE), 3),
    r_m_Min = round(min(`r/m`, na.rm = TRUE), 3),
    r_m_Max = round(max(`r/m`, na.rm = TRUE), 3),
    r_m_Range = paste0(r_m_Min, " - ", r_m_Max)
  )

write.csv(
  param_summary,
  file.path(output_dir, "Gubbins_RhoTheta_RM_Summary.csv"),
  row.names = FALSE
)

cat("\nKruskal-Wallis test for rho/theta:\n")
print(kruskal.test(`rho/theta` ~ Clade_Label, data = data))
cat("\nKruskal-Wallis test for r/m:\n")
print(kruskal.test(`r/m` ~ Clade_Label, data = data))

dunn_rho <- dunnTest(`rho/theta` ~ Clade_Label, data = data, method = "bh")
dunn_rm <- dunnTest(`r/m` ~ Clade_Label, data = data, method = "bh")

get_sig_data <- function(dunn_res_df, alpha = 0.05) {
  sig_rows <- dunn_res_df %>% filter(P.adj < alpha)
  if (nrow(sig_rows) == 0) return(NULL)
  comparisons <- lapply(strsplit(sig_rows$Comparison, " - "), trimws)
  annotations <- sapply(sig_rows$P.adj, function(p) {
    if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)
  })
  list(comparisons = comparisons, annotations = annotations, count = nrow(sig_rows))
}

sig_rho <- get_sig_data(dunn_rho$res)
sig_rm <- get_sig_data(dunn_rm$res)

clades <- levels(data$Clade_Label)
nature_colors <- c(
  "#1D70F5", "#F54242", "#D28DF0", "#F5A11D", "#20D49C", "#33A8FF",
  "#8B3AD6", "#F0E624", "#3DEDD5", "#F74DE5", "#E8828D", "#3DE340"
)
my_palette <- if (length(clades) > 12) {
  colorRampPalette(nature_colors)(length(clades))
} else {
  nature_colors[seq_along(clades)]
}

custom_theme <- theme_minimal(base_size = 14, base_family = font_choice) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, color = "black", face = "bold"),
    axis.text = element_text(color = "black", size = 13),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    plot.title = element_blank()
  )

y_max_rho <- max(data$`rho/theta`, na.rm = TRUE)
p_rho <- ggplot(data, aes(x = Clade_Label, y = `rho/theta`, fill = Clade_Label)) +
  geom_violin(trim = FALSE, alpha = 0.8, color = "black", linewidth = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA, linewidth = 0.5) +
  scale_fill_manual(values = my_palette) +
  custom_theme + labs(y = expression(rho/theta))

if (!is.null(sig_rho)) {
  y_pos_rho <- seq(from = y_max_rho * 1.05, by = y_max_rho * 0.12, length.out = sig_rho$count)
  p_rho <- p_rho + geom_signif(
    comparisons = sig_rho$comparisons,
    annotations = sig_rho$annotations,
    y_position = y_pos_rho,
    tip_length = 0.02,
    vjust = -0.5,
    textsize = 4,
    family = font_choice
  ) + scale_y_continuous(expand = expansion(mult = c(0.05, 0.12 * sig_rho$count + 0.1)))
}

y_max_rm <- max(data$`r/m`, na.rm = TRUE)
p_rm <- ggplot(data, aes(x = Clade_Label, y = `r/m`, fill = Clade_Label)) +
  geom_violin(trim = FALSE, alpha = 0.8, color = "black", linewidth = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA, linewidth = 0.5) +
  scale_fill_manual(values = my_palette) +
  custom_theme + labs(y = "r/m")

if (!is.null(sig_rm)) {
  y_pos_rm <- seq(from = y_max_rm * 1.05, by = y_max_rm * 0.12, length.out = sig_rm$count)
  p_rm <- p_rm + geom_signif(
    comparisons = sig_rm$comparisons,
    annotations = sig_rm$annotations,
    y_position = y_pos_rm,
    tip_length = 0.02,
    vjust = -0.5,
    textsize = 4,
    family = font_choice
  ) + scale_y_continuous(expand = expansion(mult = c(0.05, 0.12 * sig_rm$count + 0.1)))
}

ggsave(file.path(output_dir, "RhoTheta_Violin_Nature.png"), plot = p_rho, width = 7, height = 5, dpi = 600, type = "cairo")
ggsave(file.path(output_dir, "RtoM_Violin_Nature.png"), plot = p_rm, width = 7, height = 5, dpi = 600, type = "cairo")

cat("Recombination-parameter analysis complete. Results saved to: ", output_dir, "\n")

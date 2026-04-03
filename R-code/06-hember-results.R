library(patchwork) # For multi-panel figures
library(data.table)
library(ggplot2)
library(tidyr)
library(forcats)
# Local path to the results
res_hember_DT <- readRDS("R://ahuels/AIRCO/01_projects/019_adrc_bb_prediction_machine_learning/harf-paper/results/hember-lab/harf_arf_hember_results.rds")

# Keep iteration-level data
dt <- res_hember_DT[, .(
  Data, chunck_size, UVD.result, MMD_rbk.result,
  algorithm, evidence, iteration = 1:.N
)]
dt$chunck_size <- as.factor(dt$chunck_size)

# Pivot longer for metrics
dt_long <- pivot_longer(
  dt,
  cols = c(UVD.result, MMD_rbk.result),
  names_to = "Perf_measure",
  values_to = "Perf_value"
)

# Aggregate mean & SD **by algorithm and evidence as well**
agg <- as.data.table(dt_long)[, .(
  mean_perf = -log(mean(Perf_value)),
  sd_perf = -log(sd(Perf_value))
), by = .(Data, chunck_size, Perf_measure, algorithm, evidence)]

agg <- agg[ , `:=`(
  evidence = factor(ifelse(evidence, "Evi.", "No evi."),
                    levels = c("No evi.", "Evi.")),
  Perf_measure = fct_recode(Perf_measure,
                            UVD = "UVD.result",
                            MMD = "MMD_rbk.result")
)]


plot_chunk <- function(agg, perf_measure, data_name) {
  my_plot <- ggplot(agg[Perf_measure == perf_measure & Data == data_name], aes(x = Data, y = factor(chunck_size), fill = mean_perf)) +
    geom_tile(color = "white") +
    facet_grid(algorithm ~ evidence, scales = "free_y", switch = "x") +
    scale_fill_gradient(low = "#132B43", high = "#56B1F7", name = "Mean Perf") +
    theme_minimal(base_size = 12) +
    geom_text(aes(label = round(mean_perf, 2)), size = 3, color = "white") +
    labs(x = NULL, y = NULL, title = NULL) +
    theme(
      axis.text.x = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.spacing = unit(0.5, "lines"),
      legend.position = "none",
      legend.title = element_blank(),
      plot.margin = margin(t = 20, r = 5, b = 5, l = 5),
      panel.spacing.x = unit(0.01, "lines") # reduce spacing between columns
    ) 
  return(my_plot)
}

# =====================================
# Generate UVD plots for each dataset
# =====================================
#
lake_UVD_plot <- plot_chunk(agg = agg, "UVD", "lake") +  
  theme(strip.text.x = element_blank(),
        strip.text.y = element_blank()) + labs(y = "Chunk size")
manno_UVD_plot <- plot_chunk(agg = agg, "UVD", "manno") + labs(y = NULL) + 
  theme(strip.text.x = element_blank())
li_UVD_plot <- plot_chunk(agg = agg, "UVD", "li") + 
  theme(strip.text.y = element_blank()) + labs(y = "Chunk size") 
patel_UVD_plot <- plot_chunk(agg = agg, "UVD", "patel") + labs(y = NULL)

UVD_plot <- plot_grid(
  lake_UVD_plot,
  manno_UVD_plot,
  li_UVD_plot,
  patel_UVD_plot,
  labels = c("Lake", "Manno", "Li", "Patel"),
  ncol = 2, nrow = 2,
  align = "v", axis = "l"
)
print(UVD_plot)

library(patchwork) # For multi-panel figures
library(data.table)
library(ggplot2)
library(tidyr)
library(forcats)
# Local path to the results
res_hember_DT <- readRDS("R://ahuels/AIRCO/01_projects/019_adrc_bb_prediction_machine_learning/harf-paper/results/hember-lab/harf_arf_hember_results.rds")
res_hember_DT <- readRDS("R://ahuels/AIRCO/01_projects/019_adrc_bb_prediction_machine_learning/harf-paper/results/hember-lab/harf_arf_hember_res_inter.rds")
# Keep iteration-level data
dt <- res_hember_DT[, .(
  Data, chunck_size, UVD.result, MMD_rbk.result,
  algorithm, evidence, iteration = 1:.N, num_btwn_pcs,
  ARI_ORG, NMI_ORG, ARI_SYN, NMI_SYN, ARI_DIFF, NMI_DIFF
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
), by = .(Data, chunck_size, Perf_measure, algorithm, evidence, num_btwn_pcs)]

agg <- agg[ , `:=`(
  evidence = factor(ifelse(evidence, "Evi.", "No evi."),
                    levels = c("No evi.", "Evi.")),
  Perf_measure = fct_recode(Perf_measure,
                            UVD = "UVD.result",
                            MMD = "MMD_rbk.result")
)]


plot_chunk <- function(agg, perf_measure, data_name, n_pcs = 2) {
  my_plot <- ggplot(agg[Perf_measure == perf_measure & Data == data_name & num_btwn_pcs %in% c(NA,2)], aes(x = Data, y = factor(chunck_size), fill = mean_perf)) +
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


# =====================================
# Generate MMD plots for each dataset
# =====================================
#

lake_MMD_plot <- plot_chunk(agg = agg, "MMD", "lake") +  
  theme(strip.text.x = element_blank(),
        strip.text.y = element_blank()) + labs(y = "Chunk size")
manno_MMD_plot <- plot_chunk(agg = agg, "MMD", "manno") + labs(y = NULL) + 
  theme(strip.text.x = element_blank())
li_MMD_plot <- plot_chunk(agg = agg, "MMD", "li") + 
  theme(strip.text.y = element_blank()) + labs(y = "Chunk size") 
patel_MMD_plot <- plot_chunk(agg = agg, "MMD", "patel") + labs(y = NULL)

MMD_plot <- plot_grid(
  lake_MMD_plot,
  manno_MMD_plot,
  li_MMD_plot,
  patel_UVD_plot,
  labels = c("Lake", "Manno", "Li", "Patel"),
  ncol = 2, nrow = 2,
  align = "v", axis = "l"
)
print(MMD_plot)

# =====================================
# Generate ARI plots for each dataset
# =====================================
#
# Boxplots of ARI_ORG and ARI_SYN by algorithm and evidence and data
dt[ , `:=`(
  ARI_DIFF = ARI_ORG - ARI_SYN,
  NMI_DIFF = NMI_ORG - NMI_SYN
)]
cluster_perf_long <- pivot_longer(
  dt,
  cols = c(ARI_ORG, ARI_SYN, NMI_ORG, NMI_SYN, ARI_DIFF, NMI_DIFF),
  names_to = "cluster_metric",
  values_to = "cluster_perf"
)
cluster_perf_long <- as.data.table(cluster_perf_long)
# create a variable Data_version from cluster_metric by extracting the suffix (ORG or SYN)
cluster_perf_long[, Data_version := fcase(
  cluster_metric %in% c("ARI_ORG", "NMI_ORG"), "ORG",
  cluster_metric %in% c("ARI_SYN", "NMI_SYN"), "SYN",
  cluster_metric %in% c("ARI_DIFF", "NMI_DIFF"), "DIFF"
)]
# Replace Syn by algorithm name
cluster_perf_long[Data_version == "SYN", Data_version := algorithm]

# Extract ARI and NMI from cluster_metric
cluster_perf_long[, cluster_metric := ifelse(grepl("ARI", cluster_metric),
                                            "ARI", "NMI")] 

cluster_perf_plot <- ggplot(
  cluster_perf_long[
    num_btwn_pcs %in% c(2, NA) &
      evidence == FALSE &
      Data_version == "DIFF"
  ],
  aes(x = as.factor(chunck_size), y = cluster_perf, fill = algorithm)
) +
  geom_boxplot(outlier.alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_grid(cluster_metric ~ Data) +
  theme_bw() +
  labs(x = "Chunk size", y = "Empirical ARI and NMI") 

print(cluster_perf_plot)

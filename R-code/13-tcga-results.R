library(patchwork) # For multi-panel figures
library(data.table)
library(ggplot2)
library(tidyr)
library(forcats)
library(cowplot)
library(ggsci)
# Local path to the results
res_tcga_DT <- fread(file.path(res_tcga_data_dir, "harf_arf_tcga_res.rds"))
# res_tcga_DT <- readRDS("R://ahuels/AIRCO/01_projects/019_adrc_bb_prediction_machine_learning/harf-paper/results/tcga-tgex/harf_arf_tcga_res.rds")
# Keep iteration-level data
dt <- res_tcga_DT[, .(
  Data, chunck_size, CD, UVD.result, MMD_rbk.result,
  algorithm, evidence, iteration = 1:.N, num_btwn_pcs,
  AUC_RF_DIFF, AUC_Lasso_DIFF, time
)]
dt[algorithm == "ARF", chunck_size := 55]
dt <- dt[chunck_size %in% c(10, 55), ]
dt$chunck_size <- as.factor(dt$chunck_size)

my_theme <- theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(
      size = 18,
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      size = 12,
      hjust = 0.5,
      color = "gray40"
    ),
    
    axis.text = element_text(
      color = "black",
      size = 12
    ),
    
    strip.text = element_text(
      size = 12
    ),
    
    panel.grid.minor = element_blank(),
    
    panel.grid.major.x = element_blank(),
    
    
    legend.position = "none"
  )

# ---------------------------------------
# Keep only "No Evidence"
# ---------------------------------------
plot_dt <- copy(dt[evidence == FALSE])

plot_dt[, algorithm := factor(
  algorithm,
  levels = c("ARF", "HARF")
)]

plot_dt[, Data := factor(Data)]

# ---------------------------------------
# Convert to long format
# ---------------------------------------
plot_long <- melt(
  plot_dt,
  measure.vars = c("AUC_RF_DIFF", "AUC_Lasso_DIFF"),
  variable.name = "MODEL",
  value.name = "AUC_DIFF"
)

# Cleaner labels
plot_long[, MODEL := factor(
  MODEL,
  levels = c("AUC_RF_DIFF", "AUC_Lasso_DIFF"),
  labels = c("RF", "Lasso")
)]

# ---------------------------------------
# Combined figure
# ---------------------------------------
combined_plot <- ggplot(
  plot_long,
  aes(
    x = MODEL,
    y = AUC_DIFF,
    fill = algorithm
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray50",
    linewidth = 0.7
  ) +
  
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.65,
    alpha = 1,
    outlier.shape = 21,
    outlier.size = 1.8,
    color = "gray20"
  ) +
  
  facet_wrap(
    ~ Data,
    ncol = 2
  ) +
  
  labs(
    x = "Prediction Model",
    y = expression(Delta ~ "AUC"),
    fill = "Synthesizer"
  ) +
  
  theme_classic() +
  theme_bw() +
  scale_fill_npg(labels = c("ARF", expression(italic("h")*"-ARF")))  +
  
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  
  my_theme +
  
  theme(
    legend.position = "bottom",
    
    legend.title = element_blank(),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    
    strip.text = element_text(
      size = 13
    ),
    
    axis.text.x = element_text(
      size = 11
    )
  )

combined_plot

# Save the plot
ggsave(
  # filename = "R://ahuels/AIRCO/01_projects/019_adrc_bb_prediction_machine_learning/harf-paper/results/tcga-tgex/downstreampred.eps",
   filename = file.path(res_tcga_data_dir, "downstreampred.eps"),
  plot = combined_plot,
  width = 13, height = 12, units = "cm", dpi = 400
)



# Pivot longer for metrics
dt_long <- pivot_longer(
  dt,
  cols = c(CD, UVD.result, MMD_rbk.result),
  names_to = "Perf_measure",
  values_to = "Perf_value"
)
dt_long <- as.data.table(dt_long)
dt_long[, Perf_measure := factor(
  Perf_measure,
  levels = c("CD", "UVD.result", "MMD_rbk.result"),
  labels = c("CD", "UVD", "MMD")
)]

# =====================================
# Average performance across datasets and metric
# =====================================
#
dt_long_avg <- pivot_longer(
  dt,
  cols = c(UVD.result, MMD_rbk.result, CD, time),
  names_to = "Perf_measure",
  values_to = "Perf_value"
)
dt_long_avg <- as.data.table(dt_long_avg)[, .(
  mean_perf = mean(Perf_value, na.rm = TRUE),
  sd_perf = sd(Perf_value, na.rm = TRUE)
), by = .(Data, algorithm, evidence, Perf_measure, num_btwn_pcs)]

dt_long_avg[, Perf_measure := fcase(
  Perf_measure == "UVD.result", "UVD",
  Perf_measure == "MMD_rbk.result", "MMD",
  Perf_measure == "CD", "CD",
  Perf_measure == "time", "Time"
)]
dt_long_avg[, Data := tools::toTitleCase(Data)]
dt_long_avg[, evidence := factor(ifelse(evidence, "Evi.", "No evi."),
                                 levels = c("No evi.", "Evi."))]
perf_plot <- ggplot(
  dt_long_avg[num_btwn_pcs %in% c(2, NA) &
                Perf_measure %in% c("UVD", "MMD", "CD") & 
                evidence == "No evi."],
  aes(x = Perf_measure, y = mean_perf, fill = algorithm)
) +
  
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_perf - sd_perf, ymax = mean_perf
                    + sd_perf), width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~Data) +
  theme_classic() +
  theme_bw() +
  scale_fill_npg(labels = c("ARF", expression(italic("h")*"-ARF"))) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  labs(x = "Empirical value", y = "Empirical value") +
  scale_x_discrete(labels = c("UVD", "MMD", "CD")) +
  theme(
    strip.text = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
print(perf_plot)

## Save
ggsave(
  # filename = "R://ahuels/AIRCO/01_projects/019_adrc_bb_prediction_machine_learning/harf-paper/results/tcga-tgex/downstreampreddist.eps",
  filename = file.path(res_tcga_data_dir, "downstreampreddist.eps"),
  plot = perf_plot,
  width = 13, height = 12, units = "cm", dpi = 400
)


# Print the average time 
print(dt_long_avg[Perf_measure == "Time" & evidence == "No evi.", .(Data, algorithm, mean_perf, sd_perf)])




# ==============================================================================
# Appendix
# ==============================================================================
# Boxplot of AUC and MMD against chunk size variation
#
dt_chunk <- res_tcga_DT[, .(
  Data, chunck_size, CD, UVD.result, MMD_rbk.result,
  algorithm, evidence, iteration = 1:.N, num_btwn_pcs,
  AUC_RF_DIFF, AUC_Lasso_DIFF, time
)]
dt_chunk[algorithm == "ARF", chunck_size := 55]
dt_chunk$chunck_size <- as.factor(dt_chunk$chunck_size)

pred_perf_long <- pivot_longer(
  dt_chunk,
  cols = c(AUC_RF_DIFF, AUC_Lasso_DIFF, MMD_rbk.result),
  names_to = "prediction_metric",
  values_to = "prediction_perf"
)
pred_perf_long <- as.data.table(pred_perf_long)


# RF AUC plot
rf_auc_plot <- ggplot(
  pred_perf_long[
    num_btwn_pcs %in% c(2, NA) &
      evidence == FALSE &
      prediction_metric == "AUC_RF_DIFF"
  ],
  aes(x = as.factor(chunck_size), y = prediction_perf, fill = algorithm)
) +
  geom_boxplot(outlier.alpha = 1, outlier.size = 1) +
  # geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~ Data, ncol = 2) +
  theme_classic()+
  theme_bw() +
  scale_fill_npg(labels = c("ARF", expression(italic("h")*"-ARF"))) +
  theme(legend.position = "none", legend.title = element_blank(),
        axis.text.x = element_text(
          angle = 45,
          hjust = 1
        )) +
  labs(x = "Chunk size", y = "Empirical AUC") +
  scale_x_discrete(
    name = "Chunk size",
    labels = c("5", "10", "15", "20", "25", "30", "35", "40", "45", "50", "ARF")
  ) +
  scale_y_continuous(labels = function(x) sprintf("%.3f", x)) +
  theme(
    axis.title.x = element_blank()
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
print(rf_auc_plot)


# Lasso AUC plot
lasso_auc_plot <- ggplot(
  pred_perf_long[
    num_btwn_pcs %in% c(2, NA) &
      evidence == FALSE &
      prediction_metric == "AUC_Lasso_DIFF"
  ],
  aes(x = as.factor(chunck_size), y = prediction_perf, fill = algorithm)
) +
  geom_boxplot(outlier.alpha = 1, outlier.size = 1) +
  # geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~ Data, ncol = 2) +
  theme_classic()+
  theme_bw() +
  scale_fill_npg(labels = c("ARF", expression(italic("h")*"-ARF"))) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.text.x = element_text(
          angle = 45,
          hjust = 1
        )) +
  labs(x = "Chunk size", y = expression(Delta ~ "AUC")) +
  scale_x_discrete(
    name = "Chunk size",
    labels = c("5", "10", "15", "20", "25", "30", "35", "40", "45", "50", "ARF")
  ) +
  scale_y_continuous(labels = function(x) sprintf("%.3f", x)) +
  theme(
    axis.title.x = element_blank()
  )
print(lasso_auc_plot)


legend <- get_legend(lasso_auc_plot + theme(legend.position = "bottom"))

pred_auc_plot <- (rf_auc_plot / lasso_auc_plot) + 
  plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = 'A') & 
  theme(legend.position = "bottom",
        strip.text = element_text(size = 14),
        plot.caption = element_text(hjust = 0.5, size = 12),
  )
print(pred_auc_plot)

# Save the plot
ggsave(
  # filename = "R://ahuels/AIRCO/01_projects/019_adrc_bb_prediction_machine_learning/harf-paper/results/tcga-tgex/downstreampredchunkauc.eps",
  filename = file.path(res_tcga_data_dir, "downstreampredchunkauc.eps"),
  plot = pred_auc_plot,
  width = 13, height = 20, units = "cm", dpi = 400
)

# MMD plot
mmd_plot <- ggplot(
  pred_perf_long[
    num_btwn_pcs %in% c(2, NA) &
      evidence == FALSE &
      prediction_metric == "MMD_rbk.result"
  ],
  aes(x = as.factor(chunck_size), y = prediction_perf, fill = algorithm)
) +
  geom_boxplot(outlier.alpha = 1, outlier.size = 1) +
  # geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~ Data, ncol = 2) +
  theme_classic()+
  theme_bw() +
  scale_fill_npg(labels = c("ARF", expression(italic("h")*"-ARF"))) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.text.x = element_text(
          angle = 45,
          hjust = 1
        )) +
  labs(x = "Chunk size", y = expression(Delta ~ "MMD")) +
  scale_x_discrete(
    name = "Chunk size",
    labels = c("5", "10", "15", "20", "25", "30", "35", "40", "45", "50", "ARF")
  ) +
  scale_y_continuous(labels = function(x) sprintf("%.3f", x)) 
print(mmd_plot)

# Save the plot
ggsave(
  # filename = "R://ahuels/AIRCO/01_projects/019_adrc_bb_prediction_machine_learning/harf-paper/results/tcga-tgex/downstreampredchunkmmd.eps",
  filename = file.path(res_tcga_data_dir, "downstreampredchunkmmd.eps"),
  plot = mmd_plot,
  width = 13, height = 10, units = "cm", dpi = 400
)


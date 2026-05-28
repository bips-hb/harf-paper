library(patchwork) # For multi-panel figures
library(data.table)
library(ggplot2)
library(tidyr)
library(forcats)
library(cowplot)
library(ggsci)

ad_results_dir <- "R://ahuels/AIRCO/01_projects/019_adrc_bb_prediction_machine_learning/harf-paper/results/ad"
ad_res_pred_DT <- readRDS(file.path(ad_results_dir, "harf_arf_ad_pred.rds"))
ad_res_pred_DT$evidence <- ad_res_pred_DT$evidence.x
# Keep iteration-level data
dt <- ad_res_pred_DT[, .(
  chunck_size, CD, UVD.result, MMD_rbk.result,
  algorithm, evidence, iteration = 1:.N, num_btwn_pcs,
  AUC_RF_DIFF, AUC_Lasso_DIFF, MCC_RF_DIFF, MCC_Lasso_DIFF,
  time, prop_synth
)]
dt[algorithm == "ARF", chunck_size := 55]
dt <- dt[chunck_size %in% c(50, 55) & prop_synth == 1, ]
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

plot_long <- melt(
  ad_res_pred_DT,
  id.vars = c("algorithm"),
  measure.vars = c(
    "AUC_RF_DIFF",
    "AUC_Lasso_DIFF",
    "MCC_RF_DIFF",
    "MCC_Lasso_DIFF"
  ),
  variable.name = "Metric",
  value.name = "Value"
)

# Model labels
plot_long[, Model := fifelse(
  grepl("RF", Metric),
  "RF",
  "LASSO"
)]

# Metric labels
plot_long[, Measure := fifelse(
  grepl("AUC", Metric),
  'Delta~"ROC-AUC"',
  'Delta~"MCC"'
)]

# ---------------------------------------
# Combined figure
# ---------------------------------------

combined_plot <- ggplot(
  plot_long,
  aes(
    x = Model,
    y = Value,
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
    alpha = 0.95,
    outlier.shape = 21,
    outlier.size = 1.8,
    color = "gray20"
  ) +
  
  facet_wrap(
    ~ Measure,
    scales = "free_y",
    nrow = 1,
    labeller = label_parsed
  ) +
  
  scale_fill_npg(
    labels = c("ARF", expression(italic("h") * "-ARF"))
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  
  labs(
    x = "Prediction Model",
    y = "Performance Difference (Original - Synthetic)",
    fill = "Synthesizer"
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
    
    strip.background = element_rect(
      fill = "gray95",
      color = "black"
    ),
    
    strip.text = element_text(
      size = 13
    ),
    
    axis.text.x = element_text(
      size = 11
    )
  )

combined_plot

library(patchwork) # For multi-panel figures
library(data.table)
library(ggplot2)
library(tidyr)
library(forcats)
library(cowplot)
library(ggsci)

ad_results_dir <- "R://ahuels/AIRCO/01_projects/019_adrc_bb_prediction_machine_learning/harf-paper/results/ad"
ad_res_pred_DT <- readRDS(file.path(ad_results_dir, "harf_arf_ad_pred.rds"))
# Keep iteration-level data
dt <- ad_res_pred_DT[, .(
  chunck_size, CD, UVD.result, MMD_rbk.result,
  algorithm, evidence, iteration = 1:.N, 
  AUC_RF_ORG, AUC_Lasso_ORG, MCC_RF_ORG, MCC_Lasso_ORG,
  AUC_RF_SYN, AUC_Lasso_SYN, MCC_RF_SYN, MCC_Lasso_SYN,
  AUC_RF_DIFF, AUC_Lasso_DIFF, MCC_RF_DIFF, MCC_Lasso_DIFF,
  time, prop_synth, ApoE_ORG, ApoE_SYN
)]
dt[algorithm == "ARF", chunck_size := 55]
dt <- dt[chunck_size %in% c(15, 55) & prop_synth == 1, ]
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
  "Lasso"
)]

# Metric labels
plot_long[, Measure := fifelse(
  grepl("AUC", Metric),
  'Delta~"AUC"',
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
    alpha = 1,
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
    y = "Empirical value",
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

## ========================================================================
# Plot the ApOE OR for original vs synthetic data
## ========================================================================
##
# ---------------------------------------
# Compute Odds Ratios + 95% CI
# ---------------------------------------

# Original data
org_summary <- plot_dt[
  ,
  .(
    Mean = mean(ApoE_ORG, na.rm = TRUE),
    SD   = sd(ApoE_ORG, na.rm = TRUE),
    N    = .N
  )
]

org_summary[
  ,
  `:=`(
    OR      = exp(Mean),
    CI_low  = exp(Mean - 1.96 * SD / sqrt(N)),
    CI_high = exp(Mean + 1.96 * SD / sqrt(N))
  )
]

org_dt <- data.table(
  Group   = "Original",
  OR      = org_summary$OR,
  CI_low  = org_summary$CI_low,
  CI_high = org_summary$CI_high
)

# Synthetic data
syn_dt <- plot_dt[
  ,
  .(
    Mean = mean(ApoE_SYN, na.rm = TRUE),
    SD   = sd(ApoE_SYN, na.rm = TRUE),
    N    = .N
  ),
  by = algorithm
]

syn_dt[
  ,
  `:=`(
    OR      = exp(Mean),
    CI_low  = exp(Mean - 1.96 * SD / sqrt(N)),
    CI_high = exp(Mean + 1.96 * SD / sqrt(N))
  )
]

syn_dt[, Group := algorithm]

# Combine
plot_or <- rbindlist(list(
  org_dt,
  syn_dt[, .(Group, OR, CI_low, CI_high)]
))
# Replace HARF by h-ARF
plot_or[Group == "HARF", Group := "h-ARF"]
plot_or[, Group := factor(
  Group,
  levels = c("ARF", "h-ARF", "Original")
)]

# ---------------------------------------
# Horizontal OR plot
# ---------------------------------------
or_plot <- ggplot(
  plot_or,
  aes(
    x = Group,
    y = OR,
    fill = Group
  )
) +
  
  geom_col(
    width = 0.6
  ) +
  
  geom_errorbar(
    aes(
      ymin = CI_low,
      ymax = CI_high
    ),
    width = 0.15,
    linewidth = 0.7
  ) +
  
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "gray50"
  ) +
  
  coord_flip() +
  
  scale_fill_npg() +
  
  labs(
    y = "ApoE Odds Ratio (95% CI)",
    x = NULL
  ) +
  
  theme_bw() +
  my_theme +
  
  theme(
    legend.position = "none",
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    
    axis.text = element_text(size = 12)
  ) 
or_plot

## ========================================================================
## Plot distribution metrics
## ========================================================================
## 
metrics_long <- melt(
  plot_dt,
  measure.vars = c("CD", "UVD.result", "MMD_rbk.result"),
  variable.name = "Metric",
  value.name = "Value"
)

metrics_long[
  ,
  Metric := factor(
    Metric,
    levels = c("CD", "UVD.result", "MMD_rbk.result"),
    labels = c("CD", "UVD", "MMD")
  )
]
summary_metrics <- metrics_long[
  ,
  .(
    mean_value = mean(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    n = .N
  ),
  by = .(Metric, algorithm)
]

summary_metrics[
  ,
  `:=`(
    lower = mean_value - 1.96 * sd / sqrt(n),
    upper = mean_value + 1.96 * sd / sqrt(n)
  )
]
p_metrics <- ggplot(
  summary_metrics,
  aes(
    x = Metric,
    y = mean_value,
    fill = algorithm
  )
) +
  
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65
  ) +
  
  geom_errorbar(
    aes(
      ymin = lower,
      ymax = upper
    ),
    position = position_dodge(width = 0.75),
    width = 0.15,
    linewidth = 0.8
  ) +
  scale_fill_npg() +
  
  theme_bw() +
  my_theme +
  
  labs(
    x = "Distribution metric",
    y = "Empirical value"
  ) +
  
  my_theme +
  
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    )
  )

p_metrics

# ==============================================================================
# Combine all plots into a multi-panel figure
# ==============================================================================
combined_plot <- combined_plot  + 
  theme(legend.position = "none")
or_plot <- or_plot + theme(legend.position = "none") + guides(fill = "none")
p_metrics <- p_metrics + theme(legend.position = "none") + guides(fill = "none")

final_plot <- combined_plot /
  (p_metrics | or_plot) +
  
  plot_layout(
    heights = c(1.5, 1.5),
    guides = "collect"   # <-- ONE shared legend
  ) +
  
  plot_annotation(
    title = NULL,
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 16, hjust = 0.5, face = NULL),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
      plot.tag = element_text(size = 14, face = NULL)
    )
  ) &
  
  theme(
    legend.position = "bottom"
  )

final_plot

## Save the final figure
ggsave(
  filename = file.path(ad_results_dir, "adsyncompared.eps"),
  plot = final_plot,
  width = 16, height = 18, units = "cm", dpi = 400
)


## ========================================================================
# Prediction performance for original, h-ARF and ARF data
## ========================================================================
##

metrics <- c(
  "AUC_RF_ORG", "AUC_Lasso_ORG",
  "MCC_RF_ORG", "MCC_Lasso_ORG",
  "AUC_RF_SYN", "AUC_Lasso_SYN",
  "MCC_RF_SYN", "MCC_Lasso_SYN"
)

plot_long <- melt(
  dt,
  measure.vars = metrics,
  variable.name = "Metric",
  value.name = "Value"
)

plot_long[, c("Measure", "Model", "Type") :=
            tstrsplit(Metric, "_", fixed = TRUE)]

plot_long[, Group := fifelse(
  Type == "ORG",
  "Original",
  fifelse(algorithm == "HARF", "h-ARF", "ARF")
)]

plot_long[, Group := factor(Group, levels = c("Original", "h-ARF", "ARF"))]
# simpler fallback
# plot_long[, Group := fifelse(Type == "ORG", "Original", "Synthetic")]

plot_long[, Model := factor(Model, levels = c("RF", "Lasso"))]
plot_long[, Measure := factor(Measure, levels = c("AUC", "MCC"))]

reversed_npg_colors <- rev(pal_npg("nrc")(3))
prefperfplot <- ggplot(plot_long, aes(x = Model, y = Value, fill = Group)) +
  
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.65,
    outlier.size = 1.3
  ) +
  
  facet_wrap(~ Measure, scales = "free_y") +
  
  labs(
    x = "Prediction Model",
    y = "Empirical value",
    fill = "Data Source"
  ) +
  scale_fill_manual(values = reversed_npg_colors) +
  
  my_theme +
  
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    
    strip.text = element_text(size = 13),
    axis.text.x = element_text(size = 12),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    
    strip.background = element_rect(
      fill = "gray95",
      color = "black",
      linewidth = 0.8
    ),
    
    panel.spacing = unit(0.8, "lines")
  )
print(prefperfplot)
# save the predictive performance plot
ggsave(
  plot = prefperfplot,
    filename = file.path(ad_results_dir, "adpredperf.eps"),
  width = 16, height = 10, units = "cm", dpi = 400
)

# ==============================================================================
# Combine all plots into a multi-panel figure
# ==============================================================================
prefperfplot <- prefperfplot  + 
  theme(legend.position = "none")
or_plot <- or_plot + theme(legend.position = "none") + guides(fill = "none")
p_metrics <- p_metrics + theme(legend.position = "none") + guides(fill = "none")

pred_final_plot <- prefperfplot /
  (p_metrics | or_plot) +
  
  plot_layout(
    heights = c(1.5, 1.5),
    guides = "collect"   # <-- ONE shared legend
  ) +
  
  plot_annotation(
    title = NULL,
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 16, hjust = 0.5, face = NULL),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
      plot.tag = element_text(size = 14, face = NULL)
    )
  ) &
  
  theme(
    legend.position = "bottom"
  )

pred_final_plot

## Save the final figure
ggsave(
  filename = file.path(ad_results_dir, "adpredsyncompared.eps"),
  plot = pred_final_plot,
  width = 16, height = 18, units = "cm", dpi = 400
)


# Mean AUC and MCC for each group and predictive model, with 95% confidence intervals
summary_metrics <- plot_long[ 
  ,
  .(
    mean_value = mean(Value, na.rm = TRUE),
    lower = quantile(Value, 0.025, na.rm = TRUE),
    upper = quantile(Value, 0.975, na.rm = TRUE)
  ),
  by = .(Group, Model, Measure)
]


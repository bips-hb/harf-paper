library(patchwork) # For multi-panel figures
library(data.table)
library(ggplot2)
library(tidyr)
library(forcats)
library(cowplot)
library(ggsci)
# Local path to the results
res_tcga_DT <- readRDS("R://ahuels/AIRCO/01_projects/019_adrc_bb_prediction_machine_learning/harf-paper/results/tcga-tgex/harf_arf_tcga_res.rds")
# Keep iteration-level data
dt <- res_tcga_DT[, .(
  Data, chunck_size, CD, UVD.result, MMD_rbk.result,
  algorithm, evidence, iteration = 1:.N, num_btwn_pcs,
  AUC_RF_DIFF, AUC_Lasso_DIFF, time
)]
dt[algorithm == "ARF", chunck_size := 55]
dt <- dt[chunck_size %in% c(10, 55), ]
dt$chunck_size <- as.factor(dt$chunck_size)

# Pivot longer for metrics
dt_long <- pivot_longer(
  dt,
  cols = c(UVD.result, MMD_rbk.result),
  names_to = "Perf_measure",
  values_to = "Perf_value"
)

# =====================================
# Average performance across datasets and metric
# =====================================
#
dt_long_avg <- pivot_longer(
  dt,
  cols = c(UVD.result, MMD_rbk.result, CD, AUC_RF_DIFF, AUC_Lasso_DIFF, time),
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
  Perf_measure == "AUC_RF_DIFF", "AUC_RF_DIFF",
  Perf_measure == "AUC_Lasso_DIFF", "AUC_Lasso_DIFF",
  Perf_measure == "time", "Time"
)]
dt_long_avg[, Data := tools::toTitleCase(Data)]
dt_long_avg[, evidence := factor(ifelse(evidence, "Evi.", "No evi."),
                                 levels = c("No evi.", "Evi."))]
perf_plot <- ggplot(
  dt_long_avg[num_btwn_pcs %in% c(2, NA) &
                Perf_measure %in% c("UVD", "MMD", "CD", "AUC_RF_DIFF", "AUC_Lasso_DIFF") & 
                evidence == "No evi."],
  aes(x = Perf_measure, y = mean_perf, fill = algorithm)
) +
  # Gray background for ARI and NMI
  geom_rect(
    inherit.aes = FALSE,
    xmin = 3.5, xmax = 5.5,
    ymin = -Inf, ymax = Inf,
    fill = "honeydew",
    color = "grey50",
    linetype = "dashed",
    linewidth = 0.4
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
  scale_x_discrete(labels = c("UVD", "MMD", "CD", "AUC RF Diff", "AUC Lasso Diff")) +
  theme(
    strip.text = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
print(perf_plot)


# ggplot auc only
auc_plot <- ggplot(
  dt_long_avg[num_btwn_pcs %in% c(2, NA) &
                Perf_measure %in% c("AUC_RF_DIFF", "AUC_Lasso_DIFF") & 
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
  scale_x_discrete(labels = c("AUC RF Diff", "AUC Lasso Diff")) +
  theme(
    strip.text = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
print(auc_plot)

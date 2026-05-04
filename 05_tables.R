library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(huxtable)
library(stringr)
library(forcats)
#library(flextable)
library(magrittr)
library(forcats)
library(viridisLite)
library(patchwork)

ci_root <- file.path(
  "Simulation/model_estimates_unif"
)



# term_rate <- c("low","medium", "high")
term_rate <- c("low")
params <- c("cov1","cov2","cov3")
params_baseline <- c("shape", "rate")
models <- c("flexsurv","msm","msm_age","nhm","mipd_iter","streams","streams2")



all_ci <- expand_grid(event_rate = term_rate) %>%
  pmap_dfr(function( event_rate) {
    base_dir <- file.path(ci_root, event_rate)
    if (!dir.exists(base_dir)) return(tibble())

    subdirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
    if (!length(subdirs)) return(tibble())

    map_dfr(subdirs, function(sd) {
      # prefer explicit summary if present
      fp <- if (file.exists(file.path(sd, "ic_summary.rds"))) {
        file.path(sd, "ic_summary.rds")
      } else {
        fs <- list.files(sd, pattern = "\\.[Rr][Dd][Ss]$", full.names = TRUE)
        if (!length(fs)) return(tibble())
        fs[1]
      }

      df <- readRDS(fp)

      df %>%
        mutate(
          # use function args explicitly to avoid clobbering by any existing columns
          event_rate        = as.character(.env$event_rate),

          est_mean          = round(est_mean, 3),
          est_sd            = round(est_sd, 3),
          `rel_bias(%)`     = round(rel_bias * 100, 3),
          `coverage(%)`     = round(coverage_rate * 100, 2),
          `type_I_error(%)` = round(type_I_rate * 100, 2),
          `type_II_error(%)`= round(type_II_rate * 100, 2)
        ) %>%
        select(-rel_bias, -coverage_rate, -type_I_rate, -type_II_rate)
      # NOTE: do NOT mutate(across(everything(), as.character(.))) here.
    })
  })

save(all_ci,file = "Simulation/model_estimates_unif/overall_metrics_cleaned.RData")

all_ci <-  all_ci %>%
  filter(model %in% models)


cov_ci <- all_ci %>%
  filter(param %in% c("cov1","cov2","cov3"))

cov_ci <- cov_ci %>%
  dplyr::ungroup() %>%
  dplyr::select(
    
    model,
    trans_idx,
    event_rate,
    est_mean,
    est_sd,
    'rel_bias(%)',
    'coverage(%)',
    'type_I_error(%)',
    'type_II_error(%)'
  )%>%
  arrange(event_rate,trans_idx,model) %>% 
  mutate(
    # model = dplyr::recode(model,
    #                       "flexsurv_real"  = "flex_real",
    #                       "flexsurv_panel" = "flex_panel",
    #                       "mipd_iter_model" = "mipd",
    #                       "dropout_adaptive_tau" = "Streams",
    #                       "dropout_tau080_095" = "fixed_tresh"
    # ),
    est_mean           = suppressWarnings(as.numeric(est_mean)),
    est_sd             = suppressWarnings(as.numeric(est_sd)),
    `rel_bias(%)`      = suppressWarnings(as.numeric(`rel_bias(%)`)),
    `coverage(%)`      = suppressWarnings(as.numeric(`coverage(%)`)),
    `type_I_error(%)`  = suppressWarnings(as.numeric(`type_I_error(%)`)),
    `type_II_error(%)` = suppressWarnings(as.numeric(`type_II_error(%)`))
  )

cov_ci %>% flextable::flextable() %>% flextable::save_as_docx(path="Simulation/model_estimates/Result_covariates.docx")
base_ci <- all_ci %>%
  filter(param %in% params_baseline) %>%
  select(model, param,
         event_rate, trans_idx,
         est_mean, est_sd,
         `rel_bias(%)`, `coverage(%)`,
         `type_I_error(%)`, `type_II_error(%)`) %>%
  mutate(
    # opzionale: rinomina modelli
    # model = recode(model,
    #                "flexsurv_real"  = "flex_real",
    #                "flexsurv_panel" = "flex_panel",
    #                "mipd_iter_model" = "mipd_iter",
    #                "dropout_adaptive_tau" = "Streams",
    #                "dropout_tau080_095" = "fixed_tresh"),

    est_mean           = suppressWarnings(as.numeric(est_mean)),
    est_sd             = suppressWarnings(as.numeric(est_sd)),
    `rel_bias(%)`      = suppressWarnings(as.numeric(`rel_bias(%)`)),
    `coverage(%)`      = suppressWarnings(as.numeric(`coverage(%)`)),
    `type_I_error(%)`  = suppressWarnings(as.numeric(`type_I_error(%)`)),
    `type_II_error(%)` = suppressWarnings(as.numeric(`type_II_error(%)`))
  )



rename_coef <- function(df, trans) {
  newname <- switch(as.character(trans),
                    "1" = "beta01",
                    "2" = "beta02",
                    "3" = "beta12"
  )
  names(df)[names(df) == "est_mean"] <- newname
  df
}

by_trans <- split(cov_ci, cov_ci$trans_idx)

tables <- Map(function(df, nm) {
  rename_coef(df, nm) %>% dplyr::select(-trans_idx)
}, by_trans, names(by_trans))



build_ht_for_transition <- function(df, trans = 1, docx_file = NULL, save_rds = NULL, rds_file = NULL) {
  # --- pick beta column & label by transition ---
  beta_col <- switch(as.character(trans),
                     "1" = "beta01",
                     "2" = "beta02",
                     "3" = "beta12",
                     stop("trans must be 1, 2, or 3"))
  beta_lab <- switch(as.character(trans),
                     "1" = "beta01",
                     "2" = "beta02",
                     "3" = "beta12")

  # --- reshape to wide with the chosen beta ---
  metrics <- c(beta_col, "rel_bias(%)", "est_sd", "coverage(%)", "type_I_error(%)", "type_II_error(%)")

  # #desired <- c( "flex_real", "cox", "flex_panel", "msm", "msm_age", "nhm", "mipd", "fixed_tresh", "Streams")
  # desired <- c("flex_panel", "msm_age", "nhm", "mipd", "fixed_tresh", "Streams")

  df_wide <- df %>%
    select(scheme_id, model, event_rate,
           !!sym(beta_col), est_sd,  `rel_bias(%)`, `coverage(%)`, `type_I_error(%)`, `type_II_error(%)`) %>%
    #mutate(model = factor(model, levels = desired)) %>%
    arrange(scheme_id, model) %>%
    pivot_wider(
      names_from  = event_rate,
      values_from = all_of(metrics),
      names_glue  = "{event_rate}__{.value}"
    ) %>%
    arrange(scheme_id, model) %>%
    group_by(scheme_id) %>%
    mutate(Scenario = if_else(row_number() == 1, as.character(scheme_id), "")) %>%
    ungroup() %>%
    select(Scenario, model,
           starts_with("low__"),
           starts_with("medium__"),
           starts_with("high__"))




  # --- headers (big spanners + metric labels) ---
  rates   <- c("low","medium","high")
  n       <- length(metrics)

  hdr1 <- c("Scenario", "Model",
            rep("Low event rates",      n),
            rep("Moderate event rates", n),
            rep("High event rates",     n))

  metric_label <- function(m) dplyr::case_when(
    m == beta_col         ~ beta_lab,
    m == "rel_bias(%)"    ~ "Rel.bias(\\%)",
    m == "coverage(%)"    ~ "Cov.Rate(\\%)",
    m == "type_I_error(%)"~ "TypeI(\\%)",
    m == "type_II_error(%)" ~ "TypeII(\\%)",
    TRUE ~ m
  )
  hdr2 <- c("", "", rep(metric_label(metrics), times = length(rates)))

  mat <- rbind(hdr1, hdr2, as.matrix(df_wide))
  ht  <- as_hux(mat, add_colnames = FALSE)

  # mark header rows
  ht <- set_header_rows(ht, 1:2, TRUE)
  ht <- set_escape_contents(ht, FALSE)  # allow LaTeX in headers

  # merge big headers
  ht <- merge_cells(ht, row = 1, col = 3:(2 + n))
  ht <- merge_cells(ht, row = 1, col = (3 + n):(2 + 2*n))
  ht <- merge_cells(ht, row = 1, col = (3 + 2*n):(2 + 3*n))

  # formatting
  last <- ncol(ht)
  sep1 <- 2 + n
  sep2 <- 2 + 2*n

  ht <- set_align(ht, row = 1, col = 3:last, value = "center")
  ht <- set_valign(ht, row = 1, col = 3:last, value = "middle")
  ht <- set_bold(ht, row = 1:2, col = 1:last, value = TRUE)
  ht <- set_top_border(ht,    row = 1, col = everywhere, value = 1.2)
  ht <- set_bottom_border(ht, row = 2, col = everywhere, value = 1.2)
  ht <- set_right_border(ht,  row = everywhere, col = sep1, value = 0.6)
  ht <- set_right_border(ht,  row = everywhere, col = sep2, value = 0.6)
  ht <- set_align(ht, row = everywhere, col = 1:2,   value = "center")
  ht <- set_align(ht, row = everywhere, col = 3:last, value = "center")
  ht <- set_na_string(ht, "–")
  if ("set_all_padding" %in% getNamespaceExports("huxtable")) {
    ht <- set_all_padding(ht, 0.5)
  }

  # thick top border before each Scenario change (account for 2 header rows)
  starts <- which(df_wide$Scenario != "")
  if (length(starts) > 1) {
    rows_to_border <- starts[-1] + 2
    ht <- set_top_border(ht, row = rows_to_border, col = everywhere, value = 1.2)
  }

  if (isTRUE(save_rds)) {
    saveRDS(ht, file = rds_file)
  }

  # export to Word if requested
  if (is.null(docx_file)) {
    docx_file <- paste0("results_table_trans", trans, ".docx")
  }
  if ("quick_docx" %in% getNamespaceExports("huxtable")) {
    huxtable::quick_docx(ht, file = docx_file, open = TRUE)
  } else {
    # fallback via officer + flextable
    if (!requireNamespace("officer", quietly = TRUE) ||
        !requireNamespace("flextable", quietly = TRUE)) {
      warning("Install 'officer' and 'flextable' for DOCX export; returning huxtable only.")
      return(ht)
    }
    ft <- as_flextable(ht)
    ft <- padding(ft, padding.left = 1, padding.right = 1, padding.top = 1, padding.bottom = 1, part = "all")
    w <- c(0.7, 0.9, rep(0.55, ncol_part(ft, "body") - 2))
    ft <- width(ft, j = 1:length(w), width = w, part = "all")
    ft <- fontsize(ft, size = 7, part = "all")
    ft <- autofit(ft)


    doc <- officer::read_docx()
    doc <- officer::body_add_flextable(doc, ft)
    print(doc, target = docx_file)
    message("Saved: ", docx_file)
  }

  ht
}


ht1 <- build_ht_for_transition(tables[[1]], trans = 1, docx_file = "results_t1.docx", save_rds = TRUE, rds_file = "results_t1.rds")

ht2 <- build_ht_for_transition(tables[[2]], trans = 2, docx_file = "results_t2.docx", save_rds = TRUE, rds_file = "results_t2.rds")

ht3 <- build_ht_for_transition(tables[[3]], trans = 3, docx_file = "results_t3.docx", save_rds = TRUE, rds_file = "results_t3.rds")


saveRDS(base_ci, file = "baseline_params.rds")



## plots for baseline curves

desired <- c("flex_panel", "nhm", "mipd", "fixed_tresh", "Streams")

params_gt_all <- map_dfr(unique(cov_ci_plot$scheme_id), load_gt_for_scheme)
base_ci <- base_ci %>%
  mutate(model = factor(model, levels = desired)) %>%
  arrange(scheme_id, model)


S_gompertz <- function(t, rate, shape) {
  if (abs(shape) < 1e-10) {
    exp(-rate * t)  # exponential limit as shape -> 0
  } else {
    exp(-(rate/shape) * (exp(shape * t) - 1))
  }
}

plot_gompertz<- function(df_in,n_grid = 300,t_min = 60, t_max = 140,rate_is_log = TRUE) {


  params_wide <- df_in %>%
    select(event_rate, model, trans_idx, param, est_mean) %>%
    pivot_wider(names_from = param, values_from = est_mean) %>%
    mutate(
      rate_nat  = if (rate_is_log) exp(rate) else rate,
      shape_nat = shape
    ) %>%
    filter(is.finite(rate_nat), rate_nat > 0, is.finite(shape_nat))


  curves <- params_wide %>%
    mutate(trans_idx = as.character(trans_idx)) %>%
    group_by(event_rate, model, trans_idx) %>%
    do({
      r  <- .$rate_nat[1]
      s  <- .$shape_nat[1]
      tt <- seq(t_min, t_max, length.out = n_grid)
      tibble(t = tt, S = S_gompertz(tt, r, s))
    }) %>%
    ungroup()

  ggplot(curves, aes(x = t, y = S, color = model, group = model)) +
    geom_line(linewidth = 1) +
    facet_grid(rows = vars(event_rate),
               cols = vars(trans_idx),
               labeller = labeller(
                 trans_idx = c(`1` = "Transition 1",
                               `2` = "Transition 2",
                               `3` = "Transition 3")
               )) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_color_viridis_d(option = "D", end = 0.9) +
    theme_minimal(base_size = 13) +
    labs(
      title = "Gompertz survival by transition",
      x = "time",
      y = "S(t)",
      color = "Model"
    )
}

p <- plot_gompertz(base_ci)
print(p)

# ggsave("gompertz_t1.png", plots[["1"]], width=6, height=4, dpi=150)


#  ===============
#   Bias of cov
#  ===============

desired <- c("flex_panel", "msm_age", "nhm", "mipd", "Streams")

cov_ci <- cov_ci %>%
  filter(model %in% desired) %>%
  mutate(model = factor(model, levels = desired)) %>%
  arrange(scheme_id, model) %>%
  mutate( model = recode(model,
                         "flex_panel"  = "a",
                         "msm_age" = "b",
                         "nhm" = "c",
                         "mipd" = "d",
                         "Streams" = "e")

  )




# Facet by scheme id, event_rate fixed

event_rate_label <- "high"   # "low" | "medium" | "high"
gt_base_dir <- file.path("Simulation", "simulation_results", "ground_truth_params", event_rate_label)



load_gt_for_scheme <- function(sid) {
  f <- file.path(gt_base_dir, sprintf("params_scheme_%02d.rds", sid))
  if (!file.exists(f)) {
    stop(sprintf("Missing GT file for scheme_id=%d: %s", sid, f))
  }
  params_gt <- readRDS(f)

  as.data.frame(params_gt) %>%
    rownames_to_column("trans") %>%
    pivot_longer(cols = starts_with("cov"),
                 names_to = "param",
                 values_to = "true_value") %>%
    mutate(
      trans_idx = case_when(
        trans == "0->1" ~ 1,
        trans == "0->2" ~ 2,
        trans == "1->2" ~ 3,
        TRUE            ~ NA_real_
      ),
      scheme_id = sid
    ) %>%
    select(scheme_id, trans_idx, param, true_value)
}


params_gt_all <- map_dfr(unique(cov_ci_plot$scheme_id), load_gt_for_scheme)


cov_ci_diff <- cov_ci %>%
  inner_join(params_gt_all, by = c("scheme_id", "trans_idx")) %>%
  mutate(
    diff    = est_mean - true_value,
    ci_low  = (est_mean - 1.96 * est_sd) - true_value,
    ci_high = (est_mean + 1.96 * est_sd) - true_value
  )

cov_ci_diff<- cov_ci_diff %>%
  filter(event_rate == event_rate_label)


ggplot(cov_ci_diff, aes(x = param, y = diff, color = model, group = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.2,
                position = position_dodge(width = 0.5)) +
  facet_grid(rows = vars(scheme_id),
             cols = vars(trans_idx),
             labeller = labeller(
               trans_idx = c(`1`="Transition 1", `2`="Transition 2", `3`="Transition 3")
             )) +
  theme_minimal(base_size = 14) + scale_color_viridis_d(option = "D", end = 0.9) +
  labs(color = "Model",
       title = sprintf("Bias of the estimated exposure effect"))




# Facet by event_rate, fixed scheme_id

df_plot <- cov_ci_diff %>%
  filter(event_rate=="low")

pd <- position_dodge(width = 0.5)

for (sid in sort(unique(df_plot$scheme_id))) {
  df <- df_plot %>% filter(scheme_id == sid)

  p <- ggplot(df, aes(x = param, y = diff, color = model, group = model)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(position = pd, size = 2) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                  width = 0.2, position = pd) +
    facet_grid(
      rows = vars(event_rate),
      cols  = vars(trans_idx),
      labeller = labeller(
        trans_idx  = c(`1` = "Transition 1", `2` = "Transition 2", `3` = "Transition 3"),
        event_rate = c(high = "")   # <-- blank the “high” strip label
      )
    ) +
    theme_minimal(base_size = 14) +
    scale_color_viridis_d(option = "D", end = 0.9) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      axis.text.x   = element_blank(),   # <-- hide "cov"
      axis.ticks.x  = element_blank(),
      strip.text.y  = element_blank(),   # (optional) hide all row strip text
      strip.background = element_blank(), # (optional) remove strip background
      legend.position = "bottom"
    ) +
    labs(title = NULL, color = "Model", y = NULL, x = NULL)

  print(p)
  ggsave(sprintf("Plots/bias_low_by_scheme_%02d.png", sid),
         p, width = 8, height = 4, dpi = 300)
}


##### covergae #########

rate <-  "high"


df_cov_plot <- cov_ci %>%
  mutate(
    trans = factor(trans_idx, levels = c(1,2,3)),
    cov_color = case_when(
      `coverage(%)` >= 95 ~ "green",
      `coverage(%)` >= 80 & `coverage(%)` < 95 ~ "yellow",
      TRUE ~ "red"
    ),
    cov = `coverage(%)`
  ) %>%
  filter(event_rate == rate)

n_schemes <- dplyr::n_distinct(df_cov_plot$scheme_id)


p_cov <- ggplot(df_cov_plot, aes(x = model, y = trans, fill = cov_color)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f", cov)), size = 3) +
  facet_grid(~scheme_id) +
  scale_fill_manual(
    name = "Coverage",
    values = c("green" = "#1a9850", "yellow" = "#fee08b", "red" = "#d73027"),
    labels = c(">=95%", "80-95%", "<80%")
  ) +
  labs(x = NULL, y = NULL,
       title = "Coverage heatmap by model × transition (facets: scheme_id)",
       subtitle = "Green ≥95%, Yellow 80–95%, Red <80%") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_cov)
ggsave("heatmap_coverage.png", p_cov,
       width = max(8, 1.6 * n_schemes), height = 4, dpi = 300)

############# type I error vs power ##################


rate <- "high"  # or "medium", "high"



# --- Prepare data ---
df_avg <- cov_ci %>%
  filter(event_rate == rate) %>%
  group_by(model, trans_idx) %>%
  summarise(
    avg_typeI = mean(`type_I_error(%)`, na.rm = TRUE),
    avg_power = 100 * mean(1 - `type_II_error(%)` / 100, na.rm = TRUE),
    avg_coverage = mean(`coverage(%)`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(trans = factor(trans_idx, levels = c(3, 2, 1))) %>%  # 1 at bottom
  pivot_longer(
    cols = c(avg_typeI, avg_power, avg_coverage),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    value_color = case_when(
      metric == "avg_typeI" & value <= 5 ~ "≤ 5%",
      metric == "avg_typeI" & value > 5 & value <= 15 ~ "5–15%",
      metric == "avg_typeI" & value > 15 ~ "> 15%",
      metric == "avg_power" & value >= 80 ~ "≥ 80%",
      metric == "avg_power" & value >= 50 & value < 80 ~ "50–80%",
      metric == "avg_power" & value < 50 ~ "< 50%",
      metric == "avg_coverage" & value >= 95 ~ "≥ 95%",
      metric == "avg_coverage" & value >= 75 & value < 95 ~ "75–95%",
      metric == "avg_coverage" & value < 75 ~ "< 75%"
    ),
    metric = factor(metric,
                    levels = c("avg_typeI", "avg_power", "avg_coverage"),
                    labels = c("Type I Error (%)", "Power (%)", "Coverage (%)"))
  )

# --- Color palette ---
status_colors <- c("≤ 5%"="#1be850", "5–15%"="#fee08b", "> 15%"="firebrick",
                   "≥ 80%"="#1be850", "50–80%"="#fee08b", "< 50%"="firebrick",
                   "≥ 95%"="#1be850", "75–95%"="#fee08b", "< 75%"="firebrick")

# --- Helper function for per-metric plots ---
make_plot <- function(df, metric_name, legend_title, level_order) {
  df_metric <- df %>%
    filter(metric == metric_name) %>%
    mutate(value_color = factor(value_color, levels = level_order))

  ggplot(df_metric, aes(x = model, y = trans, fill = value_color)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.1f", value)), size = 3) +
    scale_fill_manual(
      values = status_colors[level_order],
      name = legend_title,
      drop = FALSE
    ) +
    labs(
      title = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", size = 13)
    )
}

# --- Define legend orders for each metric ---
order_typeI <- c("≤ 5%", "5–15%", "> 15%")
order_power <- c("≥ 80%", "50–80%", "< 50%")
order_cov   <- c("≥ 95%", "75–95%", "< 75%")

# --- Create each plot with its ordered legend ---

p_typeI <- make_plot(df_avg, "Type I Error (%)", "Type I Error (%)", order_typeI) + labs(x = NULL , y = NULL)
p_power <- make_plot(df_avg, "Power (%)", "Power (%)", order_power) + labs(x = NULL , y = "Transition")
p_cov   <- make_plot(df_avg, "Coverage (%)", "Coverage (%)", order_cov) + labs(x = "Model" , y = NULL)


p_final <- (p_typeI / p_power / p_cov)
p_final
ggsave(sprintf("Plots/heatmap_avg_metrics_by_trans_%s.png", rate), p_final,
       width = 8, height = 4, dpi = 300)


## Computational time



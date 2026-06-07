###############################################################################
############################# PLOTS #######################################
###############################################################################

###############################################################################
# ---- 1: DYNAMICS PLOTS ----
###############################################################################

plot_dynamics <- function(ode_result, params_vec, targets, N_h, N_c) {
  
  # Add aggregated compartments
  ode_result$S_h_total <- with(ode_result, S0_h + SA_h + S_II_h + S_III_h)
  ode_result$S_c_total <- with(ode_result, S0_c + SA_c + S_II_c + S_III_c)
  ode_result$C_h <- with(ode_result, C0_h + CA_h + C_II_h + C_III_h)
  ode_result$C_c <- with(ode_result, C0_c + CA_c + C_II_c + C_III_c)
  ode_result$I_tot_h <- with(ode_result, I_h + I_II_h + I_III_h)
  ode_result$I_tot_c <- with(ode_result, I_c + I_II_c + I_III_c)
  ode_result$N_h_sim <- with(ode_result, S0_h + SA_h + S_II_h + S_III_h + C0_h + CA_h + C_II_h + C_III_h + I_h + I_II_h + I_III_h)
  ode_result$N_c_sim <- with(ode_result, S0_c + SA_c + S_II_c + S_III_c + C0_c + CA_c + C_II_c + C_III_c + I_c + I_II_c + I_III_c)
  
  # Targets in absolute numbers
  target_C_h <- targets$portage_h * N_h
  target_C_c <- targets$portage_c * N_c
  
  # Hopital - all compartments
  p_h_all <- ode_result %>%
    dplyr::select(time, S0_h, SA_h, C0_h, CA_h, I_h, S_II_h, C_II_h, I_II_h, S_III_h, C_III_h, I_III_h) %>%
    reshape2::melt(id.vars = "time") %>%
    ggplot2::ggplot(ggplot2::aes(time, value, color = variable)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Temps (jours)", y = "Nombre d'individus",
      color = "Hopital :", title = "Dynamique du modele - Hopital (tous les compartiments)"
    )
  
  # Hopital - totals only (N, S, C, I)
  p_h_totals <- ode_result %>%
    dplyr::select(time, N_h_sim, S_h_total, C_h, I_tot_h) %>%
    reshape2::melt(id.vars = "time") %>%
    ggplot2::ggplot(ggplot2::aes(time, value, color = variable)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Temps (jours)", y = "Nombre d'individus",
      color = "Hopital :", title = "Verification des totaux hopital (N, S, C, I)"
    )
  
  # Hopital - aggregates only (S total, C total, I total, N) + target in legend
  h_agg <- ode_result %>%
    dplyr::select(time, S_h_total, C_h, I_tot_h, N_h_sim) %>%
    reshape2::melt(id.vars = "time") %>%
    dplyr::mutate(
      variable = factor(variable,
                        levels = c("S_h_total", "C_h", "I_tot_h", "N_h_sim"),
                        labels = c("S total", "C total", "I total", "N total"))
    )

  p_h_target <- ggplot2::ggplot(h_agg, ggplot2::aes(time, value, color = variable)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_hline(
      ggplot2::aes(yintercept = target_C_h, color = "Cible prevalence colonisation"),
      linetype = "dashed", linewidth = 0.8, inherit.aes = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = c("S total" = "#1F78B4",
                 "C total" = "#33A02C",
                 "I total" = "#F1C40F",
                 "N total" = "#C51BCE",
                 "Cible prevalence colonisation" = "red2"),
      breaks = c("S total", "C total", "I total", "N total", "Cible prevalence colonisation")
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes = list(linetype = c("solid", "solid", "solid", "solid", "dashed"),
                          linewidth = c(0.9, 0.9, 0.9, 0.9, 0.8))
    )) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Temps (jours)", y = "Nombre d'individus",
      color = "Hopital :", title = "Dynamique du modele - Hopital"
    )
  
  # Communaute - all compartments
  p_c_all <- ode_result %>%
    dplyr::select(time, S0_c, SA_c, C0_c, CA_c, I_c, S_II_c, C_II_c, I_II_c, S_III_c, C_III_c, I_III_c) %>%
    reshape2::melt(id.vars = "time") %>%
    ggplot2::ggplot(ggplot2::aes(time, value, color = variable)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Temps (jours)", y = "Nombre d'individus",
      color = "Communaute :", title = "Dynamique du modele - Communaute (tous les compartiments)"
    )
  
  # Communaute - totals only (N, S, C, I)
  p_c_totals <- ode_result %>%
    dplyr::select(time, N_c_sim, S_c_total, C_c, I_tot_c) %>%
    reshape2::melt(id.vars = "time") %>%
    ggplot2::ggplot(ggplot2::aes(time, value, color = variable)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Temps (jours)", y = "Nombre d'individus",
      color = "Communaute :", title = "Verification des totaux communaute (N, S, C, I)"
    )
  
  # Communaute - aggregates only (S total, C total, I total, N) + target in legend
  c_agg <- ode_result %>%
    dplyr::select(time, S_c_total, C_c, I_tot_c, N_c_sim) %>%
    reshape2::melt(id.vars = "time") %>%
    dplyr::mutate(
      variable = factor(variable,
                        levels = c("S_c_total", "C_c", "I_tot_c", "N_c_sim"),
                        labels = c("S total", "C total", "I total", "N total"))
    )

  p_c_target <- ggplot2::ggplot(c_agg, ggplot2::aes(time, value, color = variable)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_hline(
      ggplot2::aes(yintercept = target_C_c, color = "Cible prevalence colonisation"),
      linetype = "dashed", linewidth = 0.8, inherit.aes = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = c("S total" = "#1F78B4",
                 "C total" = "#33A02C",
                 "I total" = "#F1C40F",
                 "N total" = "#C51BCE",
                 "Cible prevalence colonisation" = "red2"),
      breaks = c("S total", "C total", "I total", "N total", "Cible prevalence colonisation")
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes = list(linetype = c("solid", "solid", "solid", "solid", "dashed"),
                          linewidth = c(0.9, 0.9, 0.9, 0.9, 0.8))
    )) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Temps (jours)", y = "Nombre d'individus",
      color = "Communaute :", title = "Dynamique du modele - Communaute"
    )

  return(list(hospital_all = p_h_all, hospital_totals = p_h_totals, hospital_target = p_h_target,
              community_all = p_c_all, community_totals = p_c_totals, community_target = p_c_target))
}

###############################################################################
# ---- 2: GRID SEARCH PLOTS ----
###############################################################################

plot_grid_search_all_points <- function(grid_result, param_x, param_y, title_text) {
  p <- ggplot(grid_result$grid, aes(x = .data[[param_x]], y = .data[[param_y]])) +
    geom_point(aes(color = "Valeurs testées"), alpha = 0.75, size = 1.8) +
    scale_color_manual(
      values = c("Valeurs testées" = "grey35", "Meilleure valeur" = "#2F80ED"),
      name = "Légende"
    ) +
    theme_bw() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.minor = element_line(color = "#E3E3E3"),
      panel.grid.major = element_line(color = "#D0D0D0"),
      legend.position = "right"
    ) +
    labs(
      title = title_text,
      x = param_x,
      y = param_y
    )

  if (!is.null(grid_result$best_guess)) {
    p <- p +
      geom_point(
        data = grid_result$best_guess,
        aes(x = .data[[param_x]], y = .data[[param_y]], color = "Meilleure valeur"),
        size = 3
      )
  }

  return(p)
}

plot_grid_search_beta <- function(grid_result, best_guess, targets) {
  dx <- 0.02 * diff(range(grid_result$grid$portage_h, na.rm = TRUE))
  dy <- 0.02 * diff(range(grid_result$grid$portage_c, na.rm = TRUE))
  best_label <- sprintf("Meilleure valeur\nbeta_h = %.4g\nbeta_c = %.4g",
                        best_guess$beta_h[1], best_guess$beta_c[1])
  
  ggplot(grid_result$grid, aes(x = portage_h, y = portage_c)) +
    geom_point(aes(color = "Points simulés"), alpha = 0.3, size = 0.6) +
    geom_vline(aes(xintercept = targets$portage_h, color = "Cibles", linetype = "Cibles")) +
    geom_hline(aes(yintercept = targets$portage_c, color = "Cibles", linetype = "Cibles")) +
    geom_point(data = best_guess, aes(x = portage_h, y = portage_c, color = "Meilleure valeur"), size = 3) +
    annotate("label",
             x = best_guess$portage_h[1] + dx,
             y = best_guess$portage_c[1] + dy,
             label = best_label,
             hjust = 0, vjust = 0,
             color = "#2F80ED", fill = "white", size = 3) +
    scale_color_manual(
      values = c("Points simulés" = "grey35", "Cibles" = "red2", "Meilleure valeur" = "#2F80ED"),
      name = "Légende"
    ) +
    scale_linetype_manual(values = c("Cibles" = "dashed"), guide = "none") +
    theme_bw() +
    theme(legend.position = "right") +
    labs(x = "Prévalence de portage hôpital",
         y = "Prévalence de portage communauté",
         title = "Recherche en grille (beta) : sorties simulées")
}

plot_grid_search_sigma <- function(grid_result, best_guess, targets) {
  dx <- 0.02 * diff(range(grid_result$grid$incidence_h, na.rm = TRUE))
  dy <- 0.02 * diff(range(grid_result$grid$incidence_c, na.rm = TRUE))
  best_label <- sprintf("Meilleure valeur\nsigma_h = %.4g\nsigma_c = %.4g",
                        best_guess$sigma_h[1], best_guess$sigma_c[1])
  
  ggplot(grid_result$grid, aes(x = incidence_h, y = incidence_c)) +
    geom_point(aes(color = "Points simulés"), alpha = 0.3, size = 0.6) +
    geom_vline(aes(xintercept = targets$incidence_h, color = "Cibles", linetype = "Cibles")) +
    geom_hline(aes(yintercept = targets$incidence_c, color = "Cibles", linetype = "Cibles")) +
    geom_point(data = best_guess, aes(x = incidence_h, y = incidence_c, color = "Meilleure valeur"), size = 3) +
    annotate("label",
             x = best_guess$incidence_h[1] + dx,
             y = best_guess$incidence_c[1] + dy,
             label = best_label,
             hjust = 0, vjust = 0,
             color = "#2F80ED", fill = "white", size = 3) +
    scale_color_manual(
      values = c("Points simulés" = "grey35", "Cibles" = "red2", "Meilleure valeur" = "#2F80ED"),
      name = "Légende"
    ) +
    scale_linetype_manual(values = c("Cibles" = "dashed"), guide = "none") +
    theme_bw() +
    theme(legend.position = "right") +
    labs(x = "Incidence hôpital (par jour)",
         y = "Incidence communauté (par jour)",
         title = "Recherche en grille (sigma) : sorties simulées")
}

plot_grid_search_k <- function(grid_result, best_guess, targets) {
  dx <- 0.02 * diff(range(grid_result$grid$recid_1, na.rm = TRUE))
  dy <- 0.02 * diff(range(grid_result$grid$recid_2, na.rm = TRUE))
  best_label <- sprintf("Meilleure valeur\nk_II = %.4g\nk_III = %.4g",
                        best_guess$k_II[1], best_guess$k_III[1])
  
  ggplot(grid_result$grid, aes(x = recid_1, y = recid_2)) +
    geom_point(aes(color = "Points simulés"), alpha = 0.3, size = 0.6) +
    geom_vline(aes(xintercept = targets$recid_1, color = "Cibles", linetype = "Cibles")) +
    geom_hline(aes(yintercept = targets$recid_2, color = "Cibles", linetype = "Cibles")) +
    geom_point(data = best_guess, aes(x = recid_1, y = recid_2, color = "Meilleure valeur"), size = 3) +
    annotate("label",
             x = best_guess$recid_1[1] + dx,
             y = best_guess$recid_2[1] + dy,
             label = best_label,
             hjust = 0, vjust = 0,
             color = "#2F80ED", fill = "white", size = 3) +
    scale_color_manual(
      values = c("Points simulés" = "grey35", "Cibles" = "red2", "Meilleure valeur" = "#2F80ED"),
      name = "Légende"
    ) +
    scale_linetype_manual(values = c("Cibles" = "dashed"), guide = "none") +
    theme_bw() +
    theme(legend.position = "right") +
    labs(
      x = "Taux de récidive 1",
      y = "Taux de récidive 2+",
      title = "Recherche en grille (k) : sorties simulées"
    )
}

###############################################################################
# ---- 3: TABLE PLOTS ----
###############################################################################

create_all_parameters_table <- function(params_final, N_h, N_c, alpha_eq) {
  
  sigma_A_h <- params_final["k_A"] * params_final["sigma_h"]
  sigma_A_c <- params_final["k_A"] * params_final["sigma_c"]
  sigma_II_h <- params_final["k_II"] * params_final["sigma_h"]
  sigma_II_c <- params_final["k_II"] * params_final["sigma_c"]
  sigma_III_h <- params_final["k_III"] * params_final["sigma_h"]
  sigma_III_c <- params_final["k_III"] * params_final["sigma_c"]
  
  alpha_I <- alpha_eq * params_final["w_I"]
  alpha_II <- alpha_eq * params_final["w_II"]
  alpha_III <- alpha_eq * params_final["w_III"]
  
  table_data <- data.frame(
    Parameter = c(
      "N",
      "beta",
      "sigma",
      "tau",
      "omega",
      "nu",
      "gamma",
      "epsilon",
      "p",
      "phi",
      "k_A",
      "k_II",
      "k_III",
      "sigma_A",
      "sigma_II",
      "sigma_III",
      "delta",
      "w",
      "w_I",
      "w_II",
      "w_III",
      "alpha",
      "alpha_I",
      "alpha_II",
      "alpha_III"
    ),
    Hospital = c(
      sprintf("%.0f", N_h),
      sprintf("%.6g", params_final["beta_h"]),
      sprintf("%.6g", params_final["sigma_h"]),
      sprintf("%.6g", params_final["tau_h"]),
      sprintf("%.6g", params_final["omega"]),
      sprintf("%.6g", params_final["nu"]),
      sprintf("%.6g", params_final["gamma"]),
      sprintf("%.6g", params_final["epsilon"]),
      sprintf("%.6g", params_final["p"]),
      sprintf("%.6g", params_final["phi"]),
      sprintf("%.6g", params_final["k_A"]),
      sprintf("%.6g", params_final["k_II"]),
      sprintf("%.6g", params_final["k_III"]),
      sprintf("%.6g", sigma_A_h),
      sprintf("%.6g", sigma_II_h),
      sprintf("%.6g", sigma_III_h),
      sprintf("%.6g", params_final["delta"]),
      sprintf("%.6g", params_final["w"]),
      sprintf("%.6g", params_final["w_I"]),
      sprintf("%.6g", params_final["w_II"]),
      sprintf("%.6g", params_final["w_III"]),
      sprintf("%.6g", alpha_eq),
      sprintf("%.6g", alpha_I),
      sprintf("%.6g", alpha_II),
      sprintf("%.6g", alpha_III)
    ),
    Community = c(
      sprintf("%.0f", N_c),
      sprintf("%.6g", params_final["beta_c"]),
      sprintf("%.6g", params_final["sigma_c"]),
      sprintf("%.6g", params_final["tau_c"]),
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      sprintf("%.6g", sigma_A_c),
      sprintf("%.6g", sigma_II_c),
      sprintf("%.6g", sigma_III_c),
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      ""
    ),
    Unit = c(
      "individuals",
      "day^-1",
      "day^-1",
      "day^-1",
      "day^-1",
      "-",
      "day^-1",
      "day^-1",
      "-",
      "day^-1",
      "-",
      "-",
      "-",
      "day^-1",
      "day^-1",
      "day^-1",
      "day^-1",
      "-",
      "-",
      "-",
      "-",
      "day^-1",
      "day^-1",
      "day^-1",
      "day^-1"
    )
  )
  
  table_theme <- gridExtra::ttheme_default(
    core = list(fg_params = list(fontsize = 8), bg_params = list(fill = rep(c("grey95", "white"), length.out = nrow(table_data)))),
    colhead = list(fg_params = list(fontface = "bold", fontsize = 9), bg_params = list(fill = "steelblue", col = "white"))
  )
  
  table_plot <- gridExtra::tableGrob(table_data, rows = NULL, theme = table_theme)
  title_plot <- grid::textGrob("All model parameters", gp = grid::gpar(fontsize = 14, fontface = "bold"))
  
  return(gridExtra::arrangeGrob(title_plot, table_plot, ncol = 1, heights = c(0.06, 0.94)))
}

create_calibration_table <- function(params_final, metrics, targets) {
  
  table_data <- data.frame(
    "Paramètre (vert) / cible (jaune)" = c(
      "beta_h",
      "beta_c",
      "sigma_h",
      "sigma_c",
      "k_II",
      "k_III",
      "",
      "Prévalence de colonisation hospitalière",
      "Prévalence de colonisation communautaire",
      "Incidence ICD hospitalière",
      "Incidence ICD communautaire",
      "Proportion de première récidive",
      "Proportion de récidive multiple"
    ),
    Valeur = c(
      sprintf("%.6g", params_final["beta_h"]),
      sprintf("%.6g", params_final["beta_c"]),
      sprintf("%.6g", params_final["sigma_h"]),
      sprintf("%.6g", params_final["sigma_c"]),
      sprintf("%.6g", params_final["k_II"]),
      sprintf("%.6g", params_final["k_III"]),
      "",
      sprintf("%.3f%%", metrics$portage_h * 100),
      sprintf("%.3f%%", metrics$portage_c * 100),
      sprintf("%.3f", metrics$inc_h * 100000 * 365),
      sprintf("%.3f", metrics$inc_c * 100000 * 365),
      sprintf("%.3f%%", metrics$recid_1_tot * 100),
      sprintf("%.3f%%", metrics$recid_2_tot * 100)
    ),
    Cible = c(
      rep("-", 6),
      "",
      sprintf("%.3f%%", targets$portage_h * 100),
      sprintf("%.3f%%", targets$portage_c * 100),
      sprintf("%.3f", targets$incidence_h * 100000 * 365),
      sprintf("%.3f", targets$incidence_c * 100000 * 365),
      sprintf("%.3f%%", targets$recid_1 * 100),
      sprintf("%.3f%%", targets$recid_2 * 100)
    ),
    Unité = c(
      "jour-1",
      "jour-1",
      "jour-1",
      "jour-1",
      "-",
      "-",
      "",
      "%",
      "%",
      "cas/100k/an",
      "cas/100k/an",
      "%",
      "%"
    ),
    "Erreur relative" = c(
      rep("-", 6),
      "",
      sprintf("%.2f%%", abs(metrics$errors$err_portage_h) * 100),
      sprintf("%.2f%%", abs(metrics$errors$err_portage_c) * 100),
      sprintf("%.2f%%", abs(metrics$errors$err_inc_h) * 100),
      sprintf("%.2f%%", abs(metrics$errors$err_inc_c) * 100),
      sprintf("%.2f%%", abs(metrics$errors$err_recid_1_tot) * 100),
      sprintf("%.2f%%", abs(metrics$errors$err_recid_2_tot) * 100)
    ),
    check.names = FALSE
  )
  
  table_theme <- gridExtra::ttheme_default(
    core = list(fg_params = list(fontsize = 9), bg_params = list(fill = c(rep("honeydew2", 6), "white", rep("lightyellow", 6)))),
    colhead = list(fg_params = list(fontface = "bold", fontsize = 10), bg_params = list(fill = "steelblue", col = "white"))
  )
  
  table_plot <- gridExtra::tableGrob(table_data, rows = NULL, theme = table_theme)
  title_plot <- grid::textGrob("Paramètres calibrés et cibles", gp = grid::gpar(fontsize = 14, fontface = "bold"))
  
  return(gridExtra::arrangeGrob(table_plot, top = title_plot))
}


###############################################################################
# ---- 4 : STEADY STATE CHECK TABLES ----
###############################################################################

create_initial_condition_check_table <- function(initial_condition_check) {
  table_data <- initial_condition_check[, c("initial_condition", "equilibrium_reached", "final_time", "same_equilibrium")]
  table_data$initial_condition <- gsub("^IC_", "CI ", table_data$initial_condition)
  table_data$final_time <- sprintf("%.0f jours", table_data$final_time)
  names(table_data) <- c("Condition initiale", "Équilibre atteint", "Temps final", "Même équilibre")

  table_theme <- gridExtra::ttheme_default(
    core = list(fg_params = list(fontsize = 9), bg_params = list(fill = rep(c("grey95", "white"), length.out = nrow(table_data)))),
    colhead = list(fg_params = list(fontface = "bold", fontsize = 10), bg_params = list(fill = "steelblue", col = "white"))
  )

  table_plot <- gridExtra::tableGrob(table_data, rows = NULL, theme = table_theme)
  title_plot <- grid::textGrob("Vérification de l'état stationnaire des conditions initiales", gp = grid::gpar(fontsize = 14, fontface = "bold"))

  return(gridExtra::arrangeGrob(title_plot, table_plot, ncol = 1, heights = c(0.10, 0.90)))
}

create_perturbation_check_table <- function(perturbation_check) {
  table_data <- perturbation_check[, c("perturbation_compartment", "perturbation_pct", "equilibrium_reached", "same_equilibrium")]
  table_data$perturbation_pct <- sprintf("%+.0f%%", 100 * table_data$perturbation_pct)
  names(table_data) <- c("Compartiment perturbé", "Perturbation", "Équilibre atteint", "Même équilibre")

  table_theme <- gridExtra::ttheme_default(
    core = list(fg_params = list(fontsize = 7.5), bg_params = list(fill = rep(c("grey95", "white"), length.out = nrow(table_data)))),
    colhead = list(fg_params = list(fontface = "bold", fontsize = 8), bg_params = list(fill = "steelblue", col = "white"))
  )

  table_plot <- gridExtra::tableGrob(table_data, rows = NULL, theme = table_theme)
  title_plot <- grid::textGrob("Vérification de l'état stationnaire après perturbation", gp = grid::gpar(fontsize = 14, fontface = "bold"))

  return(gridExtra::arrangeGrob(title_plot, table_plot, ncol = 1, heights = c(0.06, 0.94)))
}


###############################################################################
# ---- 5 : SCENARIO PLOTS ----
###############################################################################

# Prepare the three data frames used by hygiene or antibiotic plots
prepare_intervention_plot_data <- function(simulation_results, intervention = c("HYG", "ATB")) {
  intervention <- match.arg(intervention)
  met <- prepare_intervention_bar_metrics(simulation_results, intervention)

  df_car <- met %>%
    dplyr::select(scenario, car_c_prev, car_h_prev) %>%
    tidyr::pivot_longer(cols = c(car_c_prev, car_h_prev), names_to = "facet", values_to = "value") %>%
    dplyr::mutate(facet = factor(facet, levels = c("car_c_prev", "car_h_prev"), labels = c("Communaute", "Hopital")))

  df_inc <- met %>%
    dplyr::select(scenario, inc_c_total, inc_h_total) %>%
    tidyr::pivot_longer(cols = c(inc_c_total, inc_h_total), names_to = "facet", values_to = "value") %>%
    dplyr::mutate(facet = factor(facet, levels = c("inc_c_total", "inc_h_total"), labels = c("Communaute", "Hopital")))

  df_pr <- met %>%
    dplyr::select(scenario, inc_primo_total, inc_rec_total) %>%
    tidyr::pivot_longer(cols = c(inc_primo_total, inc_rec_total), names_to = "facet", values_to = "value") %>%
    dplyr::mutate(facet = factor(facet, levels = c("inc_primo_total", "inc_rec_total"), labels = c("Primo CDI", "Recidive CDI")))

  return(list(metrics = met, carriage = df_car, incidence = df_inc, primo_rec = df_pr))
}

# Plot bars with a dashed baseline and relative-change labels
plot_relative_bars_reference <- function(df_plot, title_text, y_label, palette_values, legend_labels) {
  df_base <- df_plot %>%
    dplyr::filter(scenario == "baseline") %>%
    dplyr::select(facet, baseline = value)

  df_scen <- df_plot %>%
    dplyr::filter(scenario != "baseline") %>%
    dplyr::left_join(df_base, by = "facet") %>%
    dplyr::mutate(
      rel_change = (value - baseline) / baseline * 100,
      rel_label = sprintf("%.0f%%", rel_change),
      scenario = factor(scenario, levels = c("faible", "moyen", "fort"))
    )

  p <- ggplot2::ggplot(df_scen, ggplot2::aes(x = scenario, y = value, fill = scenario)) +
    ggplot2::geom_col(width = 0.9, color = "black", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = rel_label), fontface = "bold", vjust = 1.3, size = 5) +
    ggplot2::geom_hline(
      data = df_base,
      ggplot2::aes(yintercept = baseline, linetype = "Scenario de reference"),
      color = "red2",
      linewidth = 0.8,
      inherit.aes = FALSE
    ) +
    ggplot2::facet_wrap(~facet, nrow = 1, scales = "free_y") +
    ggplot2::scale_fill_manual(
      values = palette_values,
      breaks = c("faible", "moyen", "fort"),
      labels = legend_labels[c("faible", "moyen", "fort")],
      name = "Scenario"
    ) +
    ggplot2::scale_linetype_manual(values = c("Scenario de reference" = "dashed"), name = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      strip.background = ggplot2::element_rect(fill = "white", color = "black"),
      panel.grid.minor = ggplot2::element_line(color = "#E3E3E3"),
      panel.grid.major = ggplot2::element_line(color = "#D0D0D0"),
      legend.position = "right",
      strip.text = ggplot2::element_text(face = "bold"),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    ) +
    ggplot2::labs(title = title_text, y = y_label)

  return(p)
}

# Create intervention plots for hygiene or antibiotic scenarios
plot_intervention_reference_panels <- function(plot_data, intervention = c("HYG", "ATB")) {
  intervention <- match.arg(intervention)
  met <- plot_data$metrics

  pal <- if (intervention == "ATB") {
    c(faible = "#C7E9C0", moyen = "#74C476", fort = "#238B45")
  } else {
    c(faible = "#D8EAFE", moyen = "#7FB3FF", fort = "#2F80ED")
  }

  lbl <- c(
    faible = sprintf("%s -%d%%", intervention, round(100 * met$reduction[met$scenario == "faible"])),
    moyen = sprintf("%s -%d%%", intervention, round(100 * met$reduction[met$scenario == "moyen"])),
    fort = sprintf("%s -%d%%", intervention, round(100 * met$reduction[met$scenario == "fort"]))
  )

  scenario_title <- if (intervention == "ATB") "Scénario Antibiotique" else "Scénario Hygiène"

  p_car <- plot_relative_bars_reference(plot_data$carriage, paste0(scenario_title, " - Diminution de la prevalence de colonisation"), "Prevalence de colonisation", pal, lbl)
  p_inc <- plot_relative_bars_reference(plot_data$incidence, paste0(scenario_title, " - Diminution de l'incidence CDI"), "Incidence CDI (/100k/an)", pal, lbl)
  p_pr <- plot_relative_bars_reference(plot_data$primo_rec, paste0(scenario_title, " - Diminution de l'incidence primo/recidive CDI"), "Incidence CDI (/100k/an)", pal, lbl)

  p_all <- p_car + p_inc + p_pr +
    patchwork::plot_layout(guides = "collect", ncol = 3) +
    patchwork::plot_annotation(title = paste0(scenario_title, " - Diminution des indicateurs CDI")) &
    ggplot2::theme(legend.position = "right")

  return(list(carriage = p_car, incidence = p_inc, primo_rec = p_pr, combined = p_all, data = met))
}


###############################################################################
# ---- 6 : VACCINATION PLOTS ----
###############################################################################

# Prepare the three data frames used by vaccination plots
prepare_vaccination_plot_data <- function(simulation_results) {
  df <- build_vacc_plot_data(simulation_results)
  ref <- compute_intervention_reference_metrics(simulation_results)

  df_car <- df %>%
    dplyr::select(VC, VE, car_c_rel, car_h_rel) %>%
    tidyr::pivot_longer(cols = c(car_c_rel, car_h_rel), names_to = "facet", values_to = "rel_change") %>%
    dplyr::mutate(facet = factor(facet, levels = c("car_c_rel", "car_h_rel"), labels = c("Communaute", "Hopital")))

  df_inc <- df %>%
    dplyr::select(VC, VE, inc_c_total_rel, inc_h_total_rel) %>%
    tidyr::pivot_longer(cols = c(inc_c_total_rel, inc_h_total_rel), names_to = "facet", values_to = "rel_change") %>%
    dplyr::mutate(facet = factor(facet, levels = c("inc_c_total_rel", "inc_h_total_rel"), labels = c("Communaute", "Hopital")))

  df_pr <- df %>%
    dplyr::select(VC, VE, inc_primo_rel, inc_rec_rel) %>%
    tidyr::pivot_longer(cols = c(inc_primo_rel, inc_rec_rel), names_to = "facet", values_to = "rel_change") %>%
    dplyr::mutate(facet = factor(facet, levels = c("inc_primo_rel", "inc_rec_rel"), labels = c("Primo CDI", "Recidive CDI")))

  ref_car <- ref %>%
    dplyr::select(reference, car_c_rel, car_h_rel) %>%
    tidyr::pivot_longer(cols = c(car_c_rel, car_h_rel), names_to = "facet", values_to = "rel_change") %>%
    dplyr::mutate(facet = factor(facet, levels = c("car_c_rel", "car_h_rel"), labels = c("Communaute", "Hopital")))

  ref_inc <- ref %>%
    dplyr::select(reference, inc_c_total_rel, inc_h_total_rel) %>%
    tidyr::pivot_longer(cols = c(inc_c_total_rel, inc_h_total_rel), names_to = "facet", values_to = "rel_change") %>%
    dplyr::mutate(facet = factor(facet, levels = c("inc_c_total_rel", "inc_h_total_rel"), labels = c("Communaute", "Hopital")))

  ref_pr <- ref %>%
    dplyr::select(reference, inc_primo_rel, inc_rec_rel) %>%
    tidyr::pivot_longer(cols = c(inc_primo_rel, inc_rec_rel), names_to = "facet", values_to = "rel_change") %>%
    dplyr::mutate(facet = factor(facet, levels = c("inc_primo_rel", "inc_rec_rel"), labels = c("Primo CDI", "Recidive CDI")))

  return(list(data = df, carriage = df_car, incidence = df_inc, primo_rec = df_pr,
              ref_carriage = ref_car, ref_incidence = ref_inc, ref_primo_rec = ref_pr))
}

# Generic vaccination line plot
plot_vacc_relative_lines <- function(df_plot, title_text, y_label, reference_lines = NULL) {
  df_plot$series <- sprintf("VE %d%%", round(100 * df_plot$VE))
  ve_labels <- sort(unique(df_plot$series))
  violet_values <- grDevices::colorRampPalette(c("#D8B4FE", "#8B5CF6", "#4C1D95"))(length(ve_labels))
  color_values <- c(stats::setNames(violet_values, ve_labels),
                    "Antibiotiques -20%" = "#238B45",
                    "Hygiene -20%" = "#2F80ED")
  linetype_values <- c(stats::setNames(rep("solid", length(ve_labels)), ve_labels),
                       "Antibiotiques -20%" = "dashed",
                       "Hygiene -20%" = "dashed")

  zero_rows <- unique(df_plot[, c("series", "facet")])
  zero_rows$VC <- 0
  zero_rows$rel_change <- 0
  df_plot <- dplyr::bind_rows(df_plot, zero_rows)
  df_plot <- df_plot[order(df_plot$series, df_plot$facet, df_plot$VC), ]

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = VC * 100, y = rel_change, color = series, linetype = series, group = series)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.7) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_color_manual(values = color_values) +
    ggplot2::scale_linetype_manual(values = linetype_values) +
    ggplot2::scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20), expand = c(0, 0)) +
    ggplot2::facet_wrap(~facet, nrow = 1, scales = "free_y") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      strip.background = ggplot2::element_rect(fill = "white", color = "black"),
      panel.grid.minor = ggplot2::element_line(color = "#E3E3E3"),
      panel.grid.major = ggplot2::element_line(color = "#D0D0D0"),
      strip.text = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::labs(title = title_text, x = "Couverture vaccinale (%)", y = y_label, color = NULL, linetype = NULL)

  if (!is.null(reference_lines)) {
    p <- p +
      ggplot2::geom_hline(
        data = reference_lines,
        ggplot2::aes(yintercept = rel_change, color = reference, linetype = reference),
        linewidth = 1,
        inherit.aes = FALSE
      )
  }

  return(p)
}

# Vaccination impact plots
plot_vacc_impact <- function(plot_data) {
  p_car <- plot_vacc_relative_lines(plot_data$carriage, "Scénario Vaccination - Diminution relative de la prevalence de colonisation", "Changement relatif (%)", plot_data$ref_carriage)
  p_inc <- plot_vacc_relative_lines(plot_data$incidence, "Scénario Vaccination - Diminution relative de l'incidence CDI", "Changement relatif (%)", plot_data$ref_incidence)
  p_pr <- plot_vacc_relative_lines(plot_data$primo_rec, "Scénario Vaccination - Diminution relative de l'incidence primo/recidive CDI", "Changement relatif (%)", plot_data$ref_primo_rec)

  p_all <- (p_car + p_inc + p_pr) +
    patchwork::plot_layout(guides = "collect", ncol = 1) +
    patchwork::plot_annotation(title = "Scénario Vaccination - Diminution relative des indicateurs CDI") &
    ggplot2::theme(legend.position = "right")

  return(list(carriage = p_car, incidence = p_inc, primo_rec = p_pr, combined = p_all, data = plot_data$data))
}


###############################################################################
# ---- 7 : PRCC PLOTS ----
###############################################################################

# Plot PRCC results
plot_prcc_results <- function(prcc_table) {
  output_labels <- c(
    portage_h = "Prévalence de colonisation hospitalière",
    portage_c = "Prévalence de colonisation communautaire",
    inc_h = "Incidence ICD hospitalière",
    inc_c = "Incidence ICD communautaire",
    inc_primo_tot = "Incidence primo-infection ICD",
    inc_rec_tot = "Incidence récidive ICD"
  )

  prcc_table$output <- factor(
    prcc_table$output,
    levels = names(output_labels),
    labels = output_labels
  )

  p <- ggplot2::ggplot(
    prcc_table,
    ggplot2::aes(x = reorder(parameter, prcc), y = prcc, fill = prcc > 0)
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~output, scales = "free_y", ncol = 2) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::scale_fill_manual(values = c("TRUE" = "#2C7FB8", "FALSE" = "#D95F0E"), guide = "none") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Analyse de sensibilité PRCC",
      subtitle = "Chaque paramètre varie de -50 % à +50 %",
      x = "Paramètre",
      y = "PRCC"
    )

  return(p)
}

# Plot PRCC results for selected output pairs
plot_prcc_pairs <- function(prcc_table) {
  output_labels <- c(
    portage_h = "Prévalence de colonisation hospitalière",
    portage_c = "Prévalence de colonisation communautaire",
    inc_h = "Incidence ICD hospitalière",
    inc_c = "Incidence ICD communautaire",
    inc_primo_tot = "Incidence primo-infection ICD",
    inc_rec_tot = "Incidence récidive ICD"
  )

  make_one_plot <- function(outputs_to_keep, title_text) {
    df_plot <- prcc_table[prcc_table$output %in% outputs_to_keep, ]
    df_plot$output <- factor(
      df_plot$output,
      levels = outputs_to_keep,
      labels = output_labels[outputs_to_keep]
    )

    p <- ggplot2::ggplot(
      df_plot,
      ggplot2::aes(x = reorder(parameter, prcc), y = prcc, fill = prcc > 0)
    ) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::facet_wrap(~output, scales = "free_y", ncol = 2) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
      ggplot2::scale_fill_manual(values = c("TRUE" = "#2C7FB8", "FALSE" = "#D95F0E"), guide = "none") +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = title_text,
        subtitle = "Chaque paramètre varie de -50 % à +50 %",
        x = "Paramètre",
        y = "PRCC"
      )

    return(p)
  }

  p_incidence <- make_one_plot(c("inc_c", "inc_h"), "PRCC - Incidence ICD")
  p_primo_rec <- make_one_plot(c("inc_primo_tot", "inc_rec_tot"), "PRCC - Incidence primo-infection et récidive ICD")
  p_carriage <- make_one_plot(c("portage_c", "portage_h"), "PRCC - Prévalence de colonisation")

  return(list(
    incidence = p_incidence,
    primo_rec = p_primo_rec,
    carriage = p_carriage
  ))
}

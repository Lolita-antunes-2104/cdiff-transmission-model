###############################################################################
############################ 3 : SCENARIOS ####################################
###############################################################################

# Compare two simulation outputs.
# If the second output comes from the vaccination model, vaccinated and
# non-vaccinated compartments are added before comparison.
compare_outputs <- function(reference_output, tested_output, comparison_name, tolerance = 1e-10) {
  reference_names <- setdiff(colnames(reference_output), "time")
  n_rows <- min(nrow(reference_output), nrow(tested_output))

  reference_values <- reference_output[seq_len(n_rows), reference_names, drop = FALSE]
  tested_values <- reference_values

  for (compartment in reference_names) {
    if (compartment %in% colnames(tested_output)) {
      tested_values[[compartment]] <- tested_output[[compartment]][seq_len(n_rows)]
    } else {
      tested_values[[compartment]] <- tested_output[[paste0(compartment, "_nv")]][seq_len(n_rows)] +
        tested_output[[paste0(compartment, "_v")]][seq_len(n_rows)]
    }
  }

  max_abs_diff <- max(abs(as.matrix(reference_values) - as.matrix(tested_values)))

  return(data.frame(
    comparison = comparison_name,
    same_outputs = max_abs_diff <= tolerance,
    max_abs_diff = max_abs_diff,
    stringsAsFactors = FALSE
  ))
}


###############################################################################
# ---- SCENARIO PLOT DATA ----
###############################################################################

# Compute plot metrics from one standard model simulation
compute_standard_plot_metrics <- function(one_result) {
  last_state <- as.list(tail(one_result$ode, 1))
  metrics <- compute_all_metrics(last_state, one_result$params)

  out <- data.frame(
    car_h_prev = as.numeric(metrics$portage_h),
    car_c_prev = as.numeric(metrics$portage_c),
    inc_h_total = as.numeric(metrics$inc_h) * 365 * 1e5,
    inc_c_total = as.numeric(metrics$inc_c) * 365 * 1e5,
    inc_primo_total = as.numeric(metrics$inc_primo_tot) * 365 * 1e5,
    inc_rec_total = as.numeric(metrics$inc_rec_tot) * 365 * 1e5
  )

  return(out)
}

# Prepare final metrics for hygiene or antibiotic scenarios
prepare_intervention_bar_metrics <- function(simulation_results, intervention = c("HYG", "ATB")) {
  intervention <- match.arg(intervention)

  scenario_list <- if (intervention == "HYG") {
    simulation_results$hygiene
  } else {
    simulation_results$antibiotics
  }

  rows <- list()

  rows$baseline <- compute_standard_plot_metrics(simulation_results$baseline)
  rows$baseline$scenario <- "baseline"
  rows$baseline$reduction <- 0

  rows$faible <- compute_standard_plot_metrics(scenario_list$low)
  rows$faible$scenario <- "faible"
  rows$faible$reduction <- scenario_list$low$reduction

  rows$moyen <- compute_standard_plot_metrics(scenario_list$medium)
  rows$moyen$scenario <- "moyen"
  rows$moyen$reduction <- scenario_list$medium$reduction

  rows$fort <- compute_standard_plot_metrics(scenario_list$strong)
  rows$fort$scenario <- "fort"
  rows$fort$reduction <- scenario_list$strong$reduction

  out <- dplyr::bind_rows(rows)
  out$scenario <- factor(out$scenario, levels = c("baseline", "faible", "moyen", "fort"))

  return(out)
}

# Compute plot metrics from one vaccination simulation
compute_vaccination_plot_metrics <- function(one_result) {
  st <- as.list(tail(one_result$ode, 1))
  VE <- as.numeric(one_result$VE)

  N_h <- with(st, S0_h_nv + SA_h_nv + S_II_h_nv + S_III_h_nv + C0_h_nv + CA_h_nv + C_II_h_nv + C_III_h_nv + I_h_nv + I_II_h_nv + I_III_h_nv +
                S0_h_v + SA_h_v + S_II_h_v + S_III_h_v + C0_h_v + CA_h_v + C_II_h_v + C_III_h_v + I_h_v + I_II_h_v + I_III_h_v)
  N_c <- with(st, S0_c_nv + SA_c_nv + S_II_c_nv + S_III_c_nv + C0_c_nv + CA_c_nv + C_II_c_nv + C_III_c_nv + I_c_nv + I_II_c_nv + I_III_c_nv +
                S0_c_v + SA_c_v + S_II_c_v + S_III_c_v + C0_c_v + CA_c_v + C_II_c_v + C_III_c_v + I_c_v + I_II_c_v + I_III_c_v)
  N_tot <- N_h + N_c

  C_h <- with(st, C0_h_nv + CA_h_nv + C_II_h_nv + C_III_h_nv + C0_h_v + CA_h_v + C_II_h_v + C_III_h_v)
  C_c <- with(st, C0_c_nv + CA_c_nv + C_II_c_nv + C_III_c_nv + C0_c_v + CA_c_v + C_II_c_v + C_III_c_v)

  sigma_h <- as.numeric(one_result$params["sigma_h"])
  sigma_c <- as.numeric(one_result$params["sigma_c"])
  k_A <- as.numeric(one_result$params["k_A"])
  k_II <- as.numeric(one_result$params["k_II"])
  k_III <- as.numeric(one_result$params["k_III"])

  sigma_h_v <- sigma_h * (1 - VE)
  sigma_c_v <- sigma_c * (1 - VE)

  inc_h_primo <- sigma_h * st$C0_h_nv + k_A * sigma_h * st$CA_h_nv + sigma_h_v * st$C0_h_v + k_A * sigma_h_v * st$CA_h_v
  inc_c_primo <- sigma_c * st$C0_c_nv + k_A * sigma_c * st$CA_c_nv + sigma_c_v * st$C0_c_v + k_A * sigma_c_v * st$CA_c_v
  inc_h_rec <- k_II * sigma_h * st$C_II_h_nv + k_III * sigma_h * st$C_III_h_nv + k_II * sigma_h_v * st$C_II_h_v + k_III * sigma_h_v * st$C_III_h_v
  inc_c_rec <- k_II * sigma_c * st$C_II_c_nv + k_III * sigma_c * st$C_III_c_nv + k_II * sigma_c_v * st$C_II_c_v + k_III * sigma_c_v * st$C_III_c_v

  out <- data.frame(
    car_h_prev = as.numeric(C_h / N_h),
    car_c_prev = as.numeric(C_c / N_c),
    inc_h_total = as.numeric((inc_h_primo + inc_h_rec) / N_tot * 365 * 1e5),
    inc_c_total = as.numeric((inc_c_primo + inc_c_rec) / N_tot * 365 * 1e5),
    inc_primo_total = as.numeric((inc_h_primo + inc_c_primo) / N_tot * 365 * 1e5),
    inc_rec_total = as.numeric((inc_h_rec + inc_c_rec) / N_tot * 365 * 1e5)
  )

  return(out)
}

# Build the complete vaccination data frame
build_vacc_plot_data <- function(simulation_results) {
  base <- compute_standard_plot_metrics(simulation_results$baseline)
  scenario_names <- names(simulation_results$vaccination)

  rows <- lapply(scenario_names, function(sc) {
    one <- simulation_results$vaccination[[sc]]
    out <- compute_vaccination_plot_metrics(one)
    out$scenario <- sc
    out$VE <- as.numeric(one$VE)
    out$VC <- as.numeric(one$VC)
    out
  })

  df <- dplyr::bind_rows(rows)

  df <- df %>%
    dplyr::mutate(
      car_h_rel = (car_h_prev - base$car_h_prev) / base$car_h_prev * 100,
      car_c_rel = (car_c_prev - base$car_c_prev) / base$car_c_prev * 100,
      inc_h_total_rel = (inc_h_total - base$inc_h_total) / base$inc_h_total * 100,
      inc_c_total_rel = (inc_c_total - base$inc_c_total) / base$inc_c_total * 100,
      inc_primo_rel = (inc_primo_total - base$inc_primo_total) / base$inc_primo_total * 100,
      inc_rec_rel = (inc_rec_total - base$inc_rec_total) / base$inc_rec_total * 100
    )

  return(df)
}

# Compute reference lines from the strongest hygiene and antibiotic scenarios
compute_intervention_reference_metrics <- function(simulation_results) {
  base <- compute_standard_plot_metrics(simulation_results$baseline)

  make_reference_row <- function(reference_name, scenario_result) {
    met <- compute_standard_plot_metrics(scenario_result)

    data.frame(
      reference = reference_name,
      car_h_rel = (met$car_h_prev - base$car_h_prev) / base$car_h_prev * 100,
      car_c_rel = (met$car_c_prev - base$car_c_prev) / base$car_c_prev * 100,
      inc_h_total_rel = (met$inc_h_total - base$inc_h_total) / base$inc_h_total * 100,
      inc_c_total_rel = (met$inc_c_total - base$inc_c_total) / base$inc_c_total * 100,
      inc_primo_rel = (met$inc_primo_total - base$inc_primo_total) / base$inc_primo_total * 100,
      inc_rec_rel = (met$inc_rec_total - base$inc_rec_total) / base$inc_rec_total * 100,
      stringsAsFactors = FALSE
    )
  }

  out <- dplyr::bind_rows(
    make_reference_row("Antibiotiques -20%", simulation_results$antibiotics$strong),
    make_reference_row("Hygiene -20%", simulation_results$hygiene$strong)
  )

  return(out)
}

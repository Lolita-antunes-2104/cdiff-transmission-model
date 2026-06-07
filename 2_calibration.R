###############################################################################
#################### 2 : CALIBRATION (GRID + OPTIMIZATION) ####################
###############################################################################

###############################################################################
# ---- 1) GRID SEARCH ----
###############################################################################

# Generic grid search function
grid_search <- function(param_names, param_ranges, target_metrics, 
                        params_base, init_cond, time_vec, n_cores = NULL) { 

  # 1) Create all parameter combinations
  grid_values <- lapply(param_ranges, function(x) seq(x[1], x[2], length.out = x[3]))
  names(grid_values) <- param_names
  param_grid <- do.call(expand.grid, grid_values)

  # 2) Choose number of cores
  if (is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }

  cat(sprintf("Nombre de combinaisons testees : %d\n", nrow(param_grid)))
  cat(sprintf("Nombre de coeurs utilises : %d\n", n_cores))

  # 3) Define what happens for one parameter combination
  one_row <- function(i) {
    param_row <- as.numeric(param_grid[i, ])
    names(param_row) <- names(param_grid)

    params_test <- params_base
    params_test[param_names] <- param_row[param_names]

    # 3.1) Run the model until equilibrium
    res <- run_model_until_equilibrium(
      params_vec = params_test,
      init_cond = init_cond, 
      time_max = max(time_vec),
      by = 1,
      equilibrium_tol = 1e-6,
      min_time_before_check = 365,
      chunk_length = 365,
      model_function = cdiff_model
    )

    last_state <- as.list(res[nrow(res), ])
    metrics <- c()

    # 3.2) Compute the target metrics for this grid search
    # beta grid search
    if ("portage_h" %in% names(target_metrics)) { 
      metrics["portage_h"] <- unname(compute_carriage_prevalence(last_state, "h"))
      metrics["portage_c"] <- unname(compute_carriage_prevalence(last_state, "c"))
    }

    # sigma grid search
    if ("incidence_h" %in% names(target_metrics)) { 
      N_tot <- compute_population_totals(last_state, "both")
      metrics["incidence_h"] <- unname(compute_CDI_incidence(last_state, params_test, "h", "total")) / N_tot
      metrics["incidence_c"] <- unname(compute_CDI_incidence(last_state, params_test, "c", "total")) / N_tot
    }

    # k grid search
    if ("recid_1" %in% names(target_metrics)) { 
      metrics["recid_1"] <- unname(compute_recurrence_prevalence(last_state, "both", "rec_1"))
      metrics["recid_2"] <- unname(compute_recurrence_prevalence(last_state, "both", "rec_2"))
    }

    # 3.3) Add equilibrium diagnostics
    metrics["equilibrium_reached"] <- as.numeric(attr(res, "equilibrium_reached"))
    metrics["max_abs_dxdt"] <- as.numeric(attr(res, "equilibrium_max_abs_dydt"))

    return(metrics)
  }

  # 4) Run all parameter combinations
  metric_list <- pbapply::pblapply(seq_len(nrow(param_grid)), one_row, cl = n_cores)
  metric_table <- as.data.frame(do.call(rbind, metric_list))

  # 5) Combine parameters, metrics, and diagnostics
  grid_result <- cbind(param_grid, metric_table)
  grid_result$equilibrium_reached <- as.logical(grid_result$equilibrium_reached)

  # 6) Compute distance to calibration targets
  distance <- rep(0, nrow(grid_result))
  for (metric_name in names(target_metrics)) {
    distance <- distance + (grid_result[[metric_name]] - target_metrics[[metric_name]])^2
  }

  grid_result$distance <- sqrt(distance)
  best_guess <- grid_result[which.min(grid_result$distance), ]

  # 7) Prepare equilibrium diagnostic summary
  n_total <- nrow(grid_result)
  n_reached <- sum(grid_result$equilibrium_reached, na.rm = TRUE)
  prop_reached <- n_reached / n_total
  if (all(is.na(grid_result$max_abs_dxdt))) {
    worst_max_abs_dxdt <- NA_real_
  } else {
    worst_max_abs_dxdt <- max(grid_result$max_abs_dxdt, na.rm = TRUE)
  }

  # 8) Return everything useful
  return(list(
    grid = grid_result,
    best_guess = best_guess,
    diagnostic = list(
      label = paste(param_names, collapse = " / "),
      n_total = n_total,
      n_reached = n_reached,
      prop_reached = prop_reached,
      best_max_abs_dxdt = best_guess$max_abs_dxdt,
      worst_max_abs_dxdt = worst_max_abs_dxdt
    )
  ))
}






###############################################################################
# ---- 2) MULTI-START OPTIMISATION ----
###############################################################################

run_optimization <- function(initial_params, target_metrics, params_base,
                            init_cond, time_vec, n_starts = 10,
                            n_cores = NULL, seed = 123) {

  # 1) Choose number of cores
  if (is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }

  # 2) Create starting points
  set.seed(seed)
  initial_log <- log(initial_params)
  start_points <- lapply(seq_len(n_starts), function(i) {
    if (i == 1) return(initial_log)
    return(initial_log + runif(length(initial_log), -0.5, 0.5))
  })

  cat(sprintf("\nNombre de demarrages : %d\n", n_starts))
  cat(sprintf("Nombre de coeurs utilises : %d\n", n_cores))

  # 3) Define one optimization start
  one_start <- function(start_log) {

    opt <- optim(
      par = start_log,
      method = "Nelder-Mead",
      control = list(maxit = 10000, reltol = 1e-8),
      fn = function(par_log) {

        params_test <- params_base
        tested_values <- exp(par_log)
        names(tested_values) <- names(initial_params)
        params_test[names(tested_values)] <- tested_values

        tryCatch({
          res <- run_model_until_equilibrium(
            params_vec = params_test,
            init_cond = init_cond,
            time_max = max(time_vec),
            by = 1,
            equilibrium_tol = 1e-6,
            min_time_before_check = 365,
            chunk_length = 365,
            model_function = cdiff_model
          )

          last_state <- as.list(res[nrow(res), ])
          if (any(is.na(unlist(last_state)))) return(1e10)
          if (any(unlist(last_state) < 0)) return(1e10)

          metrics <- compute_all_metrics(last_state, params_test, targets = target_metrics)
          return(sum(unlist(metrics$errors)^2))

        }, error = function(e) {
          return(1e10)
        })
      }
    )

    par_natural <- exp(opt$par)
    names(par_natural) <- names(initial_params)

    return(list(
      par_log = opt$par,
      par_natural = par_natural,
      value = opt$value,
      convergence = opt$convergence
    ))
  }

  # 4) Run all starts in parallel
  all_results <- pbapply::pblapply(start_points, one_start, cl = n_cores)

  # 5) Keep the best optimization result
  best_index <- which.min(sapply(all_results, function(x) x$value))
  calibrated_params <- all_results[[best_index]]$par_natural
  params_final <- params_base
  params_final[names(calibrated_params)] <- calibrated_params

  # 6) Run final calibrated model
  final_res <- run_model_until_equilibrium(
    params_vec = params_final,
    init_cond = init_cond,
    time_max = max(time_vec),
    by = 1,
    equilibrium_tol = 1e-6,
    min_time_before_check = 365,
    chunk_length = 365,
    model_function = cdiff_model
  )

  # 7) Compute final metrics and alpha
  final_state <- as.list(final_res[nrow(final_res), ])
  final_metrics <- compute_all_metrics(final_state, params_final, targets = target_metrics)
  tot_h <- compute_totals(final_state$S0_h, final_state$SA_h, final_state$S_II_h, final_state$S_III_h, final_state$C0_h, final_state$CA_h, final_state$C_II_h, final_state$C_III_h, final_state$I_h, final_state$I_II_h, final_state$I_III_h)
  tot_c <- compute_totals(final_state$S0_c, final_state$SA_c, final_state$S_II_c, final_state$S_III_c, final_state$C0_c, final_state$CA_c, final_state$C_II_c, final_state$C_III_c, final_state$I_c, final_state$I_II_c, final_state$I_III_c)
  alpha_rates <- compute_alpha_dynamic(tot_h, tot_c, final_state$I_c, final_state$I_II_c, final_state$I_III_c, params_final["delta"], params_final["w"], params_final["w_I"], params_final["w_II"], params_final["w_III"])

  # 8) Return everything useful
  return(list(
    best = all_results[[best_index]],
    all_results = all_results,
    params_final = params_final,
    metrics = final_metrics,
    ode_result = final_res,
    alpha_eq = alpha_rates$alpha
  ))
}






###############################################################################
# ---- 3) STATIONARY STATE TESTS ----
###############################################################################

# Define the model compartment names used in the tests
get_stationary_state_names <- function() {
  hospital_names <- c("S0_h", "SA_h", "C0_h", "CA_h", "I_h", "S_II_h", "C_II_h", "I_II_h", "S_III_h", "C_III_h", "I_III_h")
  community_names <- c("S0_c", "SA_c", "C0_c", "CA_c", "I_c", "S_II_c", "C_II_c", "I_II_c", "S_III_c", "C_III_c", "I_III_c")
  state_names <- c(hospital_names, community_names)

  return(list(
    hospital_names = hospital_names,
    community_names = community_names,
    state_names = state_names
  ))
}

# Check if the final calibrated run is already at equilibrium
check_final_equilibrium <- function(params_final, final_res, state_names,
                                    equilibrium_tol = 1e-6,
                                    model_function = cdiff_model) {

  cat("\n--- 1) Verification de l'etat final calibre ---\n")

  reference_state <- as.numeric(final_res[nrow(final_res), state_names])
  names(reference_state) <- state_names

  dx_reference <- model_function(0, reference_state, params_final)[[1]]
  max_abs_dxdt_reference <- max(abs(dx_reference))

  equilibrium_check <- data.frame(
    equilibrium_reached = max_abs_dxdt_reference <= equilibrium_tol,
    final_time = max(final_res$time),
    max_abs_dxdt = max_abs_dxdt_reference,
    stringsAsFactors = FALSE
  )

  cat(sprintf("Equilibrium reached: %s\n", equilibrium_check$equilibrium_reached))
  cat(sprintf("Final time: %.0f days\n", equilibrium_check$final_time))
  cat(sprintf("Final max|dX/dt|: %.3e\n", equilibrium_check$max_abs_dxdt))

  return(list(
    equilibrium_check = equilibrium_check,
    reference_state = reference_state
  ))
}

# Test if different initial conditions return to the same equilibrium
check_initial_conditions <- function(params_final, initial_conditions_list, time_vec, reference_state,
                                     state_names,
                                     equilibrium_tol = 1e-6,
                                     distance_tol = 1e-6,
                                     model_function = cdiff_model) {

  cat("\n--- 2) Test avec plusieurs conditions initiales ---\n")

  initial_condition_check <- data.frame()

  for (i in seq_along(initial_conditions_list)) {
    init_name <- names(initial_conditions_list)[i]
    init_test <- initial_conditions_list[[i]]

    cat(sprintf("%s (%d / %d)\n", init_name, i, length(initial_conditions_list)))

    result_test <- run_model_until_equilibrium(
      params_vec = params_final,
      init_cond = init_test,
      time_max = max(time_vec),
      by = 1,
      equilibrium_tol = equilibrium_tol,
      min_time_before_check = 365,
      chunk_length = 365,
      model_function = model_function
    )

    end_state <- as.numeric(result_test[nrow(result_test), state_names])
    names(end_state) <- state_names
    distance_to_reference <- sum(abs(end_state - reference_state)) / sum(abs(reference_state))

    one_row <- data.frame(
      initial_condition = init_name,
      equilibrium_reached = attr(result_test, "equilibrium_reached"),
      final_time = max(result_test$time),
      same_equilibrium = distance_to_reference <= distance_tol,
      stringsAsFactors = FALSE
    )

    initial_condition_check <- rbind(initial_condition_check, one_row)
  }

  cat(sprintf(
    "Same equilibrium from different initial conditions: %d / %d\n",
    sum(initial_condition_check$same_equilibrium),
    nrow(initial_condition_check)
  ))

  return(initial_condition_check)
}

# Test if selected perturbations return to the same equilibrium
check_perturbations <- function(params_final, perturbation_table, time_vec, reference_state,
                                hospital_names, community_names, state_names, N_h, N_c,
                                equilibrium_tol = 1e-6,
                                distance_tol = 1e-6,
                                model_function = cdiff_model) {

  cat("\n--- 3) Test avec perturbations choisies ---\n")

  perturbation_check <- data.frame()

  for (i in seq_len(nrow(perturbation_table))) {

    compartment_name <- perturbation_table$perturbation_compartment[i]
    one_perturbation_pct <- perturbation_table$perturbation_pct[i]
    cat(sprintf("Perturbation %s (%+.0f%%)\n", compartment_name, 100 * one_perturbation_pct))

    perturbed_init <- reference_state
    perturbed_init[compartment_name] <- perturbed_init[compartment_name] * (1 + one_perturbation_pct)

    if (compartment_name %in% hospital_names) {
      perturbed_init[hospital_names] <- perturbed_init[hospital_names] / sum(perturbed_init[hospital_names]) * N_h
    } else {
      perturbed_init[community_names] <- perturbed_init[community_names] / sum(perturbed_init[community_names]) * N_c
    }

    result_perturbation <- run_model_until_equilibrium(
      params_vec = params_final,
      init_cond = perturbed_init,
      time_max = max(time_vec),
      by = 1,
      equilibrium_tol = equilibrium_tol,
      min_time_before_check = 365,
      chunk_length = 365,
      model_function = model_function
    )

    end_perturbation <- as.numeric(result_perturbation[nrow(result_perturbation), state_names])
    names(end_perturbation) <- state_names
    distance_perturbation <- sum(abs(end_perturbation - reference_state)) / sum(abs(reference_state))

    one_row <- data.frame(
      perturbation_compartment = compartment_name,
      perturbation_pct = one_perturbation_pct,
      equilibrium_reached = attr(result_perturbation, "equilibrium_reached"),
      same_equilibrium = distance_perturbation <= distance_tol,
      stringsAsFactors = FALSE
    )

    perturbation_check <- rbind(perturbation_check, one_row)
  }

  cat(sprintf(
    "Same equilibrium after perturbations: %d / %d\n",
    sum(perturbation_check$same_equilibrium),
    nrow(perturbation_check)
  ))

  return(perturbation_check)
}

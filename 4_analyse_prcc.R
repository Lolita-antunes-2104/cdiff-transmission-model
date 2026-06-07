###############################################################################
######################## 4 : PRCC SENSITIVITY #################################
###############################################################################

# Steps:
# 1) Change each parameter randomly between -50% and +50%.
# 2) Run the model for each parameter set.
# 3) Compute model outputs.
# 4) Compute PRCC with sensitivity::pcc().
run_prcc_sensitivity <- function(params_final, init_cond, time_vec,
                                 n_samples = 100,
                                 seed = 123,
                                 variation = 0.50,
                                 equilibrium_tol = 1e-6,
                                 model_function = cdiff_model) {

  cat("\n=== PRCC SENSITIVITY ANALYSIS ===\n")

  ########### 1) Prepare initial condition ###########
  if (is.data.frame(init_cond)) {
    init_cond <- unlist(init_cond[1, names(init_cond) != "time"])
  }

  if ("time" %in% names(init_cond)) {
    init_cond <- init_cond[names(init_cond) != "time"]
  }

  init_names <- names(init_cond)
  init_cond <- as.numeric(init_cond)
  names(init_cond) <- init_names

  ########### 2) Choose the parameters tested ###########
  param_names <- c(
    "beta_h", "beta_c",
    "sigma_h", "sigma_c",
    "tau_h", "tau_c",
    "omega",
    "nu", "gamma", "epsilon", "p", "phi",
    "k_A", "k_II", "k_III",
    "w", "w_I", "w_II", "w_III",
    "delta",
    "alpha_const"
  )

  param_names <- param_names[param_names %in% names(params_final)]
  param_names <- param_names[as.numeric(params_final[param_names]) != 0]

  ########### 3) Create random parameter values ###########
  set.seed(seed)

  samples <- data.frame(matrix(NA_real_, nrow = n_samples, ncol = length(param_names)))
  colnames(samples) <- param_names

  for (param_name in param_names) {
    base_value <- as.numeric(params_final[param_name])
    min_value <- base_value * (1 - variation)
    max_value <- base_value * (1 + variation)
    samples[[param_name]] <- runif(n_samples, min = min_value, max = max_value)
  }

  ########### 4) Run the model and compute outputs ###########
  output_names <- c(
    "portage_h", "portage_c",
    "inc_h", "inc_c",
    "inc_primo_tot", "inc_rec_tot"
  )

  outputs <- data.frame(matrix(NA_real_, nrow = n_samples, ncol = length(output_names)))
  colnames(outputs) <- output_names

  simulation_status <- data.frame(
    simulation = seq_len(n_samples),
    success = FALSE,
    equilibrium_reached = NA,
    final_time = NA_real_
  )

  for (i in seq_len(n_samples)) {
    cat(sprintf("PRCC simulation %d / %d\n", i, n_samples))

    params_test <- params_final
    params_test[param_names] <- as.numeric(samples[i, param_names])

    res_test <- run_model_until_equilibrium(
      params_vec = params_test,
      init_cond = init_cond,
      time_max = max(time_vec),
      by = 365,
      equilibrium_tol = equilibrium_tol,
      min_time_before_check = 365,
      chunk_length = 365,
      model_function = model_function
    )

    last_state <- as.list(res_test[nrow(res_test), ])
    metrics <- compute_all_metrics(last_state, params_test)

    for (output_name in output_names) {
      outputs[i, output_name] <- as.numeric(metrics[[output_name]])
    }

    simulation_status$success[i] <- all(is.finite(as.numeric(outputs[i, ])))
    simulation_status$equilibrium_reached[i] <- attr(res_test, "equilibrium_reached")
    simulation_status$final_time[i] <- max(res_test$time)
  }

  cat(sprintf("Successful PRCC simulations: %d / %d\n", sum(simulation_status$success), n_samples))

  ########### 5) Compute PRCC with sensitivity::pcc() ###########
  prcc_rows <- list()
  pcc_results <- list()

  for (output_name in output_names) {
    ok <- complete.cases(samples) & is.finite(outputs[[output_name]])

    X <- samples[ok, , drop = FALSE]
    y <- outputs[[output_name]][ok]

    if (length(y) > length(param_names) + 2) {
      pcc_results[[output_name]] <- sensitivity::pcc(X = X, y = y, rank = TRUE)
      one_prcc <- pcc_results[[output_name]]$PRCC

      prcc_rows[[output_name]] <- data.frame(
        output = output_name,
        parameter = rownames(one_prcc),
        prcc = as.numeric(one_prcc[, "original"]),
        stringsAsFactors = FALSE
      )
    } else {
      prcc_rows[[output_name]] <- data.frame(
        output = output_name,
        parameter = param_names,
        prcc = NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }

  prcc_table <- dplyr::bind_rows(prcc_rows)
  prcc_table$abs_prcc <- abs(prcc_table$prcc)

  ########### 6) Return results ###########
  return(list(
    samples = samples,
    outputs = outputs,
    simulation_status = simulation_status,
    pcc_results = pcc_results,
    prcc = prcc_table
  ))
}

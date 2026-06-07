###############################################################################
#################################### MAIN ######################################
###############################################################################

# -----------------------------------------------------------------------------
# SCRIPT GOAL
# -----------------------------------------------------------------------------
# This script runs the full pipeline, in this order:
# 1) Define parameters of the model and targets
# 2) Check model diagnostics
# 3) Run 3 successive grid searches
# 4) Run multi-start optimization
# 5) Run model with calibrated parameters
# 6) Check equilibrium state before simulations
# 7) Run simulations:
#    - preparation for simulations
#    - precheck with neutral scenarios
#    - hospital hygiene, vaccination, and community antibiotic reduction
# 8) Run PRCC sensitivity analysis

###############################################################################
# ---- 0) SOURCE FILES ----
###############################################################################

source("0_packages.R")
source("1_modele.R")
source("2_calibration.R")
source("3_scenarios.R")
source("4_analyse_prcc.R")
source("plot.R")

# Create the results folder (if it does not already exist) to save plots, tables, and output files.
dir.create("résultats", showWarnings = FALSE, recursive = TRUE)


###############################################################################
# ---- 1) BASE PARAMETERS / CI / TIME GRID / TARGET ----
###############################################################################
# In variable names, "_h" means hospital and "_c" means community.

# Population sizes
N_h <- 100 # hospital
N_c <- 50000 # community

# Fixed parameters
nu <- 15 
gamma <- 0.013 
epsilon <- 0.07 
p <- 0.5
phi <- 0.018
k_A <- 6.67 
omega <- 0.024 
delta <- 1 / 7.5

# Weights for dynamic alpha
w <- 1
w_I <- 1.5
w_II <- 2
w_III <- 2.5

# Initial values before calibration
beta_h <- 0.06
beta_c <- 0.02
sigma_h <- 0.003
sigma_c <- 0.001
k_II <- 2
k_III <- 3


# Antibiotic exposure rates
tau_h <- -log(1 - 0.355 * 0.45) / (1 / delta)
tau_c <- -log(1 - 0.0188 * 0.35) / 1

# Full parameter vector
params <- c(
  beta_h = beta_h, beta_c = beta_c,
  sigma_h = sigma_h, sigma_c = sigma_c,
  tau_h = tau_h, tau_c = tau_c,
  omega = omega,
  nu = nu, gamma = gamma,
  epsilon = epsilon, p = p,
  phi = phi,
  k_A = k_A, k_II = k_II, k_III = k_III,
  w = w, w_I = w_I, w_II = w_II, w_III = w_III,
  delta = delta
)

# Initial conditions before calibration (arbitrary)
# 1% in each listed compartment; the rest is placed in S0.
init_prop <- c(SA_h = 0.01, C0_h = 0.01, CA_h = 0.01, I_h = 0.01, S_II_h = 0.01, C_II_h = 0.01, I_II_h = 0.01, S_III_h = 0.01, C_III_h = 0.01, I_III_h = 0.01)
init_h <- c(S0_h = 1 - sum(init_prop), init_prop) * N_h

init_prop <- c(SA_c = 0.01, C0_c = 0.01, CA_c = 0.01, I_c = 0.01,S_II_c = 0.01, C_II_c = 0.01, I_II_c = 0.01, S_III_c = 0.01, C_III_c = 0.01, I_III_c = 0.01)
init_c <- c(S0_c = 1 - sum(init_prop), init_prop) * N_c

init_cond <- c(init_h, init_c)

# Time grid: simulations can run up to 100,000 days, with one output per day.
time_vec <- seq(0, 100000, by = 1)

# Calibration targets
targets <- list(
  portage_h = 0.081, 
  portage_c = 0.015, 
  incidence_h = 14.8 / (100000 * 365), 
  incidence_c = 18.5 / (100000 * 365), 
  recid_1 = 0.25,
  recid_2 = 0.50
)





###############################################################################
# ---- 2) BASE MODEL DIAGNOSTICS ----
###############################################################################

cat("\n=== BASE MODEL DIAGNOSTICS ===\n")

# 1) Run model until equilibrium
res_eq_base <- run_model_until_equilibrium(
  params_vec = params,
  init_cond = init_cond,
  time_max = max(time_vec),
  by = 1,
  equilibrium_tol = 1e-6,
  min_time_before_check = 365,
  chunk_length = 365,
  model_function = cdiff_model
)

cat(sprintf("Equilibrium reached: %s\n", attr(res_eq_base, "equilibrium_reached")))
cat(sprintf("Final time: %.0f days\n", max(res_eq_base$time)))
cat(sprintf("Final max|dX/dt|: %.3e\n", attr(res_eq_base, "equilibrium_max_abs_dydt")))

# 2) Check population is constant (conservation)
res_eq_base$N_h <- with(res_eq_base, S0_h + SA_h + S_II_h + S_III_h + C0_h + CA_h + C_II_h + C_III_h + I_h + I_II_h + I_III_h)
res_eq_base$N_c <- with(res_eq_base, S0_c + SA_c + S_II_c + S_III_c + C0_c + CA_c + C_II_c + C_III_c + I_c + I_II_c + I_III_c)
max_N_h_diff <- max(abs(res_eq_base$N_h - N_h))
max_N_c_diff <- max(abs(res_eq_base$N_c - N_c))

cat(sprintf("Max |N_h - N_h0|: %.3e\n", max_N_h_diff))
cat(sprintf("Max |N_c - N_c0|: %.3e\n", max_N_c_diff))
cat(sprintf("Population constant: %s\n", max_N_h_diff <= 1e-6 & max_N_c_diff <= 1e-6))

# 3) Check alpha is constant over the last year
alpha_base <- sapply(seq_len(nrow(res_eq_base)), function(i) {
  state_i <- as.list(res_eq_base[i, ])
  tot_h <- compute_totals(state_i$S0_h, state_i$SA_h, state_i$S_II_h, state_i$S_III_h, state_i$C0_h, state_i$CA_h, state_i$C_II_h, state_i$C_III_h, state_i$I_h, state_i$I_II_h, state_i$I_III_h)
  tot_c <- compute_totals(state_i$S0_c, state_i$SA_c, state_i$S_II_c, state_i$S_III_c, state_i$C0_c, state_i$CA_c, state_i$C_II_c, state_i$C_III_c, state_i$I_c, state_i$I_II_c, state_i$I_III_c)
  
  compute_alpha_dynamic(tot_h, tot_c, state_i$I_c, state_i$I_II_c, state_i$I_III_c,
                        delta, w, w_I, w_II, w_III)$alpha
})

alpha_last_year <- tail(alpha_base, min(365, length(alpha_base)))
max_alpha_diff <- max(abs(alpha_last_year - tail(alpha_last_year, 1)))

cat(sprintf("Max |alpha - alpha_final| over last year: %.3e\n", max_alpha_diff))
cat(sprintf("Alpha constant over last year: %s\n", max_alpha_diff <= 1e-6))

# 4) Plot dynamics
diag_plots_base <- plot_dynamics(res_eq_base, params, targets, N_h, N_c)
print(diag_plots_base$hospital_all)
print(diag_plots_base$community_all)
print(diag_plots_base$hospital_totals)
print(diag_plots_base$community_totals)





###############################################################################
# ---- 3) GRID SEARCH (3 STEPS) ----
###############################################################################

cat("\n=== RECHERCHE EN GRILLE ===\n")

########### 1) beta_h / beta_c via carriage ########### 
# Define the tested values for beta_h and beta_c
beta_ranges <- list(
  beta_h = c(0.001, 0.5, 30),
  beta_c = c(0.001, 0.5, 30)
)

# Run the grid search:
# for each beta_h / beta_c combination, the model is run until equilibrium (run_until_equilibrium),
# then hospital and community carriage are compared with the targets.
beta_res <- grid_search(
  param_names = c("beta_h", "beta_c"),
  param_ranges = beta_ranges,
  target_metrics = list(portage_h = targets$portage_h, portage_c = targets$portage_c),
  params_base = params, # sigma and k are fixed arbitrary
  init_cond = init_cond,
  time_vec = time_vec,
  n_cores = NULL
)

# Equilibrium diagnostic
cat(sprintf("\n--- Diagnostic etat stationnaire (%s) ---\n", beta_res$diagnostic$label))
cat(sprintf("Equilibre atteint : %d / %d simulations (%.1f%%)\n", beta_res$diagnostic$n_reached, beta_res$diagnostic$n_total, 100 * beta_res$diagnostic$prop_reached))
cat(sprintf("Meilleur point retenu : max|dX/dt| = %.3e\n", beta_res$diagnostic$best_max_abs_dxdt))
cat(sprintf("Pire point testé : max|dX/dt| = %.3e\n", beta_res$diagnostic$worst_max_abs_dxdt))

# Update the parameter vector with the best beta values found (for the next grid search)
params["beta_h"] <- as.numeric(beta_res$best_guess$beta_h)
params["beta_c"] <- as.numeric(beta_res$best_guess$beta_c)

# Grid-search result
cat("\n--- MEILLEUR beta ---\n")
print(beta_res$best_guess)

p_grid_beta <- plot_grid_search_beta(beta_res, beta_res$best_guess, targets)
print(p_grid_beta)
ggplot2::ggsave("résultats/grid_search_beta.png", plot = p_grid_beta, width = 8, height = 6, dpi = 300)

p_grid_beta_points <- plot_grid_search_all_points(beta_res, "beta_h", "beta_c", "Recherche en grille (beta) : valeurs testées")
print(p_grid_beta_points)
ggplot2::ggsave("résultats/grid_search_beta_all_points.png", plot = p_grid_beta_points, width = 7, height = 6, dpi = 300)

p_grid_beta_combined <- p_grid_beta_points + p_grid_beta + patchwork::plot_layout(ncol = 2)
print(p_grid_beta_combined)
ggplot2::ggsave("résultats/grid_search_beta_combined.png", plot = p_grid_beta_combined, width = 14, height = 6, dpi = 300)




########### 2) sigma_h / sigma_c via incidence ########### 
# Define the tested values for sigma_h and sigma_c
sigma_ranges <- list(
  sigma_h = c(0.000001, 0.01, 30),
  sigma_c = c(0.000001, 0.01, 30)
)

# Run the grid search
sigma_res <- grid_search(
  param_names = c("sigma_h", "sigma_c"),
  param_ranges = sigma_ranges,
  target_metrics = list(incidence_h = targets$incidence_h, incidence_c = targets$incidence_c),
  params_base = params, # beta come from the previous grid search; and k is arbitrary
  init_cond = init_cond,
  time_vec = time_vec,
  n_cores = NULL
)

# Equilibrium diagnostic
cat(sprintf("\n--- Diagnostic etat stationnaire (%s) ---\n", sigma_res$diagnostic$label))
cat(sprintf("Equilibre atteint : %d / %d simulations (%.1f%%)\n", sigma_res$diagnostic$n_reached, sigma_res$diagnostic$n_total, 100 * sigma_res$diagnostic$prop_reached))
cat(sprintf("Meilleur point retenu : max|dX/dt| = %.3e\n", sigma_res$diagnostic$best_max_abs_dxdt))
cat(sprintf("Pire point testé : max|dX/dt| = %.3e\n", sigma_res$diagnostic$worst_max_abs_dxdt))

# Update the parameter vector with the best sigma values found (for the next grid search)
params["sigma_h"] <- as.numeric(sigma_res$best_guess$sigma_h)
params["sigma_c"] <- as.numeric(sigma_res$best_guess$sigma_c)

# Grid-search result
cat("\n--- MEILLEUR sigma ---\n")
print(sigma_res$best_guess)

p_grid_sigma <- plot_grid_search_sigma(sigma_res, sigma_res$best_guess, targets)
print(p_grid_sigma)
ggplot2::ggsave("résultats/grid_search_sigma.png", plot = p_grid_sigma, width = 8, height = 6, dpi = 300)

p_grid_sigma_points <- plot_grid_search_all_points(sigma_res, "sigma_h", "sigma_c", "Recherche en grille (sigma) : valeurs testées")
print(p_grid_sigma_points)
ggplot2::ggsave("résultats/grid_search_sigma_all_points.png", plot = p_grid_sigma_points, width = 7, height = 6, dpi = 300)

p_grid_sigma_combined <- p_grid_sigma_points + p_grid_sigma + patchwork::plot_layout(ncol = 2)
print(p_grid_sigma_combined)
ggplot2::ggsave("résultats/grid_search_sigma_combined.png", plot = p_grid_sigma_combined, width = 14, height = 6, dpi = 300)




########### 3) k_II / k_III via recurrences ########### 
# Define the tested values for k_II and k_III
k_ranges <- list(
  k_II = c(1, 1000, 30),
  k_III = c(1, 1000, 30)
)

# Run the grid search
k_res <- grid_search(
  param_names = c("k_II", "k_III"),
  param_ranges = k_ranges,
  target_metrics = list(recid_1 = targets$recid_1, recid_2 = targets$recid_2),
  params_base = params, # beta and sigma come from the previous grid search;
  init_cond = init_cond,
  time_vec = time_vec,
  n_cores = NULL
)

# Equilibrium diagnostic
cat(sprintf("\n--- Diagnostic etat stationnaire (%s) ---\n", k_res$diagnostic$label))
cat(sprintf("Equilibre atteint : %d / %d simulations (%.1f%%)\n", k_res$diagnostic$n_reached, k_res$diagnostic$n_total, 100 * k_res$diagnostic$prop_reached))
cat(sprintf("Meilleur point retenu : max|dX/dt| = %.3e\n", k_res$diagnostic$best_max_abs_dxdt))
cat(sprintf("Pire point testé : max|dX/dt| = %.3e\n", k_res$diagnostic$worst_max_abs_dxdt))

# Update the parameter vector with the best sigma values found
params["k_II"] <- as.numeric(k_res$best_guess$k_II)
params["k_III"] <- as.numeric(k_res$best_guess$k_III)

# Grid-search result
cat("\n--- MEILLEUR k ---\n")
print(k_res$best_guess)

p_grid_k <- plot_grid_search_k(k_res, k_res$best_guess, targets)
print(p_grid_k)
ggplot2::ggsave("résultats/grid_search_k.png", plot = p_grid_k, width = 8, height = 6, dpi = 300)

p_grid_k_points <- plot_grid_search_all_points(k_res, "k_II", "k_III", "Recherche en grille (k) : valeurs testées")
print(p_grid_k_points)
ggplot2::ggsave("résultats/grid_search_k_all_points.png", plot = p_grid_k_points, width = 7, height = 6, dpi = 300)

p_grid_k_combined <- p_grid_k_points + p_grid_k + patchwork::plot_layout(ncol = 2)
print(p_grid_k_combined)
ggplot2::ggsave("résultats/grid_search_k_combined.png", plot = p_grid_k_combined, width = 14, height = 6, dpi = 300)







###############################################################################
# ---- 4) MULTI-START OPTIMIZATION ----
###############################################################################

cat("\n=== MULTI-START OPTIMIZATION ===\n")

# 1) Initial values for the optimizer
# params contains all model parameters; here we keep only the 6 parameters optimized by calibration.
initial_params_optim <- c(
  beta_h = unname(params["beta_h"]),
  beta_c = unname(params["beta_c"]),
  sigma_h = unname(params["sigma_h"]),
  sigma_c = unname(params["sigma_c"]),
  k_II = unname(params["k_II"]),
  k_III = unname(params["k_III"])
)

cat("\n=== Attention : l'optimisation peut prendre plusieurs heures ===\n")

# 2) Run multi-start optimization
optimization_results <- run_optimization(
  initial_params = initial_params_optim,
  target_metrics = targets,
  params_base = params,
  init_cond = init_cond,
  time_vec = time_vec,
  n_starts = 50,
  n_cores = NULL
)

# 3) Optimization diagnostics
cat("\n=== DIAGNOSTICS DE L'OPTIMISATION ===\n")
cat("Nombre de lancements :", length(optimization_results$all_results), "\n")
cat("Valeurs de la fonction objectif :\n")
print(sapply(optimization_results$all_results, function(x) x$value))
cat("Codes de convergence :\n")
print(sapply(optimization_results$all_results, function(x) x$convergence))
cat("Meilleure valeur de la fonction objectif :", optimization_results$best$value, "\n")

# 4) Keep the best calibrated parameters
calibrated_params <- optimization_results$best$par_natural
cat("\n=== PARAMETRES CALIBRES ===\n")
print(calibrated_params)








###############################################################################
# ---- 5) RUN WITH CALIBRATED PARAMETERS ----
###############################################################################

# 1) Build the final parameter vector
final_params <- params
final_params[names(calibrated_params)] <- calibrated_params

# 2) Run the final model with calibrated parameters
final_res <- run_model_until_equilibrium(
  params_vec = final_params,
  init_cond = init_cond,
  time_max = max(time_vec),
  by = 1,
  equilibrium_tol = 1e-6,
  min_time_before_check = 365,
  chunk_length = 365,
  model_function = cdiff_model
)

# 3) Keep the final state, metrics, and alpha
final_state_eq <- tail(final_res, 1)
y_eq <- unlist(final_state_eq[1, -1], use.names = TRUE)
final_state <- as.list(final_state_eq[1, ])
final_metrics <- compute_all_metrics(final_state, final_params, targets = targets)

tot_h <- compute_totals(final_state$S0_h, final_state$SA_h, final_state$S_II_h, final_state$S_III_h, final_state$C0_h, final_state$CA_h, final_state$C_II_h, final_state$C_III_h, final_state$I_h, final_state$I_II_h, final_state$I_III_h)
tot_c <- compute_totals(final_state$S0_c, final_state$SA_c, final_state$S_II_c, final_state$S_III_c, final_state$C0_c, final_state$CA_c, final_state$C_II_c, final_state$C_III_c, final_state$I_c, final_state$I_II_c, final_state$I_III_c)
alpha_eq <- compute_alpha_dynamic(tot_h, tot_c, final_state$I_c, final_state$I_II_c, final_state$I_III_c, delta, w, w_I, w_II, w_III)$alpha

# 4) Store final outputs in one simple object for the next sections
final_results <- list(
  params_final = final_params,
  ode_result = final_res,
  metrics = final_metrics,
  alpha_eq = alpha_eq,
  optimization = optimization_results
)

# 5) Save calibration outputs
saveRDS(optimization_results, "résultats/new_calibration_results.rds")
saveRDS(final_results, "résultats/new_final_results.rds")
write.csv(data.frame(parameter = names(calibrated_params), value = as.numeric(calibrated_params)),
  "résultats/new_calibrated_parameters.csv", row.names = FALSE
)

# 6) Check population is constant (conservation)
final_res$N_h <- with(final_res, S0_h + SA_h + S_II_h + S_III_h + C0_h + CA_h + C_II_h + C_III_h + I_h + I_II_h + I_III_h)
final_res$N_c <- with(final_res, S0_c + SA_c + S_II_c + S_III_c + C0_c + CA_c + C_II_c + C_III_c + I_c + I_II_c + I_III_c)
max_N_h_diff_final <- max(abs(final_res$N_h - N_h))
max_N_c_diff_final <- max(abs(final_res$N_c - N_c))

cat(sprintf("Max |N_h - N_h0|: %.3e\n", max_N_h_diff_final))
cat(sprintf("Max |N_c - N_c0|: %.3e\n", max_N_c_diff_final))
cat(sprintf("Population constant: %s\n", max_N_h_diff_final <= 1e-6 & max_N_c_diff_final <= 1e-6))

# 7) Check alpha is constant over the last year
alpha_final <- sapply(seq_len(nrow(final_res)), function(i) {
  state_i <- as.list(final_res[i, ])
  tot_h <- compute_totals(state_i$S0_h, state_i$SA_h, state_i$S_II_h, state_i$S_III_h, state_i$C0_h, state_i$CA_h, state_i$C_II_h, state_i$C_III_h, state_i$I_h, state_i$I_II_h, state_i$I_III_h)
  tot_c <- compute_totals(state_i$S0_c, state_i$SA_c, state_i$S_II_c, state_i$S_III_c, state_i$C0_c, state_i$CA_c, state_i$C_II_c, state_i$C_III_c, state_i$I_c, state_i$I_II_c, state_i$I_III_c)

  compute_alpha_dynamic(tot_h, tot_c, state_i$I_c, state_i$I_II_c, state_i$I_III_c,
                        delta, w, w_I, w_II, w_III)$alpha
})

alpha_final_last_year <- tail(alpha_final, min(365, length(alpha_final)))
max_alpha_diff_final <- max(abs(alpha_final_last_year - tail(alpha_final_last_year, 1)))

cat(sprintf("Max |alpha - alpha_final| over last year: %.3e\n", max_alpha_diff_final))
cat(sprintf("Alpha constant over last year: %s\n", max_alpha_diff_final <= 1e-6))

# 8) Table plots
p_table_parameters <- create_all_parameters_table(final_params, N_h, N_c, alpha_eq)
p_table_calibration <- create_calibration_table(final_params, final_metrics, targets)

grid::grid.newpage()
grid::grid.draw(p_table_parameters)

grid::grid.newpage()
grid::grid.draw(p_table_calibration)

ggplot2::ggsave(filename = "résultats/table_all_parameters.png", plot = p_table_parameters, width = 11, height = 9, dpi = 300)
ggplot2::ggsave(filename = "résultats/table_calibration.png", plot = p_table_calibration, width = 11, height = 6, dpi = 300)

# 9) Dynamics after calibration
plots_dyn <- plot_dynamics(final_res, final_params, targets, N_h, N_c)
p_dyn_target <- plots_dyn$hospital_target + plots_dyn$community_target
p_dyn_all <- plots_dyn$hospital_all + plots_dyn$community_all

print(p_dyn_target)
print(p_dyn_all)

ggplot2::ggsave(filename = "résultats/dynamique_modele_agregee_hopital_communaute.png", plot = p_dyn_target, width = 14, height = 6, dpi = 300)










###############################################################################
# ---- 6) CHECK STEADY STATE BEFORE SIMULATIONS ----
###############################################################################

cat("\n=== TEST ETAT STATIONNAIRE ===\n")

########### 1) Get compartment names ###########
names_list <- get_stationary_state_names()


########### 2) Check final calibrated state ###########
final_check <- check_final_equilibrium(
  params_final = final_params,
  final_res = final_res,
  state_names = names_list$state_names,
  equilibrium_tol = 1e-6,
  model_function = cdiff_model
)


########### 3) Check different initial conditions ###########

# Define the initial conditions to be tested
init_test_values <- c(0.001, 0.002, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.02, 0.025, 0.03)
initial_conditions_list <- list()

for (i in seq_along(init_test_values)) {
  p_init <- init_test_values[i]

  init_prop <- c(SA_h = p_init, C0_h = p_init, CA_h = p_init, I_h = p_init, S_II_h = p_init, C_II_h = p_init, I_II_h = p_init, S_III_h = p_init, C_III_h = p_init, I_III_h = p_init)
  init_h <- c(S0_h = 1 - sum(init_prop), init_prop) * N_h

  init_prop <- c(SA_c = p_init, C0_c = p_init, CA_c = p_init, I_c = p_init, S_II_c = p_init, C_II_c = p_init, I_II_c = p_init, S_III_c = p_init, C_III_c = p_init, I_III_c = p_init)
  init_c <- c(S0_c = 1 - sum(init_prop), init_prop) * N_c

  initial_conditions_list[[paste0("IC_", i)]] <- c(init_h, init_c)
}

# Check different initial conditions
initial_condition_check <- check_initial_conditions(
  params_final = final_params,
  initial_conditions_list = initial_conditions_list,
  time_vec = time_vec,
  reference_state = final_check$reference_state,
  state_names = names_list$state_names,
  equilibrium_tol = 1e-6,
  distance_tol = 1e-6,
  model_function = cdiff_model
)

cat("\n--- Initial-condition test table ---\n")
cat("initial_condition: tested initial condition.\n")
cat("equilibrium_reached: TRUE if max|dX/dt| <= 1e-6.\n")
cat("final_time: final simulation time in days.\n")
cat("same_equilibrium: TRUE if distance_to_reference <= 1e-6.\n")
print(initial_condition_check[, c("initial_condition", "equilibrium_reached", "final_time", "same_equilibrium")])

p_table_initial_conditions <- create_initial_condition_check_table(initial_condition_check)
grid::grid.newpage()
grid::grid.draw(p_table_initial_conditions)
ggplot2::ggsave("résultats/table_initial_condition_check.png", plot = p_table_initial_conditions, width = 9, height = 5, dpi = 300)


########### 4) Define the perturbations to be tested ###########
perturbation_table <- data.frame(
  perturbation_compartment = c("S0_h", "S0_h", "SA_h", "SA_h", "S_II_h", "S_II_h", "S_III_h", "S_III_h", "C0_h", "C0_h", "CA_h", "CA_h", "C_II_h", "C_II_h", "C_III_h", "C_III_h", "I_h", "I_h", "I_II_h", "I_II_h", "I_III_h", "I_III_h", 
                               "S0_c", "S0_c", "SA_c", "SA_c", "S_II_c", "S_II_c", "S_III_c", "S_III_c", "C0_c", "C0_c", "CA_c", "CA_c", "C_II_c", "C_II_c", "C_III_c", "C_III_c", "I_c", "I_c", "I_II_c", "I_II_c", "I_III_c", "I_III_c"
                               ),
  perturbation_pct = c(0.01, -0.01, 0.01, -0.01, 0.01, -0.01, 0.01, -0.01, 0.01, -0.01, 0.01, -0.01, 0.01, -0.01, 0.01, -0.01, 0.01, -0.01, 0.01, -0.01, 0.01, -0.01, 
                       0.01, -0.01, 0.01, -0.01, 0.01, -0.01, 0.01, -0.01, 0.01, -0.01, 0.01, -0.01, 0.01, -0.01, 0.01, -0.01, 0.01, -0.01, 0.01, -0.01, 0.01, -0.01),
  stringsAsFactors = FALSE
)

# 6) Check +1% and -1% perturbations
perturbation_check <- check_perturbations(
  params_final = final_params,
  perturbation_table = perturbation_table,
  time_vec = time_vec,
  reference_state = final_check$reference_state,
  hospital_names = names_list$hospital_names,
  community_names = names_list$community_names,
  state_names = names_list$state_names,
  N_h = N_h,
  N_c = N_c,
  equilibrium_tol = 1e-6,
  distance_tol = 1e-6,
  model_function = cdiff_model
)

cat("\n--- Perturbation test table ---\n")
cat("perturbation_compartment: perturbed compartment.\n")
cat("perturbation_pct: perturbation size (+0.01 = +1%, -0.01 = -1%).\n")
cat("equilibrium_reached: TRUE if max|dX/dt| <= 1e-6.\n")
cat("same_equilibrium: TRUE if distance_to_reference <= 1e-6.\n")
print(perturbation_check)

p_table_perturbations <- create_perturbation_check_table(perturbation_check)
grid::grid.newpage()
grid::grid.draw(p_table_perturbations)
ggplot2::ggsave("résultats/table_perturbation_check.png", plot = p_table_perturbations, width = 9, height = 13, dpi = 300)









###############################################################################
# ---- 7) SCENARIO SIMULATIONS ----
###############################################################################

cat("\n=== SIMULATIONS ===\n")


########### 7.1) Preparation for simulations ###########

# Parameters used for all simulations
params_sim <- final_results$params_final
params_sim["alpha_const"] <- as.numeric(alpha_eq)
params_sim["use_fixed_alpha"] <- 1

# Initial condition used for non-vaccination simulations
final_state_for_sim <- as.list(final_state_eq[1, ])
init_sim <- unlist(final_state_for_sim[names(final_state_for_sim) != "time"])
init_sim <- as.numeric(init_sim)
names(init_sim) <- names(final_state_for_sim)[names(final_state_for_sim) != "time"]

# Time vector used for simulations
time_sim <- seq(0, 5 * 365, by = 1)

# Compartment order used for vaccination simulations
base_order_h <- c("S0_h", "SA_h", "C0_h", "CA_h", "I_h", "S_II_h", "C_II_h", "I_II_h", "S_III_h", "C_III_h", "I_III_h")
base_order_c <- c("S0_c", "SA_c", "C0_c", "CA_c", "I_c", "S_II_c", "C_II_c", "I_II_c", "S_III_c", "C_III_c", "I_III_c")


########### 7.2) Precheck before simulations ###########

# Baseline simulation
ode_precheck_base <- run_model(
  params_vec = params_sim,
  init_cond = init_sim,
  time_vec = time_sim,
  model_function = cdiff_model
)

# Hygiene with r = 0 must give the same output as baseline
params_hyg_r0 <- params_sim
params_hyg_r0["beta_h"] <- params_sim["beta_h"] * (1 - 0)

ode_precheck_hyg_r0 <- run_model(
  params_vec = params_hyg_r0,
  init_cond = init_sim,
  time_vec = time_sim,
  model_function = cdiff_model
)

print(compare_outputs(ode_precheck_base, ode_precheck_hyg_r0, "Hygiene r=0 vs baseline"))

# Antibiotic reduction with r = 0 must give the same output as baseline
params_atb_r0 <- params_sim
params_atb_r0["tau_c"] <- params_sim["tau_c"] * (1 - 0)

ode_precheck_atb_r0 <- run_model(
  params_vec = params_atb_r0,
  init_cond = init_sim,
  time_vec = time_sim,
  model_function = cdiff_model
)

print(compare_outputs(ode_precheck_base, ode_precheck_atb_r0, "ATB r=0 vs baseline"))

# Vaccination with VC = 0 and VE = 0 must give the same output as baseline
VC <- 0
VE <- 0

init_nv <- init_sim * (1 - VC)
init_v <- init_sim * VC
names(init_nv) <- paste0(names(init_sim), "_nv")
names(init_v) <- paste0(names(init_sim), "_v")

init_vacc_precheck <- c(
  init_nv[paste0(base_order_h, "_nv")],
  init_nv[paste0(base_order_c, "_nv")],
  init_v[paste0(base_order_h, "_v")],
  init_v[paste0(base_order_c, "_v")]
)

params_vacc_precheck <- params_sim
params_vacc_precheck["VE"] <- VE

ode_precheck_vacc <- run_model(
  params_vec = params_vacc_precheck,
  init_cond = init_vacc_precheck,
  time_vec = time_sim,
  model_function = cdiff_vacc
)

print(compare_outputs(ode_precheck_base, ode_precheck_vacc, "Vaccination VC=0, VE=0 vs baseline"))

# Vaccination with VC = 100% and VE = 0 must also give the same output as baseline
VC <- 1
VE <- 0

init_nv <- init_sim * (1 - VC)
init_v <- init_sim * VC
names(init_nv) <- paste0(names(init_sim), "_nv")
names(init_v) <- paste0(names(init_sim), "_v")

init_vacc_precheck <- c(
  init_nv[paste0(base_order_h, "_nv")],
  init_nv[paste0(base_order_c, "_nv")],
  init_v[paste0(base_order_h, "_v")],
  init_v[paste0(base_order_c, "_v")]
)

params_vacc_precheck <- params_sim
params_vacc_precheck["VE"] <- VE

ode_precheck_vacc <- run_model(
  params_vec = params_vacc_precheck,
  init_cond = init_vacc_precheck,
  time_vec = time_sim,
  model_function = cdiff_vacc
)

print(compare_outputs(ode_precheck_base, ode_precheck_vacc, "Vaccination VC=100, VE=0 vs baseline"))


########### 7.3) Scenario simulations ###########

simulation_results <- list()

# Baseline
simulation_results$baseline <- list(
  params = params_sim,
  ode = run_model(params_sim, init_sim, time_sim, model_function = cdiff_model)
)

# Hospital hygiene scenarios
simulation_results$hygiene <- list()

params_hyg_low <- params_sim
params_hyg_low["beta_h"] <- params_sim["beta_h"] * (1 - 0.05)
simulation_results$hygiene$low <- list(
  reduction = 0.05,
  params = params_hyg_low,
  ode = run_model(params_hyg_low, init_sim, time_sim, model_function = cdiff_model)
)

params_hyg_medium <- params_sim
params_hyg_medium["beta_h"] <- params_sim["beta_h"] * (1 - 0.10)
simulation_results$hygiene$medium <- list(
  reduction = 0.10,
  params = params_hyg_medium,
  ode = run_model(params_hyg_medium, init_sim, time_sim, model_function = cdiff_model)
)

params_hyg_strong <- params_sim
params_hyg_strong["beta_h"] <- params_sim["beta_h"] * (1 - 0.20)
simulation_results$hygiene$strong <- list(
  reduction = 0.20,
  params = params_hyg_strong,
  ode = run_model(params_hyg_strong, init_sim, time_sim, model_function = cdiff_model)
)

# Community antibiotic reduction scenarios
simulation_results$antibiotics <- list()

params_atb_low <- params_sim
params_atb_low["tau_c"] <- params_sim["tau_c"] * (1 - 0.05)
simulation_results$antibiotics$low <- list(
  reduction = 0.05,
  params = params_atb_low,
  ode = run_model(params_atb_low, init_sim, time_sim, model_function = cdiff_model)
)

params_atb_medium <- params_sim
params_atb_medium["tau_c"] <- params_sim["tau_c"] * (1 - 0.10)
simulation_results$antibiotics$medium <- list(
  reduction = 0.10,
  params = params_atb_medium,
  ode = run_model(params_atb_medium, init_sim, time_sim, model_function = cdiff_model)
)

params_atb_strong <- params_sim
params_atb_strong["tau_c"] <- params_sim["tau_c"] * (1 - 0.20)
simulation_results$antibiotics$strong <- list(
  reduction = 0.20,
  params = params_atb_strong,
  ode = run_model(params_atb_strong, init_sim, time_sim, model_function = cdiff_model)
)

# Vaccination scenarios
simulation_results$vaccination <- list()
vacc_VE_values <- c(0.30, 0.50, 0.70)
vacc_VC_values <- c(0.20, 0.40, 0.60, 0.80, 1)

for (VE in vacc_VE_values) {
  for (VC in vacc_VC_values) {
    scenario_name <- paste0("VE", round(100 * VE), "_VC", round(100 * VC))
    cat(sprintf("Running vaccination scenario: %s\n", scenario_name))

    init_nv <- init_sim * (1 - VC)
    init_v <- init_sim * VC
    names(init_nv) <- paste0(names(init_sim), "_nv")
    names(init_v) <- paste0(names(init_sim), "_v")

    init_vacc <- c(
      init_nv[paste0(base_order_h, "_nv")],
      init_nv[paste0(base_order_c, "_nv")],
      init_v[paste0(base_order_h, "_v")],
      init_v[paste0(base_order_c, "_v")]
    )

    params_vacc <- params_sim
    params_vacc["VE"] <- VE

    simulation_results$vaccination[[scenario_name]] <- list(
      VE = VE,
      VC = VC,
      params = params_vacc,
      init_cond = init_vacc,
      ode = run_model(params_vacc, init_vacc, time_sim, model_function = cdiff_vacc)
    )
  }
}

saveRDS(simulation_results, "résultats/new_simulation_results.rds")
cat("Simulation results saved in: résultats/new_simulation_results.rds\n")


########### 7.4) Scenario plots ###########

# Hospital hygiene plots
hygiene_plot_data <- prepare_intervention_plot_data(simulation_results, intervention = "HYG")
p_hygiene <- plot_intervention_reference_panels(hygiene_plot_data, intervention = "HYG")
print(p_hygiene$carriage)
print(p_hygiene$incidence)
print(p_hygiene$primo_rec)
print(p_hygiene$combined)
ggplot2::ggsave("résultats/plot_hygiene_carriage.png", plot = p_hygiene$carriage, width = 8, height = 5.5, dpi = 300)
ggplot2::ggsave("résultats/plot_hygiene_incidence.png", plot = p_hygiene$incidence, width = 8, height = 5.5, dpi = 300)
ggplot2::ggsave("résultats/plot_hygiene_primo_rec.png", plot = p_hygiene$primo_rec, width = 8, height = 5.5, dpi = 300)
ggplot2::ggsave("résultats/plot_hygiene_combined.png", plot = p_hygiene$combined, width = 18, height = 5.5, dpi = 300)

# Community antibiotic reduction plots
antibiotics_plot_data <- prepare_intervention_plot_data(simulation_results, intervention = "ATB")
p_antibiotics <- plot_intervention_reference_panels(antibiotics_plot_data, intervention = "ATB")
print(p_antibiotics$carriage)
print(p_antibiotics$incidence)
print(p_antibiotics$primo_rec)
print(p_antibiotics$combined)
ggplot2::ggsave("résultats/plot_antibiotics_carriage.png", plot = p_antibiotics$carriage, width = 8, height = 5.5, dpi = 300)
ggplot2::ggsave("résultats/plot_antibiotics_incidence.png", plot = p_antibiotics$incidence, width = 8, height = 5.5, dpi = 300)
ggplot2::ggsave("résultats/plot_antibiotics_primo_rec.png", plot = p_antibiotics$primo_rec, width = 8, height = 5.5, dpi = 300)
ggplot2::ggsave("résultats/plot_antibiotics_combined.png", plot = p_antibiotics$combined, width = 18, height = 5.5, dpi = 300)

# Vaccination plots
vaccination_plot_data <- prepare_vaccination_plot_data(simulation_results)
p_vaccination <- plot_vacc_impact(vaccination_plot_data)
print(p_vaccination$carriage)
print(p_vaccination$incidence)
print(p_vaccination$primo_rec)
print(p_vaccination$combined)
ggplot2::ggsave("résultats/plot_vaccination_carriage.png", plot = p_vaccination$carriage, width = 12, height = 5, dpi = 300)
ggplot2::ggsave("résultats/plot_vaccination_incidence.png", plot = p_vaccination$incidence, width = 12, height = 5, dpi = 300)
ggplot2::ggsave("résultats/plot_vaccination_primo_rec.png", plot = p_vaccination$primo_rec, width = 12, height = 5, dpi = 300)
ggplot2::ggsave("résultats/plot_vaccination_combined.png", plot = p_vaccination$combined, width = 12, height = 10, dpi = 300)










###############################################################################
# ---- 8) PRCC SENSITIVITY ANALYSIS ----
###############################################################################

# Each tested parameter is randomly changed between -50% and +50%.
prcc_results <- run_prcc_sensitivity(
  params_final = params_sim,
  init_cond = init_sim,
  time_vec = time_vec,
  n_samples = 100,
  seed = 123,
  variation = 0.50,
  equilibrium_tol = 1e-6,
  model_function = cdiff_model
)

p_prcc <- plot_prcc_results(prcc_results$prcc)
prcc_results$plot <- p_prcc
p_prcc_pairs <- plot_prcc_pairs(prcc_results$prcc)
prcc_results$pair_plots <- p_prcc_pairs

print(p_prcc)
print(p_prcc_pairs$incidence)
print(p_prcc_pairs$primo_rec)
print(p_prcc_pairs$carriage)

saveRDS(prcc_results, "résultats/new_prcc_results.rds")
write.csv(prcc_results$prcc, "résultats/new_prcc_table.csv", row.names = FALSE)
ggplot2::ggsave("résultats/prcc_plot.png", plot = p_prcc, width = 12, height = 16, dpi = 300)
ggplot2::ggsave("résultats/prcc_incidence.png", plot = p_prcc_pairs$incidence, width = 10, height = 7, dpi = 300)
ggplot2::ggsave("résultats/prcc_primo_rec.png", plot = p_prcc_pairs$primo_rec, width = 10, height = 7, dpi = 300)
ggplot2::ggsave("résultats/prcc_carriage.png", plot = p_prcc_pairs$carriage, width = 10, height = 7, dpi = 300)

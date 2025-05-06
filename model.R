# Load packages
library(deSolve)  
library(ggplot2)  
library(reshape2) 

# Define the SEIR-SEI model function 
seir_sei_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total populations
    N_h <- S_h + E_h + I_h + R_h
    N_m <- S_m + E_m + I_m
    
    # Force of infection
    lambda_h <- beta_h * I_m / N_m
    lambda_m <- beta_m * I_h / N_h
    
    # Differential equations for humans
    dS_h <- omega_h * R_h - lambda_h * S_h
    dE_h <- lambda_h * S_h - sigma_h * E_h
    dI_h <- sigma_h * E_h - gamma_h * I_h
    dR_h <- gamma_h * I_h - omega_h * R_h
    
    # Differential equations for mosquitoes with death rate
    dS_m <- -lambda_m * S_m
    dE_m <- lambda_m * S_m - sigma_m * E_m
    dI_m <- sigma_m * E_m - mu_m * I_m
    
    #rate of change
    return(list(c(dS_h, dE_h, dI_h, dR_h, dS_m, dE_m, dI_m)))
  })
}

# Set model parameters
parameters <- c(
  beta_h = 0.5,    # Transmission rate (mosquito to human)
  beta_m = 0.4,    # Transmission rate (human to mosquito)
  sigma_h = 0.167, # Human incubation rate (1/6 days)
  sigma_m = 0.143, # Mosquito incubation rate (1/7 days)
  gamma_h = 0.143, # Human recovery rate (1/7 days)
  omega_h = 0.01,  # Rate of waning immunity in humans (1/100 days)
  mu_m = 0.033     # Mosquito death rate for infectious mosquitoes (1/30 days)
)

# Set initial state
initial_state <- c(
  S_h = 9900,  # Susceptible humans
  E_h = 50,    # Exposed humans
  I_h = 50,    # Infected humans
  R_h = 0,     # Recovered humans
  S_m = 50000, # Susceptible mosquitoes
  E_m = 100,   # Exposed mosquitoes
  I_m = 100    # Infected mosquitoes
)

# 150 days
times <- seq(0, 150, by = 1)

# Solve the system of ODEs
solution <- ode(
  y = initial_state,
  times = times,
  func = seir_sei_model,
  parms = parameters,
  method = "lsoda"
)

# Convert the output to a data frame
solution_df <- as.data.frame(solution)

# Prepare data for plotting

# Human compartments
human_data <- solution_df[, c("time", "S_h", "E_h", "I_h", "R_h")]
human_melted <- melt(human_data, id.vars = "time", 
                     variable.name = "Compartment", 
                     value.name = "Population")

# Mosquito compartments
mosquito_data <- solution_df[, c("time", "S_m", "E_m", "I_m")]
mosquito_melted <- melt(mosquito_data, id.vars = "time", 
                        variable.name = "Compartment", 
                        value.name = "Population")

# Create plots

# Human population plot
human_plot <- ggplot(human_melted, 
                     aes(x = time, y = Population, color = Compartment)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("S_h" = "blue", "E_h" = "orange", 
                                "I_h" = "red", "R_h" = "green"),
                     labels = c("Susceptible", "Exposed", "Infected", "Recovered")) +
  labs(title = "Human Population",
       x = "Time (days)",
       y = "Individuals",
       color = "Status") +
  theme_minimal()

# Mosquito population plot
mosquito_plot <- ggplot(mosquito_melted,
                        aes(x = time, y = Population, color = Compartment)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("S_m" = "blue", "E_m" = "orange", "I_m" = "red"),
                     labels = c("Susceptible", "Exposed", "Infected")) +
  labs(title = "Aedes Population",
       x = "Time (days)",
       y = "Mosquitoes",
       color = "Status") +
  theme_minimal()

# Prevalence in human and mosquito populations
prevalence_df <- data.frame(
  time = solution_df$time,
  Human_prevalence = solution_df$I_h / 
    (solution_df$S_h + solution_df$E_h + solution_df$I_h + solution_df$R_h),
  Mosquito_prevalence = solution_df$I_m / 
    (solution_df$S_m + solution_df$E_m + solution_df$I_m)
)

prevalence_melted <- melt(prevalence_df, id.vars = "time", 
                          variable.name = "Population", 
                          value.name = "Prevalence")

prevalence_plot <- ggplot(prevalence_melted, 
                          aes(x = time, y = Prevalence, color = Population)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("Human_prevalence" = "darkred", 
                                "Mosquito_prevalence" = "darkblue"),
                     labels = c("Human", "Aedes")) +
  labs(title = "Disease Prevalence",
       x = "Time (days)",
       y = "Proportion infected",
       color = "Population") +
  theme_minimal() +
  ylim(0, NA)  # Start y-axis at 0

# Total mosquito population plot
solution_df$total_mosquitoes <- solution_df$S_m + solution_df$E_m + solution_df$I_m

mosquito_pop_plot <- ggplot(solution_df, aes(x = time, y = total_mosquitoes)) +
  geom_line(size = 1, color = "purple") +
  labs(title = "Total Aedes Population",
       x = "Time (days)",
       y = "Mosquitoes") +
  theme_minimal()

# Print the plots
print(human_plot)
print(mosquito_plot)
print(mosquito_pop_plot)
print(prevalence_plot)

# Arrange multiple plots using patchwork package
library(patchwork)

combined_plot <- human_plot / mosquito_plot / mosquito_pop_plot / prevalence_plot +
  plot_layout(ncol = 1, heights = c(1, 1, 1, 1)) +
  plot_annotation(title = "DENV Model Dynamics")

print(combined_plot)

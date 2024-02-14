library(groundhog)
library(lubridate)
pkgs <- c("conflicted","SimInf","ggplot2","dplyr","mlrMBO")
groundhog.library(pkgs, "2023-04-30")
conflicts_prefer(dplyr::filter)

# set the directory to find the csv in the folder data
setwd("~/Documents/GitHub/StochasticMetapop/Plague")

source("/Users/dianapli/Documents/GitHub/StochasticMetapop/Plague/Number_Infected_Parishes.R")

YSTAD_group <- read_csv("/Users/dianapli/Desktop/PythonMathematicalModeling/docs/PlagueProject/data/private/yellow_group.csv", col_types = cols())
# parish names without repetitions
parish_names <- unique(YSTAD_group$ParishName)
npatches <- length(parish_names)
patchPop <- YSTAD_group$BEF1699[1:npatches]
maxDays <- max(YSTAD_group$EndDaysPlague)

# generating the initial conditions for the model
S0 <- rep(0, npatches) 
Ilow0 <- rep(0, npatches)
Ihigh0 <- rep(0, npatches)
R0 <- rep(0, npatches)
Dcum0 <- rep(0, npatches)
Sv0 <- rep(0, npatches)
Iv0 <- rep(20, npatches)
Ilow0[1] <- 1.0
Ihigh0[1] <- 1.0
Ilow0[2] <- 1.0
Ihigh0[2] <- 1.0

for (i in 1:npatches) {
  S0[i] <- patchPop[i] - Ilow0[i] - Ihigh0[i] - R0[i]
}

for (i in 1:npatches) {
  Sv0[i] <- patchPop[i]
}

u0 <- data.frame(
  S0 = S0,
  Ilow0 = Ilow0,
  Ihigh0 = Ihigh0,
  R0 = R0,
  Dcum0 = Dcum0,
  Sv0 = Sv0,
  Iv0 = Iv0
)

# model
transitions <- c("S -> beta_v * Iv * S / (S+Ilow+Ihigh+R) -> Ilow", 
                 "S -> beta_p * Ihigh * S / (S+Ilow+Ihigh+R) -> Ihigh",# New equation
                 "Ilow -> (1 - mu) * sigma * Ilow -> R",
                 "Ilow -> mu * sigma * Ilow -> Ihigh",
                 "Ihigh -> gamma * Ihigh -> Dcum",
                 "@ -> growth_v * Sv * (1 - (Sv + Iv)/((S+Ilow+Ihigh+R) * index_ave)) -> Sv",
                 "Sv -> growth_v * Sv * (1 - (Sv + Iv)/((S+Ilow+Ihigh+R) * index_ave)) -> @", # @ 
                 "Sv -> (beta_low * Ilow + beta_high * Ihigh) * Sv / (S+Ilow+Ihigh+R) -> Iv",
                 "Iv -> gamma_v * Iv -> @")

compartments <- c("S","Ilow","Ihigh","R","Dcum","Sv","Iv")
colnames(u0) <- compartments

# EMatrix (#compartments x #events). For each column vector j we put 1 if the compartment
# participated in the event j, otherwise zero.
E <- matrix(c(0,0,0,0,0,0,1), nrow = 7, ncol = 1
            , dimnames = list(c("S", "Ilow", "Ihigh", "R", "Dcum","Sv","Iv")
                              , c("event1")))

connection_matrix <- read_csv("/Users/dianapli/Desktop/PythonMathematicalModeling/docs/PlagueProject/data/private/connection_Yellow_group.csv"
                              , col_types = cols())

connection_matrix <- as.matrix(connection_matrix)

global_parameters <- c(index_ave = 23, growth_v = 0.14, beta_v = 0.09
                       , gamma_v = 0.71, beta_low = 0.021, beta_high = 0.66
                       , beta_p = 0.24, sigma = 0.20, gamma = 0.84
                       , mu = 0.46)

# Function to calculate the error in the cumulative number of infected parishes per month between the model and the data
objectiveFunction_2 <- function(local_parameter, gdf) {
  
  # Create an external transfer event to move infected vectors from node i to
  # other nodes starting seven days after the initial plague in the node i until
  # maxDays and every seven days 
  
  list_of_events <- list()
  for(i in 1:npatches){
    initial_time <- gdf$BeginDaysPlague[i]
    infect_event <- lapply(seq(from = initial_time + 7, to = maxDays, by = 7),
                           function(t) {
                             infect <- data.frame(
                               event = "extTrans",
                               time = t,
                               node = i,
                               dest = setdiff(1:npatches, i),
                               n = 0,
                               proportion = local_parameter$Iv_prop[i] * connection_matrix[i, setdiff(1:npatches, i)],
                               select = 1,
                               shift = 0
                             )
                           })
    list_of_events[[i]] <- infect_event
  }
  events <- do.call(rbind, unlist(list_of_events, recursive = FALSE))
  
  # Defining the model
  model_subpneumonic <- mparse(transitions = transitions,
                        compartments = compartments,
                        gdata = global_parameters,
                        ldata = local_parameter,
                        u0 = u0,
                        tspan = 1:maxDays,
                        events = events,
                        E = E
  )
  result <- run(model = model_subpneumonic)
  traj_D <- trajectory(model = result, compartments = "Dcum")
  # Defining the initial date of the gdf to start counting the number of infected parishes per month
  date <- min(gdf$BeginPlaguePeriod, na.rm = TRUE)
  
  # Getting the number of infected parishes per month from the data
  infected_parishes_month <- count_infected_parishes_month(gdf,date)
  
  # Initializing the number of infected parishes per month for the model's output
  model_infected_parishes_month <- rep(0, length(infected_parishes_month))
  
  # Initializing the error between the model's output and the data
  error <- rep(0, length(infected_parishes_month))
  
  # Computing the total number of parishes in the dataframe without repetitions
  total_parishes <- length(gdf$ParishName)
  
  # Computing the number of infected parishes per month from the model's output
  for (i in 1:length(infected_parishes_month)) {
    init_days <- infected_parishes_month[i,'DaysFromInitDate']
    final_days <- infected_parishes_month[i,'DaysToEndOfMonth']
    
    for (k in 1:total_parishes) {
      for (day in init_days:final_days) {
        # Check if there are any rows that satisfy the condition
        rows <- traj_D[traj_D$node == k & traj_D$time == day, ]
        if (nrow(rows) == 0){
          break # Breaks the innermost loop when the condition is met
        }
        else if(nrow(rows) > 0){
          if(rows$Dcum >= 1){
            model_infected_parishes_month[i] <- model_infected_parishes_month[i] + 1
            break # Breaks the innermost loop when the condition is met
          }
        }
      }
    }
    error[i] <- (model_infected_parishes_month[i] - infected_parishes_month[i,'NumberInfectedParishes'])^2
  }
  
  max_error <- max(error)
  
  # Computing the error between the model's output and the data
  total_error <- sum(error) /(length(infected_parishes_month) * max_error)
  
  return(total_error)
}

# Define the objective function for mlrMBO
objectiveFunction_2_mlr <- makeSingleObjectiveFunction(
  name = "My Function",
  noisy = TRUE,
  has.simple.signature= TRUE,
  fn = function(xs){
    Iv_prop <- rep(0, npatches)
    K <- 1
    for (I in 1:npatches){
      Iv_prop[K] <- xs[I]
      K <- K + 1
    }
    # now, make a df such that the columns are the vectors beta and mu
    Iv_prop <- as.data.frame(Iv_prop)
    colnames(Iv_prop) <- c("Iv_prop")
    return (objectiveFunction_2(local_parameter = Iv_prop
                                , gdf = YSTAD_group
                                ))},
  par.set = ps,
  minimize = TRUE
)

ps<- makeParamSet(
  makeNumericVectorParam("Iv_prop", len=npatches, lower = 0, upper = 1)
)

# Generate an initial design
des = generateDesign(n = 80, par.set = ps)
des$y = apply(des, 1, objectiveFunction_2_mlr)
des <- des[order(des$y),][1:40,]

control = makeMBOControl(final.method = "best.predicted")
control = setMBOControlTermination(control, iters = 80)
control = setMBOControlInfill(control, crit = crit.eqi)

# Run the optimization
res <- mbo(fun = objectiveFunction_2_mlr, design = des, control = control, show.info = TRUE)

# Printout
best_params <- res$x
options(scipen=999)
infected_parishes <- run_infected_parishes_model(res$x, global_parameters = global_parameters, gdf = YSTAD_group, n=0)
infected_parishes
plot_infected_parishes(infected_parishes, YSTAD_group, 0)

## Run the model with the best parameters
best_params <- as.data.frame(best_params)
objectiveFunction_2(best_params, YSTAD_group)

list_of_events <- list()
for(i in 1:npatches){
  initial_time <- YSTAD_group$BeginDaysPlague[i]
  infect_event <- lapply(seq(from = initial_time + 7, to = maxDays, by = 7),
                         function(t) {
                           infect <- data.frame(
                             event = "extTrans",
                             time = t,
                             node = i,
                             dest = setdiff(1:npatches, i),
                             n = 0,
                             proportion = best_params$Iv_prop[i] * connection_matrix[i, setdiff(1:npatches, i)],
                             select = 1,
                             shift = 0
                           )
                         })
  list_of_events[[i]] <- infect_event
}
events <- do.call(rbind, unlist(list_of_events, recursive = FALSE))

# Defining the model
model_subpneumonic <- mparse(transitions = transitions,
                             compartments = compartments,
                             gdata = global_parameters,
                             ldata = best_params,
                             u0 = u0,
                             tspan = 1:maxDays,
                             events = events,
                             E = E
)
result <- run(model = model_subpneumonic)
traj_D <- trajectory(model = result, compartments = "Dcum")
ggplot(traj_D) + geom_line(aes(x = time, y = Dcum, color = factor(node)))




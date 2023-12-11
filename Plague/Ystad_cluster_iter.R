library(SimInf)
library(tools)
library(dplyr)
library(ggplot2)
library(optimx)
library(readr) # for read csv files

# set the directory to find the csv in the folder data
#setwd("~/Documents/GitHub/StochasticMetapop/Plague")

data_cluster <- read_csv("data/Ystad_group.csv")
beginDaysPlague <- data_cluster$BeginDaysPlague
endDaysPlague <- data_cluster$EndDaysPlague
parish_names <- data_cluster$ParishName
patchPop <- data_cluster$BEF1699
maxDays <- max(endDaysPlague)

# data
class_plague <- data_cluster$plague
cum_deaths <- data_cluster$VictimsNumber

# Importing the file with the gravitational term Ni*Nj/dij^2
grav_matrix <- read_csv("data/gravitational.csv")
colnames(grav_matrix) <- parish_names
weight_matrix <- as.matrix(grav_matrix)

# generating the initial conditions for the model
npatches <- length(parish_names)
S0 <- rep(0, npatches) 
E0 <- rep(0, npatches)
R0 <- rep(0, npatches)
D0 <- rep(0, npatches)
I0 <- rep(0, npatches)
I0[1] <- 1.0

for (i in 1:npatches) {
  S0[i] <- patchPop[i] - E0[i] - I0[i] - R0[i]
}

u0 <- data.frame(
  S0 = S0,
  E0 = E0,
  I0 = I0,
  R0 = R0,
  D0 = D0
)

transitions <- c(
  "S-> beta*S*I/(S+E+I+R) -> E",
  "E -> sigma*E -> I",
  "I -> gamma*(1-mu)*I -> R",
  "I -> (gamma*mu)*I -> D"
)
compartments <- c("S", "E", "I", "R", "D")
parameters <- c("beta", "mu")
colnames(u0) <- compartments
tspan <- seq(from = 1, to = maxDays + 20)

#EMatrix (#compartments x #events). For each column vector j we put 1 if the compartment
# participated in the event j, otherwise zero.
E <- matrix(
  c(0, 1, 1, 0, 0),
  nrow = 5,
  ncol = 1,
  dimnames = list(c("S", "E", "I", "R", "D"),
                  c("event1"))
)

# Adding events 

# Create an external transfer event to move exposed or infected individuals
# from node 1 (Ystad) to other nodes every seven days (tspan)

list_of_events1 <- lapply(tspan,
                          function(t) {
                            infect <- data.frame(
                              event = "extTrans",
                              time = t,
                              node = 1,
                              dest = 2:npatches,
                              n = 0,
                              proportion = 0.05,
                              select = 1,
                              shift = 0
                            )
                          })

# Create an external transfer event to move exposed or infected individuals
# from node 1 (Ystad) to other nodes starting seven days after the initial plague
# period in Ystad and every seven days until the end period in Ystad proportional to gravity term

data_init_node <- subset(data_cluster, ParishName == "YSTAD")
initial_time <- data_init_node$BeginDaysPlague[1]
final_time <- data_init_node$EndDaysPlague[1]
for (i in 2:npatches) {
  infect_event <-
    lapply(seq(from = initial_time + 7, to = final_time, by = 7),
           function(t) {
             infect <- data.frame(
               event = "extTrans",
               time = t,
               node = 1,
               dest = 2:npatches,
               n = 0,
               proportion = grav_matrix$YSTAD[i],
               select = 1,
               shift = 0
             )
           })
}

# Creating the events to use in the model
events <- do.call(rbind, list_of_events1)

# Defining the model considering some events
model <- mparse(
  transitions = transitions,
  compartments = compartments,
  gdata = c(sigma = 0.17 , gamma = 0.4),
  ldata = data.frame(beta = numeric(npatches),
                     mu = numeric(npatches)),
  u0 = u0,
  E = E,
  events = events,
  tspan = 1:maxDays + 20
)

## Optimization process

model <- mparse(
  transitions = transitions,
  compartments = compartments,
  gdata = c(sigma = 0.17 , gamma = 0.4),
  ldata = data.frame(beta = numeric(npatches),
                     mu = numeric(npatches)),
  u0 = u0,
  E = E,
  events = events,
  tspan = 1:maxDays + 20
)

# Define the objective function
objectiveFunction <- function(parameters){
  # Update local parameters of the model
  for(i in 1:npatches){
    model@ldata[1,i] <- parameters[1] #beta
    model@ldata[2,i] <- parameters[2] #mu
  }
    
  # Initialize a list to store the results of each simulation
    sim_results <- vector("list", 40)
    
  # Run the simulation k times and store the results
    for(j in 1:40){
      sim_results[[j]] <- run(model)
    }
  # Compute the average trajectory for each node
    avg_trajectory <- sapply(1:npatches, function(node){
      sapply(sim_results, function(sim){
        sim@trajectory[[node]]$D
        }) %>% rowMeans(na.rm = TRUE)
    })
  
  # Compute errors per trajectory
  errors <- sapply(1:length(gdf), function(i){
    initial_position <- beginDaysPlague[i]
    final_position <- endDaysPlague[i]
    deaths <- cum_deaths[i]
    traj_D <- trajectory(model = sol, compartments = "D")
    traj_D_node <- dplyr::filter(traj_D, node == i)
    local_error <- 0
    
    if (!is.na(initial_position)) {
      if (final_position != 0 && deaths != 0) {
        local_error <- abs(traj_D_node$D[initial_position] - 1.0)
        + abs(traj_D_node$D[final_position] - deaths[i])
      }
      else{local_error <- abs(traj_D_node$D[initial_position]-1.0)}
    }
    return(local_error)
  })
  
  total_error <- sum(unlist(errors))
  print(total_error)
  return(total_error)
  }


## Solve the optimization problem

results <- optimx::optimx(c(beta = 0.5, mu = 0.5),
                          objectiveFunction,
                          method = "L-BFGS-B",
                          lower = c(0,0),
                          upper = c(1,1),
                          show_progress())

# Print the result 
cat("error =", results$fvalues, "\n")
cat("beta =", results$par[1], "\n")
cat("mu =", results$par[2], "\n")


# ## Plot all trajectories
# for (node in 1:npatches){
#   df <- data.frame()
#   for(i in 1:length(sim_results)){
#     df_temp <- as.data.frame(sim_results[[i]]@trajectory[[node]])
#     df_temp$Time <- 1:nrow(df_temp)
#     df_temp$Simulation <- rep(i, nrow(df_temp))
#     df <- rbind(df, df_temp)
#   }
#   
#   # Add average trajectory to the data frame
#   df_avg <- data.frame(Time = 1:length(avg_trajectory[[node]]), D = avg_trajectory[[node]], Simulation = "Average")
#   df <- rbind(df, df_avg)
#   
#   # Add best parameters trajectory to the data frame
#   df_best <- data.frame(Time = 1:length(best_trajectory[[node]]), D = best_trajectory[[node]], Simulation = "Best parameters")
#   df <- rbind(df, df_best)
#   
#   # Plot
#   ggplot(df, aes(x = Time, y = D, color = factor(Simulation))) +
#     geom_line() +
#     scale_color_manual(values = c(rep("blue", length(sim_results)), "red", "green")) +
#     labs(title = paste("Node", node), x = "Time", y = "D", color = "Simulation") +
#     theme_minimal()
# }






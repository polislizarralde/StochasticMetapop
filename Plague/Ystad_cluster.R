library(SimInf)
library(tools)
library(dplyr)
library(ggplot2)
library(optimx)
library(gbutils)
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

# Fixing the file to compute the distance between polygons
#data_cluster <- (data_cluster, centroid, c("lat","long"), sep = ",", remove = TRUE, convert = TRUE)

# generating the initial conditions for the model
npatches <- 2
S0 <- rep(0, npatches) # nolint: object_name_linter.
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

# patchPop <- function(df,
#                      column_pop = 'BEF1699',
#                      column_name = 'ParishName'
#                     )
#                     {patchNames <- unique(df[[column_name]])
#                      patchPop <- sapply(patchNames, function(name)
#                                  {unique_pop <- unique(df[df[[name]] == name, column_pop])
#                                   if(length(unique_pop) > 0){
#                                     return(unique_pop[1]) # return only the first unique pop value
#                                   } else {
#                                     return(NA)
#                                   }
#                                   })
#                      return(patchPop)
# }
# pop_values <- patchPop(data_cluster, data_cluster$BEF1699, data_cluster$ParishName)
# View(pop_values)
# patchPop(data_cluster = data_cluster)

transitions <- c(
  "S-> beta*S*I/(S+E+I+R) -> E",
  "E -> sigma*E -> I",
  "I -> gamma*(1-mu)*I -> R",
  "I -> (gamma*mu)*I -> D"
)
compartments <- c("S", "E", "I", "R", "D")
parameters <- c("beta", "mu")
colnames(u0) <- compartments
tspan <- seq(from = 1, to = maxDays + 20, by = 7)

#EMatrix (#compartments x #events). For each column vector j we put 1 if the compartment
# participated in the event j, otherwise zero.
E <- matrix(
  c(0, 1, 1, 0, 0),
  nrow = 5,
  ncol = 1,
  dimnames = list(c("S", "E", "I", "R", "D"),
                  c("event1"))
)

# Adding several events 

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

# Create an external transfer event to move exposed or infected individuals from node 1 (Ystad) to other 
# nodes starting seven days after the initial plague period in each node and every seven days

# for(parish_name in data_cluster$ParishName[-1]) {
#   parish_events <- list()
#   begin_parish_name <-
#     lapply(data_cluster$ParishName, function(x) {
#       data_cluster$BeginDaysPlague[data_cluster$ParishName == x]
#     })
#   end_parish_name <-
#     lapply(data_cluster$ParishName, function(x) {
#       data_cluster$EndDaysPlague[data_cluster$ParishName == x]
#     })
#   
#   
#   }

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
 events <- do.call(rbind, infect_event)

# Create an external transfer event to move exposed or infected individuals
# from node 1 (Ystad) to other nodes starting seven days after the initial plague
# period in Ystad and every seven days until the end period in Ystad

# 
# list_of_events2 <- lapply(head(beginDaysPlague[-1], length(beginDaysPlague[-1])-2),
#                           function(t) {
#                             infect <- data.frame(
#                               event = "extTrans",
#                               time = t+7,
#                               node = 1,
#                               dest = 2:npatches,
#                               n = 0,
#                               proportion = 0.05,
#                               select = 1,
#                               shift = 0
#                             )
#                           })

# Event moving individuals from node 1 to the others proportional to the gravitational term
# for (i in 2:npatches) {
#   for (t in tspan) {
#     list_of_events2 <- lapply(beginDaysPlague[-1],
#                               function(t) {
#                                 infect <- data.frame(
#                                   event = "extTrans",
#                                   time = t,
#                                   node = 1,
#                                   dest = 2:npatches,
#                                   n = 0,
#                                   proportion = grav_matrix$YSTAD[i],
#                                   select = 1,
#                                   shift = 0
#                                 )
#                               })
#   }
# }
# 
# # Initialize list_of_events3
# list_of_events3 <- list()
# 
# # Event moving individuals from node 1 to the others proportional to the gravitational term
# for (t in seq(from = 0, to =21, by = 7)) {
#   for (i in 1:3) {
#     for (j in 1:3) {
#       if (i != j) {
#         # for (t in 0:maxDays) {
#         event <- data.frame(
#           event = "extTrans",
#           time = t,
#           node = i,
#           dest = j,
#           n = 0,
#           proportion = weight_matrix[i, j],
#           select = 1,
#           shift = 0
#         )
#         # Append the events to list_of_events3
#         list_of_events3 <- c(list_of_events3, event)
#       }
#     }
#   }
# }
 
# Creating the events to use in the model
events <- do.call(rbind, list_of_events1)

# Defining mu and beta for each node randomly 
# local_parameters <- data.frame(beta = (runif(npatches)),
#                                mu = runif(npatches))
# 
# local_parameters <- data.frame(beta = numeric(npatches),
#                                 mu = numeric(npatches))
# local_parameters$beta
# # Defining the model considering some events
# model <- mparse(
#   transitions = transitions,
#   compartments = compartments,
#   gdata = c(sigma = 0.17 , gamma = 0.4),
#   ldata = data.frame(beta = numeric(npatches),
#                      mu = numeric(npatches)),
#   u0 = u0,
#   tspan = 1:maxDays + 20
# )
# set.seed(123)
# set_num_threads(1)
# 
# result <- run(model = model)



# # Cumulative deaths
# traj_D <- trajectory(model = result, compartments = "D")
# # Show the points per node
# #View(traj_D)
# # Plot a specific trajectory
# ggplot(traj_D) + geom_line(aes(x = time, y = D, color = factor(node)))
# 
# # Infected
# traj_I <- trajectory(model = result, compartments = "I")
# # Show the points per node
# #View(traj_I)
# # Plot a specific trajectory
# ggplot(traj_I) + geom_line(aes(x = time, y = I, color = factor(node)))



########### Optimization process

model <- mparse(
  transitions = transitions,
  compartments = compartments,
  gdata = c(sigma = 0.17 , gamma = 0.4),
  ldata = matrix(c(0.5,0.3,0.5,0.3),nrow=2,ncol=2,byrow=TRUE),
  u0 = u0,
  # E = E,
  # events = events,
  tspan = 1:maxDays + 20
)

# Define the objective function

gdf <- data_cluster
gdf[c(7,8), c(12,13)] <- NA
gdf[c(3,5), c(6)] <- 0

npatches
parameters <- c(0.3,0.2,0.2,0.3)
npar<- 2
mat <- matrix(parameters,nrow=npatches,ncol=npar,byrow=TRUE)
model <- mparse(
  transitions = transitions,
  compartments = compartments,
  gdata = c(sigma = 0.17 , gamma = 0.4),
  ldata = data.frame(beta = as.numeric(mat[1,]), mu = as.numeric(mat[2,]) ),
  u0 = u0,
  # E = E,
  # events = events,
  tspan = 1:maxDays + 20
)

objectiveFunction <- function(parameters){
  mat <- matrix(parameters,nrow=npatches,ncol=npar,byrow=TRUE)
  # Update local parameters of the model
  model <- mparse(
    transitions = transitions,
    compartments = compartments,
    gdata = c(sigma = 0.17 , gamma = 0.4),
    ldata = data.frame(beta = as.numeric(mat[1,]), mu = as.numeric(mat[2,]) ),
    u0 = u0,
    # E = E,
    # events = events,
    tspan = 1:maxDays + 20
  )
  # Solve the model
  sol <- run(model)
  
  # Compute errors per trajectory
  errors <- sapply(1:npatches, function(i){
    initial_position <- beginDaysPlague[i]
    final_position <- endDaysPlague[i]
    deaths <- cum_deaths[i]
    traj_D <- trajectory(model = sol, compartments = "D")
    traj_D_node <- dplyr::filter(traj_D, node == i)
    
    if (!isNA(initial_position)) {
      if (final_position != 0 && deaths != 0) {
        # we now convert all the values to numeric to avoid problems.
        # we must take the absolute value of the difference between the
        # values of the trajectory at the initial_position minus one plus
        # the absolute value of the difference between the values of the
        # trajectory at the final_position minus the deaths.
        
        ans = ((as.double(traj_D_node$D[initial_position])-1.0))^2 + ((
          as.double(traj_D_node$D[final_position]) - deaths))^2
        print("-------")
        print(as.double(traj_D_node$D[initial_position]))
        print(as.double(traj_D_node$D[final_position]))
        print(deaths)
        print(ans)
        return(ans)
      }
      else{
        print("valores case 2")
        print(traj_D_node$D[initial_position])
        return (traj_D_node$D[initial_position]-1.0)
      }
    }
    print("valor 0 case 3")
    return(0)
  })
  print("errors:=")
  print(errors)
 
 total_error <- sum(unlist(errors))
 print("total_error")
 print(total_error)
 return(total_error)
}

# Solve the optimization problem
initial_parameters <- rep(c(0.5, 0.5), npatches)


results <- optimx::optimx(initial_parameters,
                         objectiveFunction,
                         method = "L-BFGS-B",
                         lower = rep(c(0,0), npatches),
                         upper = rep(c(1,1), npatches)
                         )


                         

# Print the result 

print(results)
# cat("error =", results$fvalues, "\n")
# cat("beta1 =", results$par[1], "\n")
# cat("mu1 =", results$par[2], "\n")
# cat("beta2 =", results$par[3], "\n")
# cat("mu2 =", results$par[4], "\n")


# Run the model with best parameters
# for(i in 1:npatches){
#   model@ldata[1,i] <- results$par[2*i-1] #beta
#   model@ldata[2,i] <- results$par[2*i]   #mu
# }
# 
# final_sol <- run(model)
# 
# # Plotting
# par(mfrow=c(ceiling(npatches/2), 2)) # Adjust layout according to number of nodes
# 
# for(i in 1:npatches) {
#   traj_D_node <- dplyr::filter(trajectory(model = final_sol, compartments = "D"), node == i)
#   plot(traj_D_node$time, traj_D_node$D, main=paste("Node", i), xlab="Time", ylab="D")
# }






library(SimInf)
library(tools)
library(dplyr)
library(ggplot2)
library(optimx)
library(gbutils)
library(readr) # for read csv files

# set the directory to find the csv in the folder data
setwd("~/Documents/GitHub/StochasticMetapop/Plague")

data_cluster <- read_csv("data/Ystad_group.csv")
beginDaysPlague <- data_cluster$BeginDaysPlague
endDaysPlague <- data_cluster$EndDaysPlague
parish_names <- data_cluster$ParishName
patchPop <- data_cluster$BEF1699
maxDays <- max(endDaysPlague)
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
tspan <- seq(from = 1, to = maxDays, by = 1)

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
# from node 1 (Ystad) to other nodes every seven days
infect_event1 <- lapply(seq(from = initial_time + 7, to = maxDays, by = 7),
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
# period in Ystad and every seven days until proportional to the gravity term

data_init_node <- subset(data_cluster, ParishName == "YSTAD")
initial_time <- data_init_node$BeginDaysPlague[1]

# Initialize an empty list
list_of_events2 <- list()
for(i in 2:npatches){
infect_event2 <-
  lapply(seq(from = initial_time + 7, to = maxDays, by = 7),
         function(t) {
               data.frame(
               event = "extTrans",
               time = t,
               node = 1,
               dest = i,
               n = 0,
               proportion = grav_matrix$YSTAD[i],
               select = 1,
               shift = 0
             )
            })
# Append the events to list_of_events2
# list_of_events2 is a list of lists of dataframes
list_of_events2[[i]] <- infect_event2
}
weight_matrix[1,5]

# Create an external transfer event to move exposed or infected individuals
# from node i to other nodes starting seven days after the initial plague in the node i
# until maxDays and every seven days with a proportion given by the respective gravitational term.

# Initialize list_of_events3
list_of_events3 <- list()
for (i in 1:npatches) {
  initial_time <- beginDaysPlague[i]
  for (j in 1:npatches) {
    if (i != j) {
      infect_event3 <-
        lapply(seq(
          from = initial_time + 7,
          to = maxDays,
          by = 7
        ),
        function(t) {
          data.frame(
            event = "extTrans",
            time = t,
            node = i,
            dest = j,
            n = 0,
            proportion = weight_matrix[i, j],
            select = 1,
            shift = 0
          )
          
        })
      # Append the dataframe to the list
      list_of_events3 <- append(list_of_events3, infect_event3)
    }
  }
}

# Creating the events to use in the model
#events <- do.call(rbind, infect_event1)
#events <- do.call(rbind, unlist(list_of_events2, recursive = FALSE))
# The unlist function flattens the list of lists into a single list
events <- do.call(rbind, list_of_events3)
events


# Defining mu and beta for each node randomly 
local_parameters <- data.frame(beta = c(0.7,0.5,0.5,0.5,0.5,0.5,0.5,0.5),
                                 mu = runif(npatches))
local_parameters

# Defining the model considering some events
model <- mparse(
  transitions = transitions,
  compartments = compartments,
  gdata = c(sigma = 0.17 , gamma = 0.4),
  ldata = local_parameters,
  E = E,
  events = events,
  u0 = u0,
  tspan = 1:maxDays + 20
)

set.seed(123)
set_num_threads(1)

result <- run(model = model)

# Cumulative deaths
traj_D <- trajectory(model = result, compartments = "D")
# Show the points per node
#View(traj_D)
# Plot a specific trajectory
ggplot(traj_D) + geom_line(aes(x = time, y = D, color = factor(node)))

# Infected
traj_I <- trajectory(model = result, compartments = "I")
# Show the points per node
#View(traj_I)
# Plot a specific trajectory
ggplot(traj_I) + geom_line(aes(x = time, y = I, color = factor(node)))


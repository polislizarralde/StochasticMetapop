library(SimInf)
library(tools)
library(dplyr)
library(ggplot2)

library(readr) # for read csv files

if (!require(tidyverse)) {
  install.packages("tidyverse")
}

# Load tidyverse
library(tidyverse)

# set the directory to find the csv in the folder data
setwd("~/Documents/GitHub/StochasticMetapop/Plague")

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

grav_matrix$YSTAD[2]

weight_matrix <- as.matrix(grav_matrix)

# Fixing the file to compute the distance between polygons
#data_cluster <- (data_cluster, centroid, c("lat","long"), sep = ",", remove = TRUE, convert = TRUE)

View(data_cluster)

# generating the initial conditions for the model
npatches <- length(parish_names)
S0 <- rep(0, npatches) # nolint: object_name_linter.
E0 <- rep(0, npatches)
R0 <- rep(0, npatches)
D0 <- rep(0, npatches)
I0 <- rep(0, npatches)
I0[1] <- 1.0

for (i in 1:npatches) {
  S0[i] <- patchPop[i] - E0[i] - I0[i] - R0[i]
}

u0 <- data_frame(
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
#
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

# Create an external transfer event to move exposed or infected individuals
# from node 1 (Ystad) to other nodes every seven days

for (t in tspan) {
  list_of_events <- lapply(beginDaysPlague[-1],
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
}

# Event moving individuals from node 1 to the others proportional to the gravitational term
for (i in 2:npatches) {
  for (t in tspan) {
    list_of_events2 <- lapply(beginDaysPlague[-1],
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
}

# Event moving individuals from node 1 to the others proportional to the gravitational term
for (i in 2:npatches) {
  for (t in tspan) {
    list_of_events2 <- lapply(beginDaysPlague[-1],
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
}

# Initialize list_of_events3
list_of_events3 <- list()

# Event moving individuals from node 1 to the others proportional to the gravitational term
for (t in seq(from = 0, to =21, by = 7)) {
  for (i in 1:3) {
    for (j in 1:3) {
      if (i != j) {
        # for (t in 0:maxDays) {
        event <- data.frame(
          event = "extTrans",
          time = t,
          node = i,
          dest = j,
          n = 0,
          proportion = weight_matrix[i, j],
          select = 1,
          shift = 0
        )
        # Append the events to list_of_events3
        list_of_events3 <- c(list_of_events3, event)
      }
    }
  }
}
 

# Combine all the dataframes into one
events <- do.call(rbind, list_of_events2)
View(events)

# Defining mu and beta for each node randomly

local_parameters <- data.frame(beta = (runif(npatches)),
                               mu = runif(npatches))

model <- mparse(
  transitions = transitions,
  compartments = compartments,
  gdata = c(sigma = 0.17 , gamma = 0.4),
  ldata = local_parameters,
  u0 = u0,
  E = E,
  events = events,
  tspan = 1:maxDays + 20
)
set.seed(123)
set_num_threads(1)
result <- run(model = model)

# Cumulative deaths
traj_D <- trajectory(model = result, compartments = "D")
# Show the points per node
View(traj_D)
# Plot a specific trajectory
ggplot(traj_D) + geom_line(aes(x = time, y = D, color = factor(node)))

# Infected
traj_I <- trajectory(model = result, compartments = "I")
# Show the points per node
View(traj_I)
# Plot a specific trajectory
ggplot(traj_I) + geom_line(aes(x = time, y = I, color = factor(node)))

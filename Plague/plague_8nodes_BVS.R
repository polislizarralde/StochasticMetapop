library(SimInf)
library(tools)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)

# scania geometry
scania <- readRDS("data/scania_geometry.RDS")

# local dynamic in each node
u0 <- read_excel("data/data_seven_nodes.xlsx")
u0 <- t(u0)

# node lookup of Ystad etc.
nodes <- c(Ystad = 128, Oja = 136, Hedeskoga = 134, Bjaresjo = 133, Bromma = 135, Stora_Herrestad = 137, Borrie = 140, Stora_Kopinge = 138)
nodes_lookup <- nodes
names(nodes_lookup) <- 1:length(nodes)

transitions <- c("S-> beta*S*I/(S+E+I+R) -> E",
                 "E -> sigma*E -> I",
                 "I -> gamma*(1-mu)*I -> R",
                 "I -> (gamma*mu)*I -> D")

compartments <- c("S", "E", "I", "R", "D")

parameters <- c("beta", "mu")

colnames(u0) <- compartments

tspan <- seq(from = 1, to = 350, by = 1) # Switched to daily.

# EMatrix (#compartments x #events). For each column vector j we put 1 if the compartment
# participated in the event j, otherwise zero.
E <- matrix(c(0, 1, 1, 0, 0), nrow = 5, ncol = 1, dimnames = list(c("S", "E", "I", "R", "D"), c("event1")))

print(E)

# Events
all_possible_events <- readRDS("data/events_distances.RDS")

# Scale the proportion back to something we want.
max(all_possible_events$proportion) # 19594 is the longest distance (in meters) between centroids of parishes
min(all_possible_events$proportion) # 1049 is the shorted distance (in meters) between centroids of parishes

# Inverse square (in km) for gravity.
all_possible_events$proportion <- 1/(all_possible_events$proportion / 1000)^2

# Not much room to scale up, but there is room to scale down.
max(all_possible_events$proportion) # 0.908
min(all_possible_events$proportion) # 0.0026

# Filter on relevant events for this cluster.
relevant_events <- dplyr::filter(all_possible_events,node %in% nodes & dest %in% nodes) %>%
  mutate(node = match(node, nodes_lookup),
         dest = match(dest, nodes_lookup))

# Expand the relevant events to include a time axis
all_relevant_events <- tidyr::crossing(relevant_events, time = tspan) %>% dplyr::relocate(time,.after = event)


# Parameters
local_parameters <- data.frame(beta = c(0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6),
                                 mu = c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7))   

gravity_scaler <- 0.01
events <- all_relevant_events
events$proportion <- events$proportion * gravity_scaler

model <- mparse(transitions = transitions,
                compartments = compartments,
                gdata = c(sigma = 0.4 , gamma = 0.1),
                ldata = local_parameters,
                u0 = u0,
                E = E,
                events = events,
                tspan = tspan)
#set.seed(123)
#set_num_threads(1)
result <- run(model = model)

# One way to plot:
ggplot(trajectory(run(model = model))) + geom_line(aes(x=time,y=I)) + facet_wrap(~node)

# One way to plot: https://www.rdocumentation.org/packages/SimInf/versions/5.1.0/topics/plot-methods
plot(result)

# One way to plot:
traj <- trajectory(result)
long_traj <- tidyr::pivot_longer(traj, cols = c(S, E, I, R, D), names_to = "type", values_to = "count")
ggplot(long_traj) + geom_line(aes(x=time,y=count,color=type)) + facet_wrap(~node)


# Cumulative deaths
traj_D <- trajectory(model=result, compartments = "D")
# Show the points per node
View(traj_D)
# Plot a specific trajectory
ggplot(traj_D) + geom_line(aes(x=time,y=D,color=factor(node)))

# Infected
traj_I <- trajectory(model=result, compartments = "I")
# Show the points per node
View(traj_I)
# Plot a specific trajectory
ggplot(traj_I) + geom_line(aes(x=time,y=I,color=factor(node)))



library(SimInf)
library(tools)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)

# scania geometry
scania <- readRDS("data/scania_geometry.RDS")

# observations

# Ystad (node 1) and Oja (node 2) start in june 1712
# Stora Kopinge (node 8) and Bjaresjo start (node 4) in july 1712
# Bromma (node 5) starts in aug 1712.
# Hedeskoga (node 3) starts in sep 1712.
# Stora Herrestad (node 6) and Borrie (node 7) unknown.

# Roughly:
# Ystad (node 1) and Oja (node 2) start at day 0, and Ystad has about 210 days of plague. Oja about 300 days.
# Stora Kopinge (node 8) and Bjaresjo (node 4) start at day 30. Stora Kopinge is done by day 240. Baresjo is done around day 60.
# Bromma (node 5) starts at day 60. Not visualized on the map anymore after its 'false start' in april.
# Hedeskoga (node 3) starts around day 90, and has about 60 days of plague.


# local dynamic in each node
u0 <- read_excel("data/data_seven_nodes_BVS.xlsx")
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

tspan <- seq(from = 1, to = 300, by = 1) # Switched to daily.

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

# Clean up a little
rm(relevant_events,all_possible_events)

# Parameters
gravity_scaler <- 0.08
events <- all_relevant_events
events$proportion <- events$proportion * gravity_scaler

set_num_threads(1)

model <- mparse(transitions = transitions,
                compartments = compartments,
                gdata = c(beta = 1.1, mu = 0.7, sigma = 0.06 , gamma = 0.5),
                u0 = u0,
                E = E,
                events = events,
                tspan = tspan)

# Manual fitting:
rseed <- sample(1:1000,1)
#set.seed(rseed)
set.seed(569)
result <- run(model = model)
traj <- trajectory(result)
ggplot(dplyr::filter(traj,node %in% c(1,2,3,4,8))) + geom_line(aes(x=time,y=I)) + facet_wrap(~node)

# Roughly what we want to observe:
# Ystad (node 1) starts at day 0 until day 210.
# Oja (node 2) start at day 0, until day 300.
# Bjaresjo (node 4) start at day 30 until day 60.
# Stora Kopinge (node 8) starts at day 30 until day 240.
# Hedeskoga (node 3) starts around day 90 until day 150.

# Manual fitting. Can we manage that without a vector? Sort of seems to work.


       
# Other ways to plot: https://www.rdocumentation.org/packages/SimInf/versions/5.1.0/topics/plot-methods
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



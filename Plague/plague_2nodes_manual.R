library(SimInf)
library(tools)
library(dplyr)
library(ggplot2)

# local dynamic in each node
u0 <- data_two_nodes
u0<- t(u0)

transitions <- c("S-> beta*S*I/(S+E+I+R) -> E",
                 "E -> sigma*E -> I",
                 "I -> gamma*(1-mu)*I -> R",
                 "I -> (gamma*mu)*I -> D"
                 )
compartments <- c("S", "E", "I", "R", "D")
parameters <- c("beta", "mu")

colnames(u0)<-compartments
print(u0)

tspan <- seq(from = 1, to = 400, by = 30)

#EMatrix (#compartments x #events). For each column vector j we put 1 if the compartment
# participated in the event j, otherwise zero.
E <- matrix(c(0, 1, 1, 0, 0), nrow = 5, ncol = 1, dimnames = list(c("S", "E", "I", "R", "D"), c("event1")))

print(E)
# Create an external transfer event to move exposed or infected individuals 
# from node 0 to node 1 at t = 60

infect <- data.frame(event = "extTrans", time = 60, node = 1, dest = 2,
                        n = 0, proportion = 0.25, select = 1, shift = 0)
events <- rbind(infect)

local_parameters <- data.frame(beta = c(0.7, 0.8),
                              mu = c(0.8, 0.4))   

print(local_parameters)


model <- mparse(transitions = transitions,
                compartments = compartments,
                gdata = c(sigma = 0.17 , gamma = 0.4),
                ldata = local_parameters,
                u0 = u0,
                E = E,
                events = events,
                tspan =tspan)
set.seed(123)
set_num_threads(1)
result <- run(model = model)
result

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

# Exposed
traj_E <- trajectory(model=result, compartments = "E")
# Show the points per node
View(traj_E)
# Plot a specific trajectory
ggplot(traj_E) + geom_line(aes(x=time,y=E,color=factor(node)))


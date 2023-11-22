library(SimInf)
library(tools)
library(dplyr)
library(ggplot2)

# local dynamic in each node
u0 <- data_seven_nodes
u0<- t(u0)

transitions <- c("S-> beta*S*I/(S+E+I+R) -> E",
                 "E -> sigma*E -> I",
                 "I -> gamma*(1-mu)*I -> R",
                 "I -> (gamma*mu)*I -> D"
)
compartments <- c("S", "E", "I", "R", "D")
parameters <- c("beta", "mu")

colnames(u0)<-compartments

tspan <- seq(from = 1, to = 350, by = 7)

#EMatrix (#compartments x #events). For each column vector j we put 1 if the compartment
# participated in the event j, otherwise zero.
E <- matrix(c(0, 1, 1, 0, 0), nrow = 5, ncol = 1, dimnames = list(c("S", "E", "I", "R", "D"), c("event1")))

print(E)
# Create an external transfer event to move exposed or infected individuals 
# from node 1 (Ystad) to node 2 (Oja) at t = 15 (15 Jun 1712)
infect_2 <- data.frame(event = "extTrans", time = 15, node = 1, dest = 2,
                     n = 0, proportion = 0.10, select = 1, shift = 0)

# Create an external transfer event to move exposed or infected individuals 
# from node 1 (Ystad) to node 3 (Hedeskoga) at t = 107 (15 Sep 1712)
infect_3 <- data.frame(event = "extTrans", time = 107, node = 1, dest = 3,
                       n = 0, proportion = 0.10, select = 1, shift = 0)

# Create an external transfer event to move exposed or infected individuals 
# from node 1 (Ystad) to node 4 (Bjaresjo) at t = 45 (15 Jul 1712)
infect_4 <- data.frame(event = "extTrans", time = 45, node = 1, dest = 4,
                       n = 0, proportion = 0.10, select = 1, shift = 0)

# Create an external transfer event to move exposed or infected individuals 
# from node 1 (Ystad) to node 5 (Bromma) at t = 76 (15 Aug 1712)
infect_5 <- data.frame(event = "extTrans", time = 76, node = 1, dest = 5,
                       n = 0, proportion = 0.10, select = 1, shift = 0)

# Create an external transfer event to move exposed or infected individuals 
# from node 1 (Ystad) to node 8 (Stora Kopinge) at t = 20 (20 Jun 1712)
infect_8 <- data.frame(event = "extTrans", time = 20, node = 1, dest = 8,
                       n = 0, proportion = 0.10, select = 1, shift = 0)

# I don't create events for move individuals from node 1 (Ystad) to node 6 and 7
# (Stora Herrestad and Borrie) because these parishes weren't affected by the plague

events <- rbind(infect_2, infect_3, infect_4, infect_5, infect_8)

local_parameters <- data.frame(beta = c(0.6, 0.5, 0.4, 0.4, 0.4, 0.1, 0.1, 0.4),
                                 mu = c(0.8, 0.8, 0.7, 0.7, 0.6, 0.2, 0.2, 0.5)
                               )   

model <- mparse(transitions = transitions,
                compartments = compartments,
                gdata = c(sigma = 0.17 , gamma = 0.4),
                ldata = local_parameters,
                u0 = u0,
                E = E,
                events = events,
                tspan = 1:400)
set.seed(123)
set_num_threads(1)
result <- run(model = model)

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



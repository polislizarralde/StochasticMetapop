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

#EMatrix
E <- matrix(c(0, 0, 1, 0, 0), nrow = 5, ncol = 1, dimnames = list(c("S", "E", "I", "R", "D"), c("1")))


# Create one enter event to introduce an infected individuals to the 2nd node at t = 50
infect <- data.frame(event = "enter", time = 50, node = 1, dest = 2,
                       n = 0, proportion = 0.25, select = 3, shift = 0)
exposed <- data.frame(event = "enter", time = 10, node = 1, dest = 2,
                      n = 0, proportion = 0.1, select = 2, shift = 0)
events <- rbind(infect, exposed)

local_parameters <- data.frame(beta = c(0.7, 0.8),
                              mu = c(0.8, 0.4))   

print(local_parameters)

# n represents the number of nodes, on this case each node corresponds to one parish
n <- 2

model <- mparse(transitions = transitions,
                compartments = compartments,
                gdata = c(sigma = 0.17 , gamma = 0.4),
                ldata = local_parameters,
                u0 = u0,
                events <- events,
                E = E,
                tspan = 1:250)


set.seed(123)
set_num_threads(1)
result <- run(model = model)

traj_D <- trajectory(model=result, compartments = "D")
# Show the points per node
View(traj_D)
# Plot a specific trajectory
ggplot(traj_D) + geom_line(aes(x=time,y=D,color=factor(node)))


traj_I <- trajectory(model=result, compartments = "I")
# Show the points per node
View(traj_I)
# Plot a specific trajectory
ggplot(traj_I) + geom_line(aes(x=time,y=I,color=factor(node)))


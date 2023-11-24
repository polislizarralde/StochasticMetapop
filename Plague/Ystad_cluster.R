library(SimInf)
library(tools)
library(dplyr)
library(ggplot2)

library(readr) # for read csv files

data_cluster <- Ystad_cluster
beginDaysPlague <- data_cluster$BeginDaysPlague
endDaysPlague <- data_cluster$EndDaysPlague

# local dynamic in each node
u0 <- init_cond_Ystad_cluster
parish_names <- colnames(u0)

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

# Create an external transfer event to move exposed or infected individuals
# from node 1 (Ystad) to other nodes every seven days

for(t in 0:213){list_of_events <- lapply(beginDaysPlague[-1],
                         function(t){
                         infect <- data.frame(event = "extTrans",
                                              time = t+7,
                                              node = 1,
                                              dest = 2:11,
                                              n = 0,
                                              proportion = 0.05,
                                              select = 1,
                                              shift = 0
                                              )
                                      }
                         )}

# Combine all the dataframes into one
events <- do.call(rbind, list_of_events)

# Number of nodes
n <- 11

# Defining mu and beta for each node randomly

local_parameters <- data.frame(beta = (runif(n)),
                               mu = runif(n)
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



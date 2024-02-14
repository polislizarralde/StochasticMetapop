library(groundhog)
library(lubridate)
pkgs <- c("conflicted","SimInf","ggplot2","dplyr","mlrMBO")
groundhog.library(pkgs, "2023-04-30")
conflicts_prefer(dplyr::filter)

# set the directory to find the csv in the folder data
setwd("~/Documents/GitHub/StochasticMetapop/Plague")

source("/Users/dianapli/Documents/GitHub/StochasticMetapop/Plague/Number_Infected_Parishes.R")

YSTAD_group <- read_csv("/Users/dianapli/Desktop/PythonMathematicalModeling/docs/PlagueProject/data/private/yellow_group.csv", col_types = cols())
# parish names without repetitions
parish_names <- unique(YSTAD_group$ParishName)
npatches <- length(parish_names)
patchPop <- YSTAD_group$BEF1699[1:npatches]
maxDays <- max(YSTAD_group$EndDaysPlague)

# generating the initial conditions for the model
S0 <- rep(0, npatches) 
Ilow0 <- rep(0, npatches)
Ihigh0 <- rep(0, npatches)
R0 <- rep(0, npatches)
Dcum0 <- rep(0, npatches)
Sv0 <- rep(0, npatches)
Iv0 <- rep(0, npatches)
Iv0[1] <- 45.0
Ilow0[1] <- 1.0
Ihigh0[1] <- 1.0

for (i in 1:npatches) {
  S0[i] <- patchPop[i] - Ilow0[i] - Ihigh0[i] - R0[i]
}

for (i in 1:npatches) {
  Sv0[i] <- patchPop[i]
}

u0 <- data.frame(
  S0 = S0,
  Ilow0 = Ilow0,
  Ihigh0 = Ihigh0,
  R0 = R0,
  Dcum0 = D0,
  Sv0 = Sv0,
  Iv0 = Iv0
)
u0

# model
transitions <- c("S -> beta_v * Iv * S / (S+Ilow+Ihigh+R) -> Ilow", 
                 "S -> beta_p * Ihigh * S / (S+Ilow+Ihigh+R) -> Ihigh",# New equation
                 "Ilow -> (1 - mu) * sigma * Ilow -> R",
                 "Ilow -> mu * sigma * Ilow -> Ihigh",
                 "Ihigh -> gamma * Ihigh -> Dcum",
                 "@ -> growth_v * Sv * (1 - (Sv + Iv)/((S+Ilow+Ihigh+R) * index_ave)) -> Sv",
                 "Sv -> growth_v * Sv * (1 - (Sv + Iv)/((S+Ilow+Ihigh+R) * index_ave)) -> @", # @ 
                 "Sv -> (beta_low * Ilow + beta_high * Ihigh) * Sv / (S+Ilow+Ihigh+R) -> Iv",
                 "Iv -> gamma_v * Iv -> @")

compartments <- c("S","Ilow","Ihigh","R","Dcum","Sv","Iv")
colnames(u0) <- compartments

# EMatrix (#compartments x #events). For each column vector j we put 1 if the compartment
# participated in the event j, otherwise zero.
E <- matrix(c(0,0,0,0,0,0,1), nrow = 7, ncol = 1
            , dimnames = list(c("S", "Ilow", "Ihigh", "R", "Dcum","Sv","Iv")
                              , c("event1")))

connection_matrix <- read_csv("/Users/dianapli/Desktop/PythonMathematicalModeling/docs/PlagueProject/data/private/connection_Yellow_group.csv"
                              , col_types = cols())

connection_matrix <- as.matrix(connection_matrix)

global_parameters <- c(index_ave = 23, growth_v = 0.14, beta_v = 0.09
                       , gamma_v = 0.71, beta_low = 0.021, beta_high = 0.66
                       , beta_p = 0.24, sigma = 0.20, gamma = 0.84
                       , mu = 0.46)

# define the local parameters per parish as a data frame
local_parameter <- data.frame(Iv_prop = rep(0.5, npatches))

gdf <- YSTAD_group
list_of_events <- list()
for(i in 1:npatches){
  initial_time <- gdf$BeginDaysPlague[i]
  infect_event <- lapply(seq(from = initial_time + 7, to = maxDays, by = 7),
                         function(t) {
                           infect <- data.frame(
                             event = "extTrans",
                             time = t,
                             node = i,
                             dest = setdiff(1:npatches, i),
                             n = 0,
                             proportion = local_parameter$Iv_prop[i] * connection_matrix[i, setdiff(1:npatches, i)],
                             select = 1,
                             shift = 0
                           )
                         })
  list_of_events[[i]] <- infect_event
}
events <- do.call(rbind, unlist(list_of_events, recursive = FALSE))

# Defining the model
model <- mparse(transitions = transitions,
                compartments = compartments,
                gdata = global_parameters,
                ldata = local_parameter,
                u0 = u0,
                tspan = 1:maxDays,
                events = events,
                E = E
)
set.seed(123)
set_num_threads(1)
result <- run(model = model)

# Cumulative deaths
traj_D <- trajectory(model = result, compartments = "Dcum")
ggplot(traj_D) + geom_line(aes(x = time, y = Dcum, color = factor(node)))

# Infected
traj_I <- trajectory(model = result, compartments = "Ilow")
# Show the points per node
View(traj_I)
# Plot a specific trajectory
ggplot(traj_I) + geom_line(aes(x = time, y = Ilow, color = factor(node)))

library(groundhog)
library(lubridate)
pkgs <- c("conflicted","SimInf","ggplot2","dplyr","mlrMBO")
groundhog.library(pkgs, "2023-04-30")
conflicts_prefer(dplyr::filter)

# set the directory to find the csv in the folder data
setwd("~/Documents/GitHub/StochasticMetapop/Plague")

source("/Users/dianapli/Documents/GitHub/StochasticMetapop/Plague/Number_Infected_Parishes.R")

infec_south_parishes <- read.csv("/Users/dianapli/Desktop/PythonMathematicalModeling/docs/PlagueProject/data/private/infectedSouthParishes.csv")
# Testing the function with a small example
# Filter rows where ParishName is 'YSTAD', 'ÖJA' or 'HEDESKOGA'
example1 <- infec_south_parishes[infec_south_parishes$ParishName %in% c('YSTAD', 'ÖJA', 'HEDESKOGA'), ]
# example1 <- infec_south_parishes[infec_south_parishes$ParishName %in% c('YSTAD'), ]
parish_names <- example1$ParishName
patchPop <- example1$BEF1699
cum_deaths <- example1$VictimsNumber
maxDays <- 303 # Value takes from python

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

# SEIRD model.
transitions <- c(
  "S-> beta*S*I/(S+E+I+R) -> E",
  "E -> sigma*E -> I",
  "I -> gamma*(1-mu)*I -> R",
  "I -> (gamma*mu)*I -> D"
)
compartments <- c("S", "E", "I", "R", "D")
parameters <- c("beta", "mu")
colnames(u0) <- compartments

# EMatrix (#compartments x #events). For each column vector j we put 1 if the compartment
# participated in the event j, otherwise zero.
E <- matrix(
  c(0, 1, 1, 0, 0),
  nrow = 5,
  ncol = 1,
  dimnames = list(c("S", "E", "I", "R", "D"),
                  c("event1"))
)

# Create an external transfer event to move exposed or infected individuals
# from node 1 (Ystad) to other nodes starting seven days after the initial plague
# period in Ystad and every seven days until proportional to the gravity term

data_init_node <- subset(example1, ParishName == "YSTAD")
initial_time <- 0
grav_matrix <- matrix(c(0.        , 0.0533247 , 0.02542587,
                        0.0533247 , 0.        , 0.00226749,
                        0.02542587, 0.00226749, 0.        ), 
                      nrow = 3, 
                      ncol = 3, 
                      byrow = TRUE)
colnames(grav_matrix) <- c('YSTAD', 'ÖJA', 'HEDESKOGA' )

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
               proportion = grav_matrix[i,"YSTAD"],
               select = 1,
               shift = 0
             )
           })
  # Append the events to list_of_events2
  # list_of_events2 is a list of lists of dataframes
  list_of_events2[[i]] <- infect_event2
}
events <- do.call(rbind, unlist(list_of_events2, recursive = FALSE))


# Function to calculate the error in the cumulative number of infected parishes per month between the model and the data
objectiveFunction_2 <- function(parameters, gdf, n) {
  
  # Defining the model
  model_SEIRD <- mparse(transitions = transitions,
                        compartments = compartments,
                        gdata = c(sigma = 0.17 , gamma = 0.4),
                        ldata = parameters,
                        u0 = u0,
                        tspan = 1:maxDays
                        # ,
                        # events = events,
                        # E = E
  )
  result <- run(model = model_SEIRD)
  traj_D <- trajectory(model = result, compartments = "D")
  
  # Defining the initial date of the gdf to start counting the number of infected parishes per month
  date <- min(gdf$BeginPlaguePeriod, na.rm = TRUE)
  
  # Getting the number of infected parishes per month from the data
  cum_infected_parishes_by_month <- count_infected_by_month(gdf,date,n)
  
  # Initializing the number of infected parishes per month for the model's output
  infected_parishes <- rep(0, length(cum_infected_parishes_by_month))
  
  # Initializing the error between the model's output and the data
  error <- rep(0, length(cum_infected_parishes_by_month))
  
  # Computing the total number of parishes in the dataframe without repetitions
  total_parishes <- length(gdf$ParishName)
  
  # Computing the number of infected parishes per month from the model's output
  for (i in 1:length(cum_infected_parishes_by_month)) {
    init_days <- cum_infected_parishes_by_month[i,'DaysFromInitDate']
    final_days <- cum_infected_parishes_by_month[i,'DaysToEndOfMonth']
    
    for (k in 1:total_parishes) {
      for (day in init_days:final_days) {
        # Check if there are any rows that satisfy the condition
        rows <- traj_D[traj_D$node == k & traj_D$time == day, ]
        if (nrow(rows) == 0){
          break # Breaks the innermost loop when the condition is met
        }
        else if(nrow(rows) > 0){
          if(rows$D >= 1){
            infected_parishes[i] <- infected_parishes[i] + 1
            break # Breaks the innermost loop when the condition is met
          }
        }
      }
    }
    error[i] <- (infected_parishes[i] - cum_infected_parishes_by_month[i,'CumInfectParishes'])^2
  }
  
  # Computing the error between the model's output and the data
  total_error <- (1/length(cum_infected_parishes_by_month)) * (1/total_parishes)^2 * sum(error)
  
  return(total_error)
}


# parameters
# objectiveFunction_2(parameters = parameters, gdf = example1, npar = 2, n = 0)
# total_parishes <- length(example1$ParishName)
# total_parishes
date <- min(example1$BeginPlaguePeriod, na.rm = TRUE)
cum_infected_parishes_by_month <- count_infected_by_month(example1,date,0)
date
cum_infected_parishes_by_month

ncopies<-20

ps<- makeParamSet(
  makeNumericVectorParam("beta", len=npatches, lower = 0, upper = 1),
  makeNumericVectorParam("mu",len=npatches, lower = 0, upper = 1)
)
ps

# Generate an initial design
des <- generateDesign(n = ncopies, par.set = ps)

#des$y = apply(des,1,objectiveFunction_2_mlr)

# Define the objective function for mlrMBO
objectiveFunction_2_mlr <- makeSingleObjectiveFunction(
  name = "My Function",
  noisy = TRUE,
  has.simple.signature= TRUE,
  fn = function(xs){
    
    beta <- rep(0, npatches)
    K <- 1
    for (I in 1:npatches){
      beta[K] <- xs[I]
      K <- K + 1
    }
    mu <- rep(0, npatches)
    K <- 1
    for (I in (npatches+1):(2*npatches)){
      mu[K] <- xs[I]
      K <- K + 1
    }
    
    # now, make a df such that the columns are the vectors beta and mu
    beta <- as.data.frame(beta)
    mu <- as.data.frame(mu)
    beta_mu <- cbind(beta, mu)
    colnames(beta_mu) <- c("beta", "mu")
    print("beta_mu")
    print(beta_mu)
    
    return (objectiveFunction_2(parameters = beta_mu, gdf = example1, n = 0))},
  par.set = ps,
  minimize = TRUE
)


control = makeMBOControl()
control = setMBOControlTermination(control, iters = 50)
control = setMBOControlInfill(control, crit = makeMBOInfillCritEI())  # Not sure if this is the default, or something Avelda added.

# Run the optimization
res <- mbo(fun = objectiveFunction_2_mlr, design = des, control = control, show.info = TRUE)

# Extract the best parameters
best_params <- res$x

# Plot the cumulative number of parishes per month from the data
#plot_cum_infected_parishes_month(example1, "JUN 1712", 0)

# Plot the cumulative number of parishes per month using the model estimation
plot_cum_infect_parishes_model(res$x, example1,0)

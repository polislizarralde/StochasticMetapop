library(SimInf)
library(tools)
library(dplyr)
library(ggplot2)
library(optimx)
library(gbutils)
library(readr) # for read csv files

# set the directory to find the csv in the folder data
#setwd("~/Documents/GitHub/StochasticMetapop/Plague")

data_cluster <- read_csv("data/Ystad_group.csv")
beginDaysPlague <- data_cluster$BeginDaysPlague
endDaysPlague <- data_cluster$EndDaysPlague
parish_names <- data_cluster$ParishName
patchPop <- data_cluster$BEF1699
maxDays <- endDaysPlague[1]

# data
class_plague <- data_cluster$plague
cum_deaths <- data_cluster$VictimsNumber

# generating the initial conditions for the model
npatches <- 1
E0 <- 0.0
R0 <- 0.0
D0 <- 0.0
I0 <- 1.0
S0 <- patchPop[1] - E0 - I0 - R0

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
tspan <- seq(from = 1, to = maxDays + 20, by = 1)


mat <- matrix(c(0.5,0.3),nrow=1,ncol=2,byrow=TRUE)
# Defining the model considering some events
model <- mparse(
  transitions = transitions,
  compartments = compartments,
  gdata = c(sigma = 0.17 , gamma = 0.4),
  ldata = data.frame(beta = as.numeric(mat[1,1]), mu = as.numeric(mat[1,2])),
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


########### Optimization process

npatches <- 1
npar <-2 # parameters to estimate
objectiveFunction <- function(parameters){
  mat <- matrix(parameters,nrow=npatches,ncol=npar,byrow=TRUE)
  # Update local parameters of the model
  model <- mparse(
    transitions = transitions,
    compartments = compartments,
    gdata = c(sigma = 0.17 , gamma = 0.4),
    ldata = data.frame(beta = as.numeric(mat[npatches,1]), mu = as.numeric(mat[npatches,2])),
    u0 = u0,
    tspan = 1:maxDays + 20
  )
  # Solve the model
  sol <- run(model)
  
  # Compute errors per trajectory
  errors <- sapply(1:npatches, function(i){
    initial_position <- beginDaysPlague[i]
    final_position <- endDaysPlague[i]
    deaths <- cum_deaths[i]
    traj_D <- trajectory(model = sol, compartments = "D")
    traj_D_node <- dplyr::filter(traj_D, node == i)
    
    if (!isNA(initial_position)) {
      if (final_position != 0 && deaths != 0) {
        
        print("valores case 1")
        print(as.double(traj_D_node$D[initial_position]))
        print(traj_D_node$D[final_position])
        print(deaths)
        rrr <- abs(traj_D_node$D[initial_position]*1.0 - 1.0) +
          abs(traj_D_node$D[final_position]*1.0 - deaths*1.0)
        print("rrr")
        print(rrr)
        return(rrr)
      }
      else{
        print("valores case 2")
        print(traj_D_node$D[initial_position])
        return (abs(traj_D_node$D[initial_position]-1.0))
      }
    }
    print("valor 0 case 3")
    return(0)
  })
  print("errors:=")
  print(errors)
  
  total_error <- sum(unlist(errors))
  return(total_error)
}

# Solve the optimization problem
initial_parameters <- rep(c(0.5, 0.5), npatches)


results <- optimx::optimx(initial_parameters,
                          objectiveFunction,
                          method = "L-BFGS-B",
                          lower = rep(c(0,0), npatches),
                          upper = rep(c(1,1), npatches)
)




# Print the result 

print(results)
# cat("error =", results$fvalues, "\n")
# cat("beta1 =", results$par[1], "\n")
# cat("mu1 =", results$par[2], "\n")
# cat("beta2 =", results$par[3], "\n")
# cat("mu2 =", results$par[4], "\n")


# Run the model with best parameters
# for(i in 1:npatches){
#   model@ldata[1,i] <- results$par[2*i-1] #beta
#   model@ldata[2,i] <- results$par[2*i]   #mu
# }
# 
# final_sol <- run(model)
# 
# # Plotting
# par(mfrow=c(ceiling(npatches/2), 2)) # Adjust layout according to number of nodes
# 
# for(i in 1:npatches) {
#   traj_D_node <- dplyr::filter(trajectory(model = final_sol, compartments = "D"), node == i)
#   plot(traj_D_node$time, traj_D_node$D, main=paste("Node", i), xlab="Time", ylab="D")
# }






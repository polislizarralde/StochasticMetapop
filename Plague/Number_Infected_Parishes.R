library(ggplot2)
library(lubridate)
library(dplyr)
library(readr) # for read csv files

# Function to count the number of parishes affected by the plague per month
count_infected_parishes_month <- function(df, date, n = 0
                                    , column_name = 'ParishName'
                                    , start_date = 'BeginPlaguePeriod'
                                    , end_date = 'EndPlaguePeriod') {
  
  # Convert your date columns to datetime format
  df[[start_date]] <- parse_date_time(df[[start_date]], orders = "my")

  df[[end_date]] <- parse_date_time(df[[end_date]], orders = "my", quiet = TRUE)

  # Replace NA with corresponding date in start_date column plus n months
  df[[end_date]][is.na(df[[end_date]])] <- df[[start_date]][is.na(df[[end_date]])] %m+% months(n)

  # Convert your date to datetime format
  date <- parse_date_time(date, orders = "my")

  # Add the converted date to a new column in df
  df$ConvertedDate <- date

  # Check if there is at least one non-NA date in end_date column
  if (any(!is.na(df[[end_date]]))) {
    # Define the range of dates
    dates <- seq(from = date, to = max(df[[end_date]], na.rm = TRUE), by = "month")
    
    # Create a unique identifier combining Parish and date ranges
    df$UniqueID <- paste0(df[[column_name]], '_', df[[start_date]], '_', df[[end_date]])

    # Create a dataframe to store the results
    results <- data.frame(Month = dates,
                          DaysFromInitDate = as.numeric(dates - min(df[[start_date]]), units="days"),
                          NumberInfectedParishes = 0,
                          CumInfectParishes = 0,
                          EndOfMonth = dates + months(1) - days(1))
    
    # Initialize an empty list to store the sets of infected parishes
    infected_parishes <- vector("list", length(dates))

    # Iterate over the dates
    for (i in 1:length(dates)) {
      date <- dates[i]
      # Count nodes where infection start date is before or on the given date
      # and either there is no end date or the end date is after the given date
      infected_nodes <- df[(df[[start_date]] <= date) & (df[[end_date]] >= date), ]
      # Store the results
      results[results$Month == date, 'NumberInfectedParishes'] <- length(unique(infected_nodes$UniqueID))
      
      # Add the set of infected parishes to the list
      infected_parishes[[i]] <- unique(infected_nodes[[column_name]])
    }
    
    # Initialize cum_infect_parishes as a vector filled with 0s and the same length as dates
    cum_infect_parishes <- rep(0, length(dates))
    
    # Check if the first set of infected parishes is non-empty
    if (length(infected_parishes[[1]]) > 0) {
      cum_infect_parishes[1] <- length(infected_parishes[[1]])
      # Defining a variable to store the union of the infected parishes
      union_infected_parishes <- unique(infected_parishes[[1]])
    } else {
      union_infected_parishes <- c()  # Empty vector in R
    }

  } else {
    stop("All values in the end_date column are NA. Cannot create a sequence of dates.")
  }
  # Add a new column to count the days from the initial date to the end of the month
  results$DaysToEndOfMonth <- as.integer(results$EndOfMonth - min(df[[start_date]]))
  
  # Add a new column with the sets of infected parishes
  results$InfectedParishes <- infected_parishes
  
  # Calculate the cumulative number of infected parishes by month using the sets
  cum_infect_parishes <- as.integer(length(dates))
  cum_infect_parishes[1] <- length(infected_parishes[[1]])
  
  for (i in 2:length(infected_parishes)) {
    if (length(intersect(infected_parishes[[i]], infected_parishes[[i-1]])) == 0) {
      cum_infect_parishes[i] <- cum_infect_parishes[i-1]
    } else {
      cum_infect_parishes[i] <- cum_infect_parishes[i-1] + length(setdiff(infected_parishes[[i]], infected_parishes[[i-1]]))
    }
  }
  
  # Add a new column with the cumulative number of infected parishes
  results$CumInfectParishes <- cum_infect_parishes  

  return(results)
}


plot_parishes_month <- function(df, date, n = 0, column_name = 'ParishName'
                                , start_date = 'BeginPlaguePeriod'
                                , end_date = 'EndPlaguePeriod') {
  results <- count_infected_parishes_month(df, date, n, column_name, start_date, end_date)
  
  # Ensure 'date' exists in 'results' and not 'df'
  if(!"Month" %in% names(results)) {
    stop("Error: 'date' column not found in 'results'")
  }
  
  results$Month <- as.Date(results$Month, format="%b %Y")
  
  ggplot(results, aes(x=Month, y=NumberInfectedParishes)) +
    geom_line(color='orange') +
    scale_x_date(labels = scales::date_format("%b %Y")) +
    labs(x='Month', y='No. infected parishes', title='South Scania') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


# Plot the results
run_infected_parishes_model <- function(best_params, global_parameters, gdf, n = 0) {
  # Create a data frame with the best parameters
  df_best_params <- as.data.frame(best_params)
  
  # Evaluate the model in the best parameters
  model <- mparse(transitions = transitions,
                 compartments = compartments,
                 gdata = global_parameters,
                 ldata = df_best_params,
                 u0 = u0,
                 tspan = 1:maxDays,
                 events = events,
                 E = E
  )
  result <- run(model = model)
  traj_D <- trajectory(model = result, compartments = "Dcum")
  
  # Defining the initial date of the gdf to start counting the number of infected parishes per month
  # date <- min(gdf$BeginPlaguePeriod, na.rm = TRUE) # This is not giving the correct date
  date <- gdf$BeginPlaguePeriod[1] # Initial. works if gdf is sorted by BeginPlaguePeriod
  
  # Getting the number of infected parishes per month from the data
  cum_infected_parishes_by_month <- count_infected_parishes_month(gdf,date,n)
  
  # Initializing the number of infected parishes per month for the model's output
  infected_parishes <- rep(0, length(cum_infected_parishes_by_month$Month))
  
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
  }

  return(infected_parishes)
}


plot_infected_parishes <- function(infected_parishes, gdf, n = 0){
  date <- gdf$BeginPlaguePeriod[1] # Initial. works if gdf is sorted by BeginPlaguePeriod
  # Getting the number of infected parishes per month from the data
  cum_infected_parishes_by_month <- count_infected_parishes_month(gdf,date,n)
  
  # Create a new dataframe to store infected parishes and date
  df <- data.frame(
    "month" = cum_infected_parishes_by_month$Month,
    "infected_parishes" = infected_parishes
  )
  # Plot the results
  df$month <- as.Date(df$month, format="%b %Y")
  
  ggplot(df, aes(x=month)) +
    # Model estimation
    geom_line(aes(y=infected_parishes), color='blue') +
    # Real data
    geom_line(aes(y=cum_infected_parishes_by_month$NumberInfectedParishes), color='orange') +
    scale_x_date(labels = scales::date_format("%b %Y")) +
    #scale_y_continuous(limits = c(0,total_parishes))+
    labs(x='Month', y='No. infected parishes', title='South Scania') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
}

# library(ggplot2)
# 
# plot_experiment <- function(iterations, best_params, gdf, n) {
#   # run iterations times the model and store the results in a list
#   # and make another function to plot the results
#   results <- vector("list", iterations)
#   colors_iterate <- rainbow(iterations)
#   
#   # Create an empty data frame to store the results
#   df <- data.frame(month = numeric(0), value = numeric(0))
#   cum_infected_parishes_by_month <- count_infected_parishes_month(gdf, date, n)
#   
#   # Ensure that the 'month' column is numeric and sorted
#   cum_infected_parishes_by_month$month <- as.numeric(cum_infected_parishes_by_month$month)
#   cum_infected_parishes_by_month <- cum_infected_parishes_by_month[order(cum_infected_parishes_by_month$month), ]
#   
#   for (i in 1:iterations) {
#     result <- run_infected_parishes_model(best_params, gdf, n)
#     results[[i]] <- result
#     df <- rbind(df, data.frame(month = seq_along(result), value = result))
#   }
#   
#   # Check if 'df' is empty before creating the plot
#   if (nrow(df) == 0) {
#     warning("No data to plot.")
#     return(NULL)
#   }
#   
#   # Create the plot
#   plot <- ggplot(df, aes(x = month, y = value, color = 'Model')) +
#     geom_line() +
#     scale_x_continuous(breaks = seq(1, nrow(df), length.out = 6),
#                        labels = seq(1, nrow(df), length.out = 6)) +
#     labs(x = 'Month', y = 'No. infected parishes', title = 'South Scania') +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))
#   
#   # Add the observed line only if there is data
#   if (!is.null(cum_infected_parishes_by_month)) {
#     plot <- plot + geom_line(data = cum_infected_parishes_by_month,
#                              aes(x = month, y = NumberInfectedParishes),
#                              color = 'orange', linetype = 'dashed')
#   }
#   return(plot)
# }


###############################
###############################
#### Count deaths per month ####
count_deaths_month <- function(gdf, column_name = 'ParishName'
                                   , begin_date = 'BeginPlaguePeriod'
                                   , victims_column = 'VictimsNumber'
                                   , end_date = 'EndPlaguePeriod'
                                   , pop_size = 'BEF1699') {

  # Create a copy of the dataframe
  gdf_copy <- gdf

  # Filter the dataframe to get only the infected parishes
  gdf_copy <- gdf_copy[!is.na(gdf_copy[[begin_date]]), ]

  # Fixing the type of the column to iterate over
  gdf_copy$new_format_BeginPlaguePeriod <- as.Date(gdf_copy$new_format_BeginPlaguePeriod, format="%Y-%m-%d")
  gdf_copy$new_format_EndPlaguePeriod <- as.Date(gdf_copy$new_format_EndPlaguePeriod, format="%Y-%m-%d")

  # Sort by the column
  gdf_copy <- gdf_copy[order(gdf_copy$new_format_BeginPlaguePeriod),]

  # adding a column with the number of days from the beginning of the plague
  gdf_copy$BeginDaysPlague <- as.integer(gdf_copy$new_format_BeginPlaguePeriod - min(gdf_copy$new_format_BeginPlaguePeriod))

  # Adding a column with the number of days to the end of the plague
  gdf_copy$EndDaysPlague <- sapply(1:nrow(gdf_copy), function(i) {
    if (!is.na(as.Date(gdf_copy$new_format_EndPlaguePeriod[i]))) {
      return(as.integer(as.Date(gdf_copy$new_format_EndPlaguePeriod[i])
                        - as.Date(gdf_copy$new_format_BeginPlaguePeriod[1])))
    } else {
      return(NA)
    }
  })
  
  # Fix the type of pop_size to integer 
  gdf_copy[[pop_size]] <- as.integer(gdf_copy[[pop_size]])
  
  # Get the gdf sorted by the end of the plague period
  gdf_copy <- gdf_copy[order(gdf_copy$new_format_EndPlaguePeriod),]
   
  # Get the unique dates
  months <- unique(gdf_copy$new_format_EndPlaguePeriod)
  days <- unique(gdf_copy$EndDaysPlague)
  
  # Filter out the NA values from months
  months <- months[!is.na(months)]
  days <- days[!is.na(days)]

  # Create a dataframe to store the results
  results <- data.frame(EndMonth = months,
                        CumDays = days,
                        NumberDeaths = 0,
                        CumDeaths = 0,
                        CumPop = 0,
                        Parishes = "")

  # Initialize total_deaths
  total_deaths <- 0
  
  # Iterate over the dates
  for (date in months) {
    
    # Subset rows where new_format_EndPlaguePeriod is equal to date
    subset_df <- gdf_copy[gdf_copy$new_format_EndPlaguePeriod == date, ]
    
    # Calculate number of deaths
    numberOfDeaths <- sum(subset_df[[victims_column]], na.rm = TRUE)
    
    # Get parishes
    parishes <- paste(gdf_copy[gdf_copy$new_format_EndPlaguePeriod <= date, column_name], collapse = ",")
    
    # Update results dataframe
    results[results$EndMonth == date, "Parishes"] <- parishes
    results[results$EndMonth == date, "NumberDeaths"] <- numberOfDeaths
    
    # Update total_deaths
    total_deaths <- total_deaths + numberOfDeaths
    results[results$EndMonth == date, "CumDeaths"] <- total_deaths
    
    # Update CumPop
    results[results$EndMonth == date, "CumPop"] <- sum(gdf_copy[gdf_copy$new_format_EndPlaguePeriod <= date, pop_size], na.rm = TRUE)
  }
  return(results)
}

####################
####################
####Plot the results
plot_cum_deaths_model <- function(best_params, gdf) {
  
  df_best_params <- as.data.frame(best_params)
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
                               proportion = df_best_params$Iv_prop[i] * connection_matrix[i, setdiff(1:npatches, i)],
                               select = 1,
                               shift = 0
                             )
                           })
    list_of_events[[i]] <- infect_event
  }
  events <- do.call(rbind, unlist(list_of_events, recursive = FALSE))
  
  # Evaluate the model in the best parameters
  model <- mparse(transitions = transitions,
                  compartments = compartments,
                  gdata = global_parameters,
                  ldata = df_best_params,
                  u0 = u0,
                  tspan = 1:maxDays,
                  events = events,
                  E = E
  )
  result <- run(model = model)
  traj_D <- trajectory(model = result, compartments = "Dcum")
  
  # Initializing the cum. number of deaths per month for the model's output
  model_cum_deaths <- rep(0, maxDays)
  
  # Computing the cum. number of deaths per month from the model's output
  for (i in 1:maxDays) {
    for (k in 1:npatches){
      model_cum_deaths[i] <- model_cum_deaths[i] + (traj_D[traj_D$node == k & traj_D$time ==i, ]$D)
    }
  }
  
  # Plot the results
  model_data <- data.frame(
    # column from 1 to maxDays
    "days" = seq(1, maxDays),
    "model_cum_deaths" = model_cum_deaths
  )
  
  # Getting the cumulative number of deaths per month from the data
  cum_deaths_month <- count_deaths_month(gdf)
  ggplot(model_data, aes(x = days, y = model_cum_deaths)) +
    geom_line(color = "blue") +
    geom_point(data = cum_deaths_month, aes(x = CumDays, y = NumberDeaths), color = "orange") +
    #scale_x_date(labels = scales::date_format("%b %Y")) +
    labs(x = "Time", y = "Cum. Deaths", title = "South Scania") 
  #+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # ggplot(df, aes(x=days)) +
  #   # Model estimation
  #   geom_line(aes(y=model_cum_deaths), color='blue') +
  #   # Real data
  #   geom_point(aes(y=cum_deaths_month$CumDeaths), color='orange') +
  #   #scale_x_date(labels = scales::date_format("%b %Y")) +
  #   #scale_y_continuous(limits = c(0,total_parishes))+
  #   labs(x='Time', y='No. deaths', title='South Scania') +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Import the csv file
# infectedParishes <- read.csv("/Users/dianapli/Desktop/PythonMathematicalModeling/docs/PlagueProject/data/private/infectedSouthParishes.csv")
# # Testing the function with a small example
# # Filter rows where ParishName is 'YSTAD', 'ÖJA' or 'HEDESKOGA'
# example1 <- infectedParishes[infectedParishes$ParishName %in% c('YSTAD', 'ÖJA', 'HEDESKOGA', 'ELJARÖD'), ]
# # Count the number of infected parishes per month
# count_infected_parishes_month(example1, "JUN 1712", 0)
# count_deaths_by_month(example1)
# 
# # Plot the results
# plot_parishes_month(example1, "JUN 1712", 0)
# plot_cum_deaths_model(c(0.1, 0.2), example1)

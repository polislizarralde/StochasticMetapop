library(ggplot2)
library(lubridate)
library(dplyr)
library(readr) # for read csv files


# Function to count the number of parishes affected by the plague per month
count_infected_parishes_month <- function(df, date, n
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
    results <- data.frame(date = dates,
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
      results[results$date == date, 'NumberInfectedParishes'] <- length(unique(infected_nodes$UniqueID))
      
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

### This doesn't work yet
# Function to count the number of deaths per month
# count_victims_by_month <- function(df, column_name = 'ParishName', begin_date = 'BeginPlaguePeriod',
#                                    victims_column = 'VictimsNumber', end_date = 'EndPlaguePeriod', pop_size = 'BEF1699') {
#   
#   # Create a copy of the dataframe
#   df_copy <- df
#   
#   # Filter the dataframe to get only the infected parishes
#   df_copy <- filter(df_copy, !is.na(df_copy[[begin_date]]))
#   
#   # Define new column names
#   new_format_BeginPlaguePeriod <- paste0(begin_date, "_new")
#   new_format_EndPlaguePeriod <- paste0(end_date, "_new")
#   
#   # Convert date columns to datetime format such that for begin_date is the first day of the month and for end_date is the last day of the month
#   df_copy[[new_format_BeginPlaguePeriod]] <- parse_date_time(df_copy[[begin_date]], orders = "my")
#   df_copy[[new_format_EndPlaguePeriod]] <- parse_date_time(df_copy[[end_date]], orders = "my") + lubridate::days(lubridate::days_in_month(df_copy[[new_format_EndPlaguePeriod]]) - lubridate::day(df_copy[[new_format_EndPlaguePeriod]]))
#   
#   # Sort df_copy by the column new_format_BeginPlaguePeriod in ascending order
#   df_copy <- df_copy[order(df_copy$new_format_BeginPlaguePeriod), ]
#   
#   # Add column with the number of days from the begining of the plague
#   df_copy$BeginDaysPlague <- as.integer(df_copy$new_format_BeginPlaguePeriod - min(df_copy$new_format_BeginPlaguePeriod))
#   
#   # Add column with the number of days from the end of the plague
#   df_copy$EndDaysPlague <- as.integer(df_copy$new_format_EndPlaguePeriod - min(df_copy$new_format_BeginPlaguePeriod))
#   
#   # Fix the type of the victims number column to integer 
#   df_copy[[victims_column]] <- as.integer(df_copy[[victims_column]])
#   
#   # Get the unique dates 
#   months <- unique(df_copy$new_format_EndPlaguePeriod)
#   
#   # Get unique end days of plague
#   unique_end_days_plague <- unique(df_copy$EndDaysPlague)
#   
#   # Ensure both vectors have the same length
#   if(length(months) > length(unique_end_days_plague)) {
#     unique_end_days_plague <- c(unique_end_days_plague, rep(NA, length(months) - length(unique_end_days_plague)))
#   } else if (length(months) < length(unique_end_days_plague)) {
#     months <- c(months, rep(NA, length(unique_end_days_plague) - length(months)))
#   }
#   
#   results <- data.frame(EndMonth = months,
#                         CumDays = unique_end_days_plague,
#                         NumberDeaths = rep(0, length(months)),
#                         CumDeaths = rep(0, length(months)),
#                         CumPop = rep(0, length(months)),
#                         Parishes = rep("", length(months)))
#   
#   # Iterate over the dates
#   total_deaths <- 0
#   
#   for (date in months) {
#     if (!is.na(date)) {
#       # fill the dataframe "results" so in the correspondant row the
#       # number of deaths is added to the column number deaths
#       numberOfDeaths <- sum(df_copy[df_copy$new_format_EndPlaguePeriod == date, victims_column])
#       parishes <- paste(df_copy[df_copy$new_format_EndPlaguePeriod <= date, column_name], collapse = ",")
#       results[results$EndMonth == date, 'Parishes'] <- parishes
#       results[results$EndMonth == date, 'NumberDeaths'] <- numberOfDeaths
#       total_deaths <- total_deaths + numberOfDeaths
#       results[results$EndMonth == date, 'CumDeaths'] <- total_deaths
#       
#       results[results$EndMonth == date, 'CumPop'] <- sum(df_copy[df_copy$new_format_EndPlaguePeriod <= date, pop_size])
#     }
#   }
#   
#   return(results)
# }

# Import the csv file
infectedParishes <- read.csv("/Users/dianapli/Desktop/PythonMathematicalModeling/docs/PlagueProject/data/private/infectedSouthParishes.csv")


# Testing the function with a small example
# Filter rows where ParishName is 'YSTAD', 'ÖJA' or 'HEDESKOGA'
example1 <- infectedParishes[infectedParishes$ParishName %in% c('YSTAD', 'ÖJA', 'HEDESKOGA'), ]

# example1$new_format_BeginPlaguePeriod <-parse_date_time(example1[["BeginPlaguePeriod"]], orders = "my")
# example1$new_format_EndPlaguePeriod <- parse_date_time(example1[["EndPlaguePeriod"]], orders = "my") 
# example1
# Call the function to count the number of deaths per month
# res_example1 <- count_deaths_month(example1)
# res_example1
# 
# Call the function
result_example1 <- count_infected_parishes_month(example1, 'JUN 1712', 0)
result_example1

# Convert 'date' column to Date class
result_example1$date <- as.Date(result_example1$date)

# Plot the infected parishes per month
ggplot(data = result_example1, aes(x = date, y = NumberInfectedParishes)) +
  geom_line(color = "blue") +
  scale_x_date(date_labels = "%b %Y", date_breaks = "4 month") +
  theme(axis.text.x = element_text(angle = 45)) +
  labs(x = "Month", y = "Number of infected parishes", title = "South Scania")

# Plot the cum number of infected parishes per month
ggplot(data = result_example1, aes(x = date, y = CumInfectParishes)) +
  geom_line(color = "orange") +
  scale_x_date(date_labels = "%b %Y", date_breaks = "4 month") +
  theme(axis.text.x = element_text(angle = 45)) +
  labs(x = "Month", y = "Cum. number of infected parishes", title = "South Scania")


# Use the count_infected_by_month() function and store the result
results <- count_infected_by_month(infectedSouthParishes, 'JAN 1711', 0
                                   , 'ParishName'
                                   , 'BeginPlaguePeriod'
                                   , 'EndPlaguePeriod')
results
# Convert 'date' column to Date class
results$date <- as.Date(results$date)

# Plot the infected parishes per month
ggplot(data = results, aes(x = date, y = NumberInfectedParishes)) +
  geom_line(color = "blue") +
  scale_x_date(date_labels = "%b %Y", date_breaks = "4 month") +
  theme(axis.text.x = element_text(angle = 45)) +
  labs(x = "Month", y = "Number of infected parishes", title = "South Scania")

# Plot the cum number of infected parishes per month
ggplot(data = results, aes(x = date, y = CumInfectParishes)) +
  geom_line(color = "orange") +
  scale_x_date(date_labels = "%b %Y", date_breaks = "4 month") +
  theme(axis.text.x = element_text(angle = 45)) +
  labs(x = "Month", y = "Cum. number of infected parishes", title = "South Scania")


### Objective Functions

# Function to calculate the error in the cumulative number of infected parishes per month between the model and the data
objectiveFunction_2 <- function(gdf = example, column_name = "ParishName", n = 0) {
  
  # Defining the initial date of the gdf to start counting the number of infected parishes per month
  date <- min(gdf$BeginPlaguePeriod, na.rm = TRUE)
  #
  # Getting the number of infected parishes per month from the data
  cum_infected_parishes_by_month <- count_infected_by_month(gdf,date,n)
  
  # Initializing the number of infected parishes per month for the model's output
  infected_parishes <- rep(0, length(cum_infected_parishes_by_month))
  
  # Initializing the error between the model's output and the data
  error <- rep(0, length(cum_infected_parishes_by_month))
  
  # Computing the total number of parishes in the dataframe without repetitions
  total_parishes <- length(gdf[[column_name]])
  
  # Computing the number of infected parishes per month from the model's output
  for (i in 1:length(cum_infected_parishes_by_month)) {
    init_days <- cum_infected_parishes_by_month[i,'DaysFromInitDate']
    final_days <- cum_infected_parishes_by_month[i,'DaysToEndOfMonth']
    
    for (k in 1:length(grouped_by_parish)) {
      for (day in init_days:final_days) {
        value_D_subset <- traj[traj$node == k & traj$time == day, ]$D
        
        # Check if the subset operation returned an empty vector
        if (length(value_D_subset) > 0) {
          value_D <- value_D_subset >= 1
          
          if (value_D) {
            infected_parishes[i] <- infected_parishes[i] + 1
            break # Breaks the innermost loop when the condition is met
          }
        }
      }
    }
    error[i] <- abs(infected_parishes[i] - cum_infected_parishes_by_month[i,'CumInfectParishes'])
  }
  
  
  # Computing the error between the model's output and the data
  total_error <- (1/length(cum_infected_parishes_by_month)) * (1/total_parishes) * sum(error)
  
  return(list(total_parishes, infected_parishes, error, total_error))
}

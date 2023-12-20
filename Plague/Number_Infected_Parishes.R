library(ggplot2)
library(lubridate)
library(dplyr)
library(readr) # for read csv files

count_infected_by_month <- function(df, date, n
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

# Import the csv file
infectedParishes <- read.csv("data/infectedSouthParishes.csv")

# Testing the function with a small example
# Filter rows where ParishName is 'YSTAD', 'ÖJA' or 'HEDESKOGA'
example1 <- infectedParishes[infectedParishes$ParishName %in% c('YSTAD', 'ÖJA', 'HEDESKOGA'), ]
example1

# Call the function
result_example1 <- count_infected_by_month(example1, 'JUN 1712', 0)
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

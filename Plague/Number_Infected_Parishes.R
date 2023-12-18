library(lubridate)
library(dplyr)
library(readr) # for read csv files

count_infected_by_month <- function(df, date, n, start_date = 'BeginPlaguePeriod', end_date = 'EndPlaguePeriod') {
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
    
    # Create a dataframe to store the results
    results <- data.frame(date = dates, NumberInfectedParishes = 0, CumulativeInfectedParishes = 0)
    
    # Iterate over the dates
    for (date in dates) {
      # Count nodes where infection start date is before or on the given date 
      # and either there is no end date or the end date is after the given date
      infected_nodes <- df[(df[[start_date]] <= date) & (df[[end_date]] >= date), ]
      # Store the results
      results[results$date == date, 'NumberInfectedParishes'] <- nrow(infected_nodes)
    }
    
    # Calculate the cumulative sum
    results$CumulativeInfectedParishes <- cumsum(results$NumberInfectedParishes)
  } else {
    stop("All values in the end_date column are NA. Cannot create a sequence of dates.")
  }
  
  return(results)
}

# Import the csv file
infectedParishes <- read.csv("data/infectedSouthParishes.csv")

# Testing the function with a small example
# Filter rows where ParishName is 'YSTAD', 'ÖJA' or 'HEDESKOGA'
example1 <- infectedParishes[infectedParishes$ParishName %in% c('YSTAD', 'ÖJA', 'HEDESKOGA'), ]


# Call the function
result <- count_infected_by_month(example1, 'JUN 1712', 0)


# Print the result
print(result)





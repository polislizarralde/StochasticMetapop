#install.packages("sf")
#install.packages("spdep")
#install.packages("lwgeom")
#install.packages("readr")
library(sf)
library(spdep)
library(lwgeom)
library(tidyr)
library(readr)

# Read in your shapefile
counties_shp <- st_read("data/ParishesScania.shp")
counties_shp <- st_set_crs(counties_shp,st_crs(3035))

# Plot your shapefile
plot(counties_shp)

# Check shapefile. 1 not valid, 437 valid.
table(st_is_valid(counties_shp))

# Clean up shapefile
counties_shp_clean <- st_make_valid(counties_shp, reason = TRUE)

# Plot again. Not that i see a difference.
plot(counties_shp_clean)

# Some names are duplicates. This is in part when shapes change over time.
# Our time period of interest is 1711-1714, so shapes that end before that period, or start after are no longer of relevance.
duplicates <- as.data.frame(table(counties_shp_clean$G_NAME)) %>% dplyr::filter(Freq > 1) %>% dplyr::pull(Var1)
duplicates <- dplyr::filter(counties_shp_clean,G_NAME %in% duplicates)

counties_shp_clean <- dplyr::filter(counties_shp_clean,GET_END_YE >= 1714)
counties_shp_clean <- dplyr::filter(counties_shp_clean,is.na(GET_START_) | GET_START_ <= 1711)

# Check for duplicates again in name.
duplicates <- as.data.frame(table(counties_shp_clean$G_NAME)) %>% dplyr::filter(Freq > 1) %>% dplyr::pull(Var1)

# Store indices of the cleaned shapefile.
counties_shp_clean$node <- 1:nrow(counties_shp_clean)

# Store the cleaned shapefile.
saveRDS(counties_shp_clean,"data/scania_geometry.RDS")


# Calculate the shared borders between all parishes.
shared_borders <- st_touches(counties_shp_clean, sparse = FALSE)

# Get shared border length
border_length <- function(i, j, counties_shp_clean) {
  if (shared_borders[i, j]) {
    intersection <- st_intersection(counties_shp_clean[i,], counties_shp_clean[j,])
    return(st_length(intersection))
  }
  return(0)
}

# Build an event dataframe, but without a time axis yet.
n <- nrow(counties_shp_clean)

# Fill the a dataframe with events. The proportion field will have to be scaled. Right now it is the border in meters, and it should be a value between 0 and 1.
events_sharedborder <- data.frame()
for (i in 1:n) {
  for (j in 1:n) {
    border <- as.numeric(border_length(i, j, counties_shp_clean))
    if (border > 0) {
        print(paste(i,j,border))
        events_sharedborder <- rbind(events_sharedborder,data.frame(event = 'extTrans',
                                                                    node = i,
                                                                    dest = j,
                                                                    n = 0,
                                                                    proportion = border,
                                                                    select=1,
                                                                    shift=0))
    }
  }
}

# Store dataframe object
saveRDS(events_sharedborder,"data/events_sharedborder.RDS")

# Calculate centroids
centroids <- st_centroid(counties_shp_clean)

# Add the node, for looking up in events_sharedborder
centroids$node <- 1:nrow(centroids)

# Visualize
plot(centroids)

# Calculate distance for parishes that share a border.
events_distance <- data.frame()
for (i in 1:n) {
  for (j in 1:n) {
    border <- as.numeric(border_length(i, j, counties_shp_clean))
    if (border > 0) {
      dist <- as.numeric(st_distance(centroids[i,],centroids[j,]))
      print(paste(i,j,dist))
      events_distance <- rbind(events_distance,data.frame(event = 'extTrans',
                                                                  node = i,
                                                                  dest = j,
                                                                  n = 0,
                                                                  proportion = dist,
                                                                  select=1,
                                                                  shift=0))
    }
  }
}

# Check. None of the distances should be 0 unless they are the same node and dest.
min(events_distance$proportion)

# Store dataframe object
saveRDS(events_distance,"data/events_distances.RDS")

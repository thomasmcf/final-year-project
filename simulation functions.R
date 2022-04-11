library(secr)
library(tidyverse)
library(survival)
library(jagsUI)
library(foreach)
library(parallel)
library(doParallel)

# Length of side of study site
SITE_LENGTH <- 200 
# Range parameter for single leopard
RANGE <- 50

# This function simulates the dropping information of one individual leopard, given period of occupation and number of droppings
# Inputs: arrived   - Time that the individual arrived in the study site
#         left      - Time that the individual left the study site
#         drop_rate - Mean number of droppings produced in one day
#                id - Individual identifier
#
# Outputs: A list of eight variables:
#          activity_centre - The geographic centre of the individuals activity
#          present         - The time interval the individual was present for
#          drop_times      - The times a dropping was created by the individual
#          surv_times      - The "survival" time of the droppings, which is the time the dropping lasts until degrading beyond genotyping
#          decay_times     - The time at which the dropping degrades beyond genotyping
#          location        - A matrix of x and y coordinates of dropping locations
#          id              - Individual identifier
#          drop_id         - A vector of unique droppign IDs
individual <- function(arrived, left, mean_life, drop_rate, id){
  # Simulate the dropping, survival, and decay times as defined above
  droppings <- rpois(1, drop_rate * (left - arrived))
  drop_times <- runif(droppings, arrived, left)
  surv_times <- rexp(droppings, 1 / mean_life)
  decay_times <- drop_times + surv_times

  # Simulate the activity centre of the individual
  # The coordinates are modeled as independent uniform random variables, centered on (0, 0)
  activity_centre <- runif(2, -SITE_LENGTH/2, SITE_LENGTH/2)
  # Simulate the location of the droppings
  # The coordinates are independent normal random variables centred on the activity centre
  location_x <- rnorm(droppings, activity_centre[1], sqrt(RANGE))
  location_y <- rnorm(droppings, activity_centre[2], sqrt(RANGE))
  
  location <- rbind(location_x, location_y)
  # Generate an ID for each dropping, used to uniquely identify them
  drop_id <- paste(id, 1:droppings, sep = "-")
  
  # Return the simulated information in a list
  return(list(activity_centre = activity_centre,
              present = c(arrived, left),
              drop_times = drop_times,
              surv_times = surv_times,
              decay_times = decay_times,
              location = location,
              id = id,
              drop_id = drop_id))
}


# This function simulates the population with a queue
# Inputs:  N - Asymptotic expected number of individuals in the population
#  mean_stay - The average length of time spent in the population before leaving
#   var_stay - The variance of the time spent in the population before leaving
#      t_max - Time to simulate the population for
#
# Outputs: A list of four elements
#          individuals - A data frame with the time interval each individual was present for
#                times - A vector with times at which the population changes
#               change - A vector with the instantaneous changes in population, corresponding with times
#           population - Each element is the population in the area as of the corresponding time in times
queue <- function(N, mean_stay, t_max, var_stay){
  # Set parameters of the gamma distribution by matching moments
  shape <- mean_stay ** 2 / var_stay
  rate <-  mean_stay / var_stay
  
  # Initialise the individuals data frame
  individuals <- data.frame(arrival = double(),
                            departure = double())
  
  # Generate the time at which the next individual moves into the area
  # The shape parameter is adjusted in order to specify the stationary expectation
  latest <- rgamma(1, shape / N, rate)
  # Simulate arrival times until the simulation period has elapsed
  while(latest < t_max){
    # Draw the new arrival's length of stay, and append this to the data frame
    stay_time <- rgamma(1, shape, rate)
    individuals <- rbind(individuals,
                         list(arrival = latest, departure = latest + stay_time))
    
    # Draw next arrival time
    latest <- latest + rgamma(1, shape / N, rate)
  }
  
  # Find how many individuals have been present in total
  n_individs <- dim(individuals)[1]
  # Convert the data frame to a matrix
  individuals <- as.matrix(individuals)
  
  # "Unroll" all times into one vector
  # The first n_individs are the arrival times, the rest are departure times
  times <- as.vector(individuals)
  # Specify the instantaneous change in population at the times in times
  # For the arrival times this is +1, and the departures it is -1
  change <- rep(c(1, -1), c(n_individs, n_individs))
  
  # Order both change and times in chronological order
  change <- change[order(times)]
  times <- sort(times)
  
  # Calculate the abundance at each time point
  population <- cumsum(change)
  
  # Return the simulated information
  list(individuals = individuals,
       times = times,
       change = change,
       population = population)
}


# This function simulates the dropping information of a population of leopards
# Inputs: N - The asymptotic expected number of individuals
#     t_max - The time at which to stop simulating the population
# Outputs: individuals - A list of individuals with the information generated by individual()
population <- function(N, t_max, mean_life = 365/2, drop_rate = 0.5, mean_stay = 365, var_stay = (56 / 1.96) ** 2){
  # Simulate occupancy history of the area
  history <- queue(N, mean_stay, t_max, var_stay)$individuals
  N_tot <- dim(history)[1]
  
  # Initialise a list to store individuals
  individuals <- vector("list", N_tot)
  
  for(i in 1:N_tot){
    present <- history[i, ]
    # Simulate information on each individual
    individuals[[i]] <- individual(present[1], present[2], mean_life, drop_rate, id = i)
  }
  
  return(individuals)
}

# Determine number of individuals present and detectable in population at time
# Inputs: population - A population object
#               time - The time at which to evaluate present and detectable
pres_det <- function(population, time){
  # Find the length of the population
  N <- length(population)
  N_present <- 0
  N_detectable <- 0
  
  # Iterate over each individual
  for(i in 1:N){
    ind <- population[[i]]
    
    # If animal j was present at time, increment the value in the presence vector
    if(ind$present[1] <= time & time <= ind$present[2]){
      N_present <- N_present + 1
    }
    
    # If any of animals j's droppings were present at time, increment the value in the detectable vector
    if(any(ind$drop_times <= time & time <= ind$decay_times)){
      N_detectable <- N_detectable + 1
    }
  }
  
  list("present" = N_present, "detectable" = N_detectable)
}

# Delete the individuals that are not present or detectable at time
# Speeds up simulation of capture histories
# Inputs: population - A population object
#               time - Delete animals not relevant at time
prune <- function(population, time){
  N <- length(population)
  
  
  keep <- NULL
  # Iterate over the population
  for(i in 1:N){
    ind <- population[[i]]
    
    # If the individual is present and/or detectable, and it to the list of ones to keep
    if(ind$present[1] <= time & time <= ind$present[2]){
      keep <- c(keep, i)
    }
    else if(any(ind$drop_times <= time & ind$decay_times >= time)){
      keep <- c(keep, i)
    }
  }
  
  # Make a new population of just the relevant individuals
  new_pop <- vector("list", length(keep))
  for(i in 1:length(keep)){
    new_pop[[i]] <- population[[keep[i]]]
  }
  
  return(new_pop)
}

# This function calculates the minimum distance from the point (a, b) to the line segment
# joining (x0, y0) and (x1, y1)
min_segment_dist <- function(a, b, x0, y0, x1, y1){
  # Calculate the parameter for the point on the segment closest to (a, b)
  t <- - ((x0 - a)*(x1 - x0) + (y0-b)*(y1-y0)) / ((x1-x0)**2 + (y1-y0)**2)
  
  # If t <= 0, (a, b) is closest to (x0, y0)
  if(t <= 0){
    x_min <- x0
    y_min <- y0
  }
  # If t >= 1, it is closest to (x1, y1)
  else if(t >= 1){
    x_min <- x1
    y_min <- y1
  }
  # If 0 < t < 1, it is closest to a point between (x0, y0) and (x1, y1)
  else{
    x_min <- x0 + (x1 - x0) * t
    y_min <- y0 + (y1 - y0) * t
  }
  
  return(list(dist = sqrt((x_min - a) ** 2 + (y_min - b) ** 2),
              x_min = x_min,
              y_min = y_min))
}

# Simulate the capture history for one occasion
#         x0 - Starting x coordinates of transects
#         y0 - Starting y coordinates of transects
#         x1 - Ending x coordinates of transects
#         y1 - Ending y coordinates of transects
# population - Simulated Population object
#       time - The time at which visit takes place
#   occasion - Which occasion is being simulated
#    session - Which session is being simulated
#         t0 - Time at which the thorough search takes place
occasion_sim <- function(x0, y0, x1, y1, population, time, occasion, session = 1, t0){
  # The number of individuals in the population object
  N <- length(population)
  transects <- length(x0)
  
  # Empty data frame to hold the capture history for secr
  capthist <- data.frame(Session = integer(),
                         AnimalID = integer(),
                         Occasion = integer(),
                         X = double(),
                         Y = double())
  
  # Empty data frame to hold capture history for survival analysis
  survhist <- data.frame(droppingID = character(),
                         time = integer(),
                         type = character(),
                         present = character())
  
  # Iterate over the population
  for(i in 1:N){
    # Extract individual i
    ind <- population[[i]]
    # Find the number of droppings produced by this individual
    droppings <- length(ind$drop_times)
    
    for(j in 1:transects){
      # Iterate over each dropping
      for(k in 1:droppings){
        # Check if the dropping was detectable at the time of the survey
        if(ind$drop_times[k] <= time){
          # Find the length between the dropping and the transect
          dist <- min_segment_dist(ind$location[1, k], ind$location[2, k], x0[j], y0[j], x1[j], y1[j])
          # Calculate the probability of detection
          prob <- HHN(dist$dist)
          # Simulate the detection/non-detection (Bernoulli)
          detect <- rbinom(1, 1, prob)
          
          # If the dropping was detected record this
          if(detect == 1){
            if(time <= ind$decay_times[k]){
              capthist <- rbind(capthist,
                                list(Session = session,
                                     AnimalID = i,
                                     Occasion = occasion,
                                     X = dist$x_min,
                                     Y = dist$y_min))
            }
            
            survhist <- rbind(survhist,
                              list(droppingID = ind$drop_id[k],
                                   time = time,
                                   type = ifelse(time <= ind$decay_times[k], "detected", "degraded"),
                                   present = ifelse(ind$drop_times[k] <= t0, TRUE, FALSE)))
          }
        }
      }
    }
  }
  
  survhist <- distinct(survhist)
  
  return(list(capthist = capthist,
              survhist = survhist))
}


# Simulate the capture history for one occasion
#         x0 - Starting x coordinates of transects
#         y0 - Starting y coordinates of transects
#         x1 - Ending x coordinates of transects
#         y1 - Ending y coordinates of transects
# population - Simulated Population object
#      times - The times at which visit takes place
#         t0 - Time at which the thorough search takes place
survey_sim <- function(x0, y0, x1, y1, population, times, t0){
  # Find the number of occasions to simulate
  n_occasions <- length(times)
  # Simulate the first one
  occasion <- occasion_sim(x0, y0, x1, y1, population, times[1], occasion = 1, t0 = t0)
  
  # Extract the capture and survival histories
  capthist <- occasion$capthist
  survhist <- occasion$survhist
  
  # Iterate over the rest of the occasions
  for(t in 2:n_occasions){
    occasion <- occasion_sim(x0, y0, x1, y1, population, times[t], occasion = t, t0 = t0)
    
    # Update the capture and survival histories
    capthist <- rbind(capthist,
                      occasion$capthist)
    
    survhist <- rbind(survhist,
                      occasion$survhist)
    
  }
  
  # This reassigns the IDs, not too important, just makes debugging easier
  captured_ids <- unique(capthist$AnimalID)
  ids <- cbind(sort(captured_ids), 1:length(captured_ids))
  
  for(i in 1:nrow(capthist)){
    old_id <- capthist[i, "AnimalID"]
    new_id <- ids[which(ids[, 1] == old_id), 2]
    
    capthist[i, "AnimalID"] <- new_id
  }
  
  # Convert to capthist and Surv objects
  capthist <- make.capthist(capthist, transect_list, fmt="XY")
  survhist <- make.surv(survhist, t0, times)
  
  return(list(capthist = capthist,
              survhist = survhist))
}

# Converts a survival history into a Surv object
# survhist - A survival history generated by survey_sim
#       t0 - Time of thorough search
#    times - Times at which occasions occur
make.surv <- function(survhist, t0, times){
  # Convert to wide format
  survhist <- survhist %>%
    pivot_wider(names_from = time, values_from = type, values_fill = "no detection")
  
  # Add columns for number of detections, and whether degradation was observed
  # Delete rows with not enough data
  survhist <- survhist %>%
    mutate(degradation_observed = apply(survhist, 1, function(row) "degraded" %in% row),
           detections = apply(survhist, 1, function(row) sum("detected" == row))) %>%
    filter(!(present & detections == 0))%>%
    filter(!(!present & !degradation_observed & (detections == 1)))
  
  droppings <- nrow(survhist)
  
  time <- numeric(droppings)
  time2 <- numeric(droppings)
  event <- numeric(droppings)
  
  # Find what information is known about the lifetime of the dropping
  for(i in 1:droppings){
    if(survhist$present[i]){
      t_last <- times[max(which(survhist[i, ] == "detected")) - 2]
      
      time[i] <- t_last - t0
      time2[i] <- NA
      
      event[i] <- 0
    }
    
    else if(survhist$degradation_observed[i]){
      t_deg <- times[min(which(survhist[i, ] == "degraded")) - 2]
      
      if(survhist$detections[i] <= 1){
        time[i] <- t_deg - t0
        time2[i] <- NA
        event[i] <- 2
      }
      
      else{
        t_init <- times[min(which(survhist[i, ] == "detected")) - 2]
        t_last <- times[max(which(survhist[i, ] == "detected")) - 2]
        
        time[i] <- t_last - t_init
        time2[i] <- t_deg - t0
        event[i] <- 3 
      }
    }
    
    else{
      t_init <- times[min(which(survhist[i, ] == "detected")) - 2]
      t_last <- times[max(which(survhist[i, ] == "detected")) - 2]
      
      time[i] <- t_last - t_init
      time2[i] <- NA
      event[i] <- 0
    }
  }
  
  # Convert to a Surv object  
  Surv(time = time, time2 = time2, event = event, type = "interval")
}

lambda <- function(d, lambda0, sigma){
  lambda0 * exp(-d**2 / (2 * sigma**2))
}

HHN <- function(d, lambda0 = 1, sigma = 1){
  1 - exp(-lambda(d, lambda0, sigma))
}
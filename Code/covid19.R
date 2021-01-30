# ===================== PACKAGES =====================
source("/Users/Antoine/Desktop/Projet Semestre 2020/CODE/Deconvolution.R")
source("/Users/Antoine/Desktop/Projet Semestre 2020/CODE/EstimatingRt.R")

library(ggplot2)
library(shinythemes)
library(openxlsx)
library(readr)
library(cowplot)
library(zoo)
library(mgcv)
library(latex2exp)
library(splines)
require(stats)
require(graphics)

# ---------------------------- THEMES -------------------------------------
colors <- c("#00AFBB", "#E7B800", "#FC4E07")
dark_theme <- FALSE
if (dark_theme){
  font_color_plot <- '#1C1E21'
  points_color <- "grey85"
  theme_for_the_plots <- theme(plot.background = element_rect(fill = font_color_plot, colour = "black"),
                               panel.background = element_rect(fill = "#1C1E21", colour = "black"),
                               panel.grid.major = element_line(colour = "grey85"),
                               panel.grid.minor = element_line(colour = "grey90"),
                               legend.position = c(.95, .95),
                               legend.justification = c("right", "top"),
                               legend.box.just = "right",
                               legend.margin = margin(6, 6, 6, 6),
                               plot.margin = unit(c(1,1,1,1), "cm"),
                               axis.text = element_text(colour ="white"),
                               # axis.line = element_line(colour = "white"),
                               axis.line = element_blank(),
                               axis.ticks = element_line(colour = "white"),
                               axis.title = element_text(colour = "white"),
                               title = element_text(colour = "white"),
                               legend.background = element_rect(fill = font_color_plot, colour = font_color_plot),
                               legend.text = element_text(colour ="white"),
                               legend.key = element_rect(colour = NA, fill = NA))
  theme_shiny <- shinytheme('slate')
}else{
  font_color_plot <- "white"
  points_color <- "black"
  theme_for_the_plots <- theme(plot.background = element_rect(fill = font_color_plot),
                               panel.background = element_rect(fill = "white", colour = "black"),
                               panel.grid.major = element_line(colour = "grey85"),
                               panel.grid.minor = element_line(colour = "grey90"),
                               legend.position = c(.95, .95),
                               legend.justification = c("right", "top"),
                               legend.box.just = "right",
                               legend.margin = margin(6, 6, 6, 6),
                               plot.margin = unit(c(1,1,1,1), "cm"),
                               legend.key = element_rect(colour = "transparent", fill = "transparent"))
  theme_shiny <- shinytheme('paper')
}


# ===================== COVID-19 =================================

# Switzerland
lockdown_events <- list("Name" = c("Store & markets closure", "Museums & nightclubs closure", "Hairdresser closure",
                                   "Restaurants and bar closure", "School closure", "Mask mandatory on public transport", 
                                   "Restaurants and bar closure", "Extension of restrictions", "Grouping restrictions"),
                        "Start" = c("2020-03-17", "2020-03-17", "2020-03-17", "2020-03-17", "2020-03-16", "2020-07-06", "2020-12-22", "2020-10-19", "2020-10-29"),
                        "End" = c("2020-05-10", "2020-06-06", "2020-04-26", "2020-05-10", "2020-05-10", "2021-01-22", "2021-01-22", "2021-01-22", "2021-01-22"),
                        "level" = c(1, 2, 3, 4, 5, 5, 4, 2, 3))
settings <- list("country" = "Switzerland", "remove_last_n_observation" = 1, "dim_deconv" = 20, "delta_cst" = 0.0001, "range_legends" = c(2.5,4),
                 "SIR_pop" = 8.57e6, "SIR_p" = 30, "SIR_q" = 10, "SIR_cst" = TRUE, "SIR_eps" = 0.8, "SIR_first_index" = 8, "Lockdowns" = TRUE)

# France
# lockdown_events <- list("Name" = c("Lockdown", "Restaurants closure", "School closure", "Few restrictions", "Lockdown"),
#                         "Start" = c("2020-03-15", "2020-03-16", "2020-03-16", "2020-07-19", "2020-10-17"),
#                         "End" = c("2020-05-10", "2020-06-14", "2020-05-10", "2020-10-16", "2020-11-04"))
# settings <- list("country" = "France", "remove_last_n_observation" = 0, "dim_deconv" = 20, "delta_cst" = 0.0001, "range_legends" = c(4,5),
#                  "SIR_pop" = 8.57e6, "SIR_p" = 30, "SIR_q" = 10, "SIR_cst" = TRUE, "SIR_eps" = 0.8, "SIR_first_index" = 20, "Lockdowns" = FALSE)

# Germany
# lockdown_events <- list("Name" = c("Mobility restrictions", "Interdiction of gathering > 2 people", "Restrictions/lockdown"),
#                         "Start" = c("2020-03-16", "2020-03-22", "2020-03-14"),
#                         "End" = c("2020-06-15", "2020-06-15", "2020-05-11"))
# settings <- list("country" = "Germany", "remove_last_n_observation" = 0, "dim_deconv" = 40, "delta_cst" = 0.0001, "range_legends" = c(6,7),
#                  "SIR_pop" = 80e6, "SIR_p" = 10, "SIR_q" = 10, "SIR_cst" = TRUE, "SIR_eps" = 0, "SIR_first_index" = 20)

# Italy
# lockdown_events <- list("Name" = c("Lockdown of 50000 people", "Full lockdown", "Non-essential market closure"),
#                         "Start" = c("2020-02-21", "2020-03-10", "2020-03-21"),
#                         "End" = c("2020-03-09", "2020-06-03", "2020-06-03"))
# settings <- list("country" = "Italy", "remove_last_n_observation" = 0, "dim_deconv" = 40, "delta_cst" = 0.0001, "range_legends" = c(6,7),
#                  "SIR_pop" = 60e6, "SIR_p" = 20, "SIR_q" = 10, "SIR_cst" = FALSE, "SIR_eps" = 0, "SIR_first_index" = 30)

# India
# lockdown_events <- list("Name" = c("Phase 1 (lockdown)", "Phase 2 (lockdown)", "Phase 3 (lockdown)", "Phase 4 (lockdown)",
#                                    "Unlock 1", "Unlock 2", "Unlock 3", "Unlock 4", "Unlock 5", "Unlock 6"), 
#                         "Start" = c("2020-03-25", "2020-04-15", "2020-05-04", "2020-05-18",
#                                     "2020-06-01", "2020-07-01", "2020-08-01", "2020-09-01", "2020-10-01", "2020-11-01"), 
#                         "End" = c("2020-04-14", "2020-05-03", "2020-05-17", "2020-05-31",
#                                   "2020-06-30", "2020-07-31", "2020-08-31", "2020-09-30", "2020-10-31", "2020-11-30"))
# settings <- list("remove_last_n_observation" = 0, "dim_deconv" = 20, "delta_cst" = 0.0001, "range_legends" = c(4.5,7), 
#                  "SIR_pop" = 1.3e9, "SIR_p" = 10, "SIR_q" = 10, "SIR_cst" = FALSE, "SIR_eps" = 0, "SIR_first_index" = 30)

# United States
# settings <- list("country" = "United States", "remove_last_n_observation" = 0, "dim_deconv" = 20, "delta_cst" = 0.0001, "range_legends" = c(4.5,7),
#                  "SIR_pop" = 1.3e9, "SIR_p" = 10, "SIR_q" = 10, "SIR_cst" = FALSE, "SIR_eps" = 0, "SIR_first_index" = 30)


# ------------------------- 0) Collect the data ------------------ 
#                                   Functions used ---------------
# Function to get the data for the corresponding country, including the daily number of infected, recovered and deaths
#' @param country country for which the data need to be collected
#' @param option 1 to get the data from bag.admin.ch (only for Switzerland), 2 for another source
#' @param plot_ TRUE to plot the collected data
ConstructDataset <- function(country, option = 1, plot_ = F){
  
  # -------------- Some functions --------------
  # Function to get the recovered cases for the corresponding country
  #' @param country country for which the data need to be collected
  getRecovered <- function(country){
    if (country == "United States"){
      country <- "US"
    }
    data_recovered <- as.data.frame(read_csv(file = "https://data.humdata.org/hxlproxy/api/data-preview.csv?url=https%3A%2F%2Fraw.githubusercontent.com%2FCSSEGISandData%2FCOVID-19%2Fmaster%2Fcsse_covid_19_data%2Fcsse_covid_19_time_series%2Ftime_series_covid19_recovered_global.csv&filename=time_series_covid19_recovered_global.csv"))
    data_recovered_matrix <- data.matrix(data_recovered[5:length(data_recovered[1,])])
    # print(data_recovered[,2])
    
    l <- length(data_recovered[1, 5:length(data_recovered[1,])])
    
    data_country <- data.frame(matrix(ncol = 3, nrow = l))
    colnames(data_country) <- c("Date", "recovered", "new_recovered")
    
    data_country$Date <- as.Date(colnames(data_recovered)[5:length(data_recovered[1,])], "%m/%d/%y")
    
    if (length(which(data_recovered[,2] == country)) <= 1){
      data_country$recovered <- data_recovered_matrix[which(data_recovered[,2] == country),]
    }else{
      data_country$recovered <- colSums(data_recovered_matrix[which(data_recovered[,2] == country),])
    } 
    data_country$new_recovered <- data_country$recovered - c(0, data_country$recovered[1:(l-1)])
    
    return(data_country)
  }
  
  # Function to collected the data from the bag.admin API
  BAGData <- function(){
    data <- read.xlsx(sep=",",startRow = 8, detectDates = TRUE,
                      "https://www.bag.admin.ch/dam/bag/de/dokumente/mt/k-und-i/aktuelle-ausbrueche-pandemien/2019-nCoV/covid-19-datengrundlage-lagebericht.xlsx.download.xlsx/200325_Datengrundlage_Grafiken_COVID-19-Bericht.xlsx")
    names(data) <- c("date","cases","casesCumul","hospitalized","hospitalizedCumul",
                     "deaths","deathsCumul")
    return(data)
  }
  
  # Function to make a data frame of the infected, deaths and recovered cases for the chosen country, from data.humdata
  #' @param country country for which the data need to be collected
  prepare_data <- function(country){
    data_confirmed <- as.data.frame(read_csv(file = "https://data.humdata.org/hxlproxy/api/data-preview.csv?url=https%3A%2F%2Fraw.githubusercontent.com%2FCSSEGISandData%2FCOVID-19%2Fmaster%2Fcsse_covid_19_data%2Fcsse_covid_19_time_series%2Ftime_series_covid19_confirmed_global.csv&filename=time_series_covid19_confirmed_global.csv"))
    data_deaths <- as.data.frame(read_csv(file = "https://data.humdata.org/hxlproxy/api/data-preview.csv?url=https%3A%2F%2Fraw.githubusercontent.com%2FCSSEGISandData%2FCOVID-19%2Fmaster%2Fcsse_covid_19_data%2Fcsse_covid_19_time_series%2Ftime_series_covid19_deaths_global.csv&filename=time_series_covid19_deaths_global.csv"))
    data_recovered <- as.data.frame(read_csv(file = "https://data.humdata.org/hxlproxy/api/data-preview.csv?url=https%3A%2F%2Fraw.githubusercontent.com%2FCSSEGISandData%2FCOVID-19%2Fmaster%2Fcsse_covid_19_data%2Fcsse_covid_19_time_series%2Ftime_series_covid19_recovered_global.csv&filename=time_series_covid19_recovered_global.csv"))
    
    data_confirmed_matrix <- data.matrix(data_confirmed[5:length(data_confirmed[1,])])
    data_deaths_matrix <- data.matrix(data_deaths[5:length(data_deaths[1,])])
    data_recovered_matrix <- data.matrix(data_recovered[5:length(data_recovered[1,])])
    
    l <- length(data_confirmed[1, 5:length(data_confirmed[1,])])
    
    data_country <- data.frame(matrix(ncol = 7, nrow = l))
    colnames(data_country) <- c("Date", "confirmed", "deaths", "recovered", "new_confirmed", "new_deaths", "new_recovered")
    
    print(colnames(data_confirmed)[5:length(data_confirmed[1,])])
    data_country$Date <- as.Date(colnames(data_confirmed)[5:length(data_confirmed[1,])], "%m/%d/%y")
    # data_country$Date <- colnames(data_confirmed)[5:length(data_confirmed[1,])]
    
    if (length(which(data_deaths[,2] == country)) <= 1){
      data_country$confirmed <- data_confirmed_matrix[which(data_confirmed[,2] == country),]
      data_country$deaths <- data_deaths_matrix[which(data_deaths[,2] == country),]
      data_country$recovered <- data_recovered_matrix[which(data_recovered[,2] == country),]
    }else{
      data_country$confirmed <- colSums(data_confirmed_matrix[which(data_confirmed[,2] == country),])
      data_country$deaths <- colSums(data_deaths_matrix[which(data_deaths[,2] == country),])
      data_country$recovered <- colSums(data_recovered_matrix[which(data_recovered[,2] == country),])
    }
    
    data_country$new_confirmed <- data_country$confirmed - c(0, data_country$confirmed[1:(l-1)])
    data_country$new_deaths <- data_country$deaths - c(0, data_country$deaths[1:(l-1)])
    data_country$new_recovered <- data_country$recovered - c(0, data_country$recovered[1:(l-1)])
    
    # par(mfrow = c(2,3))
    # plot(data_country$Date, data_country$confirmed, type = 'l')
    # plot(data_country$Date, data_country$deaths, type = 'l')
    # plot(data_country$Date, data_country$recovered, type = 'l')
    # plot(data_country$Date, data_country$new_confirmed, type = 'l')
    # plot(data_country$Date, data_country$new_deaths, type = 'l')
    # plot(data_country$Date, data_country$new_recovered, type = 'l')
    
    d1 <- which(data_country$confirmed!=0)
    d2 <- which(data_country$deaths!=0)
    d3 <- which(data_country$recovered!=0)
    date_first_not_0 <- intersect(d1, d2)[1]
    
    return(list(data_country, date_first_not_0))
  }
  
  # Function to replace the NA and negative values by 0 values
  removeNA <- function(x){
    x[which(is.na(x))] <- 0
    x[which(x < 0)] <- 0
    return(x)
  }
  
  # Function to plot the collected data
  #' @param DAT dataframe of the data
  plot_data <- function(DAT, obs, title_, points_color = "black", moving_average_window = 7, smooth = FALSE){
    plot_of_data <- ggplot(data = DAT, aes(x = Date, y = obs)) + geom_line(color = points_color) + 
      labs(title=title_, x="Time", y="count") + 
      theme_for_the_plots
    if (smooth){
      plot_of_data <- plot_of_data + geom_line(aes(y = rollapply(obs, 
                                                                 width=moving_average_window, FUN=function(x) mean(x, na.rm=TRUE), 
                                                                 by=1, by.column=TRUE, partial=TRUE, fill=NA, align="center")), 
                                               color = colors[3], size=0.8)
    }
    return(plot_of_data)
  }
  # --------------------------------------------
  
  # Collect data depending on the option
  if (option == 1){
    data <- BAGData()
    data_cases <- data$cases
    data_dates <- data$date
    data_deaths <- data$deaths
    data_dates <- data$date
    
  }else if (option == 2){
    data <- as.data.frame(read_csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv"))
    # Select data 
    data_c <- data[which(data[,"location"] == country), ]
    
    # Some preprocessing, remove first zeros
    first_value = 1
    data_c$new_cases <- removeNA(data_c$new_cases)
    while(data_c$new_cases[first_value] < 1){ first_value = first_value+1 }
    data_c <- data_c[first_value:length(data_c$date),]
    data_c$new_cases[which(data_c$new_cases < 1)] <- 1
    
    data_cases <- data_c$new_cases
    data_deaths <- data_c$new_deaths
    data_dates <- data_c$date
    
  }else if (option == 3){
    data <- prepare_data(country)
    
    return(data)
  }
  
  data_rec <- getRecovered(country)
  inter_dates <- as.Date(intersect(data_dates, data_rec$Date))
  data_deaths <- removeNA(data_deaths)
  data_cases <- removeNA(data_cases)
  data_rec$new_recovered <- removeNA(data_rec$new_recovered)
  
  # Put all data into a data frame
  data_country <- data.frame(matrix(ncol = 7, nrow = length(inter_dates)))
  colnames(data_country) <- c("Date", "confirmed", "deaths", "recovered", "new_confirmed", "new_deaths", "new_recovered")
  data_country$Date <- inter_dates
  data_country$new_confirmed <- data_cases[which(data_dates %in% inter_dates)]
  data_country$new_deaths <- data_deaths[which(data_dates %in% inter_dates)]
  data_country$new_recovered <- data_rec$new_recovered[which(data_rec$Date%in%inter_dates)]
  data_country$confirmed <- cumsum(data_country$new_confirmed)
  data_country$deaths <- cumsum(data_country$new_deaths)
  data_country$recovered <- cumsum(data_country$new_recovered)
  
  if (plot_){
    plot_new_confirmed <- plot_data(data_country, data_country$new_confirmed, paste("New confirmed cases (",country, ")" ))
    plot_new_deaths <- plot_data(data_country, data_country$new_deaths, paste("New deaths (",country, ")" ))
    plot_new_recovered <- plot_data(data_country, data_country$new_recovered, paste("New recovered cases (",country, ")" ))
    
    plot_Total_confirmed <- plot_data(data_country, data_country$confirmed, paste("Total confirmed cases (",country, ")" ))
    plot_Total_deaths <- plot_data(data_country, data_country$deaths, paste("Total deaths (",country, ")" ))
    plot_Total_recovered <- plot_data(data_country, data_country$recovered, paste("Total recovered cases (",country, ")" ))
    
    # plot(plot_new_confirmed)
    plot_general <- plot_grid(plot_new_confirmed, plot_new_deaths, plot_new_recovered, 
                              plot_Total_confirmed, plot_Total_deaths, plot_Total_recovered, nrow = 2)
    plot(plot_general)
    
  }
  
  return(data_country)
}

# Function to replace the NA and negative values by 0 values
removeNA <- function(x){
  x[which(is.na(x))] <- 0
  x[which(x < 0)] <- 0
  return(x)
}

# Function that redistributes the total cases for sunday, saturday and monday equally on these 3 days
#' @param cases vector of new cases
#' @param dates vector of dates, of class Dates
reDistributeMondayCases <- function(cases, dates){
  weekdays_ <- weekdays(dates)
  indices_ <- c(1:length(dates))[weekdays_ %in% c("Monday")]
  cases_reworked <- cases
  for (j in 1:length(indices_)){
    index_set_week <- c((indices_[j]-2):indices_[j])
    index_set_week <- index_set_week[which(index_set_week >= 1)]
    mean_weekend_and_monday <- mean(cases[index_set_week])
    cases_reworked[index_set_week] <- mean_weekend_and_monday
  }
  return(cases_reworked)
}


#                                   Collect data -----------------

selected_country <- settings$country
if (selected_country == "Switzerland"){ option <- 1 }else{ option <- 2 }

swiss_data <- ConstructDataset(selected_country, option, T)
data_cases <- swiss_data$new_confirmed
data_dates <- swiss_data$Date
# data_cases <- reDistributeMondayCases(data_cases, data_dates)

# factor(weekdays(swiss_data$Date))

indices <- c(1:(length(data_cases)-settings$remove_last_n_observation))

# List of measures, for the plots 
space <- (settings$range_legends[2] - settings$range_legends[1])/(length(lockdown_events$Name)-1)
lockdown_data <- data.frame(x=as.Date(lockdown_events$Start), y=seq(settings$range_legends[1], settings$range_legends[2], space), 
                            vx=as.Date(lockdown_events$End), vy=seq(settings$range_legends[1], settings$range_legends[2], space),
                            name = lockdown_events$Name, ymin = rep(0, length(lockdown_events$Start)), ymax = rep(max(settings$range_legends), length(lockdown_events$Start)))
if (!is.null(lockdown_events$level)){
  lockdown_data$y <- settings$range_legends[1] +  (settings$range_legends[2] - settings$range_legends[1]) / length(unique(lockdown_events$level)) * lockdown_events$level
  lockdown_data$vy <- lockdown_data$y
}

# to replicate plots of the report:
indices <- c(1:303)

#                                   Some plots and extensions -----------------
# First, we extend the data: the 3 curves of infected, deceased and recovered individuals
# We will predict (extrapolate) for the next 20 days
prediction_window <- 20

# Function to plot the result
plot_country <- function(extension_results, indices_, times = NULL, extra = NULL, ...){
  if (is.null(times)){
    times <- extension_results$DATA$time
  }
  
  plot_final <- ggplot(data = NULL, aes(x = times))
  if (!is.null(extra)){
    plot_final <- plot_final + geom_line(aes(y = extra), alpha = 0.8, color = "gray")
  }
  plot_final <- plot_final + geom_vline(xintercept = times[length(indices_)], color = "red") + 
    geom_ribbon(aes(ymin = extension_results$extension_results$lower, 
                    ymax = extension_results$extension_results$upper, 
                    fill = "Confidence interval extension", colour = NA), alpha = 0.3, fill = colors[3]) +
    geom_line(aes(y = extension_results$extension_results$fitted, colour = "extension"), size = 0.7) + 
    geom_segment(aes(x = times[length(indices_)], y = 0, xend = times[length(extension_results$DATA$time)], yend = 0), 
                 color = "#e63946", arrow = arrow(length = unit(0.2, "cm"))) + 
    geom_text(data=NULL, aes(x= times[length(indices_)] + (times[length(extension_results$DATA$time)] - times[length(indices)])/2, y = 0, 
                             label = "Prediction"), colour = "#e63946", vjust = 1.5, size = 3) + 
    
    scale_fill_manual("", 
                      breaks = c("Confidence interval extension", "deconvolution", "Simlutaneous confidence bands"),
                      values = c("#e63946", colors[2], colors[2])) +
    scale_colour_manual("", 
                        breaks = c('observed incidence', "cases", "extension", "deconvolution"),
                        values = c("black", "#1d3557", "#364958",  colors[2])) +
    
    labs(...) + 
    theme_for_the_plots + 
    guides(color=guide_legend(override.aes=list(fill=NA))) + theme(legend.position = 'bottom') + 
    coord_cartesian(ylim = c(0, 1.2*max(extension_results$DATA$cases, na.rm = TRUE)))
  
  if (!is.null(times)){
    plot_final <- plot_final + scale_x_date(date_breaks = "1 month", date_labels = "%b")
  }
  
  return(plot_final)
}

# Infection events:
extensionSwiss <- extension_prediction(data_cases[indices], extension = list('len' = prediction_window, 'family' = nb(link = 'log'), 
                                                                             'keep_original' = TRUE, 'min_weights' = 0.1), 
                                       fm = "cases ~ s(time, k = 100, bs = \'tp\', m = c(0, 1)) + s(day, k = 5)")
# Report: Figure 21a
infected_plot <- plot_country(extensionSwiss, indices, title="New infected", x="time", y="number cases", 
                              times = seq.Date(data_dates[1], as.Date(data_dates[length(indices)]) + prediction_window, by = "1 day"))
infected_plot

# Death curve:
extensionSwiss_deaths <- extension_prediction(swiss_data$new_deaths[indices], extension = list('len' = prediction_window, 'family' = nb(link = 'log'), 
                                                                             'keep_original' = TRUE, 'min_weights' = 0.1),
                                              fm = "cases ~ s(time, k = 20, bs = \'tp\', m = c(1, 0)) + s(day, k = 5)",
                                              weights = c(rep(1, length(swiss_data$new_deaths[indices])-1), 0))
# Report: Figure 21b
deaths_plot <- plot_country(extensionSwiss_deaths, indices, title="New deaths", x="time", y="number cases",
                            times = seq.Date(data_dates[1], as.Date(data_dates[length(indices)]) + prediction_window, by = "1 day"))
deaths_plot

# Curve of recovery cases:
recovered_smooth <- rollapply(swiss_data$new_recovered[indices], width=7, FUN=function(x) mean(x, na.rm=TRUE), by=1,
                            by.column=TRUE, partial=TRUE, fill=NA, align="center")
extensionSwiss_recovered <- extension_prediction(recovered_smooth, 
                                                 extension = list('len' = prediction_window, 'family' = nb(link = 'log'), 'keep_original' = TRUE, 'min_weights' = 0.1),
                                                 weights = c(rep(1, 290), rep(0, length(indices) - 290)),
                                                 fm = "cases ~ s(time, k = 20, bs = \'tp\', m = c(0, 0)) + s(day, k = 5)")
# Report: Figure 21c
recovered_plot <- plot_country(extensionSwiss_recovered, indices, title="New recovered, smoothed", x="time", y="number cases",
                               times = seq.Date(data_dates[1], as.Date(data_dates[length(indices)]) + prediction_window, by = "1 day"),
                               extra = c(swiss_data$new_recovered[indices], rep(NaN, prediction_window)))

recovered_plot <- recovered_plot + 
  geom_line(aes(y = c(swiss_data$new_recovered[indices], rep(NaN, prediction_window))), alpha = 0.8, color = "gray")

plot_general <- plot_grid(infected_plot, deaths_plot, recovered_plot, nrow = 3)
plot(plot_general)

# Results:
# the curves are quite noisy, and so the best deconvolution method will be the one using splines.
# The prediction for the next steps are ok (for the recovered cases, we took the total number of cases, otherwise the fit is not good).
# The log link works much better here.

# Model check (good)
par(mfrow = c(2, 2))
gam.check(extensionSwiss$fit)
summary(extensionSwiss$fit)



# ------------------------- 1) Deconvolution ---------------------  
# Incubation period distribution:
n_w_conv <- 20
shape_w_conv <- 10.2
scale_w_conv <- 0.5
w_conv<- discrete_gamma(n = n_w_conv, normalize = TRUE, shape = shape_w_conv, scale = scale_w_conv)

#                                   Original data ----------------
# Uncomment to see difference:
# smoothing_dim <- 1; smoothing_lambda <- 7;   # --> select high dimension basis, compensate by chosing large lambda for the penalization 
smoothing_dim <- 7; smoothing_lambda <- 7;   # --> small order, small lambda  
# smoothing_dim <- 1; smoothing_lambda <- 1;   # high dimension basis, low lambda  

# Now we search for the best basis dimension, as well as the best penalization based on the chosen dimension.
dimensionBasis <- findBestDimension(data_cases[indices], w_conv,
                                    family = quasipoisson(link = 'log'), smoothing = smoothing_dim, 
                                    extension = list('len' = 40, 'family' = nb(link = 'log'), 'keep_original' = TRUE, 'min_weights' = 0.1),
                                    range_values = c(10:150), diagnostic_plot = TRUE,
                                    fm = "cases ~ s(time, k = 20, bs = \'tp\', m = c(1, 0)) + s(day, k = 5)")
# dimensionBasis <- 38
lambda_delta <- findBestLambda(data_cases[indices], dimensionBasis, w_conv,
                               range_values = 10**(20:-10), family = quasipoisson(link = 'log'), smoothing = smoothing_lambda,
                               extension = list('len' = 40, 'family' = nb(link = 'log'), 'keep_original' = TRUE, 'min_weights' = 0.1), 
                               diagnostic_plot = TRUE, 
                               fm = "cases ~ s(time, k = 20, bs = \'tp\', m = c(1, 0)) + s(day, k = 5)")

delta <- construct_delta(dimensionBasis + 4)
delta <- lambda_delta * delta

# We use deconvolution with the extensions that we have computed before. 
results_IWLS_original <- DeconvolutionGAM(data_cases[indices], w_conv, dimensionBasis, plot_ = TRUE, delta = delta,
                                 get_derivatives = FALSE, family = quasipoisson(link = 'log'),
                                 extension = list('len' = 40, 'family' = nb(link = 'log'), 'keep_original' = TRUE, 'min_weights' = 0.1),
                                 fm = "cases ~ s(time, k = 20, bs = \'tp\', m = c(1, 0)) + s(day, k = 5)")
ConfidenceBands_original <- SimultaneousIntervals(results_IWLS_original$SplineBasis, plot_ = TRUE,
                                         coefficients_model = results_IWLS_original$fit$coefficients,
                                         CovMatrix = vcov(results_IWLS_original$fit),
                                         sample_splines = TRUE, ilink_fct = family(results_IWLS_original$fit)$linkinv)

ggplot(data = NULL, aes(x = data_dates[indices], y = data_cases[indices], colour = "Observed new cases")) + geom_line() + 
  geom_ribbon(aes(ymin = ConfidenceBands_original$Point_wise_CI_lower[indices], ymax = ConfidenceBands_original$Point_wise_CI_upper[indices], fill = "Point-wise confidence bands", colour = NA), alpha = 0.3, fill = colors[3]) +
  geom_ribbon(aes(ymin = ConfidenceBands_original$Simultaneous_lower[indices], ymax = ConfidenceBands_original$Simultaneous_upper[indices], fill = "Simlutaneous confidence bands", colour = NA), alpha = 0.3, fill = colors[3]) +
  geom_line(aes(y = results_IWLS_original$fitted_values[indices], colour = 'Deconvolution')) + 
  scale_colour_manual("", 
                      breaks = c("incidence", 'Observed new cases', 'Deconvolution', "Point-wise confidence bands", "Simlutaneous confidence bands"),
                      values = c("black", "black", colors[3], colors[3], colors[3])) +
  scale_fill_manual("", 
                    breaks = c("Point-wise confidence bands", "Simlutaneous confidence bands"),
                    values = c(colors[3], colors[3])) +
  labs(title=paste("Deconvolution, ", settings$country), x="Time", y="number of new cases", subtitle = latex2exp::TeX(paste("$q = ", dimensionBasis, "$ and $\\lambda = ", lambda_delta, "$"))) + 
  theme_for_the_plots + 
  # coord_cartesian(ylim = c(0, 1.1 * max(data_cases))) + 
  guides(color=guide_legend(override.aes=list(fill=NA))) + 
  theme(legend.position="bottom") + 
  # scale_y_log10() + 
  # scale_y_continuous(trans='log10') + 
  scale_x_date(date_breaks = "1 month", date_labels = "%b")    # or : date_labels = "%Y (%b)" 


ggplot(data = NULL, aes(x = data_dates[indices], y = (results_IWLS_original$fitted_values[indices] - data_cases[indices])/data_cases[indices])) + 
  geom_point() + 
  geom_hline(yintercept = 0, colour = "red") + 
  labs(title=latex2exp::TeX(paste("Ratio", "$ (\\hat{I}_t - I_t)/I_t $")), x="Time", y="ratio") + 
  theme_for_the_plots + 
  scale_x_date(date_breaks = "1 month", date_labels = "%b")    # or : date_labels = "%Y (%b)"


# First wave:
print((data_cases[c(1:150)] - results_IWLS_original$fitted_values[c(1:150)])[which(data_cases[c(1:150)] - results_IWLS_original$fitted_values[c(1:150)] == min(data_cases[c(1:150)] - results_IWLS_original$fitted_values[c(1:150)]))])
# around more than 600 new cases difference

# All time high:
ggplot(data = NULL, aes(x = data_dates[indices], y = data_cases[indices] - results_IWLS_original$fitted_values[indices])) + geom_point() + 
  geom_vline(xintercept = data_dates[indices][which(data_cases[indices] - results_IWLS_original$fitted_values[indices] == min(data_cases[indices] - results_IWLS_original$fitted_values[indices]))], color = "red") + 
  labs(title=paste("Difference between observed incidence curve and deconvolved one, for ", settings$country), x="Time", y="Difference in number of new cases") + 
  theme_for_the_plots +
  guides(color=guide_legend(override.aes=list(fill=NA))) + theme(legend.position="bottom") + 
  scale_x_date(date_breaks = "1 month", date_labels = "%b")    # or : date_labels = "%Y (%b)"
# More than 5000 new cases difference, on the 2020-10-25


#                                   Adding 1's at the beginning -----------------
# Here, we add some 1's before the data, to increase the domain definition for the spline basis and avoid wiggly splines at the beginning of the region of interest
extension_first_length <- 20
indices_extension_start <- c(c(1:extension_first_length), extension_first_length + indices)
data_cases_extension_start <- c(rep(1, extension_first_length), data_cases)
new_dates <- seq.Date(data_dates[1] - extension_first_length, length=length(indices) + extension_first_length, by='1 day' )
original_indices_translated <- extension_first_length + indices

smoothing_dim <- 7; smoothing_lambda <- 7;   # --> small order, small lambda   # Report: Figure 23a
# smoothing_dim <- 1; smoothing_lambda <- 7;   # --> select high dimension basis, compensate by chosing large lambda for the penalization   # Report: Figure 23b
# smoothing_dim <- 1; smoothing_lambda <- 1;   # high dimension basis, low lambda    # Report: Figure 23c


# Now we search for the best basis dimension, as well as the best penalization based on the chosen dimension.
dimensionBasis <- findBestDimension(data_cases_extension_start[indices_extension_start], w_conv,
                                    family = quasipoisson(link = 'log'), smoothing = smoothing_dim, 
                                    extension = list('len' = 40, 'family' = nb(link = 'log'), 'keep_original' = TRUE, 'min_weights' = 0.1),
                                    range_values = c(30:150), diagnostic_plot = FALSE,
                                    fm = "cases ~ s(time, k = 50, bs = \'tp\', m = c(0, 0)) + s(day, k = 5)")
lambda_delta <- findBestLambda(data_cases_extension_start[indices_extension_start], dimensionBasis, w_conv,
                               range_values = 10**(20:-10), family = quasipoisson(link = 'log'), smoothing = smoothing_lambda,
                               extension = list('len' = 40, 'family' = nb(link = 'log'), 'keep_original' = TRUE, 'min_weights' = 0.1), 
                               diagnostic_plot = TRUE, 
                               fm = "cases ~ s(time, k = 50, bs = \'tp\', m = c(0, 0)) + s(day, k = 5)")

delta <- construct_delta(dimensionBasis + 4)
delta <- lambda_delta * delta

length(data_cases_extension_start[indices_extension_start])/dimensionBasis

# We use deconvolution with the extensions that we have computed before. 
results_IWLS_extended <- DeconvolutionGAM(data_cases_extension_start[indices_extension_start], w_conv, dimensionBasis, plot_ = FALSE, delta = delta,
                                 get_derivatives = FALSE, family = quasipoisson(link = 'log'),
                                 extension = list('len' = 40, 'family' = nb(link = 'log'), 'keep_original' = TRUE, 'min_weights' = 0.1),
                                 fm = "cases ~ s(time, k = 100, bs = \'tp\', m = c(0, 1)) + s(day, k = 5)")
ConfidenceBands_extended <- SimultaneousIntervals(results_IWLS_extended$SplineBasis, plot_ = TRUE,
                                         coefficients_model = results_IWLS_extended$fit$coefficients,
                                         CovMatrix = vcov(results_IWLS_extended$fit),
                                         sample_splines = TRUE, ilink_fct = family(results_IWLS_extended$fit)$linkinv, N_keep = 100)

plot_deconvolution <- ggplot(data = NULL, aes(x = new_dates[original_indices_translated], y = data_cases_extension_start[original_indices_translated], colour = "Observed new cases")) + geom_line() + 
  geom_ribbon(aes(ymin = ConfidenceBands_extended$Point_wise_CI_lower[original_indices_translated], ymax = ConfidenceBands_extended$Point_wise_CI_upper[original_indices_translated], fill = "Point-wise confidence bands", colour = NA), alpha = 0.3, fill = colors[3]) +
  geom_ribbon(aes(ymin = ConfidenceBands_extended$Simultaneous_lower[original_indices_translated], ymax = ConfidenceBands_extended$Simultaneous_upper[original_indices_translated], fill = "Simlutaneous confidence bands", colour = NA), alpha = 0.3, fill = colors[3]) +
  geom_line(aes(y = results_IWLS_extended$fitted_values[original_indices_translated], colour = 'Deconvolution')) + 
  scale_colour_manual("", 
                      breaks = c("incidence", 'Observed new cases', 'Deconvolution', "Point-wise confidence bands", "Simlutaneous confidence bands"),
                      values = c("black", "black", colors[3], colors[3], colors[3])) +
  scale_fill_manual("", 
                    breaks = c("Point-wise confidence bands", "Simlutaneous confidence bands"),
                    values = c(colors[3], colors[3])) +
  labs(title=paste("Deconvolution, ", settings$country), x="Time", y="number of new cases", subtitle = latex2exp::TeX(paste("$q = ", dimensionBasis, "$ and $\\lambda = ", lambda_delta, "$"))) + 
  theme_for_the_plots + 
  coord_cartesian(ylim = c(0, 1.1 * max(data_cases)), xlim = c(data_dates[1], as.Date("2021-01-22"))) + 
  guides(color=guide_legend(override.aes=list(fill=NA))) + 
  # theme(legend.position="bottom") + 
  scale_x_date(date_breaks = "1 month", date_labels = "%b")    # or : date_labels = "%Y (%b)"

plot_deconvolution   # Report: Figure 15


# Here are some plots to see the difference between the true and observed incidence cases at time t:
ggplot(data = NULL, aes(x = new_dates[original_indices_translated], y = (results_IWLS_extended$fitted_values[original_indices_translated] - data_cases_extension_start[original_indices_translated])/data_cases_extension_start[original_indices_translated])) + 
  geom_point() + 
  geom_hline(yintercept = 0, colour = "red") + 
  labs(title=latex2exp::TeX(paste("Ratio", "$ (\\hat{I}_t - I_t)/I_t $")), x="Time", y="ratio") + 
  theme_for_the_plots + 
  coord_cartesian(ylim = c(-1, 10)) + 
  scale_x_date(date_breaks = "1 month", date_labels = "%b")    # or : date_labels = "%Y (%b)"

ggplot(data = NULL, aes(x = new_dates[original_indices_translated], y = (results_IWLS_extended$fitted_values[original_indices_translated] - data_cases_extension_start[original_indices_translated]))) + 
  geom_point() + 
  geom_hline(yintercept = 0, colour = "red") + 
  labs(title=latex2exp::TeX(paste("Difference", "$ (\\hat{I}_t - I_t) $")), x="Time", y="difference") + 
  theme_for_the_plots + 
  #coord_cartesian(ylim = c(-1, 10)) + 
  scale_x_date(date_breaks = "1 month", date_labels = "%b")    # or : date_labels = "%Y (%b)"


# And here is a plot showing also samples of spline functions for the deconvolution:
ConfidenceBands_extended$plot + coord_cartesian(ylim = c(0, 11000), xlim = c(20, 323)) 


#                                   Intermediate conclusion ------------

# Estimates for R_t are smaller at the beginning when 1's are included before the first index of the observations, and look more sensible. 
# In the following, we will thus use the results obtained in the previous last section:

Swiss_deconvolution <- results_IWLS_extended$fitted_values[original_indices_translated]
results_deconvolution <- results_IWLS_extended
condidence_bands <- ConfidenceBands_extended
indices_used <- original_indices_translated

# ------------------------- 2) Reproductive number estimation with GAM ---------------------  
# W_gam <- getGenerationInterval(4.2, 4.21, 30, TRUE)

# From Tapiwa et al. (2020)
W_gam <- discrete_normal()$prob_vector

#                                   Without deconvolution, raw data ------------ 
# Offset values:
offsetValues_data_init <- convolve(c(0* 1:(length(c(0,W_gam))-1), data_cases[indices]), rev(c(0,W_gam)), type = "filter")
offsetValues_data_init[which(offsetValues_data_init <= 1)] <- 1


# We start with the poisson model, log link. We set the spline basis dimension to 70, corresponding roughly to one parameter every 4 days (as of 26.12.20)
which_model <- poisson(link = "log")
DAT_data_init <- data.frame(resp = data_cases[indices], offs = offsetValues_data_init, x = c(1:length(data_cases[indices])))
fit_data_init <- fitModel(DAT_data_init, weights_ = NULL, which_model = which_model,  fm = 'resp ~ s(x, k = 70)', ignore_first_max = 5, 
                          deriv_plots = c(0,1), detect_flat = TRUE)

ggplot(data = NULL, aes(x = data_dates[indices])) +
  geom_ribbon(aes(ymin = fit_data_init$lower, ymax = fit_data_init$upper, fill = "95% CI", colour = "white"), alpha = 0.3, fill = colors[3]) +
  # geom_line(aes(y = data_init/offsetValues_data_init, colour = "Cori et al. estimate"), size = 1) + 
  geom_line(aes(y = fit_data_init$fitted.values, colour = "Estimated reproductive number"), size = 1) + 
  labs(title="Reproductive number", x="Time", y=latex2exp::TeX("$R_t$")) +
  scale_colour_manual("", 
                      breaks = c("True reproductive number", "Estimated reproductive number", "Cori et al. estimate", "95% CI"),
                      values = c("#1d3557", colors[3], colors[1], colors[3])) +
  scale_fill_manual("", breaks = c("95% CI"), values = colors[3]) + 
  geom_hline(yintercept = 1, col = "red") +
  coord_cartesian(ylim = c(0, 3)) + 
  theme_for_the_plots + 
  guides(color=guide_legend(override.aes=list(fill=NA))) + 
  scale_x_date(date_breaks = "1 month", date_labels = "%b")

# Results:
# Estimation very variable, small confidence bands, noisy

# With negative binomial
which_model <- nb(link = "log")
DAT_data_init <- data.frame(resp = data_cases[indices], offs = offsetValues_data_init, x = c(1:length(data_cases[indices])))
fit_data_init <- fitModel(DAT_data_init, weights_ = NULL, which_model = which_model,  fm = 'resp ~ s(x, k = 100)', ignore_first_max = 5, 
                          deriv_plots = c(0,1), detect_flat = TRUE)

ggplot(data = NULL, aes(x = data_dates[indices])) +
  geom_ribbon(aes(ymin = fit_data_init$lower, ymax = fit_data_init$upper, fill = "95% CI", colour = "white"), alpha = 0.3, fill = colors[3]) +
  # geom_line(aes(y = data_init/offsetValues_data_init, colour = "Cori et al. estimate"), size = 1) + 
  geom_line(aes(y = fit_data_init$fitted.values, colour = "Estimated reproductive number"), size = 1) + 
  labs(title="Reproductive number", x="Time", y=latex2exp::TeX("$R_t$")) +
  scale_colour_manual("", 
                      breaks = c("True reproductive number", "Estimated reproductive number", "Cori et al. estimate", "95% CI"),
                      values = c("#1d3557", colors[3], colors[1], colors[3])) +
  scale_fill_manual("", breaks = c("95% CI"), values = colors[3]) + 
  geom_hline(yintercept = 1, col = "red") +
  coord_cartesian(ylim = c(0, 3)) + 
  theme_for_the_plots + 
  guides(color=guide_legend(override.aes=list(fill=NA))) + 
  scale_x_date(date_breaks = "1 month", date_labels = "%b")

# Results:
# Estimates are smooth, with more realistic confidence intervals. The reproductive number went below the critical threshold of 1 from the end
# of March, and went above 1 from June to November. 
# Of course, here we assume that observed incidence cases are representative of the true one, which is false. We will compare these results
# with the ones obtained with a deconvolution step before.
# This method can however deal with noisy data without preprocessing the inputs, which is better than Cori et al. method in this case.
# As we will see later, the estimates shown in this figure are quite close to those obtained with deconvolution, except for the shift in time:
# overall, there is the same specific periods with the increase of Rt in June and in October, of the same magnitude.


#                                   With deconvolution ------------ 
offsetValues_IWLS <- convolve(c(0* 1:(length(W_gam)-1), Swiss_deconvolution), rev(W_gam), type = "filter")
offsetValues_IWLS[which(offsetValues_IWLS <= 0)] <- 1

# Weights for the fit
weights <- abs(condidence_bands$Point_wise_CI_upper[indices_used] - condidence_bands$Point_wise_CI_lower[indices_used])/4
weights <- weights**2
weights <- weights /sqrt(results_deconvolution$fitted_values[indices_used])
weights <- 1/weights
# weights <- capWeights(weights, 0.05)
weights <- weights/mean(weights)

# Plot for the weights:
ggplot(data = NULL, aes(x = data_dates[indices], y = weights)) + geom_point() + 
  theme_for_the_plots + labs(title = 'Weights', x = "time", y = "weight") + scale_x_date(date_breaks = "1 month", date_labels = "%b")


which_model <- nb(link = "log")
DAT_IWLS <- data.frame(resp = Swiss_deconvolution, offs = offsetValues_IWLS, x = c(1:length(results_deconvolution$fitted_values[indices_used])))
fit_IWLS <- fitModel(DAT_IWLS, which_model = which_model, weights_ = weights, fm = 'resp ~ s(x, k = 70)', ignore_first_max = 10,
                     deriv_plots = c(0,1), detect_flat = TRUE, axis_x = data_dates[indices])

# summary(fit_IWLS$fit)
# plot(fit_IWLS$derivPlots$plot)

# Report: Figure 16 (top)
plot_estim_r <- fit_IWLS$derivPlots$list_plots[[1]] + 
  labs(title = "", y = latex2exp::TeX(paste("Estimates of", "$R_t$")), x = "time") +
  coord_cartesian(ylim = c(0, 4.2)) +
  geom_segment(data=lockdown_data, mapping=aes(x=x, y=y, xend=vx, yend=vy), color="black") +
  geom_vline(data=lockdown_data,aes(xintercept=x), color = 'black', linetype = "dashed", alpha = 0.8) +
  geom_vline(data=lockdown_data,aes(xintercept=vx), color = 'black', linetype = "dashed", alpha = 0.8) +
  
  geom_point(data=lockdown_data, mapping=aes(x=x, y=y), size=2, color="black") +
  geom_point(data=lockdown_data, mapping=aes(x=vx, y=vy), size=2, color="black") + 
  geom_vline(data = NULL, aes(xintercept = as.Date("2020-06-22")), color = "red",linetype = "dashed", alpha = 0.8) + 
  geom_text(data=lockdown_data, aes(x=as.Date(lockdown_events$Start) + floor((as.Date(lockdown_events$End)-as.Date(lockdown_events$Start))/2), 
                                    y= y, label = name), colour = "black", vjust = -0.5, size = 4)
plot_estim_r

# Report: Figure 16 (bottom)
plot_estim_r_deriv <- fit_IWLS$derivPlots$list_plots[[2]] + 
  labs(title = "", y = latex2exp::TeX(paste("Estimated first derivative of", "$R_t$")), x = "time") +
  coord_cartesian(ylim = c(-0.2, 0.2)) +
  # geom_segment(data=lockdown_data, mapping=aes(x=x, y=y/20, xend=vx, yend=vy/20), color="black") +
  geom_vline(data=lockdown_data,aes(xintercept=x), color = 'black', linetype = "dashed", alpha = 0.8) +
  geom_vline(data=lockdown_data,aes(xintercept=vx), color = 'black', linetype = "dashed", alpha = 0.8) + 
  
  geom_point(data=lockdown_data, mapping=aes(x=x, y=y/20), size=0.002, color="black") +
  geom_point(data=lockdown_data, mapping=aes(x=vx, y=vy/20), size=0.002, color="black") +
  geom_vline(data = NULL, aes(xintercept = as.Date("2020-06-22")), color = "red",linetype = "dashed", alpha = 0.8)
  # geom_text(data=lockdown_data, aes(x=as.Date(lockdown_events$Start) + floor((as.Date(lockdown_events$End)-as.Date(lockdown_events$Start))/2), y= y/20, label = name), colour = "black", vjust = -0.5, size = 4)
plot_estim_r_deriv

# Report: Figure 16
plot_grid(plot_estim_r + theme(plot.margin = unit(c(0.5, 1, 0, 1), "cm")), 
          plot_estim_r_deriv + theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm")), 
          nrow = 2, align = 'hv')


# Final plot
plot_deconvolution <- plot_deconvolution + 
  labs(title = "", y = "Number of new cases", x = "time", subtitle = "") +
  geom_vline(data=lockdown_data,aes(xintercept=x), color = 'black', linetype = "dashed", alpha = 0.8) +
  geom_vline(data=lockdown_data,aes(xintercept=vx), color = 'black', linetype = "dashed", alpha = 0.8) + 
  
  geom_point(data=lockdown_data, mapping=aes(x=x, y=y/20), size=0.002, color="black") +
  geom_point(data=lockdown_data, mapping=aes(x=vx, y=vy/20), size=0.002, color="black") +
  geom_vline(data = NULL, aes(xintercept = as.Date("2020-06-22")), color = "red",linetype = "dashed", alpha = 0.8)
  
title <- ggdraw() + draw_label("COVID-19 in Switzerland", fontface='bold') + theme(plot.margin = unit(c(0.1, 1, -0.5, 1), "cm"))
plot_grid(title,
          plot_deconvolution + theme(plot.margin = unit(c(-0.5, 1, -0.5, 1), "cm")),
          plot_estim_r + theme(plot.margin = unit(c(-0.5, 1, -0.5, 1), "cm")), 
          plot_estim_r_deriv + theme(plot.margin = unit(c(-0.5, 1, 0.5, 1), "cm")), 
          nrow = 4, align = 'v', rel_heights = c(0.05, 1/2, 1/4, 1/4))




# Just a simpler plot:
plot_R_t <- ggplot(data = NULL, aes(x = data_dates[indices])) +
  geom_ribbon(aes(ymin = fit_IWLS$lower, ymax = fit_IWLS$upper, fill = "95% CI", colour = "white"), alpha = 0.3, fill = colors[3]) +
  geom_line(aes(y = fit_IWLS$fitted.values, colour = "Estimated reproductive number"), size = 1) + 
  # geom_line(aes(y = results_IWLS$fitted_values[indices]/offsetValues_IWLS, colour = "Cori et al. estimate"), size = 1) + 
  labs(title=paste("Estimation of the reproductive number for ", settings$country), x="Time", y=latex2exp::TeX("$R_t$")) +
  scale_colour_manual("", 
                      breaks = c("True reproductive number", "Estimated reproductive number", "Cori et al. estimate", "95% CI"),
                      values = c("#1d3557", colors[3], colors[1], colors[3])) +
  scale_fill_manual("", breaks = c("95% CI"), values = colors[3]) + 
  geom_hline(yintercept = 1, col = "red") +
  coord_cartesian(ylim = c(0, settings$range_legends[2]+0.5)) + 
  theme_for_the_plots + 
  guides(color=guide_legend(override.aes=list(fill=NA))) + 
  scale_x_date(date_breaks = "1 month", date_labels = "%b") + theme(legend.position = 'bottom')+ 
  geom_segment(data=lockdown_data, mapping=aes(x=x, y=y, xend=vx, yend=vy), color="black") +
  geom_vline(data=lockdown_data,aes(xintercept=x), color = 'black', linetype = "dashed", alpha = 0.8) +
  geom_vline(data=lockdown_data,aes(xintercept=vx), color = 'black', linetype = "dashed", alpha = 0.8) +
  
  geom_point(data=lockdown_data, mapping=aes(x=x, y=y), size=2, color="black") +
  geom_point(data=lockdown_data, mapping=aes(x=vx, y=vy), size=2, color="black") + 
  geom_text(data=lockdown_data, aes(x=as.Date(lockdown_events$Start) + floor((as.Date(lockdown_events$End)-as.Date(lockdown_events$Start))/2), 
                                    y= y, label = name), colour = "black", vjust = -0.5, size = 4)
plot(plot_R_t)


# ------------------------- 3) Reproductive number estimation with SIR model -----------------

# First, we get the data:
swiss_data_considered <- swiss_data[indices,]

# we move the deaths counts into the recovered one, as the are treated the sum way in the function
swiss_data_considered$deaths <- 0
swiss_data_considered$new_recovered <- swiss_data_considered$new_recovered + swiss_data_considered$new_deaths
swiss_data_considered$new_deaths <- 0

#                                   Without deconvolution, raw data ------------ 
# One can try to see the estimation without deconvolving the data:
result_gam_SIR <- gamSIR(swiss_data_considered, 100, 100, 8e6, family = poisson(link = 'identity'), 
                         first_ind = NULL, splineSimulations = FALSE, debug = TRUE, basic_r = FALSE)

ggplot(data = NULL, aes(x = swiss_data_considered$Date[result_gam_SIR$indices])) +
  geom_ribbon(aes(ymin = result_gam_SIR$lower, ymax = result_gam_SIR$upper, fill = "95% CI", colour = "white"), alpha = 0.3, fill = colors[3]) +
  # geom_path(data = result_gam_SIR$fits_splines, mapping = aes(y = values, x = age, group = ind), alpha = 0.4, colour = "grey20") +
  # geom_line(aes(y = df$true_rt[1:(length(df$true_rt)-1)][result_gam_SIR$indices], colour = "True reproductive number"), size = 1) + 
  geom_line(aes(y = result_gam_SIR$Reproductive_number, colour = "Estimated reproductive number"), size = 1) + 
  labs(title="Reproductive number", x="Time", y=latex2exp::TeX("$R_t$")) +
  scale_colour_manual("", 
                      breaks = c("True reproductive number", "Estimated reproductive number", "Cori et al. estimate", "95% CI"),
                      values = c("#1d3557", colors[3], colors[1], colors[3])) +
  scale_fill_manual("", breaks = c("95% CI"), values = colors[3]) + 
  geom_hline(yintercept = 1, col = "red") +
  coord_cartesian(ylim = c(0, 3)) + 
  theme_for_the_plots + 
  guides(color=guide_legend(override.aes=list(fill=NA)))

# Results: really bad...


result_gam_SIRExp <- gamSIRExp(swiss_data_considered, 100, 100, 8e6, first_ind = NULL, family = poisson(link = 'log'), 
                               splineSimulations = FALSE, debug = TRUE, basic_r = FALSE)

ggplot(data = NULL, aes(x = swiss_data_considered$Date[result_gam_SIRExp$indices])) +
  geom_ribbon(aes(ymin = result_gam_SIRExp$Reproductive_number$lower, ymax = result_gam_SIRExp$Reproductive_number$upper, fill = "95% CI", colour = "white"), alpha = 0.3, fill = colors[3]) +
  # geom_path(data = result_gam_SIRExp$fits_splines, mapping = aes(y = values, x = age, group = ind), alpha = 0.4, colour = "grey20") +
  # geom_line(aes(y = df$true_rt[1:(length(df$true_rt)-1)][result_gam_SIRExp$indices], colour = "True reproductive number"), size = 1) + 
  geom_line(aes(y = result_gam_SIRExp$Reproductive_number$fitted, colour = "Estimated reproductive number"), size = 1) + 
  labs(title="Reproductive number", x="Time", y=latex2exp::TeX("$R_t$")) +
  scale_colour_manual("", 
                      breaks = c("True reproductive number", "Estimated reproductive number", "Cori et al. estimate", "95% CI"),
                      values = c("#1d3557", colors[3], colors[1], colors[3])) +
  scale_fill_manual("", breaks = c("95% CI"), values = colors[3]) + 
  geom_hline(yintercept = 1, col = "red") +
  coord_cartesian(ylim = c(0, 3)) + 
  theme_for_the_plots + 
  guides(color=guide_legend(override.aes=list(fill=NA)))

# Results: same as before: bad

#                                   With deconvolution ------------ 
#                                            I) Deconvolution of the two curves ---------------
# Unlike the previous method, we need to make the deconvolution of two curves: the infected and the recovered (+deaths) ones.
# This implies that we need to have an estimate for the recovery interval. For now we take it as the incubation period, but after
# we will see how to change these two estimation

# Incubation period distribution:
n_w_conv <- 20
shape_w_conv <- 10
scale_w_conv <- 0.5
w_conv<- discrete_gamma(n = n_w_conv, normalize = TRUE, shape = shape_w_conv, scale = scale_w_conv)

# Recovered + death curves:
recovered_smooth <- rollapply(swiss_data_considered$new_recovered[indices], width=7, FUN=function(x) mean(x, na.rm=TRUE), by=1,
                              by.column=TRUE, partial=TRUE, fill=NA, align="center")
extensionSwiss_recovered <- extension_prediction(recovered_smooth, 
                                                 extension = list('len' = 20, 'family' = nb(link = 'log'), 'keep_original' = TRUE, 'min_weights' = 0.1),
                                                 fm = NULL,
                                                 weights = c(rep(1, 290), rep(0, length(indices) - 290)))
recovered_plot <- plot_country(extensionSwiss_recovered, indices, title="total recovered", x="time", y="number cases", 
                               times = seq.Date(data_dates[1], as.Date(data_dates[length(indices)]) + 20, by = "1 day"))
plot(recovered_plot)

# Now we search for the best basis dimension, as well as the best penalization based on the chosen dimension.
dimensionBasis <- findBestDimension(recovered_smooth, w_conv,
                                    family = poisson(link = 'log'), smoothing = 7, 
                                    extension = list('len' = 20, 'family' = nb(link = 'log'), 'keep_original' = TRUE, 'min_weights' = 0.1),
                                    range_values = c(30:100), diagnostic_plot = TRUE, 
                                    weights_extension = c(rep(1, 290), rep(0, length(indices) - 290)))
lambda_delta <- findBestLambda(recovered_smooth, dimensionBasis, w_conv,
                               range_values = 10**(20:-10), family = poisson(link = 'log'), smoothing = 7,
                               extension = list('len' = 20, 'family' = nb(link = 'log'), 'keep_original' = TRUE, 'min_weights' = 0.1), 
                               diagnostic_plot = TRUE,
                               weights_extension = c(rep(1, 290), rep(0, length(indices) - 290)))

delta <- construct_delta(dimensionBasis + 4)
delta <- lambda_delta * delta

# We use deconvolution with the extensions that we have computed before. 
results_IWLS_recovered <- DeconvolutionGAM(recovered_smooth, w_conv, dimensionBasis, plot_ = TRUE, delta = delta,
                                           get_derivatives = FALSE, family = poisson(link = 'log'),
                                           extension = list('len' = 20, 'family' = nb(link = 'log'), 'keep_original' = FALSE, 'min_weights' = 0.1),
                                           weights_extension = c(rep(1, 290), rep(0, length(indices) - 290)))
ConfidenceBands <- SimultaneousIntervals(results_IWLS_recovered$SplineBasis, plot_ = TRUE,
                                         coefficients_model = results_IWLS_recovered$fit$coefficients,
                                         CovMatrix = vcov(results_IWLS_recovered$fit),
                                         sample_splines = TRUE, ilink_fct = family(results_IWLS_recovered$fit)$linkinv)

plot_recovered <- ggplot(data = NULL, aes(x = swiss_data_considered$Date, y = swiss_data_considered$new_recovered, colour = "Observed new cases")) + geom_line() + 
  geom_ribbon(aes(ymin = ConfidenceBands$Point_wise_CI_lower[indices], ymax = ConfidenceBands$Point_wise_CI_upper[indices], fill = "Point-wise confidence bands", colour = NA), alpha = 0.3, fill = colors[3]) +
  geom_ribbon(aes(ymin = ConfidenceBands$Simultaneous_lower[indices], ymax = ConfidenceBands$Simultaneous_upper[indices], fill = "Simlutaneous confidence bands", colour = NA), alpha = 0.3, fill = colors[3]) +
  geom_line(aes(y = results_IWLS_recovered$fitted[indices], colour = 'Deconvolution')) + 
  scale_colour_manual("", 
                      breaks = c("incidence", 'Observed new cases', 'Deconvolution', "Point-wise confidence bands", "Simlutaneous confidence bands"),
                      values = c("black", "black", colors[3], colors[3], colors[3])) +
  scale_fill_manual("", 
                    breaks = c("Point-wise confidence bands", "Simlutaneous confidence bands"),
                    values = c(colors[3], colors[3])) +
  labs(title=paste("Deconvolution, ", settings$country), x="Time", y="number of new cases") + 
  theme_for_the_plots + coord_cartesian(ylim = c(0, 1.1 * max(data_cases))) + guides(color=guide_legend(override.aes=list(fill=NA))) + theme(legend.position="bottom") + 
  scale_x_date(date_breaks = "1 month", date_labels = "%b")    # or : date_labels = "%Y (%b)"

plot_recovered
# it is hard to say of the fitted values are realistic or not, given the poor data quality/ high variability of the recovery curve 


# Update the dataframe:
swiss_data_considered$new_confirmed <- Swiss_deconvolution
swiss_data_considered$confirmed <- cumsum(swiss_data_considered$new_confirmed)

swiss_data_considered$new_recovered <- results_IWLS_recovered$fitted[indices]
swiss_data_considered$recovered <- cumsum(swiss_data_considered$new_recovered)


#                                           II) SIR model ---------------
result_gam_SIR <- gamSIR(swiss_data_considered, 100, 100, 8e6, first_ind = NULL, family = quasipoisson(link = 'identity'), 
                               splineSimulations = FALSE, debug = TRUE, basic_r = FALSE, log_offset = FALSE)

ggplot(data = NULL, aes(x = swiss_data_considered$Date[result_gam_SIR$indices])) +
  geom_ribbon(aes(ymin = result_gam_SIR$lower, ymax = result_gam_SIR$upper, fill = "95% CI", colour = "white"), alpha = 0.3, fill = colors[3]) +
  # geom_path(data = result_gam_SIR$fits_splines, mapping = aes(y = values, x = age, group = ind), alpha = 0.4, colour = "grey20") +
  # geom_line(aes(y = df$true_rt[1:(length(df$true_rt)-1)][result_gam_SIR$indices], colour = "True reproductive number"), size = 1) + 
  geom_line(aes(y = result_gam_SIR$Reproductive_number, colour = "Estimated reproductive number"), size = 1) + 
  labs(title="Reproductive number", x="Time", y=latex2exp::TeX("$R_t$")) +
  scale_colour_manual("", 
                      breaks = c("True reproductive number", "Estimated reproductive number", "Cori et al. estimate", "95% CI"),
                      values = c("#1d3557", colors[3], colors[1], colors[3])) +
  scale_fill_manual("", breaks = c("95% CI"), values = colors[3]) + 
  geom_hline(yintercept = 1, col = "red") +
  coord_cartesian(ylim = c(0, 10)) + 
  theme_for_the_plots + 
  guides(color=guide_legend(override.aes=list(fill=NA)))


# Some experiments with the weights ....
# weights = abs(results_IWLS_confirmed$point_wise_lower[indices] - results_IWLS_confirmed$point_wise_upper[indices])/4
weights <- abs(condidence_bands$Point_wise_CI_upper[indices_used] - condidence_bands$Point_wise_CI_lower[indices_used])/4
weights <- weights**2
# weights <- weights/results_IWLS_confirmed$fitted_values[indices]    # because of quasipoisson -> variance no more equal to the mean
weights <- 1/weights
weights <- capWeights(weights, 0.01)
weights <- weights/mean(weights)
# plot(weights)

weights1 <- abs(condidence_bands$Point_wise_CI_upper[indices_used] - condidence_bands$Point_wise_CI_lower[indices_used])/4
# weights1 <- abs(results_IWLS_confirmed$point_wise_lower[indices] - results_IWLS_confirmed$point_wise_upper[indices])/4
weights2 <- abs(results_IWLS_recovered$point_wise_lower[indices] - results_IWLS_recovered$point_wise_upper[indices])/4
weights1 <- capWeights(weights1, 0.01)
weights2 <- capWeights(weights2, 0.01)
weights1 <- weights1**2
weights2 <- weights2**2
# weights1 <- weights1/results_IWLS_confirmed$fitted_values[indices]
# weights2 <- weights2/results_IWLS_recovered$fitted_values[indices]
weights1 <- 1/weights1
weights2 <- 1/weights2
weights1 <- weights1/mean(weights1)
weights2 <- weights2/mean(weights2)
weights_I <- weights1 + weights2

# weights for estimating beta (which depends on curve I (infection) and R (recovered), and there product, hence here the attempt to compute the variance for each term (and product terms))
V <- abs(condidence_bands$Point_wise_CI_upper[indices_used] - condidence_bands$Point_wise_CI_lower[indices_used])/4
V <- V**2
Vr <- abs(results_IWLS_recovered$point_wise_lower[indices] - results_IWLS_recovered$point_wise_upper[indices])/4
Vr <- Vr**2
weights_for_beta <- V - (3*V - V**2 - Swiss_deconvolution**4 - 2*V*condidence_bands$Point_wise_CI_lower[indices_used]**2)/(8e6)**2
weights_for_beta <- weights_for_beta - (V * Vr + V*results_IWLS_recovered$fitted_values[indices]**2 + Vr*Swiss_deconvolution**2)/(8e6)**2
weights_for_beta[which(weights_for_beta <= 0)] <- min(weights_for_beta[which(weights_for_beta > 0)])
weights_for_beta <- (1/weights_for_beta)
weights_for_beta <- weights_for_beta/mean(weights_for_beta)


result_gam_SIRExp <- gamSIRExp(swiss_data_considered, 100, 100, 8e6, first_ind = NULL, family = poisson(link = 'log'), 
                               splineSimulations = FALSE, debug = TRUE, basic_r = FALSE,
                               weights_for_r = weights1[1:(length(weights1)-1)], weights_for_i = weights_for_beta[1:(length(weights_for_beta)-1)])

ggplot(data = NULL, aes(x = swiss_data_considered$Date[result_gam_SIRExp$indices])) +
  geom_ribbon(aes(ymin = result_gam_SIRExp$Reproductive_number$lower, ymax = result_gam_SIRExp$Reproductive_number$upper, fill = "95% CI", colour = "white"), alpha = 0.3, fill = colors[3]) +
  # geom_path(data = result_gam_SIRExp$fits_splines, mapping = aes(y = values, x = age, group = ind), alpha = 0.4, colour = "grey20") +
  # geom_line(aes(y = df$true_rt[1:(length(df$true_rt)-1)][result_gam_SIRExp$indices], colour = "True reproductive number"), size = 1) + 
  geom_line(aes(y = result_gam_SIRExp$Reproductive_number$fitted, colour = "Estimated reproductive number"), size = 1) + 
  labs(title="Reproductive number", x="Time", y=latex2exp::TeX("$R_t$")) +
  scale_colour_manual("", 
                      breaks = c("True reproductive number", "Estimated reproductive number", "Cori et al. estimate", "95% CI"),
                      values = c("#1d3557", colors[3], colors[1], colors[3])) +
  scale_fill_manual("", breaks = c("95% CI"), values = colors[3]) + 
  geom_hline(yintercept = 1, col = "red") +
  coord_cartesian(ylim = c(0, 10)) + 
  theme_for_the_plots + 
  guides(color=guide_legend(override.aes=list(fill=NA))) + 
  scale_x_date(date_breaks = "1 month", date_labels = "%b")

# Results:
# High and low values of the reproductive number seem to be amplified in comparison with the GAM method.
# But the estimates are quite smooth


#                                          III) SIR model, with different delay mean ---------------

# Here, we assume that the incubation period is known (we will not change it here), but we will vary the mean 
# for the delay of recovery.

n_w_conv <- 20
shape_w_conv <- 10.2
scale_w_conv <- 0.5
w_conv<- discrete_gamma(n = n_w_conv, normalize = TRUE, shape = shape_w_conv, scale = scale_w_conv)

# initialise the data
swiss_data_considered <- swiss_data[indices,]
swiss_data_considered$deaths <- 0
swiss_data_considered$new_recovered <- swiss_data_considered$new_recovered + swiss_data_considered$new_deaths
swiss_data_considered$new_deaths <- 0

init_data_new_recovered <- swiss_data_considered$new_recovered[indices]
recovered_smooth <- rollapply(swiss_data_considered$new_recovered[indices], width=7, FUN=function(x) mean(x, na.rm=TRUE), by=1,
                              by.column=TRUE, partial=TRUE, fill=NA, align="center")

swiss_data_considered$new_confirmed <- Swiss_deconvolution
swiss_data_considered$confirmed <- cumsum(swiss_data_considered$new_confirmed)

# Now, we change the shape parameter:
w_conv_R_deconv_shape <- 2*c(1:20)
w_conv_R_deconv_scale <- 0.5+numeric(length(w_conv_R_deconv_shape))

l_score <- c()

best_score <- Inf
best_res <- NULL

pb <- txtProgressBar(min = 0, max = length(w_conv_R_deconv_shape), style = 3)
for (j in c(1:length(w_conv_R_deconv_shape))){
  w_conv_R_deconv <- discrete_gamma(n = 20, normalize = TRUE, shape = w_conv_R_deconv_shape[j], scale = w_conv_R_deconv_scale[j])
  
  # DECONVOLUTION:
  results_IWLS_recovered <- DeconvolutionGAM(recovered_smooth, w_conv_R_deconv, dimensionBasis, plot_ = FALSE, delta = delta,
                                             get_derivatives = FALSE, family = quasipoisson(link = 'log'),
                                             extension = list('len' = 20, 'family' = nb(link = 'log'), 'keep_original' = FALSE, 'min_weights' = 0.1),
                                             weights_extension = c(rep(1, 290), rep(0, length(indices) - 290)))
  
  swiss_data_considered$new_recovered <- results_IWLS_recovered$fitted_values[indices]
  swiss_data_considered$recovered <- cumsum(swiss_data_considered$new_recovered)
  
  # SIR model
  # result_gam_SIRExp <- gamSIRExp(swiss_data_considered, 100, 100, 8e6, first_ind = NULL, family = poisson(link = 'log'), 
  #                                splineSimulations = FALSE, debug = FALSE, basic_r = FALSE)
  weights1 <- abs(condidence_bands$Point_wise_CI_upper[indices_used] - condidence_bands$Point_wise_CI_lower[indices_used])/4
  weights1 <- capWeights(weights1, 0.01)
  weights1 <- weights1**2
  weights1 <- 1/weights1
  weights1 <- weights1/mean(weights1)
  
  V <- abs(condidence_bands$Point_wise_CI_upper[indices_used] - condidence_bands$Point_wise_CI_lower[indices_used])/4
  V <- V**2
  Vr <- abs(results_IWLS_recovered$point_wise_lower[indices] - results_IWLS_recovered$point_wise_upper[indices])/4
  Vr <- Vr**2
  weights_for_beta <- V - (3*V - V**2 - Swiss_deconvolution**4 - 2*V*condidence_bands$Point_wise_CI_lower[indices_used]**2)/(8e6)**2
  weights_for_beta <- weights_for_beta - (V * Vr + V*results_IWLS_recovered$fitted_values[indices]**2 + Vr*Swiss_deconvolution**2)/(8e6)**2
  weights_for_beta[which(weights_for_beta <= 0)] <- min(weights_for_beta[which(weights_for_beta > 0)])
  weights_for_beta <- (1/weights_for_beta)
  weights_for_beta <- weights_for_beta/mean(weights_for_beta)
  
  result_gam_SIRExp <- gamSIRExp(swiss_data_considered, 100, 100, 8e6, first_ind = 2, family = poisson(link = 'log'), 
                                 splineSimulations = FALSE, debug = FALSE, basic_r = FALSE,
                                 weights_for_r = weights1[1:(length(weights1)-1)], weights_for_i = weights_for_beta[1:(length(weights_for_beta)-1)])

  
  l_score[j] <- mean((result_gam_SIRExp$Reproductive_number$fitted - fit_IWLS$fitted.values[result_gam_SIRExp$indices])**2)
  if (best_score > l_score[j]){
    best_score <- l_score[j]
    best_res <- result_gam_SIRExp
  }
  setTxtProgressBar(pb, j)
}
close(pb)

sprintf("Best mean: %s", w_conv_R_deconv_shape[which(l_score == min(l_score))] * w_conv_R_deconv_scale[which(l_score == min(l_score))])
# plot(discrete_gamma(n = n_w_conv, normalize = TRUE, shape = w_conv_R_deconv_shape[which(l_score == min(l_score))], scale = w_conv_R_deconv_scale[which(l_score == min(l_score))]))


plot_R_t <- ggplot(data = NULL, aes(x = swiss_data_considered$Date[best_res$indices])) +
  geom_ribbon(aes(ymin = best_res$Reproductive_number$lower, 
                  ymax = best_res$Reproductive_number$upper, fill = "95% CI, SIR", colour = "white"), alpha = 0.3) +
  geom_ribbon(aes(ymin = fit_IWLS$lower[best_res$indices], 
                  ymax = fit_IWLS$upper[best_res$indices], fill = "95% CI, GAM", colour = "white"), alpha = 0.3) +
  # geom_path(data = best_res$fits_splines, mapping = aes(y = values, x = age, group = ind), alpha = 0.4, colour = "grey20") +
  # geom_line(aes(y = df$true_rt[1:(length(df$true_rt)-1)][best_res$indices], colour = "True reproductive number"), size = 1) + 
  geom_line(aes(y = best_res$Reproductive_number$fitted, colour = "Estimated reproductive number, SIR"), size = 1) + 
  geom_line(aes(y = fit_IWLS$fitted.values[best_res$indices], colour = "Estimated reproductive number, GAM"), size = 1) + 
  labs(title=paste("Estimated reproductive number for ", settings$country), x="Time", y=latex2exp::TeX("$R_t$")) +
  scale_colour_manual("", 
                      breaks = c("True reproductive number", "Estimated reproductive number, SIR", "Estimated reproductive number, GAM", "95% CI"),
                      values = c("#1d3557", colors[3], "#457b9d", colors[3])) +
  scale_fill_manual("", breaks = c("95% CI, GAM", "95% CI, SIR"), 
                    values = c("#457b9d", colors[3])) + 
  geom_hline(yintercept = 1, col = "red") +
  coord_cartesian(ylim = c(0, 4.5)) + 
  theme_for_the_plots + 
  guides(color=guide_legend(override.aes=list(fill=NA))) + 
  scale_x_date(date_breaks = "1 month", date_labels = "%b") + theme(legend.position = "bottom")

plot_R_t <- plot_R_t + geom_segment(data=lockdown_data, mapping=aes(x=x, y=y, xend=vx, yend=vy), color="black") +
  geom_vline(data=lockdown_data,aes(xintercept=x), color = 'black', linetype = "dashed", alpha = 0.8) +
  geom_vline(data=lockdown_data,aes(xintercept=vx), color = 'black', linetype = "dashed", alpha = 0.8) +
  
  geom_point(data=lockdown_data, mapping=aes(x=x, y=y), size=2, color="black") +
  geom_point(data=lockdown_data, mapping=aes(x=vx, y=vy), size=2, color="black") + 
  geom_text(data=lockdown_data, aes(x=as.Date(lockdown_events$Start) + floor((as.Date(lockdown_events$End)-as.Date(lockdown_events$Start))/2), 
                                    y= y, label = name), colour = "black", vjust = -0.5, size = 4)

# Report: Figure 17
plot(plot_R_t)



# Transmission rate:
transmission_rate_plot <- ggplot(data = NULL, aes(x = swiss_data_considered$Date[best_res$indices])) +
  geom_ribbon(aes(ymin = best_res$transmission_rate$lower, 
                  ymax = best_res$transmission_rate$upper, fill = "95% CI, SIR", colour = "white"), alpha = 0.3) +
  # geom_path(data = best_res$fits_splines, mapping = aes(y = values, x = age, group = ind), alpha = 0.4, colour = "grey20") +
  # geom_line(aes(y = df$true_rt[1:(length(df$true_rt)-1)][best_res$indices], colour = "True reproductive number"), size = 1) + 
  geom_line(aes(y = best_res$transmission_rate$fitted, colour = "Estimated transmission rate"), size = 1) + 
  labs(title=paste("Estimated transmission rate for ", settings$country), x="Time", y=latex2exp::TeX("$\\beta_t$")) +
  scale_colour_manual("", 
                      breaks = c("True reproductive number", "Estimated transmission rate", "Estimated reproductive number, GAM", "95% CI"),
                      values = c("#1d3557", colors[3], "#457b9d", colors[3])) +
  scale_fill_manual("", breaks = c("95% CI, GAM", "95% CI, SIR"), 
                    values = c("#457b9d", colors[3])) + 
  # geom_hline(yintercept = 1, col = "red") +
  # coord_cartesian(ylim = c(0, 4.5)) + 
  theme_for_the_plots + 
  guides(color=guide_legend(override.aes=list(fill=NA))) + 
  scale_x_date(date_breaks = "1 month", date_labels = "%b") + theme(legend.position = "bottom") 

# removal rate:
removal_rate_plot <- ggplot(data = NULL, aes(x = swiss_data_considered$Date[best_res$indices])) +
  geom_ribbon(aes(ymin = best_res$removal_rate$lower, 
                  ymax = best_res$removal_rate$upper, fill = "95% CI, SIR", colour = "white"), alpha = 0.3) +
  # geom_path(data = best_res$fits_splines, mapping = aes(y = values, x = age, group = ind), alpha = 0.4, colour = "grey20") +
  # geom_line(aes(y = df$true_rt[1:(length(df$true_rt)-1)][best_res$indices], colour = "True reproductive number"), size = 1) + 
  geom_line(aes(y = best_res$removal_rate$fitted, colour = "Estimated removal rate"), size = 1) + 
  labs(title=paste("Estimated removal rate for ", settings$country), x="Time", y=latex2exp::TeX("$\\gamma_t$")) +
  scale_colour_manual("", 
                      breaks = c("True reproductive number", "Estimated removal rate", "Estimated reproductive number, GAM", "95% CI"),
                      values = c("#1d3557", colors[3], "#457b9d", colors[3])) +
  scale_fill_manual("", breaks = c("95% CI, GAM", "95% CI, SIR"), 
                    values = c("#457b9d", colors[3])) + 
  # geom_hline(yintercept = 1, col = "red") +
  # coord_cartesian(ylim = c(0, 4.5)) + 
  theme_for_the_plots + 
  guides(color=guide_legend(override.aes=list(fill=NA))) + 
  scale_x_date(date_breaks = "1 month", date_labels = "%b") + theme(legend.position = "bottom") 


plot_grid(transmission_rate_plot + theme(plot.margin = unit(c(0.5, 1, 0, 1), "cm")), 
          removal_rate_plot + theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm")), 
          nrow = 2, align = 'hv')


# Without constraints, but one step estimation:
w_conv_R_deconv <- discrete_gamma(n = 20, normalize = TRUE, shape = w_conv_R_deconv_shape[which(l_score == min(l_score))], scale = w_conv_R_deconv_scale[which(l_score == min(l_score))])
results_IWLS_recovered <- DeconvolutionGAM(recovered_smooth, w_conv_R_deconv, dimensionBasis, plot_ = FALSE, delta = delta,
                                           get_derivatives = FALSE, family = quasipoisson(link = 'log'),
                                           extension = list('len' = 20, 'family' = nb(link = 'log'), 'keep_original' = FALSE, 'min_weights' = 0.1),
                                           weights_extension = c(rep(1, 290), rep(0, length(indices) - 290)))
swiss_data_considered$new_recovered <- results_IWLS_recovered$fitted_values[indices]#[(length(w_conv_R_deconv)+1):length(results_swiss_recovered)]
swiss_data_considered$recovered <- cumsum(swiss_data_considered$new_recovered)

result_gam_SIR <- gamSIR(swiss_data_considered, 100, 100, 8e6, family = quasipoisson(link = 'identity'), 
                         first_ind = NULL, splineSimulations = FALSE, debug = TRUE, basic_r = FALSE)

ggplot(data = NULL, aes(x = swiss_data_considered$Date[result_gam_SIR$indices])) +
  geom_ribbon(aes(ymin = result_gam_SIR$lower, ymax = result_gam_SIR$upper, fill = "95% CI", colour = "white"), alpha = 0.3, fill = colors[3]) +
  geom_line(aes(y = result_gam_SIR$Reproductive_number, colour = "Estimated reproductive number"), size = 1) + 
  labs(title=paste("Estimation of the reproductive number for ", settings$country), x="Time", y=latex2exp::TeX("$R_t$")) +
  scale_colour_manual("", 
                      breaks = c("True reproductive number", "Estimated reproductive number", "Cori et al. estimate", "95% CI"),
                      values = c("#1d3557", colors[3], colors[1], colors[3])) +
  scale_fill_manual("", breaks = c("95% CI"), values = colors[3]) + 
  geom_hline(yintercept = 1, col = "red") +
  coord_cartesian(ylim = c(0, 3)) + 
  theme_for_the_plots + 
  guides(color=guide_legend(override.aes=list(fill=NA)))





# Now that we have the shape, we can try to change also the scale parameter:

# Now, we change the shape parameter:
w_conv_R_deconv_scale <- seq(0.1, 2, 0.1)
w_conv_R_deconv_shape <- w_conv_R_deconv_shape[which(l_score == min(l_score))]+numeric(length(w_conv_R_deconv_scale))

l_score_s <- c()

best_score_s <- Inf
best_res_s <- NULL

pb <- txtProgressBar(min = 0, max = length(w_conv_R_deconv_shape), style = 3)
for (j in c(1:length(w_conv_R_deconv_shape))){
  w_conv_R_deconv <- discrete_gamma(n = 20, normalize = TRUE, shape = w_conv_R_deconv_shape[j], scale = w_conv_R_deconv_scale[j])
  
  # DECONVOLUTION:
  results_IWLS_recovered <- DeconvolutionGAM(recovered_smooth, w_conv_R_deconv, dimensionBasis, plot_ = FALSE, delta = delta,
                                             get_derivatives = FALSE, family = quasipoisson(link = 'log'),
                                             extension = list('len' = 20, 'family' = nb(link = 'log'), 'keep_original' = FALSE, 'min_weights' = 0.1),
                                             weights_extension = c(rep(1, 290), rep(0, length(indices) - 290)))
  
  swiss_data_considered$new_recovered <- results_IWLS_recovered$fitted_values[indices]
  swiss_data_considered$recovered <- cumsum(swiss_data_considered$new_recovered)
  
  # SIR model
  weights1 <- abs(condidence_bands$Point_wise_CI_upper[indices_used] - condidence_bands$Point_wise_CI_lower[indices_used])/4
  weights1 <- capWeights(weights1, 0.01)
  weights1 <- weights1**2
  weights1 <- 1/weights1
  weights1 <- weights1/mean(weights1)
  
  V <- abs(condidence_bands$Point_wise_CI_upper[indices_used] - condidence_bands$Point_wise_CI_lower[indices_used])/4
  V <- V**2
  Vr <- abs(results_IWLS_recovered$point_wise_lower[indices] - results_IWLS_recovered$point_wise_upper[indices])/4
  Vr <- Vr**2
  weights_for_beta <- V - (3*V - V**2 - Swiss_deconvolution**4 - 2*V*condidence_bands$Point_wise_CI_lower[indices_used]**2)/(8e6)**2
  weights_for_beta <- weights_for_beta - (V * Vr + V*results_IWLS_recovered$fitted_values[indices]**2 + Vr*Swiss_deconvolution**2)/(8e6)**2
  weights_for_beta[which(weights_for_beta <= 0)] <- min(weights_for_beta[which(weights_for_beta > 0)])
  weights_for_beta <- (1/weights_for_beta)
  weights_for_beta <- weights_for_beta/mean(weights_for_beta)
  
  result_gam_SIRExp <- gamSIRExp(swiss_data_considered, 100, 100, 8e6, first_ind = 2, family = poisson(link = 'log'), 
                                 splineSimulations = FALSE, debug = FALSE, basic_r = FALSE,
                                 weights_for_r = weights1[1:(length(weights1)-1)], weights_for_i = weights_for_beta[1:(length(weights_for_beta)-1)])

  l_score_s[j] <- mean((result_gam_SIRExp$Reproductive_number$fitted - fit_IWLS$fitted.values[result_gam_SIRExp$indices])**2)
  if (best_score > l_score_s[j]){
    best_score_s <- l_score_s[j]
    best_res_s <- result_gam_SIRExp
  }
  setTxtProgressBar(pb, j)
}
close(pb)

sprintf("Best mean: %s", w_conv_R_deconv_shape[which(l_score_s == min(l_score_s))] * w_conv_R_deconv_scale[which(l_score_s == min(l_score_s))])
# plot(discrete_gamma(n = n_w_conv, normalize = TRUE, shape = w_conv_R_deconv_shape[which(l_score == min(l_score))], scale = w_conv_R_deconv_scale[which(l_score == min(l_score))]))





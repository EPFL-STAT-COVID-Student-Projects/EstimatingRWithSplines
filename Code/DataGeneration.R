# ===================== PACKAGES =====================

library(zoo)   
library(ggplot2)
library(shinythemes)
library(mgcv)  
require(stats)
library(splines)
library(cowplot)

# ---------------------------- THEMES ---------
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

# ===================== Generated data =====================
source("/Users/Antoine/Desktop/Projet Semestre 2020/CODE/simulation.R")
source("/Users/Antoine/Desktop/Projet Semestre 2020/CODE/Deconvolution.R")


# ------------------------- 1) Synthetic data and reproductive number ----------------------------------
# Generate data
sim_list <- simulation()
df <- sim_list$sim_df

# This is the set of data that we will use, i.e. the incidence curve. However, there is some delay between an infection event and the time 
# it is reported. This will be modelled by convolving the data.
data_init <- df$incidence[2:(length(df$incidence)-1)]  # first, we remove the NA

# We assume that the infectivity profile is a discretized gamma, with mean 5 days.
n_w_conv <- 20
shape_w_conv <- 10
scale_w_conv <- 0.5
w_conv<- discrete_gamma(n = n_w_conv, normalize = TRUE, shape = shape_w_conv, scale = scale_w_conv)

# Then we convolve the data. This is what we will observed and use for deconvolving the data and estimate R (after)
data_convolved <- convolveData(w_conv, data_init)

# plot(cumsum(data_init), type = 'l')
# lines(cumsum(data_convolved), col = "blue")

# DATASET
# Report: Figure 1a
ggplot(data = NULL, aes(x = 1:length(data_init))) + 
  annotate(geom = "rect", xmin = 60, xmax = 67, ymin = -Inf, ymax = Inf, fill = "#2a9d8f", alpha = 0.2)+
  annotate(geom = "rect", xmin = 90, xmax = 97, ymin = -Inf, ymax = Inf, fill = "#e63946", alpha = 0.2)+
  geom_line(aes(y = data_init, colour = 'Incidence curve'), size = 0.7) + 
  geom_line(aes(y = data_convolved, colour = 'observed incidence'), size = 0.7) + 
  scale_colour_manual("", 
                      breaks = c("Incidence curve", "S", "E", "I", "R", 'observed incidence'),
                      values = c("#1d3557", "gray", colors[1], colors[2], colors[3], '#2a9d8f')) +
  labs(title="Generated SIR data", x="Time", y="Number of new cases") + 
  theme_for_the_plots + 
  xlim(0, 250)

# Report: Figure 1b
# TRUE Reproductive number R_0 and R_t
ggplot(data = NULL, aes(x = 1:length(df$true_rt[2:(length(df$incidence)-1)]))) + #geom_line(size = 1) +
  annotate(geom = "rect", xmin = 60, xmax = 67, ymin = -Inf, ymax = Inf, fill = "#2a9d8f", alpha = 0.2)+
  annotate(geom = "rect", xmin = 90, xmax = 97, ymin = -Inf, ymax = Inf, fill = "#e63946", alpha = 0.2)+
  # geom_segment(aes(x = 60, y = 2.2, xend = 67, yend = 2.2), arrow = arrow(length = unit(0.2, "cm"))) + 
  # geom_segment(aes(x = 90, y = 2.2, xend = 97, yend = 2.2), arrow = arrow(length = unit(0.2, "cm"))) + 
  geom_segment(aes(x = 60, y = 2.1, xend = 90, yend = 2.1), color = "#e63946", arrow = arrow(length = unit(0.2, "cm"))) + 
  geom_segment(aes(x = 90, y = 2.1, xend = 60, yend = 2.1), color = "#e63946", arrow = arrow(length = unit(0.2, "cm"))) + 
  geom_text(data=NULL, aes(x=75, y= 2.1, label = "30 days"), colour = "#e63946", vjust = -0.5, size = 3) + 

  # geom_vline(xintercept = 60, color = "#a8dadc", linetype = "twodash", size = 0.7)+
  # geom_vline(xintercept = 67, color = "#a8dadc", linetype = "twodash", size = 0.7)+
  # geom_vline(xintercept = 90, color = "#f1faee", linetype = "twodash", size = 0.7)+
  # geom_vline(xintercept = 97, color = "#f1faee", linetype = "twodash", size = 0.7) + 
  geom_hline(yintercept = 1, color = "red") + 
  geom_line(aes(y = df$true_r0[2:(length(df$incidence)-1)], colour = "R0"), size = 1) + 
  geom_line(aes(y = df$true_rt[2:(length(df$incidence)-1)], colour = "R_t"), size = 1) + 
  scale_colour_manual("", 
                      breaks = c("incidence", "S", "E", "I", "R", 'observed incidence', "R_t", "R0"),
                      values = c("black", "gray", colors[1], colors[2], colors[3], colors[2], "#1d3557", "#e63946")) +
  labs(title="Reproductive number for generated SIR data", x="time", y="Reproductive number") + 
  theme_for_the_plots  +
  xlim(0, 250)


# ------------------------- 2) Adding weekend effects ----------------------------------

# Takes the matrix of delay distribution (constructed with the function 'construct_B'), and add weekend effects
weekend_effects <- function(B, d_factor = 0.6){
  # data_conv_with_weekend_effects <- data_conv
  # Times <- c(1:length(data_conv))
  # for (j in Times){
  #   if (j %% 7 == 5){
  #     data_conv_with_weekend_effects[j] <- 0.6*data_conv[j]
  #     data_conv_with_weekend_effects[j+1] <- data_conv_with_weekend_effects[j+1] + 0.2*data_conv[j]
  #     data_conv_with_weekend_effects[j+2] <- data_conv_with_weekend_effects[j+2] + 0.2*data_conv[j]
  #   }else if(j %% 7 == 6){
  #     data_conv_with_weekend_effects[j] <- 0.6*data_conv[j]
  #     data_conv_with_weekend_effects[j+1] <- data_conv_with_weekend_effects[j+1] + 0.4*data_conv[j]
  #   }
  # }
  for (j in c(1:dim(B)[1])){
    if (j %% 7 == 5){
      B[j,] <- d_factor*B[j,]
    }else if (j %% 7 == 6){
      B[j,] <- 0.4*B[j,]
    }
  }
  for (j in c(1:dim(B)[2])){
    B[,j] <- B[,j]/sum(B[,j])
  }
  return(B)
}

B_weeks <- construct_B(length(data_init), w_conv)
B_weeks <- weekend_effects(B_weeks, 0.6)

data_convolved_weekend_effects <- B_weeks %*% data_init

data_convolved_weekend_effects_smoothed <- rollapply(data_convolved_weekend_effects, width=7, FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                         by.column=TRUE, partial=TRUE, fill=NA, align="center")

# Report: Figure 1c
ggplot(data = NULL, aes(x = 1:length(data_convolved_weekend_effects))) + 
  # annotate(geom = "rect", xmin = 60, xmax = 67, ymin = -Inf, ymax = Inf, fill = "#2a9d8f", alpha = 0.1)+
  # annotate(geom = "rect", xmin = 90, xmax = 97, ymin = -Inf, ymax = Inf, fill = "#e63946", alpha = 0.1)+
  geom_line(aes(y = data_init, colour = 'Incidence curve'), size = 0.7) + 
  geom_line(aes(y = data_convolved, colour = 'observed incidence'), size = 0.7) + 
  geom_line(aes(y = data_convolved_weekend_effects, colour = 'observed incidence, with weekend effects'), size = 0.5, alpha = 0.5) + 
  geom_line(aes(y = data_convolved_weekend_effects_smoothed, colour = 'observed incidence, with weekend effects, smoothed'), size = 0.7, alpha = 1, linetype = 2) + 
  geom_point(aes(x = which(c(1:dim(B_weeks)[1]) %% 7 == 6), 
                 y = data_convolved_weekend_effects[which(c(1:dim(B_weeks)[1]) %% 7 == 6)],
                 colour = "Sundays"), size = 0.5) + 
  scale_colour_manual("", 
                      breaks = c("Incidence curve", "S", "E", "I", "R", 'observed incidence', 
                                 'observed incidence, with weekend effects', 'observed incidence, with weekend effects, smoothed', "Sundays"),
                      values = c("#1d3557", "gray", colors[1], colors[2], colors[3], '#d90429', "#335c67", "#fdc500", "#1d3557"),
                      guide = guide_legend(override.aes = list(
                        linetype = c(rep("solid", 4), "blank"),
                        shape = c(rep(NA, 4), 16)))) +
  labs(title="Generated SIR data", x="Time", y="Number of new cases") + 
  theme_for_the_plots + 
  xlim(0, 250)


indices_original_synthetic <- c(1:250)  # to remove increase in the end (for the weekend effects)
data_convolved <- data_convolved[indices_original_synthetic]
data_convolved_weekend_effects <- data_convolved_weekend_effects[indices_original_synthetic]
data_convolved_weekend_effects_smoothed <- data_convolved_weekend_effects_smoothed[indices_original_synthetic]
data_init <- data_init[indices_original_synthetic]

# ===================== I) Deconvolution methods =====================
# For the deconvolution step, we compare 3 methods. We asssume that we know the true infectivity profile


# ------------------------- 1) Roll back -------------------------

# We call the function rollBack
data_country_roll_back <- rollBack(data_convolved, w_conv)
data_country_roll_back <- data_country_roll_back[(n_w_conv+1):length(data_country_roll_back)]

results_roll_back <- rollBackCI(data_convolved, w_conv, 100)

# Report: Figure 3a
ggplot(data = NULL, aes(x = 1:length(data_init), y = data_init, colour = "incidence")) + geom_point() + 
  geom_line(aes(y = data_convolved, colour = 'observed incidence')) + 
  geom_ribbon(aes(ymin = results_roll_back$q2.5, ymax = results_roll_back$q97.5, fill = "Confidence interval roll back", colour = NA), alpha = 0.3, fill = colors[3]) +
  geom_line(aes(y = results_roll_back$mean, colour = 'Roll back mean')) + 
  scale_colour_manual("", 
                      breaks = c("incidence", 'observed incidence', 'Roll back mean', "Confidence interval roll back"),
                      values = c("black", colors[2], colors[3], colors[3])) +
  scale_fill_manual("", 
                      breaks = c("Confidence interval roll back"),
                      values = c(colors[3])) +
  labs(title="Roll back deconvolution on generated SIR data", x="Time", y="number of new cases") + 
  theme_for_the_plots + guides(color=guide_legend(override.aes=list(fill=NA)))

# Results:
# (+) Easy to apply, fast to compute, fast to optain confidence bands
# (-) tends to convolving again the data, but in the other direction. This implies that the spikes will be flatten
# (-) the end is heavily under estimating the number of cases
# (-) true incidence curve not often in the confidence bands...

# solution for the end:
fit_splines <- smooth.spline(c(1:length(data_convolved)), data_convolved)
next_values <- predict(fit_splines, x = c((length(data_convolved)+1):(length(data_convolved)+20)))

data_convolved_extended <- c(data_convolved, next_values$y)
  
results_roll_back_ext <- rollBackCI(data_convolved_extended, w_conv, 100)
data_init_extended <- c(data_init, rep(NA, 20))

# Report: Figure 3b
ggplot(data = NULL, aes(x = 1:length(data_init_extended), y = data_init_extended, colour = "incidence")) + geom_point() + 
  geom_line(aes(x = 1:length(data_convolved_extended), y = data_convolved_extended, colour = 'observed incidence (extended)')) + 
  geom_ribbon(aes(x = 1:length(data_convolved_extended), ymin = results_roll_back_ext$q2.5, ymax = results_roll_back_ext$q97.5, fill = "Confidence interval roll back (extended)", colour = NA), alpha = 0.3, fill = colors[3]) +
  geom_line(aes(x = 1:length(data_convolved_extended), y = results_roll_back_ext$mean, colour = 'Roll back mean (extended)')) + 
  geom_vline(xintercept = length(data_init), color = "red") + 
  scale_colour_manual("", 
                      breaks = c("incidence", 'observed incidence (extended)', 'Roll back mean (extended)', "Confidence interval roll back (extended)"),
                      values = c("black", colors[2], colors[3], colors[3])) +
  scale_fill_manual("", 
                    breaks = c("Confidence interval roll back (extended)"),
                    values = c(colors[3])) +
  labs(title="Roll back deconvolution on generated SIR data, with extension", x="Time", y="number of new cases") + 
  theme_for_the_plots + 
  guides(color=guide_legend(override.aes=list(fill=NA)))

# Results:
# This solves the problem for the end part, although this might not be a solution for every types of incidence curve, as this one is relatively smooth at the end.

# - - - - - - - - - - - - - - - - - Noisy curve ----------------------------------------------
data_country_roll_back <- rollBack(data_convolved_weekend_effects, w_conv)
data_country_roll_back <- data_country_roll_back[(n_w_conv+1):length(data_country_roll_back)]

results_roll_back <- rollBackCI(data_convolved_weekend_effects, w_conv, 100)
# Report: Figure 6a
ggplot(data = NULL, aes(x = 1:length(data_init), y = data_init, colour = "incidence")) + geom_point() + 
  geom_line(aes(y = data_convolved_weekend_effects, colour = 'observed incidence')) + 
  geom_ribbon(aes(ymin = results_roll_back$q2.5, ymax = results_roll_back$q97.5, fill = "Confidence interval roll back", colour = NA), alpha = 0.3, fill = colors[3]) +
  geom_line(aes(y = results_roll_back$mean, colour = 'Roll back mean')) + 
  scale_colour_manual("", 
                      breaks = c("incidence", 'observed incidence', 'Roll back mean', "Confidence interval roll back"),
                      values = c("black", colors[2], colors[3], colors[3])) +
  scale_fill_manual("", 
                    breaks = c("Confidence interval roll back"),
                    values = c(colors[3])) +
  labs(title="Roll back deconvolution on generated SIR data", x="Time", y="number of new cases") + 
  theme_for_the_plots + guides(color=guide_legend(override.aes=list(fill=NA)))

# Results:
# Weekly pattern is still very visible, confidence bands too narrow, and drop near the end of the curve

# Smoothed noisy curve first
data_country_roll_back <- rollBack(data_convolved_weekend_effects_smoothed, w_conv)
data_country_roll_back <- data_country_roll_back[(n_w_conv+1):length(data_country_roll_back)]

results_roll_back <- rollBackCI(data_convolved_weekend_effects_smoothed, w_conv, 100)
ggplot(data = NULL, aes(x = 1:length(data_init), y = data_init, colour = "incidence")) + geom_point() + 
  geom_line(aes(y = data_convolved_weekend_effects_smoothed, colour = 'observed incidence')) + 
  geom_ribbon(aes(ymin = results_roll_back$q2.5, ymax = results_roll_back$q97.5, fill = "Confidence interval roll back", colour = NA), alpha = 0.3, fill = colors[3]) +
  geom_line(aes(y = results_roll_back$mean, colour = 'Roll back mean')) + 
  scale_colour_manual("", 
                      breaks = c("incidence", 'observed incidence', 'Roll back mean', "Confidence interval roll back"),
                      values = c("black", colors[2], colors[3], colors[3])) +
  scale_fill_manual("", 
                    breaks = c("Confidence interval roll back"),
                    values = c(colors[3])) +
  labs(title="Roll back deconvolution on generated SIR data", x="Time", y="number of new cases") + 
  theme_for_the_plots + guides(color=guide_legend(override.aes=list(fill=NA)))

# Results:
# (+) No more weekly pattern (since it was remove during the smoothing)
# (-) Large underestimation of the number of cases between time 60 and 70. Also, transitions are smoother
# Basically, we are smoothing a convolution...


# For the drop at the end:
data_country_roll_back_extension <- extension_prediction(data_convolved_weekend_effects, 
                                                         extension = list('len' = 20, 'family' = nb(link = 'log'), 
                                                                          'keep_original' = TRUE, 'min_weights' = 0.1))

data_country_roll_back <- rollBack(data_country_roll_back_extension$extension_results$fitted, w_conv)
data_country_roll_back <- data_country_roll_back[(n_w_conv+1):length(data_country_roll_back)]

results_roll_back <- rollBackCI(data_country_roll_back_extension$extension_results$fitted, w_conv, 100)

ggplot(data = NULL, aes(x = 1:length(data_country_roll_back_extension$DATA$cases), y = c(data_init, rep(NaN, 20)), colour = "incidence")) + geom_point() + 
  geom_vline(xintercept = length(data_convolved_weekend_effects), color = "red") + 
  geom_line(aes(y = data_country_roll_back_extension$extension_results$fitted, colour = 'observed incidence')) + 
  geom_ribbon(aes(ymin = results_roll_back$q2.5, ymax = results_roll_back$q97.5, fill = "Confidence interval roll back", colour = NA), alpha = 0.3, fill = colors[3]) +
  geom_line(aes(y = results_roll_back$mean, colour = 'Roll back mean')) + 
  scale_colour_manual("", 
                      breaks = c("incidence", 'observed incidence', 'Roll back mean', "Confidence interval roll back"),
                      values = c("black", colors[2], colors[3], colors[3])) +
  scale_fill_manual("", 
                    breaks = c("Confidence interval roll back"),
                    values = c(colors[3])) +
  labs(title="Roll back deconvolution on generated SIR data", x="Time", y="number of new cases") + 
  theme_for_the_plots + guides(color=guide_legend(override.aes=list(fill=NA)))

# Results:
# (+) Better for the end
# (-) here, we do not account for the uncertainty of the prediction. To do this, one can use bootstrapping techniques, such as 
# sampling multiple epidemic extended curves (using the confidence intervals/ parameters estimates with their covariance), then use roll back
# deconvolution for each curve. We don't do this here, as this will not solve the problem for the underestimation of large values/ or the 
# weekly pattern still here.

# - - - - - - - - - - - - - - - - - Sensitivity analysis ----------------------------------------------
# What if we change the delay distribution ?
shapes_values <- c(2, 5, 10, 20, 40)
shapes <- c("Mean: 1", "Mean: 2.5", "Mean: 5", "Mean: 10", "Mean: 20")
colors_comparison <- c("#d56062", "#f37748", "#ecc30b", "#84bcda", "#067bc2")

results_roll_back <- list()
for (j in c(1:length(shapes_values))){
  w_conv_misspecified <- discrete_gamma(n = 50, normalize = TRUE, shape = shapes_values[j], scale = 0.5)
  results_roll_back[[j]] <- rollBackCI(data_convolved, w_conv_misspecified, 100)
}

plot_roll_back_sensitivity <- ggplot(data = NULL, aes(x = 1:length(data_init), y = data_init, colour = "incidence")) + geom_point() + 
  geom_line(aes(y = data_convolved, colour = 'observed incidence')) + 
  
  geom_ribbon(aes(ymin = results_roll_back[[1]]$q2.5, ymax = results_roll_back[[1]]$q97.5, fill = shapes[1], colour = NA), alpha = 0.3, fill = colors_comparison[1]) +
  geom_line(aes(y = results_roll_back[[1]]$mean, colour = shapes[1])) +
  geom_ribbon(aes(ymin = results_roll_back[[2]]$q2.5, ymax = results_roll_back[[2]]$q97.5, fill = shapes[2], colour = NA), alpha = 0.3, fill = colors_comparison[2]) +
  geom_line(aes(y = results_roll_back[[2]]$mean, colour = shapes[2])) +
  geom_ribbon(aes(ymin = results_roll_back[[3]]$q2.5, ymax = results_roll_back[[3]]$q97.5, fill = shapes[3], colour = NA), alpha = 0.3, fill = colors_comparison[3]) +
  geom_line(aes(y = results_roll_back[[3]]$mean, colour = shapes[3])) +
  geom_ribbon(aes(ymin = results_roll_back[[4]]$q2.5, ymax = results_roll_back[[4]]$q97.5, fill = shapes[4], colour = NA), alpha = 0.3, fill = colors_comparison[4]) +
  geom_line(aes(y = results_roll_back[[4]]$mean, colour = shapes[4])) +
  geom_ribbon(aes(ymin = results_roll_back[[5]]$q2.5, ymax = results_roll_back[[5]]$q97.5, fill = shapes[5], colour = NA), alpha = 0.3, fill = colors_comparison[5]) +
  geom_line(aes(y = results_roll_back[[5]]$mean, colour = shapes[5])) +
  
  scale_colour_manual("", 
                      breaks = c("incidence", 'observed incidence', 'Roll back mean', "Confidence interval roll back", shapes),
                      values = c("black", colors[2], colors[3], colors[3], colors_comparison)) +
  scale_fill_manual("", 
                    breaks = c("Confidence interval roll back", shapes),
                    values = c(colors[3], colors_comparison)) +
  labs(title="Roll back deconvolution on generated SIR data, for several deconvolution shape parameters", x="Time", y="number of new cases") + 
  theme_for_the_plots + guides(color=guide_legend(override.aes=list(fill=NA)))

# Report: Figure 3c
plot(plot_roll_back_sensitivity)

# And with weekend effects ?
# What if we change the delay distribution ?
results_roll_back <- list()
for (j in c(1:length(shapes_values))){
  w_conv_misspecified <- discrete_gamma(n = 50, normalize = TRUE, shape = shapes_values[j], scale = 0.5)
  results_roll_back[[j]] <- rollBackCI(data_convolved_weekend_effects, w_conv_misspecified, 100)
}

plot_roll_back_sensitivity <- ggplot(data = NULL, aes(x = 1:length(data_init), y = data_init, colour = "incidence")) + geom_point() + 
  geom_line(aes(y = data_convolved_weekend_effects, colour = 'observed incidence')) + 
  
  geom_ribbon(aes(ymin = results_roll_back[[1]]$q2.5, ymax = results_roll_back[[1]]$q97.5, fill = shapes[1], colour = NA), alpha = 0.3, fill = colors_comparison[1]) +
  geom_line(aes(y = results_roll_back[[1]]$mean, colour = shapes[1])) +
  geom_ribbon(aes(ymin = results_roll_back[[2]]$q2.5, ymax = results_roll_back[[2]]$q97.5, fill = shapes[2], colour = NA), alpha = 0.3, fill = colors_comparison[2]) +
  geom_line(aes(y = results_roll_back[[2]]$mean, colour = shapes[2])) +
  geom_ribbon(aes(ymin = results_roll_back[[3]]$q2.5, ymax = results_roll_back[[3]]$q97.5, fill = shapes[3], colour = NA), alpha = 0.3, fill = colors_comparison[3]) +
  geom_line(aes(y = results_roll_back[[3]]$mean, colour = shapes[3])) +
  geom_ribbon(aes(ymin = results_roll_back[[4]]$q2.5, ymax = results_roll_back[[4]]$q97.5, fill = shapes[4], colour = NA), alpha = 0.3, fill = colors_comparison[4]) +
  geom_line(aes(y = results_roll_back[[4]]$mean, colour = shapes[4])) +
  geom_ribbon(aes(ymin = results_roll_back[[5]]$q2.5, ymax = results_roll_back[[5]]$q97.5, fill = shapes[5], colour = NA), alpha = 0.3, fill = colors_comparison[5]) +
  geom_line(aes(y = results_roll_back[[5]]$mean, colour = shapes[5])) +
  
  scale_colour_manual("", 
                      breaks = c("incidence", 'observed incidence', 'Roll back mean', "Confidence interval roll back", shapes),
                      values = c("black", colors[2], colors[3], colors[3], colors_comparison)) +
  scale_fill_manual("", 
                    breaks = c("Confidence interval roll back", shapes),
                    values = c(colors[3], colors_comparison)) +
  labs(title="Roll back deconvolution on generated SIR data, for several deconvolution shape parameters", x="Time", y="number of new cases") + 
  theme_for_the_plots + guides(color=guide_legend(override.aes=list(fill=NA)))

plot(plot_roll_back_sensitivity)

# Some curves are smoother than others, maybe due to weekends ?
curvature_score <- c()
shapes_values <- 2*c(1:20)
results_roll_back <- list()

for (j in c(1:length(shapes_values))){
  w_conv_misspecified <- discrete_gamma(n = 20, normalize = TRUE, shape = shapes_values[j], scale = 0.5)
  results_roll_back[[j]] <- rollBackCI(data_convolved, w_conv_misspecified, 10)
  
  center_curve <- results_roll_back[[j]]$mean[2:(length(results_roll_back[[j]]$mean) - 1)]
  
  derivs_second <- results_roll_back[[j]]$mean[3:(length(results_roll_back[[j]]$mean))] - 2*center_curve + results_roll_back[[j]]$mean[1:(length(results_roll_back[[j]]$mean) - 2)]
  ind1 <- c(2:length(derivs_second))
  ind2 <- c(1:(length(derivs_second)-1))
  # curvature_score <- c(curvature_score, sum(abs(derivs_second**2)))
  curvature_score <- c(curvature_score, sum((derivs_second[ind2]**2 + derivs_second[ind1]*derivs_second[ind2] + derivs_second[ind1]**2))/6)
}
plot(shapes_values/2, curvature_score)    # no particular pattern, although deconvolutions with mean 8 and 14 give the smoothest curves


# ------------------------- 2) Optimization with Poisson response -------------------------
# We move on to the next deconvolution method.

results_optim <- optimPoissonDeconvolution(data_convolved, w_conv, smooth_observations = FALSE, maxIter = 10, method = "CG", data_start = NULL)
# results_optim$result
# 
# hess <- log_lik_hessian(results$fitted_values, data_convolved, w_conv)
# print(diag(solve(hess[1:(dim(hess)[1]-1), 1:(dim(hess)[1]-1)])))

fitted_values_optim <- results_optim$fitted_values
fitted_values_optim <- c(0, fitted_values_optim[1:(length(fitted_values_optim)-1)])

ggplot(data = NULL, aes(x = 1:length(data_init), y = data_init, colour = "incidence")) + geom_point() + 
  geom_line(aes(y = data_convolved, colour = 'observed incidence')) + 
  # geom_ribbon(aes(ymin = results_roll_back$q2.5, ymax = results_roll_back$q97.5, fill = "Confidence interval", colour = "Confidence interval"), alpha = 0.3, fill = colors[3]) +
  geom_line(aes(y = fitted_values_optim, colour = 'Optimization')) + 
  scale_colour_manual("", 
                      breaks = c("incidence", 'observed incidence', 'Optimization', "Confidence interval"),
                      values = c("black", colors[2], colors[3], colors[3])) +
  scale_fill_manual("", 
                    breaks = c("Confidence interval"),
                    values = c(colors[3])) +
  labs(title="Deconvolution with optimization on generated SIR data", x="Time", y="number of new cases") + 
  theme_for_the_plots

# lines(rollapply(results$fitted_values, width=3, FUN=function(x) mean(x, na.rm=TRUE), by=1, 
#                 by.column=TRUE, partial=TRUE, fill=NA, align="center"), col = "green")

# Results:
# (+) Gives very good results, in term of fitted values especially: good approximation for high values for which we had problems with roll back
# (+) No negative values
# (-) Carefull not to apply to much smoothing before applying the optimization, otherwise the spike will not be estimated properly (see after)
# (-) at the beginning: some disturbations, and at the end, since the gradient is close to 0 (by the specification of the infectivity profile), it is 
# quite close to the initial values, which are in this case given by the observed data
# (-) Hessian not always invertible, problems to get the CI. Also, it depends on the optimization method (CG, BFGS, SANN, L-BFGS, etc)


results_optim_bfgs <- optimPoissonDeconvolution(data_convolved, w_conv, smooth_observations = FALSE, maxIter = 20, method = "BFGS", data_start = NULL)

fitted_values_optim_bfgs <- results_optim_bfgs$fitted_values
fitted_values_optim_bfgs <- c(0, fitted_values_optim_bfgs[1:(length(fitted_values_optim_bfgs)-1)])

ggplot(data = NULL, aes(x = 1:length(data_init), y = data_init, colour = "incidence")) + geom_point() + 
  geom_line(aes(y = data_convolved, colour = 'observed incidence')) + 
  # geom_ribbon(aes(ymin = results_roll_back$q2.5, ymax = results_roll_back$q97.5, fill = "Confidence interval", colour = "Confidence interval"), alpha = 0.3, fill = colors[3]) +
  geom_line(aes(y = fitted_values_optim_bfgs, colour = 'Optimization')) + 
  scale_colour_manual("", 
                      breaks = c("incidence", 'observed incidence', 'Optimization', "Confidence interval"),
                      values = c("black", colors[2], colors[3], colors[3])) +
  scale_fill_manual("", 
                    breaks = c("Confidence interval"),
                    values = c(colors[3])) +
  labs(title="Deconvolution with optimization on generated SIR data", x="Time", y="number of new cases") + 
  theme_for_the_plots


# Plot with both curves, conjugate gradient and BFGS:
# Report: Figure 4
ggplot(data = NULL, aes(x = 1:length(data_init), y = data_init, colour = "true incidence")) + geom_point() + 
  geom_line(aes(y = data_convolved, colour = 'observed incidence')) + 
  # geom_ribbon(aes(ymin = results_roll_back$q2.5, ymax = results_roll_back$q97.5, fill = "Confidence interval", colour = "Confidence interval"), alpha = 0.3, fill = colors[3]) +
  geom_line(aes(y = fitted_values_optim_bfgs, colour = 'BFGS')) + 
  geom_line(aes(y = fitted_values_optim, colour = 'Conjugate gradient')) + 
  scale_colour_manual("", 
                      breaks = c("true incidence", 'observed incidence', 'BFGS', "Confidence interval", 'Conjugate gradient'),
                      values = c("black", colors[2], colors[1], colors[3], colors[3])) +
  scale_fill_manual("", 
                    breaks = c("Confidence interval"),
                    values = c(colors[3])) +
  labs(title="Deconvolution with optimization on synthetic SIR data", x="Time", y="number of new cases") + 
  theme_for_the_plots

# - - - - - - - - - - - - - - - - - Noisy curve ----------------------------------------------
results_optim <- optimPoissonDeconvolution(data_convolved_weekend_effects, w_conv, 
                                           smooth_observations = FALSE, maxIter = 10, method = "CG", data_start = NULL)
fitted_values_optim <- results_optim$fitted_values
fitted_values_optim <- c(0, fitted_values_optim[1:(length(fitted_values_optim)-1)])

# Report: Figure 6b
ggplot(data = NULL, aes(x = 1:length(data_init), y = data_init, colour = "incidence")) + geom_point() + 
  geom_line(aes(y = data_convolved_weekend_effects, colour = 'observed incidence')) + 
  # geom_ribbon(aes(ymin = results_roll_back$q2.5, ymax = results_roll_back$q97.5, fill = "Confidence interval", colour = "Confidence interval"), alpha = 0.3, fill = colors[3]) +
  geom_line(aes(y = fitted_values_optim, colour = 'Optimization')) + 
  scale_colour_manual("", 
                      breaks = c("incidence", 'observed incidence', 'Optimization', "Confidence interval"),
                      values = c("black", colors[2], colors[3], colors[3])) +
  scale_fill_manual("", 
                    breaks = c("Confidence interval"),
                    values = c(colors[3])) +
  labs(title="Deconvolution with optimization on generated SIR data", x="Time", y="number of new cases") + 
  theme_for_the_plots

# Results:
# Fitted values are highly variable, with strong weekly pattern.
# We can check what gives the convolution with this deconvolved data:
re_convolved_results <- convolveData(w_conv, fitted_values_optim)
ggplot(data = NULL, aes(x = 1:length(data_init), y = data_init, colour = "incidence")) + geom_point() + 
  geom_line(aes(y = data_convolved_weekend_effects, colour = 'observed incidence')) + 
  # geom_ribbon(aes(ymin = results_roll_back$q2.5, ymax = results_roll_back$q97.5, fill = "Confidence interval", colour = "Confidence interval"), alpha = 0.3, fill = colors[3]) +
  geom_line(aes(y = re_convolved_results, colour = 'Optimization')) + 
  scale_colour_manual("", 
                      breaks = c("incidence", 'observed incidence', 'Optimization', "Confidence interval"),
                      values = c("black", colors[2], colors[3], colors[3])) +
  scale_fill_manual("", 
                    breaks = c("Confidence interval"),
                    values = c(colors[3])) +
  labs(title="Deconvolution with optimization on generated SIR data", x="Time", y="number of new cases") + 
  theme_for_the_plots

# The method catches the weekly trend (even though the dfit is not perfect). Also, the optimization can take some time, even with few steps. 
# This is mostly due to the number of parameters (in fact one for each day), and so adding a constraint migh be a good idea to solve this issue (see next deconvolution method).

# Also, note that the last values are estimated not as badly as for roll back (but this mainly depends on the initial values given to the optimization).

# One last check: with smoothed data (from noisy data):
results_optim <- optimPoissonDeconvolution(data_convolved_weekend_effects_smoothed, w_conv, 
                                           smooth_observations = FALSE, maxIter = 10, method = "CG", data_start = NULL)
fitted_values_optim <- results_optim$fitted_values
fitted_values_optim <- c(0, fitted_values_optim[1:(length(fitted_values_optim)-1)])

ggplot(data = NULL, aes(x = 1:length(data_init), y = data_init, colour = "incidence")) + geom_point() + 
  geom_line(aes(y = data_convolved_weekend_effects_smoothed, colour = 'observed incidence')) + 
  # geom_ribbon(aes(ymin = results_roll_back$q2.5, ymax = results_roll_back$q97.5, fill = "Confidence interval", colour = "Confidence interval"), alpha = 0.3, fill = colors[3]) +
  geom_line(aes(y = fitted_values_optim, colour = 'Optimization')) + 
  scale_colour_manual("", 
                      breaks = c("incidence", 'observed incidence', 'Optimization', "Confidence interval"),
                      values = c("black", colors[2], colors[3], colors[3])) +
  scale_fill_manual("", 
                    breaks = c("Confidence interval"),
                    values = c(colors[3])) +
  labs(title="Deconvolution with optimization on generated SIR data", x="Time", y="number of new cases") + 
  theme_for_the_plots

# The fit is not very good



# ------------------------- 3) IWLS with Poisson response and smoothing splines -------------------------

# Function for plot of the resulting deconvolved curve
#' @param results_fit list returned by the function 'DeconvolutionGAM'
#' @param confidence_intervals_results list returned by the function 'SimultaneousIntervals', or NULL if unknown
#' @param data_initial true incidence to be plotted
#' @param sim_fits TRUE to plot the simulated splines (from 'SimultaneousIntervals')
#' @param title title of the plot
plot_deconvolution_gam <- function(results_fit, confidence_intervals_results, data_initial, sim_fits = FALSE, title = "Deconvolution with IWLS algorithm and smoothing splines on generated SIR data"){
  pl <- ggplot(data = NULL, aes(x = 1:length(data_initial), y = data_initial, colour = "incidence")) + geom_point() + 
    geom_line(aes(y = results_fit$DATA_extension$cases, colour = 'observed incidence'))
  
  if (!is.null(confidence_intervals_results)){
    pl <- pl + geom_ribbon(aes(ymin = confidence_intervals_results$Point_wise_CI_lower, 
                               ymax = confidence_intervals_results$Point_wise_CI_upper, 
                               fill = "Confidence interval", colour = NA), alpha = 0.3, fill = colors[3]) +
      geom_ribbon(aes(ymin = confidence_intervals_results$Simultaneous_lower, 
                      ymax = confidence_intervals_results$Simultaneous_upper, 
                      fill = "Confidence interval", colour = NA), alpha = 0.3, fill = colors[3])
    if (sim_fits){
      pl <- pl + 
        geom_path(data = confidence_intervals_results$fits, mapping = aes(y = values, x = age, group = ind), alpha = 0.4, colour = "grey20")
    }
  }
  pl <- pl + 
    geom_line(aes(y = results_fit$fitted, colour = 'IWLS splines')) + 
    scale_colour_manual("", 
                        breaks = c("incidence", 'observed incidence', 'IWLS splines', "Confidence interval"),
                        values = c("black", colors[2], colors[3], colors[3])) +
    scale_fill_manual("", 
                      breaks = c("Confidence interval"),
                      values = c(colors[3])) +
    labs(title=title, x="Time", y="number of new cases") + 
    theme_for_the_plots + guides(color=guide_legend(override.aes=list(fill=NA))) + #+ scale_y_log10()
    coord_cartesian(ylim = c(0, 1.18 * max(data_initial, na.rm = T))) 
  
  return(pl)
}

# - - - - - - - - - - - - - - - - - Poisson identity link -----------------------------------------------
# Taking notes of the issues of the two previous methods, we now turn too an approach that does not include too many parameters to estimate,
# by making the deconvolved curve a smooth function of the time
# We start by setting the number of knots to 70, which implies that we have a new parameter every 3.5 days (so two parameters per week).
nb_knots <- 70

# First, we use a poisson response with identity link function, no extension, and the correct delay specification
results_IWLS <- DeconvolutionGAM(data_convolved, w_conv, nb_knots, plot_ = FALSE, get_derivatives = FALSE, family = poisson(link = 'identity'),
                                 extension = NULL)
# Compute the point-wise and simultaneous confidence bands
ConfidenceBands <- SimultaneousIntervals(results_IWLS$SplineBasis, plot_ = FALSE,
                                         coefficients_model = results_IWLS$fit$coefficients,
                                         CovMatrix = vcov(results_IWLS$fit),
                                         sample_splines = TRUE, ilink_fct = family(results_IWLS$fit)$linkinv)

# Plot the results
# Report: Figure 5a
plot_deconvolution_gam(results_IWLS, ConfidenceBands, data_init, FALSE)

# Diagnostic plots:
par(mfrow = c(2,2))
gam.check(results_IWLS$fit)

# Results:
# (+) Fast computations, get the CI, good fit
# (-) At the end, very large confidence intervals, and at the start, wavy behaving
# (+) one can control also the estimates by adding a custom penalization matrix


# - - - - - - - - - - - - - - - - - Poisson log link -----------------------------------------------
# One can also change the link function (so that values are always positive, in particular the confidence intervals)
results_IWLS <- DeconvolutionGAM(data_convolved, w_conv, nb_knots, plot_ = FALSE, get_derivatives = FALSE, family = poisson(link = 'log'),
                                 extension = NULL)
# Compute the point-wise and simultaneous confidence bands
ConfidenceBands <- SimultaneousIntervals(results_IWLS$SplineBasis, plot_ = FALSE,
                                         coefficients_model = results_IWLS$fit$coefficients,
                                         CovMatrix = vcov(results_IWLS$fit),
                                         sample_splines = TRUE, ilink_fct = family(results_IWLS$fit)$linkinv)

# Plot the results
# Report: Figure 5b
plot_deconvolution_gam(results_IWLS, ConfidenceBands, data_init, FALSE)


# - - - - - - - - - - - - - - - - - Poisson identity link with extension -----------------------------------------------
# With extensions:
results_IWLS_final_no_weekend_effect <- DeconvolutionGAM(data_convolved, w_conv, nb_knots, plot_ = FALSE, get_derivatives = FALSE, family = poisson(link = 'log'),
                                 extension = list('len' = 20, 'family' = quasipoisson(link = 'identity'), 'keep_original' = TRUE, 'min_weights' = 0.1))

ConfidenceBands <- SimultaneousIntervals(results_IWLS_final_no_weekend_effect$SplineBasis, plot_ = TRUE,
                                         coefficients_model = results_IWLS_final_no_weekend_effect$fit$coefficients,
                                         CovMatrix = vcov(results_IWLS_final_no_weekend_effect$fit),
                                         sample_splines = TRUE, ilink_fct = family(results_IWLS_final_no_weekend_effect$fit)$linkinv)

# Report: Figure 5c
plot_deconvolution_gam(results_IWLS_final_no_weekend_effect, ConfidenceBands, c(data_init, rep(NaN, 20)), FALSE)

# Results: 
# The deconvolution is much better, even we taken into account the uncertainty from the prediction. Note that for around 250, we still have the 
# confidence bands that are larger, integrating the uncertainty from the prediction


# We can add some penalization (spoil: here the penalization is too strong)
delta <- construct_delta(nb_knots + 4)
delta <- 0.01 * delta
results_IWLS <- DeconvolutionGAM(data_convolved, w_conv, nb_knots, plot_ = FALSE, delta = delta, get_derivatives = FALSE, family = poisson(link = 'identity'),
                                 extension = list('len' = 20, 'family' = quasipoisson(link = 'identity'), 'keep_original' = TRUE, 'min_weights' = 0.1))

plot_deconvolution_gam(results_IWLS, NULL, c(data_init, rep(NaN, 20)), FALSE)


# Extension with ones before the start of the epidemic:
data_convolved_extended_before <- c(rep(1, 20), data_convolved)
results_IWLS_final_no_weekend_effect <- DeconvolutionGAM(data_convolved_extended_before, w_conv, nb_knots+5, plot_ = FALSE, get_derivatives = FALSE, family = poisson(link = 'log'),
                                                         extension = list('len' = 20, 'family' = quasipoisson(link = 'identity'), 'keep_original' = TRUE, 'min_weights' = 0.1))
ConfidenceBands <- SimultaneousIntervals(results_IWLS_final_no_weekend_effect$SplineBasis, plot_ = TRUE,
                                         coefficients_model = results_IWLS_final_no_weekend_effect$fit$coefficients,
                                         CovMatrix = vcov(results_IWLS_final_no_weekend_effect$fit),
                                         sample_splines = TRUE, ilink_fct = family(results_IWLS_final_no_weekend_effect$fit)$linkinv)

plot_deconvolution_gam(results_IWLS_final_no_weekend_effect, ConfidenceBands, c(rep(NaN, 20), data_init, rep(NaN, 20)), FALSE)
# more stable at the beginning


# - - - - - - - - - - - - - - - - - Problem: choice for the number of knots -----------------------------------------------
# The choice for the number of knots is not trivial
deg_splines <- c(10:200)
matrix_data <- matrix(nrow = length(deg_splines), ncol = length(data_convolved))
for (j in 1:length(deg_splines)){
  print(j)
  results_IWLS <- DeconvolutionGAM(data_convolved, w_conv, deg_splines[j], plot_ = FALSE, delta = NULL, get_derivatives = FALSE, 
                                   family = poisson(link = 'identity'),
                                   extension = NULL)
  matrix_data[j, ] <- results_IWLS$fitted_values
}
matrix_data[which(matrix_data < 0)] <- 0
matrix_data[which(matrix_data > 10000)] <- 10000

# matrix_data_flatten <- as.vector(t(matrix_data))
# time_flatten <- rep(c(1:length(data_convolved)), length(deg_splines))
# cat_flatten <- rep(deg_splines, each = length(data_convolved))
# dd <- cbind(matrix_data_flatten, time_flatten, cat_flatten)
# dd <- data.frame(dd)
# colnames(dd) <- c("deconvolved", "time", "degree")
# plot_ly(data.frame(dd), x = ~time, y = ~degree, z = ~deconvolved, split = ~degree, type = "scatter3d", mode = "lines", zaxis = list(range = c(0, 9000)))

install.packages("plotly")
library(plotly)

# Report: Figure 18 Appendix A
plot_ly(z = ~matrix_data) %>% add_surface() %>% layout(
  title = "Estimated true incidence curve for different basis dimension",
  scene = list(
    xaxis = list(title = "time"),
    yaxis = list(title = "basis dimension"),
    zaxis = list(title = "incidence cases")
  ))


# - - - - - - - - - - - - - - - - - Noisy curve ----------------------------------------------
# We start with the same setting as before, with log link function, and extension
nb_knots <- 70
results_IWLS <- DeconvolutionGAM(data_convolved_weekend_effects, w_conv, nb_knots, plot_ = FALSE, 
                                 get_derivatives = FALSE, family = poisson(link = 'log'),
                                 extension = list('len' = 20, 'family' = quasipoisson(link = 'identity'), 'keep_original' = TRUE, 'min_weights' = 0.1))
# Compute the point-wise and simultaneous confidence bands
ConfidenceBands <- SimultaneousIntervals(results_IWLS$SplineBasis, plot_ = TRUE,
                                         coefficients_model = results_IWLS$fit$coefficients,
                                         CovMatrix = vcov(results_IWLS$fit),
                                         sample_splines = TRUE, ilink_fct = family(results_IWLS$fit)$linkinv)
plot_deconvolution_gam(results_IWLS, ConfidenceBands, c(data_init, rep(NaN, 20)), FALSE)

# Results:
# Well... the fit is not very good. The splines begin to fit the weekly pattern. 

# The first thing we can try is to add penalization term
delta <- construct_delta(nb_knots + 4)
delta <- 1000 * delta    # 1000 because we are using log link
results_IWLS <- DeconvolutionGAM(data_convolved_weekend_effects, w_conv, nb_knots, plot_ = FALSE, delta = delta,
                                 get_derivatives = FALSE, family = poisson(link = 'log'),
                                 extension = list('len' = 20, 'family' = nb(link = 'log'), 'keep_original' = TRUE, 'min_weights' = 0.1))
# Compute the point-wise and simultaneous confidence bands
ConfidenceBands <- SimultaneousIntervals(results_IWLS$SplineBasis, plot_ = TRUE,
                                         coefficients_model = results_IWLS$fit$coefficients,
                                         CovMatrix = vcov(results_IWLS$fit),
                                         sample_splines = TRUE, ilink_fct = family(results_IWLS$fit)$linkinv)
plot_deconvolution_gam(results_IWLS, ConfidenceBands, c(data_init, rep(NaN, 20)), FALSE)

# Results: 
# Much better fit, but we chose the penalization constant arbitrarily: too high and the curve is too smooth, too low and the weekly pattern remains.
# One must search for an automatic way of dealing with this choice


# - - - - - - - - - - - - - - - - - Automatic selection + Noisy curve ----------------------------------------------
# First, we check if the extension is OK:
extension_sir <- extension_prediction(data_convolved_weekend_effects,
                                      extension = list('len' = 40, 'family' = nb(link = 'log'),
                                                       'keep_original' = FALSE, 'min_weights' = 0.1),
                                      fm = "cases ~ s(time, k = 20, bs = \'tp\', m = c(0,1)) + s(day, k = 3)",
                                      weights = NULL)
par(mfrow = c(1,1))
plot(extension_sir$extension_results$fitted, type = 'l', xlab = "time", ylab = "Incidence curve", main = "Observed incidence curve, with extension") 
abline(v = length(data_convolved_weekend_effects), col = "red")


# Now we search for the best basis dimension
dimensionBasis <- findBestDimension(data_convolved_weekend_effects, w_conv,
                                    family = quasipoisson(link = 'log'), smoothing = 7, 
                                    extension = list('len' = 40, 'family' = nb(link = 'log'), 'keep_original' = TRUE, 'min_weights' = 0.1),
                                    range_values = c(30:100), diagnostic_plot = TRUE,
                                    fm = "cases ~ s(time, k = 20, bs = \'tp\', m = c(0, 1)) + s(day, k = 5)")
# Then, based on the basis dimension found, we search for the best regularization
lambda_delta <- findBestLambda(data_convolved_weekend_effects, dimensionBasis, w_conv,
                               range_values = 10**(20:-10), family = quasipoisson(link = 'log'), smoothing = 7,
                               extension = list('len' = 40, 'family' = nb(link = 'log'), 'keep_original' = TRUE, 'min_weights' = 0.1), 
                               fm = "cases ~ s(time, k = 20, bs = \'tp\', m = c(0, 1)) + s(day, k = 5)", diagnostic_plot = TRUE)

delta <- construct_delta(dimensionBasis + 4)
delta <- lambda_delta * delta

results_IWLS <- DeconvolutionGAM(data_convolved_weekend_effects, w_conv, dimensionBasis, plot_ = TRUE, delta = delta,
                                 get_derivatives = FALSE, family = quasipoisson(link = 'log'),
                                 extension = list('len' = 40, 'family' = nb(link = 'log'), 'keep_original' = TRUE, 'min_weights' = 0.1),
                                 fm = "cases ~ s(time, k = 20, bs = \'tp\', m = c(0, 1)) + s(day, k = 5)")
ConfidenceBands <- SimultaneousIntervals(results_IWLS$SplineBasis, plot_ = TRUE,
                                         coefficients_model = results_IWLS$fit$coefficients,
                                         CovMatrix = vcov(results_IWLS$fit),
                                         sample_splines = TRUE, ilink_fct = family(results_IWLS$fit)$linkinv)
# Report: Figure 7a
plot_deconvolution_gam(results_IWLS, ConfidenceBands, c(data_init, rep(NaN, 40)), FALSE)

# Results: 
# The fit is not bad: the confidence intervals are larger than when we didn't have weekend effects. The large peak is still a bit underestimated, 
# but are still between the confidence bands. For the last values, the confidence bands are larger, as expected. Also, there is no weekly pattern
# in the deconvolve data, which is a good news.

# One can try also without regularization term:
results_IWLS_final_with_weekend_effect <- DeconvolutionGAM(data_convolved_weekend_effects, w_conv, dimensionBasis, plot_ = FALSE, delta = NULL,
                                 get_derivatives = FALSE, family = quasipoisson(link = 'log'),
                                 extension = list('len' = 40, 'family' = nb(link = 'log'), 'keep_original' = TRUE, 'min_weights' = 0.1),
                                 fm = "cases ~ s(time, k = 20, bs = \'tp\', m = c(0, 1)) + s(day, k = 5)")
ConfidenceBands <- SimultaneousIntervals(results_IWLS_final_with_weekend_effect$SplineBasis, plot_ = TRUE,
                                         coefficients_model = results_IWLS_final_with_weekend_effect$fit$coefficients,
                                         CovMatrix = vcov(results_IWLS_final_with_weekend_effect$fit),
                                         sample_splines = TRUE, ilink_fct = family(results_IWLS_final_with_weekend_effect$fit)$linkinv)
# Report: Figure 7b
plot_deconvolution_gam(results_IWLS_final_with_weekend_effect, ConfidenceBands, c(data_init, rep(NaN, 40)), FALSE)

# Results:
# Confidence bands are (much) larger, especially for the first and last values. This time, all the transition periods are correctly estimated
# including the peak between 60 and 70. However, the last values begin to exhibit wavy pattern, as well as the confidence bands.

# ===================== II) Estimating R =====================
# In this part, we will be estimating the reproductive number, based on the true incidence curve, but also on the estimated one.
# For this, we use two approaches: one based on a generalized additive model and the generative distribution, and one based on a SIR model
# with smoothed transmission and removal rates.

# _-**-_-**-_-**-_-**-_-**-_ 1) GAM -------------------------
source("/Users/Antoine/Desktop/Projet Semestre 2020/CODE/EstimatingRt.R")

true_reproductive_number <- df$true_rt[2:(length(df$true_rt)-1)][indices_original_synthetic]
true_r0 <- df$true_r0[2:(length(df$true_rt)-1)][indices_original_synthetic]



# RECALL of the true (basic) reproductive numbers: 
ggplot(data = NULL, aes(x = 1:length(true_reproductive_number))) + #geom_line(size = 1) +
  annotate(geom = "rect", xmin = 60, xmax = 67, ymin = -Inf, ymax = Inf, fill = "#2a9d8f", alpha = 0.2)+
  annotate(geom = "rect", xmin = 90, xmax = 97, ymin = -Inf, ymax = Inf, fill = "#e63946", alpha = 0.2)+
  geom_segment(aes(x = 60, y = 2.1, xend = 90, yend = 2.1), color = "#e63946", arrow = arrow(length = unit(0.2, "cm"))) + 
  geom_segment(aes(x = 90, y = 2.1, xend = 60, yend = 2.1), color = "#e63946", arrow = arrow(length = unit(0.2, "cm"))) + 
  geom_text(data=NULL, aes(x=75, y= 2.1, label = "30 days"), colour = "#e63946", vjust = -0.5, size = 3) + 
  geom_hline(yintercept = 1, color = "red") + 
  geom_line(aes(y = true_r0, colour = "R0"), size = 1) + 
  geom_line(aes(y = true_reproductive_number, colour = "R_t"), size = 1) + 
  scale_colour_manual("", 
                      breaks = c("incidence", "S", "E", "I", "R", 'observed incidence', "R_t", "R0"),
                      values = c("black", "gray", colors[1], colors[2], colors[3], colors[2], "#1d3557", "#e63946")) +
  labs(title="Reproductive number for generated SIR data", x="time", y="Reproductive number") + 
  theme_for_the_plots

# Generation interval, based on sum of two poisson distribution, with mean 4
W_gam <- getGenerationInterval(7, 32, 50)

# - - - - - - - - - - - - - - - - - First, without deconvolution (using true events) - - - - - - - - - - - - - - - - -----------
indices <- c(1:length(data_convolved))

# First, we start be specifying the offset values
offsetValues_data_init <- convolve(c(0* 1:(length(c(0,W_gam))-1), data_init), rev(c(0,W_gam)), type = "filter")
offsetValues_data_init[which(offsetValues_data_init <= 1)] <- 1


which_model <- poisson(link = "log")
DAT_data_init <- data.frame(resp = data_init, offs = offsetValues_data_init, x = c(1:length(data_init)))
fit_data_init <- fitModel(DAT_data_init, weights_ = NULL, which_model = which_model,  fm = 'resp ~ s(x, k = 70)', ignore_first_max = 5, 
                          deriv_plots = c(0,1), detect_flat = TRUE, x_axis = indices)

plot_r <- ggplot(data = NULL, aes(x = c(1:length(fit_data_init$fitted.values)))) +
  annotate(geom = "rect", xmin = 60, xmax = 67, ymin = -Inf, ymax = Inf, fill = "#2a9d8f", alpha = 0.2)+
  annotate(geom = "rect", xmin = 90, xmax = 97, ymin = -Inf, ymax = Inf, fill = "#e63946", alpha = 0.2)+
  geom_ribbon(aes(ymin = fit_data_init$lower, ymax = fit_data_init$upper, fill = "95% CI", colour = "white"), alpha = 0.3, fill = colors[3]) +
  geom_line(aes(y = true_reproductive_number, colour = "True reproductive number"), size = 1) + 
  geom_line(aes(y = data_init/offsetValues_data_init, colour = "Cori et al. estimate"), size = 1) + 
  geom_line(aes(y = fit_data_init$fitted.values, colour = "Estimated reproductive number"), size = 1) + 
  labs(title="Estimated reproductive number, based on the true infection events", x="Time", y=latex2exp::TeX("$R_t$")) +
  scale_colour_manual("", 
                      breaks = c("True reproductive number", "Estimated reproductive number", "Cori et al. estimate", "95% CI"),
                      values = c("#1d3557", colors[3], colors[1], colors[3])) +
  scale_fill_manual("", breaks = c("95% CI"), values = colors[3]) + 
  geom_hline(yintercept = 1, col = "red") +
  coord_cartesian(ylim = c(0, 3)) + 
  theme_for_the_plots + 
  guides(color=guide_legend(override.aes=list(fill=NA)))

# Results: 
# The estimated values seem to correctly fit the true values. At the start, the method is heavily overestimating the reproductive number (like the Cori estimate).

# Here is a plot showing the derivatives for the fitted values, with their confidence intervals. The 'critical' periods, where R_t is probably equal to 1,
# or the first derivative is equal to 0 (R_t constant, or very slightly changing).
plot(fit_data_init$derivPlots$plot)

# The method correctly catches the two transitions.


plot_r <- plot_r +
  geom_vline(xintercept = 60, color = 'black', linetype = "dashed", alpha = 0.8)+
  geom_vline(xintercept = 67, color = 'black', linetype = "dashed", alpha = 0.8)+
  geom_vline(xintercept = 90, color = 'black', linetype = "dashed", alpha = 0.8)+
  geom_vline(xintercept = 97, color = 'black', linetype = "dashed", alpha = 0.8)

plot_estim_r_deriv <- fit_data_init$derivPlots$list_plots[[2]] + 
  coord_cartesian(ylim = c(-0.3, 0.2)) +
  geom_vline(xintercept = 60, color = 'black', linetype = "dashed", alpha = 0.8)+
  geom_vline(xintercept = 67, color = 'black', linetype = "dashed", alpha = 0.8)+
  geom_vline(xintercept = 90, color = 'black', linetype = "dashed", alpha = 0.8)+
  geom_vline(xintercept = 97, color = 'black', linetype = "dashed", alpha = 0.8)

# Report: Figure 8
plot_grid(plot_r + theme(plot.margin = unit(c(0.5, 1, 0, 1), "cm")), 
          plot_estim_r_deriv + theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm")), 
          nrow = 2, align = 'hv')



#                                             Sensitivity analysis -----------------
W_gam <- getGenerationInterval(7, 32, 50)

w_mean <- c(4, 6, 7, 9, 11)
w_variance <- c(32, 32, 32, 32, 32)
parameter_changed <- c("Mean: 5", "Mean: 7", "Mean: 8", "Mean: 10", "Mean: 12")
# w_mean <- c(7, 7, 7, 7, 7)
# w_variance <- c(8, 16, 32, 64, 128)
# parameter_changed <- c("Variance: 8", "Variance: 16", "Variance: 32", "Variance: 64", "Variance: 128")
colors_comparison <- c("#ffbe0b", "#fb5607", "#ff006e", "#8338ec", "#3a86ff")
# colors_comparison <- c("#00AFBB", "#E7B800", "#FC4E07", "red", "blue")

plot_comparison <- ggplot(data = NULL, aes(x = c(1:(length(true_reproductive_number))))) + 
  annotate(geom = "rect", xmin = 60, xmax = 67, ymin = -Inf, ymax = Inf, fill = "#2a9d8f", alpha = 0.2)+
  annotate(geom = "rect", xmin = 90, xmax = 97, ymin = -Inf, ymax = Inf, fill = "#e63946", alpha = 0.2)+
  geom_hline(yintercept = 1, col = "red") +
  geom_line(aes(y = true_reproductive_number, colour = "True reproductive number"), size = 1) 
l_res <- list()

for (j in c(1:length(parameter_changed))){
  print(parameter_changed[j])
  
  W_gam_c <- getGenerationInterval(w_mean[j], w_variance[j], 50)
  offsetValues_data_init <- convolve(c(0* 1:(length(c(0,W_gam_c))-1), data_init), rev(c(0,W_gam_c)), type = "filter")
  
  offsetValues_data_init[which(offsetValues_data_init <= 1)] <- 1
  which_model <- poisson
  DAT_data_init <- data.frame(resp = data_init, offs = offsetValues_data_init, x = c(1:length(data_init)))
  fit_data_init <- fitModel(DAT_data_init, weights_ = NULL, which_model = poisson,  fm = 'resp ~ s(x, k = 70)', ignore_first_max = 5, derivatives = FALSE)
  
  l_res[[j]] <- fit_data_init
}

plot_comparison <- plot_comparison + 
  geom_line(aes(x = c(1:length(l_res[[1]]$fitted.values)), y = l_res[[1]]$fitted.values, color = parameter_changed[1])) + 
  geom_ribbon(aes(x = c(1:length(l_res[[1]]$upper)), 
                  ymin = l_res[[1]]$lower, 
                  ymax = l_res[[1]]$upper), alpha = 0.2, fill = colors_comparison[1]) + 
  geom_line(aes(x = c(1:length(l_res[[2]]$fitted.values)), y = l_res[[2]]$fitted.values, color = parameter_changed[2])) + 
  geom_ribbon(aes(x = c(1:length(l_res[[2]]$upper)), 
                  ymin = l_res[[2]]$lower, 
                  ymax = l_res[[2]]$upper), alpha = 0.2, fill = colors_comparison[2]) + 
  geom_line(aes(x = c(1:length(l_res[[3]]$fitted.values)), y = l_res[[3]]$fitted.values, color = parameter_changed[3])) + 
  geom_ribbon(aes(x = c(1:length(l_res[[3]]$upper)), 
                  ymin = l_res[[3]]$lower, 
                  ymax = l_res[[3]]$upper), alpha = 0.2, fill = colors_comparison[3]) + 
  geom_line(aes(x = c(1:length(l_res[[4]]$fitted.values)), y = l_res[[4]]$fitted.values, color = parameter_changed[4])) + 
  geom_ribbon(aes(x = c(1:length(l_res[[4]]$upper)), 
                  ymin = l_res[[4]]$lower, 
                  ymax = l_res[[4]]$upper), alpha = 0.2, fill = colors_comparison[4]) + 
  geom_line(aes(x = c(1:length(l_res[[5]]$fitted.values)), y = l_res[[5]]$fitted.values, color = parameter_changed[5])) + 
  geom_ribbon(aes(x = c(1:length(l_res[[5]]$upper)), 
                  ymin = l_res[[5]]$lower, 
                  ymax = l_res[[5]]$upper), alpha = 0.2, fill = colors_comparison[5])

plot_comparison <- plot_comparison + 
  # geom_line(aes(x = c(1:length(true_reproductive_number)), y = true_reproductive_number)) + 
  scale_colour_manual("", 
                      breaks = c(parameter_changed, "True reproductive number"),
                      values = c(colors_comparison, "#1d3557") ) +
  # scale_fill_manual("Mean", 
  #                   breaks = c(parameter_changed, "True reproductive number"),
  #                   values = c(colors_comparison, "#1d3557") ) + 
  theme_for_the_plots + 
  labs(title="Estimation of the reproductive number for different mean of the generation interval", x="Time", y=latex2exp::TeX("$R_t$")) +
  coord_cartesian(ylim = c(0, 4))
  

# Report: Figure 11a and 11b
plot(plot_comparison)


# - - - - - - - - - - - - - - - - - Roll back - - - - - - - - - - - - - - - ----
offsetValues_roll_back <- convolve(c(0* 1:(length(c(0,W_gam))-1), results_roll_back_ext$mean[indices]), rev(c(0,W_gam)), type = "filter")
DAT_roll_back <- data.frame(resp = results_roll_back_ext$mean[indices], offs = offsetValues_roll_back, x = c(1:length(results_roll_back_ext$mean[indices])))
DAT_roll_back$resp[which(DAT_roll_back$resp < 1)] <- 1
DAT_roll_back$offs[which(DAT_roll_back$offs < 1)] <- 1

weights = abs(results_roll_back_ext$q2.5[indices] - results_roll_back_ext$q97.5[indices])/4
weights[which(weights == 0)] <- 1
weights <- 1/weights**2
weights[(length(weights)-5):length(weights)] <- 0
weights <- weights/mean(weights)

which_model <- poisson(link = 'log')
fit_roll_back <- fitModel(DAT_roll_back, which_model = which_model, weights_ = weights, fm = 'resp ~ s(x, k = 70)', ignore_first_max = 5,
                          deriv_plots = c(0,1), detect_flat = TRUE)

ggplot(data = NULL, aes(x = c(1:length(fit_roll_back$fitted.values)))) +
  annotate(geom = "rect", xmin = 60, xmax = 67, ymin = -Inf, ymax = Inf, fill = "#2a9d8f", alpha = 0.2)+
  annotate(geom = "rect", xmin = 90, xmax = 97, ymin = -Inf, ymax = Inf, fill = "#e63946", alpha = 0.2)+
  geom_ribbon(aes(ymin = fit_roll_back$lower, ymax = fit_roll_back$upper, fill = "95% CI", colour = "white"), alpha = 0.3, fill = colors[3]) +
  geom_line(aes(y = true_reproductive_number, colour = "True reproductive number"), size = 1) + 
  geom_line(aes(y = results_roll_back_ext$mean[indices]/offsetValues_roll_back, colour = "Cori et al. estimate"), size = 1) + 
  geom_line(aes(y = fit_roll_back$fitted.values, colour = "Estimated reproductive number"), size = 1) + 
  labs(title="Reproductive number", x="Time", y=latex2exp::TeX("$R_t$")) +
  scale_colour_manual("", 
                      breaks = c("True reproductive number", "Estimated reproductive number", "Cori et al. estimate", "95% CI"),
                      values = c("#1d3557", colors[3], colors[1], colors[3])) +
  scale_fill_manual("", breaks = c("95% CI"), values = colors[3]) + 
  geom_hline(yintercept = 1, col = "red") +
  coord_cartesian(ylim = c(0, 3)) + 
  theme_for_the_plots + 
  guides(color=guide_legend(override.aes=list(fill=NA)))

# Results: 
# Because we are doing a convolution, but in the opposite time direction, the transitions are not correctly captured. 
# Also, because we have a drop in the estimated incidence cases for the last values, the R_t also drops, which is really problematic for
# a monitoring the evolution of R_t on a real-time basis.
plot(fit_roll_back$derivPlots$plot)


# - - - - - - - - - - - - - - - - - Optimization with Poisson response - - - - - ----
offsetValues_optim <- convolve(c(0* 1:(length(c(0,W_gam))-1), fitted_values_optim), rev(c(0,W_gam)), type = "filter")
offsetValues_optim[which(offsetValues_optim <= 1)] <- 1

DAT_optim <- data.frame(resp = fitted_values_optim, offs = offsetValues_optim, x = c(1:length(fitted_values_optim)))
which_model <- poisson(link = 'log')
fit_optim <- fitModel(DAT_optim, which_model = which_model, weights_ = NULL, fm = 'resp ~ s(x, k = 70)', ignore_first_max = 5,
                          deriv_plots = c(0,1), detect_flat = TRUE)

ggplot(data = NULL, aes(x = c(1:length(fit_optim$fitted.values)))) +
  annotate(geom = "rect", xmin = 60, xmax = 67, ymin = -Inf, ymax = Inf, fill = "#2a9d8f", alpha = 0.2)+
  annotate(geom = "rect", xmin = 90, xmax = 97, ymin = -Inf, ymax = Inf, fill = "#e63946", alpha = 0.2)+
  geom_ribbon(aes(ymin = fit_optim$lower, ymax = fit_optim$upper, fill = "95% CI", colour = "white"), alpha = 0.3, fill = colors[3]) +
  geom_line(aes(y = true_reproductive_number, colour = "True reproductive number"), size = 1) + 
  geom_line(aes(y = fitted_values_optim/offsetValues_optim, colour = "Cori et al. estimate"), size = 1) + 
  geom_line(aes(y = fit_optim$fitted.values, colour = "Estimated reproductive number"), size = 1) + 
  labs(title="Reproductive number", x="Time", y=latex2exp::TeX("$R_t$")) +
  scale_colour_manual("", 
                      breaks = c("True reproductive number", "Estimated reproductive number", "Cori et al. estimate", "95% CI"),
                      values = c("#1d3557", colors[3], colors[1], colors[3])) +
  scale_fill_manual("", breaks = c("95% CI"), values = colors[3]) + 
  geom_hline(yintercept = 1, col = "red") +
  coord_cartesian(ylim = c(0, 3)) + 
  theme_for_the_plots + 
  guides(color=guide_legend(override.aes=list(fill=NA)))

# Results:
# Not very good, mainly due to the bad quality of the deconvolution based on the optimization


# - - - - - - - - - - - - - - - - - IWLS with Poisson response and smoothing splines- - ----
#                                             Without weekend effects ----------------------
offsetValues_IWLS <- convolve(c(0* 1:(length(c(0,W_gam))-1), results_IWLS_final_no_weekend_effect$fitted[indices]), rev(c(0,W_gam)), type = "filter")
offsetValues_IWLS[which(offsetValues_IWLS <= 0)] <- 1

weights = abs(results_IWLS_final_no_weekend_effect$point_wise_lower[indices] - results_IWLS_final_no_weekend_effect$point_wise_upper[indices])/4
weights <- weights**2
# weights <- weights/results_IWLS_final_no_weekend_effect$fitted_values[indices]
weights <- 1/weights
weights <- capWeights(weights, 0.2)
weights <- weights/mean(weights)
par(mfrow = c(1,1))
plot(weights)

which_model <- nb(link = "log")
DAT_IWLS <- data.frame(resp = results_IWLS_final_no_weekend_effect$fitted_values[indices], offs = offsetValues_IWLS, x = c(1:length(results_IWLS_final_no_weekend_effect$fitted_values[indices])))
fit_IWLS <- fitModel(DAT_IWLS, which_model = which_model, weights_ = weights, fm = 'resp ~ s(x, k = 70)', ignore_first_max = 5,
                      deriv_plots = c(0,1), detect_flat = TRUE, x_axis = indices)
summary(fit_IWLS$fit)

plot_r <- ggplot(data = NULL, aes(x = c(1:length(fit_IWLS$fitted.values)))) +
  annotate(geom = "rect", xmin = 60, xmax = 67, ymin = -Inf, ymax = Inf, fill = "#2a9d8f", alpha = 0.2)+
  annotate(geom = "rect", xmin = 90, xmax = 97, ymin = -Inf, ymax = Inf, fill = "#e63946", alpha = 0.2)+
  geom_ribbon(aes(ymin = fit_IWLS$lower, ymax = fit_IWLS$upper, fill = "95% CI", colour = "white"), alpha = 0.3, fill = colors[3]) +
  geom_line(aes(y = true_reproductive_number, colour = "True reproductive number"), size = 1) + 
  geom_line(aes(y = results_IWLS_final_no_weekend_effect$fitted[indices]/offsetValues_IWLS, colour = "Cori et al. estimate"), size = 1) + 
  geom_line(aes(y = fit_IWLS$fitted.values, colour = "Estimated reproductive number"), size = 1) + 
  labs(title="Estimated reproductive number from observed incidence events without weekly pattern", x="Time", y=latex2exp::TeX("$R_t$")) +
  scale_colour_manual("", 
                      breaks = c("True reproductive number", "Estimated reproductive number", "Cori et al. estimate", "95% CI"),
                      values = c("#1d3557", colors[3], colors[1], colors[3])) +
  scale_fill_manual("", breaks = c("95% CI"), values = colors[3]) + 
  geom_hline(yintercept = 1, col = "red") +
  coord_cartesian(ylim = c(0, 3)) + 
  theme_for_the_plots + 
  guides(color=guide_legend(override.aes=list(fill=NA)))

# Results: 
# We can see clearly the confidence bands than become larger in the end, when the uncertainties for the last values are larger.
# The estimates are good, and, unlike the cori estimates that depends heavily on the deconvolution, here we add an extra smoothing step,
# which, for the first values, becomes quite useful, preventing from large oscillations.



plot_r <- plot_r +
  geom_vline(xintercept = 60, color = 'black', linetype = "dashed", alpha = 0.8)+
  geom_vline(xintercept = 67, color = 'black', linetype = "dashed", alpha = 0.8)+
  geom_vline(xintercept = 90, color = 'black', linetype = "dashed", alpha = 0.8)+
  geom_vline(xintercept = 97, color = 'black', linetype = "dashed", alpha = 0.8)

plot_estim_r_deriv <- fit_IWLS$derivPlots$list_plots[[2]] + 
  coord_cartesian(ylim = c(-0.3, 0.2)) +
  geom_vline(xintercept = 60, color = 'black', linetype = "dashed", alpha = 0.8)+
  geom_vline(xintercept = 67, color = 'black', linetype = "dashed", alpha = 0.8)+
  geom_vline(xintercept = 90, color = 'black', linetype = "dashed", alpha = 0.8)+
  geom_vline(xintercept = 97, color = 'black', linetype = "dashed", alpha = 0.8)

plot_grid(plot_r + theme(plot.margin = unit(c(0.5, 1, 0, 1), "cm")), 
          plot_estim_r_deriv + theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm")), 
          nrow = 2, align = 'hv')


#                                             With weekend effects ----------------------
offsetValues_IWLS <- convolve(c(0* 1:(length(c(0,W_gam))-1), results_IWLS_final_with_weekend_effect$fitted_values[indices]), rev(c(0,W_gam)), type = "filter")
offsetValues_IWLS[which(offsetValues_IWLS <= 0)] <- 1


weights = abs(results_IWLS_final_with_weekend_effect$point_wise_lower[indices] - results_IWLS_final_with_weekend_effect$point_wise_upper[indices])/4
weights <- weights**2
weights <- weights/sqrt(results_IWLS_final_with_weekend_effect$fitted_values[indices])
weights <- 1/weights
weights <- capWeights(weights, 0.01)
weights <- weights/mean(weights)
par(mfrow = c(1,1))
plot(weights)

which_model <- nb(link = "log")
DAT_IWLS <- data.frame(resp = results_IWLS_final_with_weekend_effect$fitted_values[indices], offs = offsetValues_IWLS, x = c(1:length(results_IWLS_final_with_weekend_effect$fitted_values[indices])))
fit_IWLS_no_delta <- fitModel(DAT_IWLS, which_model = which_model, weights_ = weights, fm = 'resp ~ s(x, k = 70)', ignore_first_max = 5,
                     deriv_plots = c(0,1), detect_flat = TRUE, x_axis = indices)

plot_r_weekend <- ggplot(data = NULL, aes(x = c(1:length(fit_IWLS_no_delta$fitted.values)))) +
  annotate(geom = "rect", xmin = 60, xmax = 67, ymin = -Inf, ymax = Inf, fill = "#2a9d8f", alpha = 0.2)+
  annotate(geom = "rect", xmin = 90, xmax = 97, ymin = -Inf, ymax = Inf, fill = "#e63946", alpha = 0.2)+
  geom_ribbon(aes(ymin = fit_IWLS_no_delta$lower, ymax = fit_IWLS_no_delta$upper, fill = "95% CI", colour = "white"), alpha = 0.3, fill = colors[3]) +
  geom_line(aes(y = true_reproductive_number, colour = "True reproductive number"), size = 1) + 
  geom_line(aes(y = results_IWLS_final_with_weekend_effect$fitted[indices]/offsetValues_IWLS, colour = "Cori et al. estimate"), size = 1) + 
  geom_line(aes(y = fit_IWLS_no_delta$fitted.values, colour = "Estimated reproductive number"), size = 1) + 
  labs(title="Reproductive number", x="Time", y=latex2exp::TeX("$R_t$")) +
  scale_colour_manual("", 
                      breaks = c("True reproductive number", "Estimated reproductive number", "Cori et al. estimate", "95% CI"),
                      values = c("#1d3557", colors[3], colors[1], colors[3])) +
  scale_fill_manual("", breaks = c("95% CI"), values = colors[3]) + 
  geom_hline(yintercept = 1, col = "red") +
  coord_cartesian(ylim = c(0, 3)) + 
  theme_for_the_plots + 
  guides(color=guide_legend(override.aes=list(fill=NA)))
plot_r_weekend
# plot(fit_IWLS_no_delta$derivPlots$list_plots[[2]])



# With delta:
offsetValues_IWLS <- convolve(c(0* 1:(length(c(0,W_gam))-1), results_IWLS$fitted_values[indices]), rev(c(0,W_gam)), type = "filter")
offsetValues_IWLS[which(offsetValues_IWLS <= 0)] <- 1


# weights = abs(results_IWLS$point_wise_lower[indices] - results_IWLS$point_wise_upper[indices])/4
weights <- family(results_IWLS$fit)$linkinv(results_IWLS$se[indices]**2)
# weights <- weights**2
# weights <- weights/sqrt(results_IWLS$fitted_values[indices])
weights <- 1/weights
weights <- capWeights(weights, 0.01)
weights <- weights/mean(weights)
# weights[c(1:50)] <- 0.1
par(mfrow = c(1,1))
plot(weights)

which_model <- nb(link = "log")
DAT_IWLS <- data.frame(resp = results_IWLS$fitted_values[indices], offs = offsetValues_IWLS, x = c(1:length(results_IWLS_final_with_weekend_effect$fitted_values[indices])))
fit_IWLS_delta <- fitModel(DAT_IWLS, which_model = which_model, weights_ = weights, fm = 'resp ~ s(x, k = 70)', ignore_first_max = 5,
                     deriv_plots = c(0,1), detect_flat = TRUE, x_axis = indices)

ggplot(data = NULL, aes(x = c(1:length(fit_IWLS_delta$fitted.values)))) +
  annotate(geom = "rect", xmin = 60, xmax = 67, ymin = -Inf, ymax = Inf, fill = "#2a9d8f", alpha = 0.2)+
  annotate(geom = "rect", xmin = 90, xmax = 97, ymin = -Inf, ymax = Inf, fill = "#e63946", alpha = 0.2)+
  geom_ribbon(aes(ymin = fit_IWLS_delta$lower, ymax = fit_IWLS_delta$upper, fill = "95% CI", colour = "white"), alpha = 0.3, fill = colors[3]) +
  geom_line(aes(y = true_reproductive_number, colour = "True reproductive number"), size = 1) + 
  geom_line(aes(y = results_IWLS$fitted[indices]/offsetValues_IWLS, colour = "Cori et al. estimate"), size = 1) + 
  geom_line(aes(y = fit_IWLS_delta$fitted.values, colour = "Estimated reproductive number"), size = 1) + 
  labs(title="Reproductive number", x="Time", y=latex2exp::TeX("$R_t$")) +
  scale_colour_manual("", 
                      breaks = c("True reproductive number", "Estimated reproductive number", "Cori et al. estimate", "95% CI"),
                      values = c("#1d3557", colors[3], colors[1], colors[3])) +
  scale_fill_manual("", breaks = c("95% CI"), values = colors[3]) + 
  geom_hline(yintercept = 1, col = "red") +
  coord_cartesian(ylim = c(0, 3)) + 
  theme_for_the_plots + 
  guides(color=guide_legend(override.aes=list(fill=NA)))

# plot(fit_IWLS_delta$derivPlots$plot)


# Comparison plot delta/ no delta:
ggplot(data = NULL, aes(x = c(1:length(fit_IWLS_delta$fitted.values)))) +
  annotate(geom = "rect", xmin = 60, xmax = 67, ymin = -Inf, ymax = Inf, fill = "#2a9d8f", alpha = 0.2)+
  annotate(geom = "rect", xmin = 90, xmax = 97, ymin = -Inf, ymax = Inf, fill = "#e63946", alpha = 0.2)+
  geom_ribbon(aes(ymin = fit_IWLS_no_delta$lower, ymax = fit_IWLS_no_delta$upper), alpha = 0.3, fill = colors[3]) +
  geom_ribbon(aes(ymin = fit_IWLS_delta$lower, ymax = fit_IWLS_delta$upper), alpha = 0.3, fill = colors[1]) +
  geom_line(aes(y = true_reproductive_number, colour = "True reproductive number"), size = 1) + 
  geom_line(aes(y = fit_IWLS_no_delta$fitted.values, colour = "Estimated reproductive number, no penalization"), size = 1) + 
  geom_line(aes(y = fit_IWLS_delta$fitted.values, colour = "Estimated reproductive number, with penalization"), size = 1) + 
  labs(title="Estimated reproductive number, based on estimated infection events with weekly pattern", x="Time", y=latex2exp::TeX("$R_t$")) +
  scale_colour_manual("", 
                      breaks = c("True reproductive number", 
                                 "Estimated reproductive number, no penalization", 
                                 "Estimated reproductive number, with penalization"),
                      values = c("#1d3557", colors[3], colors[1])) +
  geom_hline(yintercept = 1, col = "red") +
  coord_cartesian(ylim = c(0, 3)) + 
  theme_for_the_plots + 
  guides(color=guide_legend(override.aes=list(fill=NA))) #+ theme(legend.position = 'bottom')



# - - - - - - - - - - - - - - - - - Adding lockdown effect --------------------------------
# In this part, we explore the effects of adding lockdown term in the estiamation of R. In particular, we add a constant term 
# in the linear term of the generalized linear model.
b_lockdown <- 1*(c(1:length(offsetValues_IWLS)) > 60)
b_lockdown[90:length(b_lockdown)] <- 0
event_name <- "lockdown"

DAT_IWLS_lockdown <- DAT_IWLS
DAT_IWLS_lockdown <- cbind(DAT_IWLS_lockdown, b_lockdown)
colnames(DAT_IWLS_lockdown) <- c(colnames(DAT_IWLS_lockdown)[1:3], event_name)

fit_IWLS_lockdown <- fitModel(DAT_IWLS_lockdown, which_model = nb(link = 'log'), weights_ = weights, fm = 'resp ~ s(x, k = 70) + lockdown', ignore_first_max = 5,
                     deriv_plots = c(0,1), detect_flat = TRUE)

# summary(fit_IWLS_lockdown$fit)

ggplot(data = NULL, aes(x = c(1:length(fit_IWLS_lockdown$fitted.values)))) +
  annotate(geom = "rect", xmin = 60, xmax = 67, ymin = -Inf, ymax = Inf, fill = "#2a9d8f", alpha = 0.2)+
  annotate(geom = "rect", xmin = 90, xmax = 97, ymin = -Inf, ymax = Inf, fill = "#e63946", alpha = 0.2)+
  geom_ribbon(aes(ymin = fit_IWLS_lockdown$lower, ymax = fit_IWLS_lockdown$upper, fill = "95% CI", colour = "white"), alpha = 0.3, fill = colors[3]) +
  geom_line(aes(y = true_reproductive_number, colour = "True reproductive number"), size = 1) + 
  # geom_line(aes(y = results_IWLS_final_with_weekend_effect$fitted[indices]/offsetValues_IWLS, colour = "Cori et al. estimate"), size = 1) + 
  geom_line(aes(y = fit_IWLS_lockdown$fitted.values, colour = "Estimated reproductive number"), size = 1) + 
  labs(title="Reproductive number", x="Time", y=latex2exp::TeX("$R_t$")) +
  scale_colour_manual("", 
                      breaks = c("True reproductive number", "Estimated reproductive number", "Cori et al. estimate", "95% CI"),
                      values = c("#1d3557", colors[3], colors[1], colors[3])) +
  scale_fill_manual("", breaks = c("95% CI"), values = colors[3]) + 
  geom_hline(yintercept = 1, col = "red") +
  coord_cartesian(ylim = c(0, 3)) + 
  theme_for_the_plots + 
  guides(color=guide_legend(override.aes=list(fill=NA)))

# Results:
# The effects of the lockdown are significant only when the dimension for the spline basis is low enough, otherwise it is completely absorbed into the 
# spline function. 


# We can try two add two effects: one for the lockdown, and one for after the lockdown:
b_lockdown <- 1*(c(1:length(offsetValues_IWLS)) >= 60)
b_lockdown[90:length(b_lockdown)] <- 0
b_after_lockdown <- 1*(c(1:length(offsetValues_IWLS)) >= 90)

DAT_IWLS_lockdown <- DAT_IWLS
DAT_IWLS_lockdown <- cbind(DAT_IWLS_lockdown, b_lockdown, b_after_lockdown)
colnames(DAT_IWLS_lockdown) <- c(colnames(DAT_IWLS_lockdown)[1:3], "lockdown", "after_lockdown")

fit_IWLS_lockdown <- fitModel(DAT_IWLS_lockdown, which_model = quasipoisson(link = 'log'), weights_ = weights, fm = 'resp ~ s(x, k = 70) + lockdown + after_lockdown', ignore_first_max = 5,
                              deriv_plots = c(0,1), detect_flat = TRUE)

# summary(fit_IWLS_lockdown$fit)

# Report: Figure 12
ggplot(data = NULL, aes(x = c(1:length(fit_IWLS_lockdown$fitted.values)))) +
  annotate(geom = "rect", xmin = 60, xmax = 67, ymin = -Inf, ymax = Inf, fill = "#2a9d8f", alpha = 0.2)+
  annotate(geom = "rect", xmin = 90, xmax = 97, ymin = -Inf, ymax = Inf, fill = "#e63946", alpha = 0.2)+
  geom_ribbon(aes(ymin = fit_IWLS_lockdown$lower, ymax = fit_IWLS_lockdown$upper, fill = "95% CI", colour = "white"), alpha = 0.3, fill = colors[3]) +
  geom_line(aes(y = true_reproductive_number, colour = "True reproductive number"), size = 1) + 
  # geom_line(aes(y = results_IWLS_final_with_weekend_effect$fitted[indices]/offsetValues_IWLS, colour = "Cori et al. estimate"), size = 1) + 
  geom_line(aes(y = fit_IWLS_lockdown$fitted.values, colour = "Estimated reproductive number"), size = 1) + 
  labs(title="Estimated reproductive number, with basis dimension 70 and lockdown effects", x="Time", y=latex2exp::TeX("$R_t$")) +
  scale_colour_manual("", 
                      breaks = c("True reproductive number", "Estimated reproductive number", "Cori et al. estimate", "95% CI"),
                      values = c("#1d3557", colors[3], colors[1], colors[3])) +
  scale_fill_manual("", breaks = c("95% CI"), values = colors[3]) + 
  geom_hline(yintercept = 1, col = "red") +
  coord_cartesian(ylim = c(0, 3)) + 
  theme_for_the_plots + 
  guides(color=guide_legend(override.aes=list(fill=NA)))

# Results:
# The estimated parameters for the lockdown and after-lockdown highly depends on the dimension for the spline basis, which is problematic 
# for the situations where we want to quantify the effects of restrictive measures...


# _-**-_-**-_-**-_-**-_-**-_ 2) SIR model -------------------------

# Data preparation
data_SIR <- data.frame(matrix(ncol = 7, nrow = length(df$incidence)))
colnames(data_SIR) <- c("Date", "confirmed", "deaths", "recovered", "new_confirmed", "new_deaths", "new_recovered")
data_SIR$Date <- c(1:length(df$incidence))
data_SIR$confirmed <- df$E + df$I
data_SIR$deaths <- 0
data_SIR$recovered <- df$R
data_SIR$confirmed <- data_SIR$confirmed + data_SIR$recovered

data_SIR$new_confirmed <- c(0, diff(df$E + df$I))
data_SIR$new_deaths <- 0
data_SIR$new_recovered <- c(0, diff(data_SIR$recovered))

data_SIR <- data_SIR[indices_original_synthetic,]


# - - - - - - - - - - - - - - - - - First, without deconvolution (using true events) ----------------------
#                                             Without constrained rates, simultaneous estimation ----------
# Estimate R
result_gam_SIR <- gamSIR(data_SIR, 100, 100, 2e6, family = poisson(link = 'identity'), 
                         first_ind = NULL, splineSimulations = TRUE, debug = TRUE, basic_r = FALSE)

ggplot(data = NULL, aes(x = result_gam_SIR$indices)) +
  annotate(geom = "rect", xmin = 60, xmax = 67, ymin = -Inf, ymax = Inf, fill = "#2a9d8f", alpha = 0.2)+
  annotate(geom = "rect", xmin = 90, xmax = 97, ymin = -Inf, ymax = Inf, fill = "#e63946", alpha = 0.2)+
  geom_ribbon(aes(ymin = result_gam_SIR$lower, ymax = result_gam_SIR$upper, fill = "95% CI", colour = "white"), alpha = 0.3, fill = colors[3]) +
  geom_path(data = result_gam_SIR$fits_splines, mapping = aes(y = values, x = age, group = ind), alpha = 0.4, colour = "grey20") +
  geom_line(aes(y = df$true_rt[2:(length(df$true_rt))][result_gam_SIR$indices], colour = "True reproductive number"), size = 1) + 
  geom_line(aes(y = result_gam_SIR$Reproductive_number, colour = "Estimated reproductive number"), size = 1) + 
  labs(title="Estimates of the reproductive number", x="Time", y=latex2exp::TeX("$R_t$")) +
  scale_colour_manual("", 
                      breaks = c("True reproductive number", "Estimated reproductive number", "Cori et al. estimate", "95% CI"),
                      values = c("#f9c80e", colors[3], colors[1], colors[3])) +
  scale_fill_manual("", breaks = c("95% CI"), values = colors[3]) + 
  geom_hline(yintercept = 1, col = "red") +
  coord_cartesian(ylim = c(0, 3)) + 
  theme_for_the_plots + 
  guides(color=guide_legend(override.aes=list(fill=NA)))

# Results:
# The fit is good, bu there are large instabilities (sample splines), and large confidence bands at the beginning and the end.
# The main reason is that the rates are not constrained to be positive, which can lead to explosion in the reproductive number, as 
# R_t = transmission_rate_t / removal_rate_t.

# Here, we use a quasipoisson model
result_gam_SIR <- gamSIR(data_SIR, 100, 100, 2e6, family = quasipoisson(link = 'identity'), 
                         first_ind = NULL, splineSimulations = TRUE, debug = TRUE, basic_r = FALSE)

ggplot(data = NULL, aes(x = result_gam_SIR$indices)) +
  annotate(geom = "rect", xmin = 60, xmax = 67, ymin = -Inf, ymax = Inf, fill = "#2a9d8f", alpha = 0.2)+
  annotate(geom = "rect", xmin = 90, xmax = 97, ymin = -Inf, ymax = Inf, fill = "#e63946", alpha = 0.2)+
  geom_ribbon(aes(ymin = result_gam_SIR$lower, ymax = result_gam_SIR$upper, fill = "95% CI", colour = "white"), alpha = 0.3, fill = colors[3]) +
  geom_path(data = result_gam_SIR$fits_splines, mapping = aes(y = values, x = age, group = ind), alpha = 0.4, colour = "grey20") +
  geom_line(aes(y = df$true_rt[1:(length(df$true_rt)-1)][result_gam_SIR$indices], colour = "True reproductive number"), size = 1) + 
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

# Results: 
# Much better !

#                                             With constrained rates, two step estimation -----------------
# Here, we will add constraints to the removal rates, but we will need to perform the estimation in two steps,
# accounting for the estimates uncertainties for the first one in the estimation of the second. We estimate the 
# removal rate first (since removals depends only on the block 'infected'), and then the transmission rate, which depends
# on the block 'suceptible' and 'removal')
result_gam_SIRExp <- gamSIRExp(data_SIR, 100, 100, 2e6, first_ind = NULL, family = poisson(link = 'log'), 
                               splineSimulations = TRUE, debug = TRUE, basic_r = FALSE)

# Reproduction number
# Report: Figure 13
ggplot(data = NULL, aes(x = result_gam_SIR$indices)) +
  annotate(geom = "rect", xmin = 60, xmax = 67, ymin = -Inf, ymax = Inf, fill = "#2a9d8f", alpha = 0.2)+
  annotate(geom = "rect", xmin = 90, xmax = 97, ymin = -Inf, ymax = Inf, fill = "#e63946", alpha = 0.2)+
  geom_ribbon(aes(ymin = result_gam_SIRExp$Reproductive_number$lower, ymax = result_gam_SIRExp$Reproductive_number$upper, 
                  fill = "95% CI", colour = "white"), alpha = 0.3, fill = colors[3]) +
  # geom_path(data = result_gam_SIR$fits_splines, mapping = aes(y = values, x = age, group = ind), alpha = 0.4, colour = "grey20") +
  geom_line(aes(y = df$true_rt[1:(length(df$true_rt)-1)][result_gam_SIR$indices], colour = "True reproductive number"), size = 1) + 
  geom_line(aes(y = result_gam_SIRExp$Reproductive_number$fitted, colour = "Estimated reproductive number"), size = 1) + 
  labs(title="Reproductive number", x="Time", y=latex2exp::TeX("$R_t$")) +
  scale_colour_manual("", 
                      breaks = c("True reproductive number", "Estimated reproductive number", "Cori et al. estimate", "95% CI"),
                      values = c("#1d3557", colors[3], colors[1], colors[3])) +
  scale_fill_manual("", breaks = c("95% CI"), values = colors[3]) + 
  geom_hline(yintercept = 1, col = "red") +
  coord_cartesian(ylim = c(0, 3)) + 
  theme_for_the_plots + 
  guides(color=guide_legend(override.aes=list(fill=NA))) + theme(legend.position="bottom")

# Removal rate:
p_removal <- ggplot(data = NULL, aes(x = result_gam_SIR$indices)) +
  annotate(geom = "rect", xmin = 60, xmax = 67, ymin = -Inf, ymax = Inf, fill = "#2a9d8f", alpha = 0.2)+
  annotate(geom = "rect", xmin = 90, xmax = 97, ymin = -Inf, ymax = Inf, fill = "#e63946", alpha = 0.2)+
  geom_ribbon(aes(ymin = result_gam_SIRExp$removal_rate$lower, ymax = result_gam_SIRExp$removal_rate$upper, 
                  fill = "95% CI", colour = "white"), alpha = 0.3, fill = colors[3]) +
  geom_line(aes(y = result_gam_SIRExp$removal_rate$fitted, colour = "Estimated removal rate"), size = 1) + 
  labs(title="Estimation of the removal rate", x="Time", y="removal rate") +
  scale_colour_manual("", 
                      breaks = c("True reproductive number", "Estimated removal rate", "Cori et al. estimate", "95% CI"),
                      values = c("#1d3557", colors[3], colors[1], colors[3])) +
  scale_fill_manual("", breaks = c("95% CI"), values = colors[3]) + 
  #geom_hline(yintercept = 1, col = "red") +
  # coord_cartesian(ylim = c(0, 3)) + 
  theme_for_the_plots + 
  guides(color=guide_legend(override.aes=list(fill=NA))) + theme(legend.position="bottom")

# Transmission rate:
p_transmission <- ggplot(data = NULL, aes(x = result_gam_SIR$indices)) +
  annotate(geom = "rect", xmin = 60, xmax = 67, ymin = -Inf, ymax = Inf, fill = "#2a9d8f", alpha = 0.2)+
  annotate(geom = "rect", xmin = 90, xmax = 97, ymin = -Inf, ymax = Inf, fill = "#e63946", alpha = 0.2)+
  geom_ribbon(aes(ymin = result_gam_SIRExp$transmission_rate$lower, ymax = result_gam_SIRExp$transmission_rate$upper, 
                  fill = "95% CI", colour = "white"), alpha = 0.3, fill = colors[3]) +
  geom_line(aes(y = result_gam_SIRExp$transmission_rate$fitted, colour = "Estimated transmission rate"), size = 1) + 
  labs(title="Estimation of the transmission rate", x="Time", y="transmission rate") +
  scale_colour_manual("", 
                      breaks = c("True reproductive number", "Estimated transmission rate", "Cori et al. estimate", "95% CI"),
                      values = c("#1d3557", colors[3], colors[1], colors[3])) +
  scale_fill_manual("", breaks = c("95% CI"), values = colors[3]) + 
  #geom_hline(yintercept = 1, col = "red") +
  # coord_cartesian(ylim = c(0, 3)) + 
  theme_for_the_plots + 
  guides(color=guide_legend(override.aes=list(fill=NA))) + theme(legend.position="bottom")


plot_grid(p_removal, p_transmission, nrow = 2, align = 'hv')


# With the sample splines : 
ggplot(data = NULL, aes(x = result_gam_SIRExp$indices)) +
  annotate(geom = "rect", xmin = 60, xmax = 67, ymin = -Inf, ymax = Inf, fill = "#2a9d8f", alpha = 0.2)+
  annotate(geom = "rect", xmin = 90, xmax = 97, ymin = -Inf, ymax = Inf, fill = "#e63946", alpha = 0.2)+
  geom_ribbon(aes(ymin = result_gam_SIRExp$Reproductive_number$lower, ymax = result_gam_SIRExp$Reproductive_number$upper, 
                  fill = "95% CI", colour = "white"), alpha = 0.3, fill = colors[3]) +
  # geom_path(data = result_gam_SIR$fits_splines, mapping = aes(y = values, x = age, group = ind), alpha = 0.4, colour = "grey20") +
  geom_line(aes(y = df$true_rt[1:(length(df$true_rt)-1)][result_gam_SIRExp$indices], colour = "True reproductive number"), size = 1) + 
  geom_line(aes(y = result_gam_SIRExp$Reproductive_number$fitted, colour = "Estimated reproductive number"), size = 1) + 
  geom_path(data = result_gam_SIRExp$sim2, mapping = aes(y = values, x = age, group = ind), alpha = 0.4, colour = "grey20") + 
  labs(title="Reproductive number", x="Time", y=latex2exp::TeX("$R_t$")) +
  scale_colour_manual("", 
                      breaks = c("True reproductive number", "Estimated reproductive number", "Cori et al. estimate", "95% CI"),
                      values = c("#1d3557", colors[3], colors[1], colors[3])) +
  scale_fill_manual("", breaks = c("95% CI"), values = colors[3]) + 
  geom_hline(yintercept = 1, col = "red") +
  coord_cartesian(ylim = c(0, 3)) + 
  theme_for_the_plots + 
  guides(color=guide_legend(override.aes=list(fill=NA))) + theme(legend.position="bottom")




# - - - - - - - - - - - - - - - - - With deconvolution, sensitivity analysis ----------------------------------------------------

data_considered <- data_SIR

n_w_conv <- 20
shape_w_conv <- 10
scale_w_conv <- 0.5
w_conv<- discrete_gamma(n = n_w_conv, normalize = TRUE, shape = shape_w_conv, scale = scale_w_conv)

# We convolve both curves:
data_confirmed_convolved <- convolveData(w_conv, data_considered$confirmed - data_considered$recovered)
data_recovered_convolved <- convolveData(w_conv, c(0, diff(data_considered$recovered)))

data_recovered_convolved[which(data_recovered_convolved <= 0)] <- 1e-10

w_conv_I_deconv_shape <- c(10, 10, 10, 10, 10)
w_conv_I_deconv_scale <- c(0.5, 0.5, 0.5, 0.5, 0.5)

w_conv_R_deconv_shape <- 2*c(2, 5, 6, 10, 15)
w_conv_R_deconv_scale <- c(0.5, 0.5, 0.5, 0.5, 0.5)
parameter_changed <- c("Mean: 2", "Mean: 5", "Mean: 6","Mean: 10", "Mean: 15")
# colors_comparison <- c("#00AFBB", "#E7B800", "#FC4E07", "red")
# colors_comparison <- c("#ffbe0b", "#3a86ff", "#ff006e", "#fb5607", "#8338ec")
colors_comparison <- c("#d56062", "#f37748", "#ecc30b", "#84bcda", "#067bc2")

l_res <- list()

for (j in c(1:length(parameter_changed))){
  print(parameter_changed[j])
  
  # Next, we want to see the effects of having the incidence curves and recovered curves not 'synchronised', that is, when one of the two was not deconvolved properly...
  data_SIR_bis <- data_considered
  
  n_w_conv_I <- 20
  w_conv_I_deconv <- discrete_gamma(n = 20, normalize = TRUE, shape = w_conv_I_deconv_shape[j], scale = w_conv_I_deconv_scale[j])
  w_conv_R_deconv <- discrete_gamma(n = 20, normalize = TRUE, shape = w_conv_R_deconv_shape[j], scale = w_conv_R_deconv_scale[j])
  
  # DECONVOLUTION:
  nb_knots <- 100
  results_IWLS_I <- DeconvolutionGAM(data_confirmed_convolved, w_conv_I_deconv, nb_knots, plot_ = FALSE, get_derivatives = FALSE, family = quasipoisson(link = 'log'),
                                     extension = list('len' = 20, 'family' = nb(link = 'log'), 'keep_original' = TRUE, 'min_weights' = 0.1),
                                     fm = "cases ~ s(time, k = 20, bs = \'tp\', m = c(0, 1)) + s(day, k = 5)")
  
  results_IWLS_R <- DeconvolutionGAM(data_recovered_convolved, w_conv_R_deconv, nb_knots, plot_ = FALSE, get_derivatives = FALSE, family = quasipoisson(link = 'log'),
                                     extension = list('len' = 20, 'family' = nb(link = 'log'), 'keep_original' = TRUE, 'min_weights' = 0.1),
                                     fm = "cases ~ s(time, k = 20, bs = \'tp\', m = c(0, 1)) + s(day, k = 5)")
  

  data_SIR_bis$confirmed <- results_IWLS_I$fitted_values[indices_original_synthetic]
  data_SIR_bis$deaths <- 0
  data_SIR_bis$recovered <- cumsum(results_IWLS_R$fitted_values[indices_original_synthetic])
  data_SIR_bis$confirmed <- data_SIR_bis$confirmed + data_SIR_bis$recovered
  
  data_SIR_bis$new_confirmed <- c(0, diff(results_IWLS_I$fitted_values[indices_original_synthetic]))
  data_SIR_bis$new_deaths <- 0
  data_SIR_bis$new_recovered <- c(0, diff(data_SIR_bis$recovered))
  
  # Then, estimate R
  R_SIR_miss <- gamSIRExp(data_SIR_bis, 100, 100, 2e6, first_ind = NULL, family = poisson(link = 'log'), 
                          splineSimulations = FALSE, debug = FALSE, basic_r = FALSE)
  l_res[[j]] <- R_SIR_miss
}

plot_comparison <- ggplot(data = NULL, aes(x = c(1:(length(df$time)-1)))) + 
  geom_line(aes(x = c(1:length(df$true_rt[1:(length(df$true_rt)-1)][result_gam_SIR$indices])), y = df$true_rt[1:(length(df$true_rt)-1)][result_gam_SIR$indices]))
plot_comparison <- plot_comparison + 
  annotate(geom = "rect", xmin = 60, xmax = 67, ymin = -Inf, ymax = Inf, fill = "#2a9d8f", alpha = 0.2)+
  annotate(geom = "rect", xmin = 90, xmax = 97, ymin = -Inf, ymax = Inf, fill = "#e63946", alpha = 0.2)+
  geom_line(aes(x = c(1:length(l_res[[1]]$Reproductive_number$fitted)), y = l_res[[1]]$Reproductive_number$fitted, color = parameter_changed[1])) + 
  geom_ribbon(aes(x = c(1:length(l_res[[1]]$Reproductive_number$lower)), 
                  ymin = l_res[[1]]$Reproductive_number$lower, 
                  ymax = l_res[[1]]$Reproductive_number$upper, 
                  fill = parameter_changed[1]), alpha = 0.2) + 
  geom_line(aes(x = c(1:length(l_res[[2]]$Reproductive_number$fitted)), y = l_res[[2]]$Reproductive_number$fitted, color = parameter_changed[2]), size = 1) + 
  geom_ribbon(aes(x = c(1:length(l_res[[2]]$Reproductive_number$lower)), 
                  ymin = l_res[[2]]$Reproductive_number$lower, 
                  ymax = l_res[[2]]$Reproductive_number$upper, 
                  fill = parameter_changed[2]), alpha = 0.2) + 
  geom_line(aes(x = c(1:length(l_res[[3]]$Reproductive_number$fitted)), y = l_res[[3]]$Reproductive_number$fitted, color = parameter_changed[3])) + 
  geom_ribbon(aes(x = c(1:length(l_res[[3]]$Reproductive_number$lower)), 
                  ymin = l_res[[3]]$Reproductive_number$lower, 
                  ymax = l_res[[3]]$Reproductive_number$upper, 
                  fill = parameter_changed[3]), alpha = 0.2) + 
  geom_line(aes(x = c(1:length(l_res[[4]]$Reproductive_number$fitted)), y = l_res[[4]]$Reproductive_number$fitted, color = parameter_changed[4])) + 
  geom_ribbon(aes(x = c(1:length(l_res[[4]]$Reproductive_number$lower)), 
                  ymin = l_res[[4]]$Reproductive_number$lower, 
                  ymax = l_res[[4]]$Reproductive_number$upper, 
                  fill = parameter_changed[4]), alpha = 0.2) + 
  geom_line(aes(x = c(1:length(l_res[[5]]$Reproductive_number$fitted)), y = l_res[[5]]$Reproductive_number$fitted, color = parameter_changed[5])) +
  geom_ribbon(aes(x = c(1:length(l_res[[5]]$Reproductive_number$lower)),
                  ymin = l_res[[5]]$Reproductive_number$lower,
                  ymax = l_res[[5]]$Reproductive_number$upper,
                  fill = parameter_changed[5]), alpha = 0.2)

plot_comparison <- plot_comparison + 
  geom_hline(yintercept = 1, col = "red") + 
  scale_colour_manual("", 
                      breaks = parameter_changed,
                      values = colors_comparison) +
  scale_fill_manual("", breaks = parameter_changed, values = colors_comparison) + 
  theme_for_the_plots + 
  labs(title = "Estimation of the reproduction number with SIR model, using different mean of delay recovery", x = "Time", y = latex2exp::TeX("$R_t$")) + 
  coord_cartesian(ylim = c(0, 4))

# Report: Figure 14
plot(plot_comparison)












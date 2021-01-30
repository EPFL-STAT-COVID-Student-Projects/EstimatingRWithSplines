# ===================== PACKAGES =====================
library(mgcv)
library(gridExtra)

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

# -------------------------------- GAM model --------------------------------
#' Approximation of the discretized gamma interval with a negative binomial distribution, so that the expectation and variance match the specified one
#' @param mean mean of the distribution
#' @param variance variance of the distribution
#' @param n maximum number of bins, ie. get the distribution between 0 and n
#' @param check TRUE to print the empirical mean and variance 
#' @note The variance must be greater than the mean, since pr/(1-p) < pr/(1-p)**2. 
getGenerationInterval <- function(mean, variance, n, check = FALSE){
  p <- mean/variance
  m <- mean*p/(1-p)
  W_gam <- dnbinom(c(0:n), m, p)
  if (check){
    print(sum(c(0:(length(W_gam)-1)) * W_gam))
    print(sum(c(0:(length(W_gam)-1))**2 * W_gam) - sum(c(0:(length(W_gam)-1)) * W_gam)**2)
  }
  return(W_gam)
}
# W_gam <- getGenerationInterval(7, 8, 50, TRUE)


#' Function to estimate the reproductive number, assuming it is a smooth function of the time.
#' @param DATA_ dataframe with columns 'resp' (response variable), 'x' (the time, what the smooth function is based on) and 'offs' (the offset)
#' @param weights weights for the gam object. Usually the inverse of the variance of the response variable.
#' @param which_model poisson, quasipoisson, nb(link = 'log'), poisson(link = 'identity'), etc
#' @param formula formula oject (see lm, gam), by default 'as.formula('resp ~ s(x)')'
#' @return List of the fitted values, the upper and lower confidence bands and the gam object from the fit
#' @examples
#' fitModel(DAT_roll_back, weights = NULL, which_model = quasipoisson,  fm = as.formula('resp ~ s(x)'))
#' fitModel(DAT_roll_back, weights = weights, which_model = which_model,  fm = as.formula('resp ~ s(x, k = 290)'))
fitModel <- function(DATA_, which_model = poisson, fm = 'resp ~ s(x)', weights_ = NULL, derivatives = TRUE, ...){
  fm <- as.formula(fm)
  fit <- gam(fm, family = which_model, offset = log(offs), data = DATA_, weights = weights_)
  ilink <- family(fit)$linkinv
  fitted_values <- predict(fit,newdata=DATA_,type="link",se.fit = TRUE) # type = link returns the linear predictor, without using the link function
  upr <- ilink(fitted_values$fit + 2*fitted_values$se.fit)
  lwr <- ilink(fitted_values$fit - 2*fitted_values$se.fit)
  
  if (derivatives){
    plot_deriv <- FirstAndSecondDerivativesGAM(DATA_, fit, 1e-7, sample_splines = TRUE, ilink_fct = ilink, ...)
  }else{
    plot_deriv <- NULL
  }

  return(list("fitted.values" = ilink(fitted_values$fit), "upper" = upr, "lower" = lwr, "fit" = fit, "derivPlots" = plot_deriv, "se" = fitted_values$se.fit))
}

#' Function to bound the weights that are above (resp. below) the (1-q)-quantile (resp. q-quantile)
#' @param weights_ vector of weights that need to be bounded
#' @param q quantile
capWeights <- function(weights_, q){
  # quantile_q <- quantile(weights_, q)
  quantiles_1_q <- quantile(weights_, 1-q)
  weights_[which(weights_>quantiles_1_q)] <- quantiles_1_q
  # weights_[which(weights_<quantile_q)] <- quantile_q
  return(weights_)
}

#' Divide (element-wise) all weights by some quantity
readjust_weights <- function(weights, observations){
  weights_ <- weights/observations
  weights_[which(is.na(weights_))] <- 0
  return(weights_)
}

#' Function to compute/estimate the first and second derivatives of the reproductive number
#' @param DATA_ dataframe with columns 'resp' (response variable), 'x' (the time, what the smooth function is based on) and 'offs' (the offset)
#' @param fit model returned fy the function 'fitModel'
#' @param eps small epsilon for the finite difference approximation
#' @param deriv_plots option for computing the derivatives: if deriv_plots = c(0, 1, 2), then all derivatives up to order 2 are estimated
#' @param ilink_fct link function from the model
#' @param ... other options
FirstAndSecondDerivativesGAM <- function(DATA_, fit, eps = 1e-7, deriv_plots = c(0, 1, 2), ilink_fct = NULL, ...){
  # Finite difference matrices
  DATA_plus <- DATA_
  DATA_plus$x <- DATA_plus$x + eps  # data increased by a small eps
  DATA_minus <- DATA_
  DATA_minus$x <- DATA_minus$x - eps  # data decreased by a small eps

  # for each data (i.e. time points), we estimate R based on the fit
  pred <- predict(fit,newdata=DATA_,type="lpmatrix")
  pred2 <- predict(fit,newdata=DATA_plus,type="lpmatrix")
  pred3 <- predict(fit,newdata=DATA_minus,type="lpmatrix")
  pred_deriv <- (pred2 - pred) / eps    # first derivative
  pred_deriv2_second <- (pred2 + pred3 - 2*pred)  / eps^2   # second derivative

  # ----------------------- PLOTS, with point-wise / simultaneous CI
  list_of_plots <- list()
  simultaneousIntervals_results <- list()
  if (0 %in% deriv_plots){
    simultaneousIntervals_results[[length(simultaneousIntervals_results)+1]] <- SimultaneousIntervalsGAM(pred, fit$coefficients, vcov(fit), plot_ = TRUE, intercept_plot = 1,
                                                                                                      ilink_fct = ilink_fct, title = "Fitted values", x = "time", y = "fitted values",...)
    list_of_plots[[length(list_of_plots)+1]] <- simultaneousIntervals_results[[length(simultaneousIntervals_results)]]$plot
  }
  if (1 %in% deriv_plots){
    simultaneousIntervals_results[[length(simultaneousIntervals_results)+1]] <- SimultaneousIntervalsGAM(pred_deriv, fit$coefficients, vcov(fit), plot_ = TRUE, intercept_plot = 0,
                                                                                                      title = "First derivative", x = "time", y = "first derivative", ...)
    list_of_plots[[length(list_of_plots)+1]] <- simultaneousIntervals_results[[length(simultaneousIntervals_results)]]$plot
  }
  if (2 %in% deriv_plots){
    simultaneousIntervals_results[[length(simultaneousIntervals_results)+1]] <- SimultaneousIntervalsGAM(pred_deriv2_second, fit$coefficients, vcov(fit), plot_ = TRUE, intercept_plot = 0,
                                                                                                      title = "Second derivative", x = "time", y = "second derivative", ...)
    list_of_plots[[length(list_of_plots)+1]] <- simultaneousIntervals_results[[length(simultaneousIntervals_results)]]$plot
  }
  # Remark: for the first and second derivatives, we are not providing the inverse link function, as the finite difference method operates outside the link function
  # (ie. the difference should be [ilink(fitted_linear + eps) - ilink(fitted_linear)] / eps). Since the inverse link is exponential, we can not simplify it further 
  # (we can do the finite difference after the inverse link, but we would not be able to compute the confidence bands effectively). Since the derivative are close to 
  # 0 most of the time, we can simply approximate it with the linear term difference

  return(list("plot" = plot_grid(plotlist = list_of_plots, nrow = length(unique(deriv_plots)), align = 'hv'), "list_plots" = list_of_plots, "simultaneousIntervals" = simultaneousIntervals_results))
}


# (SAME AS IN Deconvolution.R)
#' Function to get the point-wise and simultaneous confidence bands from the gam object
#' @param SplineBasis basis fonction evaluated on the data
#' @param CovMatrix covariance matrix of the estimated parameters from the fitted model
#' @param N number of samples for computing the criterion for the simultaneous CI
#' @section ------------------ PLOT OPTIONS ------------------
#' @param plot_ TRUE to return the plot of the confidence bands
#' @param intercept_plot intercept for the horizontal line (set to NULL if you do not want to plot it)
#' @param sample_splines TRUE to also plot the splines generated by sampling N times the coefficients (from a multivariate normal distribution, with covariance CovMatrix)
#' @param x_axis the values in the x axis (NULL if not specified)
#' @param N_keep number of samples plotted
#' @references For more informations, see Ruppert et al. (2003), and the blog post https://fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/
SimultaneousIntervalsGAM <- function(SplineBasis, coefficients_model, CovMatrix, N = 10000 ,plot_ = FALSE, intercept_plot = NULL,
                                  sample_splines = FALSE, axis_x = NULL, ilink_fct = NULL, ignore_first_max = NULL, detect_flat = FALSE, N_keep = 30,...){
  # Function to sample from a multivariate gaussian distribution
  rmvn <- function(n, mu, sig) {
    L <- mroot(sig)
    m <- ncol(L)
    t(mu + L %*% matrix(rnorm(m*n), m, n))
  }

  # get the fitted values, standard errors and point-wise CI
  result_gam <- SplineBasis %*% coefficients_model
  se.fit <- sqrt(diag(SplineBasis %*% CovMatrix %*% t(SplineBasis)))
  q_0.025 <- result_gam - 2*se.fit
  q_0.975 <- result_gam + 2*se.fit

  # Sample new parameters from their estimates and the covariance matrix
  BUdiff <- rmvn(N, mu = rep(0, nrow(CovMatrix)), sig = CovMatrix)
  simDev <- SplineBasis %*% t(BUdiff)
  absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))  # divide by se.fit

  # then compute simultaneous CI
  masd <- apply(absDev, 2L, max, na.rm = TRUE)
  crit <- quantile(masd, prob = 0.95, type = 8, na.rm = TRUE)
  q_0.025_crit <- result_gam - crit*se.fit
  q_0.975_crit <- result_gam + crit*se.fit

  if (is.null(axis_x)){ axis_x <- c(1:length(se.fit)) }

  # if sample_splines, plot samples
  if (sample_splines){
    sims <- rmvn(N, mu = as.numeric(coefficients_model), sig = CovMatrix)
    fits <- SplineBasis %*% t(sims)
    nrnd <- N_keep
    rnd <- sample(N, nrnd)
    stackFits <- stack(as.data.frame(fits[, rnd]))
    stackFits <- transform(stackFits, age = rep(axis_x, length(rnd)))
  }else{
    stackFits <- NULL
  }

  # if a link function is provided, we apply it to the fitted values and CIs
  if (!is.null(ilink_fct)){
    result_gam <- ilink_fct(result_gam)
    q_0.025 <- ilink_fct(q_0.025)
    q_0.975 <- ilink_fct(q_0.975)
    q_0.025_crit <- ilink_fct(q_0.025_crit)
    q_0.975_crit <- ilink_fct(q_0.975_crit)
    if (sample_splines){
      stackFits$values <- ilink_fct(stackFits$values)
    }
  }

  # if we want to detect time ntervals for which CI exclude some particular threshold
  if (detect_flat && !is.null(intercept_plot)){
    detect_flat_points_wise <- compareWithThreshold(q_0.025, q_0.975, intercept_plot)
    detect_flat_simultaneous_wise <- compareWithThreshold(q_0.025_crit, q_0.975_crit, intercept_plot)
  }

  # -------- PLOT
  if (plot_){
    if (!is.null(ignore_first_max)){
      indices_ylim <- c(ignore_first_max:length(result_gam))
    }else{
      indices_ylim <- c(1:length(result_gam))
    }
    min_val <- min(result_gam[indices_ylim]) - 1.2*crit*se.fit[indices_ylim][which(result_gam[indices_ylim] == min(result_gam[indices_ylim]))[1]]
    max_val <- max(result_gam[indices_ylim]) + 1.2*crit*se.fit[indices_ylim][which(result_gam[indices_ylim] == max(result_gam[indices_ylim]))[1]]
    # max_val <- min(max_val, 10)
    p <- ggplot(NULL, aes(x = axis_x))

    if (detect_flat && !is.null(intercept_plot)){
      for (j in c(1:length(detect_flat_simultaneous_wise$start))){
        p <- p + annotate(geom = "rect", xmin = axis_x[detect_flat_simultaneous_wise$start[j]], xmax = axis_x[detect_flat_simultaneous_wise$end[j]],
                          ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.2)
      }
      for (j in c(1:length(detect_flat_points_wise$start))){
        p <- p + annotate(geom = "rect", xmin = axis_x[detect_flat_points_wise$start[j]], xmax = axis_x[detect_flat_points_wise$end[j]],
                          ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.2)
      }
    }
    p <- p +
      geom_hline(yintercept = intercept_plot, col = "red") +
      geom_ribbon(aes(x = axis_x, ymin = q_0.025, ymax = q_0.975), alpha = 0.2, fill = "red") +
      geom_ribbon(aes(x = axis_x, ymin = q_0.025_crit, ymax = q_0.975_crit), alpha = 0.2, fill = "red") +
      theme_for_the_plots +
      labs(...) +
      coord_cartesian(ylim = c(min_val, max_val))

    if (inherits(axis_x, "Date")){
      p <- p + scale_x_date(date_breaks = "1 month", date_labels = "%b")
    }

    if (sample_splines){
      p <- p + geom_path(data = stackFits, mapping = aes(y = values, x = age, group = ind), alpha = 0.4, colour = "grey20")
    }
    p <- p + geom_line(aes(y = result_gam), colour = '#f4e04d')   # #e09f3e
  }else{
    p <- NULL
  }
  return(list("Point_wise_CI_upper" = q_0.975, "Point_wise_CI_lower" = q_0.025, "Simultaneous_upper"=q_0.975_crit, "Simultaneous_lower"=q_0.025_crit, "fits" = stackFits, "plot" = p))
}

# -------------------------------- SIR model --------------------------------

# Function that estimates the reproductive number based on the SIR model.
# The transmission and removal rates are not constrained (must use 'identity' link function), but 
# are estimated in one step.
#' @param DATA dataframe. Must contain curves for total recovered, deaths and confirmed cases
#' @param q1 dimension basis for transmission rate
#' @param q2 dimension basis for removal rate
#' @param family family used for the model (poisson, quasipoisson, negative binomial, ...), with 'identity' link
#' @param first_ind index to start analysis (avoid having 0 cases). If NULL, computed automatically.
#' @param splineSimulations TRUE to sample splines from fitted values
#' @param debug TRUE to print computation steps
#' @param basic_r TRUE to return the basic reproduction number R0, FALSE to return the basic reproduction number multiplied by the proportion of susceptible cases (R_t)
gamSIR <- function(DATA, q1, q2, tot_pop, family = poisson(link = 'identity'), first_ind = NULL, 
                   splineSimulations = FALSE, debug = FALSE, basic_r = FALSE, log_offset = FALSE){
  # Function to sample from a multivariate gaussian distribution
  rmvn <- function(n, mu, sig) { 
    L <- mroot(sig)
    m <- ncol(L)
    t(mu + L %*% matrix(rnorm(m*n), m, n))
  }
  
  # if 'first_ind' is NULL, then we compute the first index, i.e. the index for which the data is not 0
  if (is.null(first_ind)){
    indr <- which(DATA$recovered + DATA$deaths <= 0)[1]
    indi <- which(DATA$confirmed - DATA$recovered - DATA$deaths <= 0)[1]
    if (is.na(indr)){indr <- 0}
    if (is.na(indi)){indi <- 0}
    
    if (debug){ cat(sprintf("First index computed! incidence : %s, recovered : %s", indi+1, indr+1), '\n') }
    first_ind <- max(indr+1, indi+1)
  }
  sub_indices <- c(first_ind:length(DATA$Date))

  # Time interval
  Time <- c(1:(length(DATA$Date[sub_indices])-1))
  T_len <- length(Time)
  
  # Sequence of knots for the spline basis
  knotsQ1 <- seq(from=min(Time),to=max(Time),length=q1+2)[2:(q1+1)]
  knotsQ2 <- seq(from=min(Time),to=max(Time),length=q2+2)[2:(q2+1)]
  
  # Spline basis
  SplineBasisQ1 <- bs(Time ,knots = knotsQ1, degree = 3)
  SplineBasisQ2 <- bs(Time ,knots = knotsQ2, degree = 3)
  
  r <- DATA$recovered[sub_indices] + DATA$deaths[sub_indices]
  i <- DATA$confirmed[sub_indices] - r
  i[which(i <= 0 )] <- 1
  r[which(r <= 0 )] <- 1
  
  # Build the matrix for the model
  if (debug){ cat("Building the specification matrix ... ") }
  BMat <- matrix(nrow = 2*T_len, ncol = dim(SplineBasisQ1)[2] + dim(SplineBasisQ2)[2])
  SplineBasisQ1_1 <- SplineBasisQ1 * replicate(dim(SplineBasisQ1)[2], i[1:(length(i)-1)])
  SplineBasisQ1_2 <- SplineBasisQ1 * replicate(dim(SplineBasisQ1)[2], i[1:(length(i)-1)]) * replicate(dim(SplineBasisQ1)[2], i[1:(length(i)-1)])/tot_pop
  SplineBasisQ1_3 <- SplineBasisQ1 * replicate(dim(SplineBasisQ1)[2], i[1:(length(i)-1)]) * replicate(dim(SplineBasisQ1)[2], r[1:(length(r)-1)])/tot_pop
  
  BMat[1:T_len, 1:dim(SplineBasisQ1)[2]] <- SplineBasisQ1_1 - SplineBasisQ1_2 - SplineBasisQ1_3
  BMat[1:T_len, (dim(SplineBasisQ1)[2] + 1):(dim(SplineBasisQ1)[2] + dim(SplineBasisQ2)[2])] <- - SplineBasisQ2 * replicate(dim(SplineBasisQ2)[2], i[1:(length(i)-1)])
  BMat[(T_len+1):(2*T_len), 1:dim(SplineBasisQ1)[2]] <- 0
  BMat[(T_len+1):(2*T_len), (dim(SplineBasisQ1)[2] + 1):(dim(SplineBasisQ1)[2] + dim(SplineBasisQ2)[2])] <- SplineBasisQ2 * replicate(dim(SplineBasisQ2)[2], i[1:(length(i)-1)])
  if (debug){ cat(sprintf("done. Dimension: %s x %s", dim(BMat)[1], dim(BMat)[2]), '\n')}
  
  # Prepare data for the model
  DAT_gam <- as.data.frame(BMat)
  xnam <- paste("x", 1:dim(BMat)[2], sep="")
  colnames(DAT_gam) <- xnam
  
  DAT_gam$resp <- c(i[2:length(i)], r[2:length(r)])
  DAT_gam$offset <- c(i[1:(length(i)-1)], r[1:(length(r)-1)])
  if (log_offset){
    DAT_gam$offset <- log(DAT_gam$offset)
  }
  
  # Fit the model
  if (debug){ cat("Fitting GAM ... "); start_time <- Sys.time() }
  
  fmla <- as.formula(paste("resp ~ ", paste(xnam, collapse= "+"), " - 1"))
  fit <- gam(fmla, family = family, data = DAT_gam, weights = NULL, offset = offset)
  ilink <- family(fit)$linkinv
  
  if (debug){ end_time <- Sys.time(); cat(sprintf("done. Time ellapsed: %.3fs", end_time - start_time), "\n")}
  
  # Compute the transmission and removal rates
  beta_t <- ilink(SplineBasisQ1 %*% fit$coefficients[1:dim(SplineBasisQ1)[2]])
  gamma_t <- ilink(SplineBasisQ2 %*% fit$coefficients[(dim(SplineBasisQ1)[2] + 1):(dim(SplineBasisQ1)[2] + dim(SplineBasisQ2)[2])]) 
  R0 <- beta_t/gamma_t
  
  # Sample splines from fitted values
  if (splineSimulations){
    if (debug){ start_time <- Sys.time(); cat("Sampling splines ... ") }
    
    sims <- rmvn(10000, mu = fit$coefficients , sig = vcov(fit))   # sample coefficients
    fits_beta <- SplineBasisQ1 %*% t(sims[,1:dim(SplineBasisQ1)[2]] )  # compute beta based on these coefficients
    fits_gamma <- SplineBasisQ2 %*% t(sims[,(dim(SplineBasisQ1)[2] + 1):(dim(SplineBasisQ1)[2] + dim(SplineBasisQ2)[2])] )  # compute gamma based on these coefficients
    if (!basic_r){ fits_beta <- fits_beta * (tot_pop - i[1:(length(i)-1)] - r[1:(length(i)-1)])/tot_pop}
    R0_fits <- fits_beta/fits_gamma # compute the BASIC reproductive number
    
    nrnd <- 30  # select some of them, for the visualisation
    rnd <- sample(10000, nrnd)
    stackFits <- stack(as.data.frame(R0_fits[, rnd]))
    stackFits <- transform(stackFits, age = rep(sub_indices[1:(length(sub_indices)-1)], length(rnd)))
    
    if (debug){ end_time <- Sys.time(); cat(sprintf("done. Time ellapsed: %.3fs", end_time - start_time), "\n")}
  }else{ stackFits <- NULL }
  
  # Confidence intervals
  inf_b <- sup_b <- numeric((dim(BMat)[1])/2)
  for (j in 1:(dim(BMat)[1]/2)){
    der <- as.matrix(c(SplineBasisQ1[j,]/gamma_t[j], -(R0[j]/gamma_t[j]) * SplineBasisQ2[j,]))
    va <- t(der) %*% vcov(fit) %*% der
    se <- sqrt(va)
    inf_b[j] <- R0[j] - 2*se
    sup_b[j] <- R0[j] + 2*se
  }
  
  # Compute the reproductive number
  if (!basic_r){
    beta_t <- beta_t * (tot_pop - i[1:(length(i)-1)] - r[1:(length(i)-1)])/tot_pop
    R0 <- beta_t/gamma_t
    inf_b <- inf_b * (tot_pop - i[1:(length(i)-1)] - r[1:(length(i)-1)])/tot_pop
    sup_b <- sup_b * (tot_pop - i[1:(length(i)-1)] - r[1:(length(i)-1)])/tot_pop
  }
  
  # p <- ggplot(NULL, aes(x = Time)) +
  #   geom_line(aes(y = R0)) +
  #   geom_hline(yintercept = 1, col = "red") +
  #   geom_ribbon(aes(ymin = inf_b, ymax = sup_b, fill = "95% CI", colour = "white"), alpha = 0.5, fill = colors[3]) +
  #   theme_for_the_plots +
  #   geom_path(data = stackFits, mapping = aes(y = values, x = age, group = ind), alpha = 0.4, colour = "grey20") +
  #   ylim(0, 5) +
  #   scale_colour_manual("",
  #                       breaks = c("True reproductive number", "Estimated reproductive number", "Cori et al. estimate", "95% CI"),
  #                       values = c("black", colors[3], colors[1], colors[3])) +
  #   scale_fill_manual("", breaks = c("95% CI"), values = colors[3])
  # 
  # plot(p)
  
  return(list("Reproductive_number"=R0, "lower" = inf_b, "upper"=sup_b, 
              "transmission_rate" = beta_t, "removal_rate" = gamma_t, 
              "fits_splines" = stackFits, "indices" = sub_indices[1:(length(sub_indices)-1)]))
}


# Function that estimates the reproductive number based on the SIR model. Positive transmission and removal rates, but estimation in two steps
#' @param DATA dataframe. Must contain curves for total recovered, deaths and confirmed cases
#' @param q1 dimension basis for transmission rate
#' @param q2 dimension basis for removal rate
#' @param family family used for the model (poisson, quasipoisson, negative binomial, ...), with 'identity' link
#' @param first_ind index to start analysis (avoid having 0 cases). If NULL, computed automatically.
#' @param splineSimulations TRUE to sample splines from fitted values
#' @param debug TRUE to print computation steps
#' @param basic_r TRUE to return the basic reproduction number R0, FALSE to return the basic reproduction number multiplied by the proportion of susceptible cases (R_t)
#' @param penalty_gamma penalty for removal rate
#' @param penalty_beta penalty for transmission rate
#' @param weights_for_r weights for estimating gamma
#' @param weights_for_i weights for estimating beta
gamSIRExp <- function(DATA, q1, q2, tot_pop, family = poisson(link = "log"), first_ind = NULL, 
                      splineSimulations = FALSE, debug = TRUE, basic_r = FALSE, 
                      penalty_gamma = NULL, penalty_beta = NULL, weights_for_r = NULL, weights_for_i = NULL,
                      log_offset = TRUE){
  # Function to sample from a multivariate gaussian distribution
  rmvn <- function(n, mu, sig) { 
    L <- mroot(sig)
    m <- ncol(L)
    t(mu + L %*% matrix(rnorm(m*n), m, n))
  }
  
  # if 'first_ind' is NULL, then we compute the first index, i.e. the index for which the data is not 0
  if (is.null(first_ind)){
    indr <- which(DATA$recovered + DATA$deaths <= 0)[1]
    indi <- which(DATA$confirmed - DATA$recovered - DATA$deaths <= 0)[1]
    if (is.na(indr)){indr <- 0}
    if (is.na(indi)){indi <- 0}
    
    if (debug){ cat(sprintf("First index computed! incidence : %s, recovered : %s", indi+1, indr+1), '\n') }
    first_ind <- max(indr+1, indi+1)
  }
  sub_indices <- c(first_ind:length(DATA$Date))
  
  # Time interval
  Time <- c(1:(length(DATA$Date[sub_indices])-1))
  T_len <- length(Time)

  r <- DATA$recovered[sub_indices] + DATA$deaths[sub_indices]
  i <- DATA$confirmed[sub_indices] - r
  i[which(i <= 0 )] <- 1
  r[which(r <= 0 )] <- 1
  
  
  # ----- First part: estimating the gamma_t --------
  if (debug){ cat("1) Estimate the removal rate ...") }
  
  diff_r <-  r[2:length(r)] - r[1:(length(r)-1)]
  diff_r[which(diff_r <= 0)] <- 1
  offset_r <- i[1:(length(i)-1)]

  DAT_gam_r <- data.frame(resp = diff_r, 
                          offs = offset_r,
                          x = c(1:length(diff_r)))
  if (log_offset){
    DAT_gam_r$offs <- log(DAT_gam_r$offs)
  }
  fm <- as.formula(paste("resp ~ s(x, k = ",q1, ")"))
  fit_r <- gam(fm, family = family, offset = offs, data = DAT_gam_r, H = penalty_gamma, weights = weights_for_r[sub_indices[c(1:(length(sub_indices)-1))]])

  ilink_r <- family(fit_r)$linkinv
  fitted_values_r <- predict(fit_r,newdata=DAT_gam_r,type="link",se.fit = TRUE) # type = link returns the linear predictor, without using the link function
  upr_r <- ilink_r(fitted_values_r$fit + 2*fitted_values_r$se.fit)
  lwr_r <- ilink_r(fitted_values_r$fit - 2*fitted_values_r$se.fit)
  
  pred_r <- predict(fit_r,newdata=DAT_gam_r,type="lpmatrix")
  simultaneous_r <- SimultaneousIntervals(pred_r, fit_r$coefficients, vcov(fit_r), plot_ = TRUE, intercept_plot = 0,
                                            title = "Fitted values for gamma_t", x = "time", y = "fitted values", ilink_fct = ilink_r, 
                                            x_axis = c(1:length(ilink_r(fitted_values_r$fit))), sample_splines = TRUE)
  
  weights_r <- ilink_r(fitted_values_r$se.fit)
  # weights_r[which(weights_r == 0)] <- 1
  weights_r <- 1/weights_r**2
  weights_r <- weights_r/ilink_r(fitted_values_r$fit)
  weights_r <- capWeights(weights_r, 0.005)
  weights_r <- weights_r/mean(weights_r)
  # can then use these weights for estimating the beta
  
  removal_rate <- list("fitted" = ilink_r(fitted_values_r$fit), "lower" = lwr_r, "upper" = upr_r, "intervals" = simultaneous_r)
  if (debug){ cat("done. \n") }
  
  # ----- Second part: estimating the beta_t, given the estimate of gamma_t --------
  if (debug){ cat("2) Estimate the transmission rate ... ") }
  
  diff_i <-  i[2:length(r)] - i[1:(length(r)-1)] + ilink_r(fitted_values_r$fit)*i[1:(length(r)-1)]
  offset_i <- i[1:(length(i)-1)] - (i[1:(length(i)-1)]**2)/tot_pop - i[1:(length(i)-1)]*r[1:(length(i)-1)]/tot_pop
  diff_i[which(diff_i <= 0)] <- 1
  
  DAT_gam_i <- data.frame(resp = diff_i, 
                          offs = offset_i,
                          x = c(1:length(diff_i)))
  if (log_offset){
    DAT_gam_i$offs <- log(DAT_gam_i$offs)
  }
  fm <- as.formula(paste("resp ~ s(x, k = ",q2, ")"))
  fit_i <- gam(fm, family = family, offset = offs, data = DAT_gam_i, H = penalty_beta, weights = weights_for_i[sub_indices[c(1:(length(sub_indices)-1))]])
  ilink_i <- family(fit_i)$linkinv
  fitted_values_i <- predict(fit_i,newdata=DAT_gam_i,type="link",se.fit = TRUE) # type = link returns the linear predictor, without using the link function
  upr_i <- ilink_i(fitted_values_i$fit + 2*fitted_values_i$se.fit)
  lwr_i <- ilink_i(fitted_values_i$fit - 2*fitted_values_i$se.fit)
  pred_i <- predict(fit_i,newdata=DAT_gam_i,type="lpmatrix")
  simultaneous_i <- SimultaneousIntervals(pred_i, fit_i$coefficients, vcov(fit_i), plot_ = TRUE, intercept_plot = 0,
                                          title = "Fitted values for beta_t", x = "time", y = "fitted values", ilink_fct = ilink_i, 
                                          x_axis = c(1:length(ilink_i(fitted_values_i$fit))), sample_splines = TRUE)
  
  transmission_rate <- list("fitted" = ilink_i(fitted_values_i$fit), "lower" = lwr_i, "upper" = upr_i, "intervals" = simultaneous_i)
  if (debug){ cat("done. \n") }
  
  # ------- Then the reproduction number -----
  # final_se <- sqrt((fitted_values_i$se.fit**2)/(fitted_values_r$fit**2) + (fitted_values_i$fit**2)*(fitted_values_r$se.fit**2)/(fitted_values_r$fit**4))
  final_se <- sqrt(fitted_values_i$se.fit**2 + fitted_values_r$se.fit**2)
  R_t <- ilink_i(fitted_values_i$fit)/ilink_r(fitted_values_r$fit)
  R_lwr <- ilink_i(fitted_values_i$fit - fitted_values_r$fit - 2*final_se)
  R_upr <- ilink_i(fitted_values_i$fit - fitted_values_r$fit + 2*final_se)

  
  # ------------------ Spline simulations, for the confidence intervals
  if (splineSimulations){
    if (debug){ cat("3) Sampling splines for transmission and removal rates ... ", '\n') }
    # Simulations : 
    N_sim <- 100
    N_for_each_sim <- 30
    
    if (TRUE){
      # 1) Estimate gamma_t (removal rate), then sample, then estimate beta_t (transmission rate)
      R_ts <- matrix(ncol = N_sim, nrow = length(fitted_values_r$fit))
      sims <- rmvn(N_sim, mu = fit_r$coefficients , sig = vcov(fit_r))
      fits_gamma <- pred_r %*% t(sims)
      
      ALLsimulations <- data.frame(values = 0, ind = 0, age = 0)
      FITS <- matrix(nrow = length(fitted_values_r$fit), ncol = 1)
      
      pb <- txtProgressBar(min = 0, max = N_sim, style = 3)
      for (j in c(1:N_sim)){
        diff_i_sim <-  i[2:length(r)] - i[1:(length(r)-1)] + ilink_r(fits_gamma[,1])*i[1:(length(r)-1)]
        DAT_gam_i_sim <- data.frame(resp = diff_i_sim,
                                    offs = offset_i,
                                    x = c(1:length(diff_i_sim)))
        fm <- as.formula(paste("resp ~ s(x, k = ",q2, ")"))
        fit_i_sim <- gam(fm, family = poisson, offset = offs, data = DAT_gam_i_sim)
        fitted_values_i_sim <- predict(fit_i_sim,newdata=DAT_gam_i_sim,type="link",se.fit = FALSE)
        R_ts[,j] <- ilink_i(fitted_values_i_sim - fits_gamma[,j])
        
        # Part 3...
        pred_i_sim <- predict(fit_i_sim,newdata=DAT_gam_i_sim,type="lpmatrix")
        sims <- rmvn(N_for_each_sim, mu = fit_i_sim$coefficients, sig = vcov(fit_i_sim))
        fits <- pred_i_sim %*% t(sims)
        
        stackFits <- stack(as.data.frame(fits))
        stackFits <- transform(stackFits, age = rep(sub_indices[1:(length(sub_indices)-1)], N_for_each_sim))
        stackFits$values <- ilink_i(stackFits$values)
        stackFits$values <- stackFits$values/rep(ilink_r(fits_gamma[,1]), N_for_each_sim)
        
        if (!basic_r){ stackFits$values <- stackFits$values * rep((tot_pop - i[1:(length(i)-1)] - r[1:(length(i)-1)])/tot_pop, N_for_each_sim)}
        
        stackFits$ind <- paste(stackFits$ind, "_", j)
        
        ALLsimulations <- rbind(ALLsimulations, stackFits)
        
        fits <- fits - replicate(N_for_each_sim, fits_gamma[,1])
        
        FITS <- cbind(FITS, fits)
        setTxtProgressBar(pb, j)
      }
      close(pb)
      ALLsimulations <- ALLsimulations[2:length(ALLsimulations[,1]),]
      FITS <- FITS[, 2:dim(FITS)[2]]
      
      absDev <- abs(sweep(FITS, 1, final_se, FUN = "/"))
      masd <- apply(absDev, 2L, max, na.rm = TRUE)
      crit <- quantile(masd, prob = 0.95, type = 8, na.rm = TRUE)
      
      # if (FALSE){
      #   plot_R_sim_3 <- ggplot(data = NULL, aes(x = c(1:length(df$true_r0[3:(length(df$true_r0)-1)])), y = df$true_r0[3:(length(df$true_r0)-1)])) + geom_line() + 
      #     geom_hline(yintercept = 1, colour = "red") + 
      #     geom_ribbon(aes(x = c(1:length(fitted_values_i$fit)), ymin = ilink_i(fitted_values_i$fit - fitted_values_r$fit - 2*final_se), 
      #                     ymax = ilink_i(fitted_values_i$fit - fitted_values_r$fit + 2*final_se)), alpha = 0.2, fill = "red") + 
      #     geom_line(aes(y = ilink_i(fitted_values_i$fit - fitted_values_r$fit)), colour = colors[1]) + 
      #     theme_for_the_plots + 
      #     labs(title = "Estimation of the reproduction number with SIR model", x = "Time", y = "Rt") + 
      #     ylim(0, 3) + 
      #     geom_path(data = ALLsimulations, mapping = aes(y = values, x = age, group = ind), alpha = 0.4, colour = "grey20")
      # }
      # plot(plot_R_sim_1)
      
      sim3 <- ALLsimulations
      
      
      nrnd <- 30
      rnd <- sample(N_sim, nrnd)
      stackFits <- stack(as.data.frame(R_ts[, rnd]))
      stackFits <- transform(stackFits, age = rep(sub_indices[1:(length(sub_indices)-1)], length(rnd)))
      
      if (!basic_r){ stackFits$values <- stackFits$values * rep((tot_pop - i[1:(length(i)-1)] - r[1:(length(i)-1)])/tot_pop, N_for_each_sim)}
      
      # if (FALSE){
      #   plot_R_sim_1 <- ggplot(data = NULL, aes(x = c(1:length(df$true_r0[3:(length(df$true_r0)-1)])), y = df$true_r0[3:(length(df$true_r0)-1)])) + geom_line() + 
      #     geom_hline(yintercept = 1, colour = "red") + 
      #     geom_ribbon(aes(x = c(1:length(fitted_values_i$fit)), ymin = ilink_i(fitted_values_i$fit - fitted_values_r$fit - 2*final_se), 
      #                     ymax = ilink_i(fitted_values_i$fit - fitted_values_r$fit + 2*final_se)), alpha = 0.2, fill = "red") + 
      #     geom_line(aes(y = ilink_i(fitted_values_i$fit)/ilink_r(fitted_values_r$fit)), colour = colors[1]) + 
      #     theme_for_the_plots + 
      #     labs(title = "Estimation of the reproduction number with SIR model", x = "Time", y = "Rt") + 
      #     ylim(0, 3) + geom_path(data = stackFits, mapping = aes(y = values, x = age, group = ind), alpha = 0.4, colour = "grey20")
      #   # plot(plot_R_sim_1)
      # }
      
      sim1 <- stackFits
      
    }
    
    
    # 2) Estimate gamma_t, them beta_t, then sample, then plot them all
    sims <- rmvn(N_sim, mu = fit_i$coefficients, sig = vcov(fit_i))
    fits <- pred_i %*% t(sims)
    nrnd <- 100
    rnd <- sample(N_sim, nrnd)
    stackFits <- stack(as.data.frame(fits[, rnd]))
    stackFits <- transform(stackFits, age = rep(sub_indices[1:(length(sub_indices)-1)], length(rnd)))
    stackFits$values <- ilink_i(stackFits$values)
    stackFits$values <- stackFits$values/rep(ilink_r(fitted_values_r$fit), length(rnd))
    
    # if (FALSE){
    #   plot_R_sim_2 <- ggplot(data = NULL, aes(x = c(1:length(df$true_r0[3:(length(df$true_r0)-1)])), y = df$true_r0[3:(length(df$true_r0)-1)])) + geom_line() + 
    #     geom_hline(yintercept = 1, colour = "red") + 
    #     geom_ribbon(aes(x = c(1:length(fitted_values_i$fit)), ymin = ilink_i(fitted_values_i$fit - fitted_values_r$fit - 2*final_se), 
    #                     ymax = ilink_i(fitted_values_i$fit - fitted_values_r$fit + 2*final_se)), alpha = 0.2, fill = "red") + 
    #     geom_line(aes(y = ilink_i(fitted_values_i$fit)/ilink_r(fitted_values_r$fit)), colour = colors[1]) + 
    #     theme_for_the_plots + 
    #     labs(title = "Estimation of the reproduction number with SIR model", x = "Time", y = "Rt") + 
    #     ylim(0, 3) + geom_path(data = stackFits, mapping = aes(y = values, x = age, group = ind), alpha = 0.4, colour = "grey20")
    #   
    #   plot(plot_grid(plot_R_sim_1, plot_R_sim_2, plot_R_sim_3, nrow = 3, align = 'hv'))
    # }
    if (!basic_r){ stackFits$values <- stackFits$values * rep((tot_pop - i[1:(length(i)-1)] - r[1:(length(i)-1)])/tot_pop, N_for_each_sim)}
    
    sim2 <- stackFits
    
  }else{ sim1 <- sim2 <- sim3 <- NULL}
  
  
  # sims <- rmvn(10000, mu = fit_i$coefficients , sig = vcov(fit_i))
  # fits_beta <- SplineBasisQ1 %*% t(sims[,1:dim(SplineBasisQ1)[2]] )
  # print(dim(fits_beta))
  # fits_gamma <- SplineBasisQ2 %*% t(sims[,(dim(SplineBasisQ1)[2] + 1):(dim(SplineBasisQ1)[2] + dim(SplineBasisQ2)[2])] )
  # R0_fits <- fits_beta/fits_gamma
  # # R0_fits <- fits_beta/rep.col(as.numeric(gamma_t), 10000)
  # # print(dim(rep.col(as.numeric(gamma_t), 10000)))
  # 
  # nrnd <- 30
  # rnd <- sample(10000, nrnd)
  # stackFits <- stack(as.data.frame(R0_fits[, rnd]))
  # stackFits <- transform(stackFits, age = rep(Time, length(rnd)))
  
  # return(list("R"=R_t, "R_lwr"=R_lwr, "R_upr" = R_upr))
  
  if (!basic_r){
    R_t <- R_t * (tot_pop - i[1:(length(i)-1)] - r[1:(length(i)-1)])/tot_pop
    R_lwr <- R_lwr * (tot_pop - i[1:(length(i)-1)] - r[1:(length(i)-1)])/tot_pop
    R_upr <- R_upr * (tot_pop - i[1:(length(i)-1)] - r[1:(length(i)-1)])/tot_pop
  }
  reproduction_number <- list("fitted" = R_t, "lower" = R_lwr, "upper" = R_upr)
  
  return(list("Reproductive_number"=reproduction_number, 
              "removal_rate" = removal_rate, 
              "transmission_rate" = transmission_rate, 
              "indices" = sub_indices[1:(length(sub_indices)-1)],
              "sim1" = sim1,"sim2" = sim2,"sim3" = sim3))
}

# --------------------------------
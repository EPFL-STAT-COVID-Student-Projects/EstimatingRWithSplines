# ------------ COVID functions and packages ------------
# Get data
library(openxlsx)
library(readr)
library(zoo)

# Stats
library(mgcv)
library(MASS)
library(splines)
require(stats); 

# Plots
library(ggplot2)
library(plotly)
library(cowplot)
library(gridExtra)
library(latex2exp)
require(graphics)
library(scales)

# Shiny packages
library(shiny)
library(shinyjs)
library(shinythemes)
library(shinyWidgets)
library(shinyBS)
library(shinyalert)
library(shinydashboard)
library(dashboardthemes)
library(shinydashboardPlus)
library(shinycssloaders)
library(waiter)

# not tested here yet:
library(rintrojs)
# library(egg)
library(shinydisconnect)
library(xtable)

# ============================ Deconvolution methods ============================
library(shinythemes)
library(splines)

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
# -------------------------------- 0) Auxiliary functions --------------------------------

# Function to compute the vector of probabilities from a (finite) discrete gamma distribution
#' @param n number of intervals 
#' @param normalize TRUE to normalize the probabilities
#' @param ... parameters for the gamma distriution (shape, scale, ...)
#' @return the discretized gamma distribution, with finite support of size n
discrete_gamma <- function(n = 20, normalize = TRUE, ...){
  discr_gamma <- 1:n
  for (i in 1:n){
    discr_gamma[i] <- pgamma(i, ...) - pgamma(i-1, ...)
  }
  if (normalize){
    discr_gamma <- discr_gamma/sum(discr_gamma)
  }
  return(discr_gamma)
}

# Function to compute the vector of probabilities from a (finite) discrete normal distribution
#' @param n number of intervals 
#' @param mean_ mean for the normal distribution
#' @param se for the normal distribution
#' @return the discretized normal distribution, with finite support of size n (j = 0, ..., n-1)
discrete_normal <- function(n, mean_ = 6.2, se_ = 1.72, check = FALSE, ...){
  sample_values <-round(rnorm(100000, mean_, se_), 0)
  sample_values <- sample_values[which(sample_values >= 1)]
  
  if (check){
    print(mean(sample_values))
    print(var(sample_values))
    plot(table(sample_values)/length(sample_values))
  }
  prob_vector <- as.numeric(table(sample_values)/length(sample_values))
  
  return(list("prob_vector" = as.numeric(table(sample_values)/length(sample_values)), "len" = length(prob_vector)))
}

# Function to convolve the data, based on the delay profile w
#' @param w_x delay profile (distribution, for example a discretized gamma)
#' @param data vector of data to be convolved
#' @return convolved data, of same length as the initial vector of data
convolveData <- function(w_x, data){
  data_extended <- c(rep(0, length(w_x)-1), data)
  data_convolved <- convolve(data_extended, rev(w_x), type = "filter")
  return(data_convolved)
}

# Compute the sum of two discretized gamma distribution with finite support
#' @param shape_1 shape for the first discretized gamma distribution
#' @param scale_1 scale for the first discretized gamma distribution
#' @param shape_2 shape for the second discretized gamma distribution
#' @param scale_2 scale for the second discretized gamma distribution
#' @param n_ half the size of the support of each of the discretized gamma
#' @param normalize TRUE to normalize the final distribution
sumDiscretizedGamma <- function(shape_1, scale_1, shape_2, scale_2, n_, normalize = TRUE){
  discr_gamma <- c(1:2*n_)
  d1 <- discrete_gamma(2*n_, normalize = TRUE, shape = shape_1, scale = scale_1)
  d2 <- discrete_gamma(2*n_, normalize = TRUE, shape = shape_2, scale = scale_2)
  
  for (k in 1:2*n_){
    discr_gamma[k] <- 0
    for (i in 0:(k-1)){
      discr_gamma[k] <- discr_gamma[k] + d1[k-i+1]*d2[i+1]
    }
  }
  
  if (normalize){
    discr_gamma <- discr_gamma/sum(discr_gamma, na.rm = TRUE)
  }
  
  return(discr_gamma)
}

# Compute the approximation of the sum of two gamma distribution, and then discretized it
#' @param shape_1 shape for the first discretized gamma distribution
#' @param scale_1 scale for the first discretized gamma distribution
#' @param shape_2 shape for the second discretized gamma distribution
#' @param scale_2 scale for the second discretized gamma distribution
#' @param n_ half the size of the support of each of the discretized gamma
#' @param normalize TRUE to normalize the final distribution
#' @references https://stats.stackexchange.com/questions/72479/generic-sum-of-gamma-random-variables
sumDiscretizedGammaApprox <- function(shape_1, scale_1, shape_2, scale_2, n_, normalize = TRUE){
  k_sum <- (shape_1*scale_1 + shape_2*scale_2)**2
  k_sum <- k_sum/((shape_1*scale_1**2 + shape_2*scale_2**2))
  theta_sum <- (shape_1*scale_1 + shape_2*scale_2)/k_sum
  return(discrete_gamma(n_, normalize = TRUE, shape = k_sum, scale = theta_sum))
}

# Construct the delta matrix that is used for the penalization for the deconvolution
#' @param n size of the square matrix
#' @return matrix of size size n x n, with 2s on the diagonal, and -1s on the upper and lower diagonals
construct_delta <- function(n){
  delta <- matrix(0L, nrow = n, ncol = n)
  diag(delta) <- 2
  delta[1,1] <- 1; delta[dim(delta)[1], dim(delta)[1]] <- 1
  sdiag(delta, 1) <- -1
  sdiag(delta, -1) <- -1
  return(delta)
}

# Function to predict the next values, based on a generalized additive model, with smooth trend for the days, and for the time
#' @param observed_data observed incidence cases
#' @param extension list of extensions parameters: ex: list('len' = 20, 'family' = quasipoisson(link = 'identity'), 'keep_original' = TRUE, 'min_weights' = 0.1),
#' where 'len' is number of steps in the future to be predicted, 'family' is the family for the model for prediction, 'keep_original' a boolean specifying 
#' if we keep the real data for deconvolution (or all prediction, even known values), and 'min_weights' a minimum value for weights of predicted values (otherwise could be 0 sometimes)
#' @return list with 'DATA' the dataset used for the fitting, 'extension_results' a list for the fitted values, and upper and lower CI, and 'weights' the weights associated with each values (for the whole dataset, included known data)
extension_prediction <- function(observed_data, extension = list('len' = 10, 'family' = nb(link = 'log'), 'keep_original' = TRUE, 'min_weights' = 0.1), fm = NULL, weights = NULL, ...){
  if (is.null(fm)){
    fm <- as.formula("cases ~ s(time, k = 50, bs = \'tp\', m = c(3, 0)) + s(day, k = 3)")
  }else{
    fm <- as.formula(fm)
  }
  indices <- c(1:length(observed_data))
  
  DAT_extrapol <- data.frame('cases' = observed_data, "time" = indices, 'day' = (indices%%7))
  DAT_extrapol_plus <- data.frame('cases' = c(observed_data, rep(NaN, extension$len)), 
                                  "time" = c(1:(length(indices)+extension$len)), 
                                  'day' = (c(1:(length(indices)+extension$len))%%7))
  
  extension_model <- gam(fm, method = "REML", data = DAT_extrapol, 
                         knots = list(time = c(-10, 0, indices[length(indices)]+1, indices[length(indices)]+extension$len+10)), 
                         family = extension$family, weights = weights)
  ilink <- family(extension_model)$linkinv
  pred_glm <- predict(extension_model, newdata = DAT_extrapol_plus, se.fit = TRUE)
  
  pred <- ilink(pred_glm$fit)
  lower <- ilink(pred_glm$fit - 2*pred_glm$se.fit) 
  upper <- ilink(pred_glm$fit + 2*pred_glm$se.fit)
  
  # observed_data <- extension_results$fitted
  weights_ <- 1/pred_glm$se.fit**2
  weights_ <- weights_/mean(weights_)
  
  if(extension$keep_original){
    pred[indices] <- observed_data
    lower[indices] <- observed_data
    upper[indices] <- observed_data
    
    weights_ <- weights_/min(max(weights_), 1.5)
    weights_[indices] <- 1
    if (!is.null(extension$min_weights)){
      weights_[which(weights_ < extension$min_weights)] <- extension$min_weights
    }
  }
  extension_results <- list("fitted" = pred, 
                            "upper" = upper, 
                            "lower" = lower)
  
  
  return(list("DATA" = DAT_extrapol_plus, "extension_results" = extension_results, "weights" = weights_, "fit" = extension_model))
}



# -------------------------------- 1) Roll back ---------------------------------------

#' Function to get the true infection events, by deconvolving the infection events
#' @param observed_data vector of observed events
#' @param prob_vector delay profile 
#' @return the deconvolved data
rollBack <- function(observed_data, prob_vector){
  n <- length(prob_vector)
  I <- 0*c(1:n, observed_data)
  for (i in length(observed_data):1){
    size <- round(max(0,observed_data[i]))
    if (is.na(size)){
      size <- 0
    }
    a <- sample(c(1:n), size = size, replace = TRUE, prob = prob_vector)
    b <- matrix(table(c(1:n, a)))
    for (j in 1:n){
      I[i+n-j] = I[i+n-j] +  b[j,1] - 1
    }
  }
  return(I)
}

#' Function to get the confidence intervals based on the roll back deconvolution method
#' @param observed_data vector of observed events
#' @param prob_vector delay profile 
#' @param n number of samples to compute the confidence intervals
#' @return matrix of the mean deconvolution, and the 0.025 and 0.975 quantiles
rollBackCI <- function(observed_data, prob_vector, n){
  roll_back_matrix = matrix(nrow = n, ncol = length(observed_data))
  for (j in 1:n){
    print(j)
    data_country_roll_back <- rollBack(observed_data, prob_vector)
    roll_back_matrix[j,] <- data_country_roll_back[(length(prob_vector)+1):length(data_country_roll_back)]
  }
  q_0.025 <- apply(roll_back_matrix , 2 , quantile , probs = 0.025 , na.rm = TRUE)
  q_0.975 <- apply(roll_back_matrix , 2 , quantile , probs = 0.975 , na.rm = TRUE)
  mean_roll_back <- apply(roll_back_matrix , 2 , quantile , probs = 0.5 , na.rm = TRUE)
  return(list("matrix" = roll_back_matrix, "q2.5" = q_0.025, "q97.5" = q_0.975, "mean" = mean_roll_back))
}




# -------------------------------- 2) Optimization ------------------------------------

# Function to compute the log likelihood for the deconvolution with poisson response
#' @param x true infection events, to be estimated
#' @param y observed incidence cases
#' @param W delay profile
#' @param exp_ TRUE use the convention x = exp(lambda), and lambda are now the parameter to be estimated
#' @param eps regularizartion parameter
#' @param use_delta TRUE to use the delta matrix as regulariztion, FALSE otherwise
log_lik <- function(x, y, W, exp_ = FALSE, eps = 0, use_delta = FALSE){
  W <- c(rev(W), 0)
  if (use_delta){
    delta <- construct_delta(length(y))
  }
  if (exp_){
    r <- convolve(c(1+numeric(length(W)-1), x), W, type = "filter")   # 0* 1:(length(W)-1)
    if (eps > 0){
      if (use_delta){
        term_control <- eps * t(x) %*% delta %*% x
      }else{
        term_control <- eps*sum(x**2)/length(x)
      }
    }
  }else{
    r <- convolve(c(1+numeric(length(W)-1), exp(x)), W, type = "filter")   # 0* 1:(length(W)-1)
    if (eps > 0){
      if (use_delta){
        term_control <- eps * t(exp(x)) %*% delta %*% exp(x)
      }else{
        term_control <- eps*sum(exp(x)**2)/length(x)
      }
    }
  }
  f <- - sum(y * log(r) - r)
  if (eps > 0){
    f <- f + term_control
  }
  return(f)
}

# Function to compute the derivative for the log likelihood (without regularization term)
#' @param x true infection events, to be estimated
#' @param y observed incidence cases
#' @param W delay profile
#' @note here x = exp(lambda), and we evaluate the derivative based on x (not on lambda !) 
log_lik_derivative <- function(x, y, W){
  W_new <- c(rev(W), 0)
  # W_rev = rev(W)
  r <- convolve(c(0* 1:(length(W_new)-1), x), W_new, type = "filter")
  derivative <- numeric(length(x))
  for (i in seq_along(derivative)){
    for (j in 1:length(W)){
      if (i+j <= length(y)){
        derivative[i] = derivative[i] + (y[j+i]/r[j+i] - 1)*W[j]
      }
    }
  }
  return(derivative)
}

# Function to compute the hessian for the log likelihood (without regularization term)
#' @param x true infection events, to be estimated
#' @param y observed incidence cases
#' @param W delay profile
#' @note here x = exp(lambda), and we evaluate the derivative based on x (not on lambda !) 
log_lik_hessian <- function(x, y, W){
  W_new <- c(rev(W), 0)
  # W_rev = rev(W)
  r <- convolve(c(0* 1:(length(W_new)-1), x), W_new, type = "filter")
  size_W <- length(W)
  
  hessian <- matrix(0, length(x), length(x))
  for (i in 1:length(x)){
    for (j in 1:i){
      if (i - j < size_W){
        for (k in 1:(size_W - (i-j))){
          if (i+k <= length(x)){
            hessian[j,i] <- hessian[j,i] - y[i+k]*W[k]*W[i-j+k]/(r[i+k]**2)
          }
        }
      }
    }
  }
  hessian <- hessian + t(hessian) - diag(diag(hessian))
  return(hessian)
}

# Function to compute the deconvolution based on a poisson response
#' @param observed_data vector of observed incidence cases
#' @param prob_vector delay profile
#' @param smooth_observations integer to smooth the obserations with a (centered) moving average of 'smooth_observations' days
#' @param maxIter maximum number of iterations for the optimization algorithm
#' @param method method used for the optimization: CG, BFGS, SANN, etc
#' @param data_start start parameters for the likelihood 
#' @param hessian TRUE to estimate the hessian matrix
#' @return a list of the fitted values and the object returned by the optimization method
optimPoissonDeconvolution <- function(observed_data, prob_vector, smooth_observations = 1, maxIter = 10, method = "CG", data_start = NULL, hessian = TRUE){
  if (smooth_observations > 1){
    # Smooth the observations (for computing the 'true' infection events by optimization)
    observed_data <- rollapply(observed_data, width=smooth_observations, FUN=function(x) mean(x, na.rm=TRUE), by=1, 
                               by.column=TRUE, partial=TRUE, fill=NA, align="center")
  }
  if (is.null(data_start)){
    data_start <- log(observed_data)
  }
  control <- list()
  control[['maxit']] <- maxIter
  optim_results <- optim(data_start, log_lik, 
                         y = observed_data, W = prob_vector, 
                         method = method, control = control, hessian = hessian)
  values_opt <- exp(optim_results$par)
  
  return(list("fitted_values" = values_opt, "result" = optim_results))
}


# -------------------------------- 3) IWLS --------------------------

#' Function to perform the deconvolution with the IWLS algorithm based on a poisson response, and splines regression
#' @param y observed incidence cases
#' @param beta_dim number of parameters
#' @param Beta_Beta_tilde specification matrix, to be multiplied with beta
#' @param T_len length of the data
#' @param lambda parameter for the regularization
#' @param beta_0 initial values for the parameters (if NULL, set to 0s, with length T_len)
#' @param tol tolerance, set to 1e-5
#' @param max_iter maximum number of iterations
#' @param delta regularization matrix
#' @return list of the estimates and the covariance matrix for the estimates
IWLS_algorithm <- function(y, beta_dim, Beta_Beta_tilde, T_len, lambda, beta_0 = NULL, tol = 1e-5, max_iter = 10000, delta = NULL){
  if (is.null(beta_0)){
    beta <- as.matrix(numeric(beta_dim))
  }else{
    beta <- beta_0
  }
  checkCondition <- TRUE
  iter = 1
  while(checkCondition){
    fitted <- Beta_Beta_tilde %*% beta  
    fitted[which(fitted <= 0)] <- 1
    W <- fitted  
    if (is.null(delta)){
      delta <- diag(lambda*(1+numeric(beta_dim)))     
    }else{
      delta <- lambda * delta
    }
    z <- fitted + (y - fitted) * 1/fitted 
    inter_B <- as.numeric(W) * Beta_Beta_tilde
    to_be_inv <- t(Beta_Beta_tilde) %*% inter_B + delta    
    beta_new <- solve(to_be_inv, t(inter_B) %*% z)
    # print(max(abs(beta_new - beta)))
    if (max(abs(beta_new - beta)) < tol){
      checkCondition <- FALSE
    }
    beta <- beta_new
    iter <- iter + 1
    if (iter > max_iter){
      checkCondition <- FALSE
    }
    if (iter %% 1000 == 0){
      print(iter/max_iter)
    }
  }
  cov_matrix <- solve(to_be_inv)
  # print(beta)
  # print(sqrt(diag(cov_matrix)))
  
  return(list("estimates"=beta, "cov"=cov_matrix))
}


# Function to construct the matrix of delay profiles
#' @param n dimension of the matrix
#' @param W_ delay profile
#' @return the matrix of delay profiles, with size n x n
construct_B <- function(n, W_){
  B <- matrix(0L, nrow = n, ncol = n)
  W_rev <- rev(W_)
  
  for (j in 1:n){
    if (j < length(W_)){
      B[j, 1:j] <- rev(W_[1:j])
    }else{
      B[j, (j-length(W_)+1):(j)] <- W_rev
    }
  }
  return(B)
}


#' Function to deconvolve the data with the IWLSE algorithm, using the gam functions
#' @param observed_data observed incidence cases
#' @param prob_vector delay profile
#' @param dim_spline_basis dimension for the spline basis
#' @param delta regularization matrix (/!\ of size (dim_spline_basis+4)x(dim_spline_basis+4))
#' @param weights_ weights for the IWLS algorithm
#' @param family family for the generalized additive model (with link function)
#' @param extension list of extensions parameters: ex: list('len' = 20, 'family' = quasipoisson(link = 'identity'), 'keep_original' = TRUE, 'min_weights' = 0.1),
#' where 'len' is number of steps in the future to be predicted, 'family' is the family for the model for prediction, 'keep_original' a boolean specifying 
#' if we keep the real data for deconvolution (or all prediction, even known values), and 'min_weights' a minimum value for weights of predicted values (otherwise could be 0 sometimes)
#' @param fitted_convolved TRUE to get the convolve fitted values, FALSE to get the deconvolved ones
#' @param get_derivatives TRUE to plot the first and second derivatives, FALSE otherwise
#'@section ------------------ PLOT OPTIONS ------------------
#' @param plot_ TRUE to plot the deconvolved data
DeconvolutionGAM <- function(observed_data, prob_vector, dim_spline_basis = 10, delta = NULL, weights_ = NULL, family = poisson(link = 'identity'),
                             extension = NULL, fitted_convolved = FALSE, get_derivatives = FALSE, plot_ = FALSE, weights_extension = NULL, ...){
  # Some checks
  if (is.null(dim_spline_basis) && is.null(delta)){stop("A dimension basis must be provided !")}
  if (!is.null(dim_spline_basis) && !is.null(delta) && dim_spline_basis != (dim(delta)[1]-4)){stop("Dimension basis do not match !")}
  
  indices <- c(1:length(observed_data))
  
  # To deal with drop effects near the end, we can add some extension
  # the uncertainties for the prediction are taken into account for the deconvolution
  if (!is.null(extension) && extension$len > 0){
    ext <- extension_prediction(observed_data, extension = extension, weights = weights_extension, ...)
    extension_results <- ext$extension_results
    weights_ <- ext$weights
    DAT_extrapol_plus <- ext$DATA
    observed_data <- extension_results$fitted
  }else{
    DAT_extrapol_plus <- data.frame('cases' = observed_data, "time" = indices, 'day' = (indices%%7))
    extension_results <- NULL
  }
  
  # Construct the spline basis
  Time <- c(1:length(observed_data))
  T_len <- length(Time)
  knots1 <- seq(from=min(Time),to=max(Time),length=dim_spline_basis+2)[2:(dim_spline_basis+1)]
  SplineBasis <- bs(Time, knots  = knots1, degree = 3) #, Boundary.knots = c(-20, T_len + 50))
  
  B <- construct_B(T_len, prob_vector)
  Beta_Beta_tilde <- B %*% SplineBasis
  
  
  # -- TRAINING 
  DAT <- as.data.frame(Beta_Beta_tilde)
  xnam <- paste("x", 1:dim(Beta_Beta_tilde)[2], sep="")
  colnames(DAT) <- xnam
  DAT$resp <- abs(observed_data)
  
  fmla <- as.formula(paste("resp ~ ", paste(xnam, collapse= "+")))
  
  fit <- gam(fmla, family = family, data = DAT, H = delta, weights = weights_)
  ilink <- family(fit)$linkinv
  V <- vcov(fit)
  
  
  if (fitted_convolved){
    pred_matrix <- predict(fit,newdata=DAT,type="lpmatrix")
    linear_pred <- pred_matrix %*% as.matrix(fit$coefficients) 
    se <- sqrt(diag(pred_matrix %*% V %*% t(pred_matrix)))
    q_0.025 <- ilink(linear_pred - 2*se)
    q_0.975 <- ilink(linear_pred + 2*se)
    BASIS <- pred_matrix
  }else{
    BASIS <- cbind(1+numeric(length(observed_data)), SplineBasis)
    linear_pred <- BASIS %*% as.matrix(fit$coefficients)    # get the deconvolved ones
    se <- sqrt(diag(BASIS %*% V %*% t(BASIS)))
    q_0.025 <- ilink(linear_pred - 2*se)
    q_0.975 <- ilink(linear_pred + 2*se)
  }
  fitted_values <- ilink(linear_pred)
  
  # function used for computing the derivatives (transform the basis, here just add column of 1s)
  pred_mat_fct <- function(basis){
    return(cbind(1+numeric(length(observed_data)), basis))
  }
  
  if (get_derivatives){
    res_derivatives <- FirstAndSecondDerivatives(as.matrix(fit$coefficients), 1e-5,
                                                 length(observed_data), 5*length(observed_data), knots1, 
                                                 V, 
                                                 plot_ = TRUE, sample_splines = TRUE, fct_pred_matrix = pred_mat_fct)
    first_deriv <- res_derivatives$first_deriv
    second_deriv <- res_derivatives$second_deriv
    curvature <- res_derivatives$curvature
  }else{
    first_deriv <- NULL
    second_deriv <- NULL
    curvature <- NULL
  }
  
  # -------- PLOT
  if (plot_){
    plot_final <- ggplot(data = NULL, aes(x = DAT_extrapol_plus$time)) + 
      geom_line(aes(y = DAT_extrapol_plus$cases, colour = "cases"), size = 0.7) + 
      
      geom_ribbon(aes(ymin = q_0.025, 
                      ymax = q_0.975, 
                      fill = "deconvolution", colour = NA), alpha = 0.3, fill = colors[2]) +
      geom_line(aes(y = fitted_values, colour = "deconvolution"), size = 0.7)
    
    if (!is.null(extension)){
      plot_final <- plot_final + 
        geom_vline(xintercept = length(indices), color = "red") + 
        geom_ribbon(aes(ymin = extension_results$lower, 
                        ymax = extension_results$upper, 
                        fill = "Confidence interval extension", colour = NA), alpha = 0.3, fill = colors[3]) +
        geom_line(aes(y = extension_results$fitted, colour = "predicted cases"), size = 0.7) + 
        geom_segment(aes(x = length(indices), y = 0, xend = length(DAT_extrapol_plus$time), yend = 0), 
                     color = "#e63946", arrow = arrow(length = unit(0.2, "cm"))) + 
        geom_text(data=NULL, aes(x= length(indices) + (length(DAT_extrapol_plus$time) - length(indices))/2, y = 0, 
                                 label = "Prediction"), colour = "#e63946", vjust = 1.5, size = 3)
    }
    
    plot_final <- plot_final + 
      scale_fill_manual("", 
                        breaks = c("Confidence interval extension", "deconvolution", "Simlutaneous confidence bands"),
                        values = c("#e63946", colors[2], colors[2])) +
      scale_colour_manual("", 
                          breaks = c("T", "S", "E", "I", "R", 'observed incidence', "cases", "predicted cases", "deconvolution"),
                          values = c("black", "gray", colors[1], colors[2], colors[3], colors[2], "#1d3557", "#e63946",  colors[2])) +
      
      labs(title="Cases", x="time", y="number cases") + 
      theme_for_the_plots + 
      guides(color=guide_legend(override.aes=list(fill=NA))) + theme(legend.position = 'bottom') + 
      coord_cartesian(ylim = c(0, 1.2*max(DAT_extrapol_plus$cases, na.rm = TRUE)))
    
    plot(plot_final)
  }else{
    plot_final <- NULL
  }
  
  return(list('fitted_values' = fitted_values, "point_wise_lower" = q_0.025, "point_wise_upper" = q_0.975, 
              "SplineBasis" = BASIS, "fit" = fit, "cov" = V, "extension_results" = extension_results, 
              "plot" = plot_final, "DATA_extension" = DAT_extrapol_plus, "se" = se))
}

#' Function that compute time intervals for which the confidence bands do not include the threshold value
#' @param lower_band vector for the lower band
#' @param upper_band vector for the upper band
#' @param threshold threshold value
#' @return list containing the vector of indices for the start and the end of each time interval
compareWithThreshold <- function(lower_band, upper_band, threshold){
  start_periods <- end_periods <- c()
  first_ind <- 1
  search_new_zone <- TRUE
  for (j in c(1:length(lower_band))){
    if (search_new_zone){
      if (lower_band[j] < threshold && upper_band[j] > threshold){
        search_new_zone <- FALSE
        start_periods <- c(start_periods, j)
      }
    }else{
      if (lower_band[j] > threshold || upper_band[j] < threshold){
        search_new_zone <- TRUE
        end_periods <- c(end_periods, j)
      }
      if (j == length(lower_band)){
        search_new_zone <- TRUE
        end_periods <- c(end_periods, j)
      }
    }
  }
  return(list("start" = start_periods, "end" = end_periods))
}


#' Function to get the point-wise and simultaneous confidence bands from the gam object
#' @param SplineBasis basis fonction evaluated on the data
#' @param CovMatrix covariance matrix of the estimated parameters from the fitted model
#' @param N number of samples for computing the criterion for the simultaneous CI
#' @section ------------------ PLOT OPTIONS ------------------
#' @param plot_ TRUE to return the plot of the confidence bands
#' @param intercept_plot intercept for the horizontal line (set to NULL if you do not want to plot it)
#' @param sample_splines TRUE to also plot the splines generated by sampling N times the coefficients (from a multivariate normal distribution, with covariance CovMatrix)
#' @param axis_x the values in the x axis (NULL if not specified)
#' @param N_keep number of samples plotted
#' @references For more informations, see Ruppert et al. (2003), and the blog post https://fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/
SimultaneousIntervals <- function(SplineBasis, coefficients_model, CovMatrix, N = 10000 ,plot_ = FALSE, intercept_plot = NULL, 
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


#' Function to compute the first and second derivatives for the fitted values of the model, using finite difference
#' @param coefficients_model coefficients of the model
#' @param eps step for the finite difference
#' @param max_time maximum time (for the x's)
#' @param length_x number of timesteps in the model (usually equivalent to max_time)
#' @param knots number of knots for the splines
#' @param CovMatrix covariance matrix of the estimated parameters from the fitted model
#' @param CI TRUE to get the simultaneous confidence intervals
#' @param curvature_score TRUE to get the curvature score 
#' @param fct_pred_matrix function provided to transform the basis (for example, adding column of 1 at the beginning)
#' @section ------------------ PLOT OPTIONS ------------------
#' @param plot_ TRUE to plot the results
FirstAndSecondDerivatives <- function(coefficients_model, eps, max_time, length_x, knots, CovMatrix, plot_ = TRUE, CI = TRUE, curvature_score = FALSE, fct_pred_matrix = NULL, ...){
  
  # Finite difference matrices:
  # 1) compute time intervals
  Time <- seq(1, max_time, length.out = length_x)
  Time_plus_eps <- Time + eps
  Time_minus_eps <- Time - eps
  # 2) compute spline basis for each time interval
  SplineBasis <- bs(Time ,knots = knots, degree = 3)#, Boundary.knots = c(1, max_time))
  SplineBasis_plus_eps <- bs(Time_plus_eps ,knots = knots, degree = 3)#, Boundary.knots = c(1, max_time))
  SplineBasis_minus_eps <- bs(Time_minus_eps ,knots = knots, degree = 3)#, Boundary.knots = c(1, max_time))
  
  # if fct_pred_matrix provided, apply it (ex: add column of 1s)
  if (!is.null(fct_pred_matrix)){
    SplineBasis <- fct_pred_matrix(SplineBasis)
    SplineBasis_plus_eps <- fct_pred_matrix(SplineBasis_plus_eps)
    SplineBasis_minus_eps <- fct_pred_matrix(SplineBasis_minus_eps)
  }
  result_gam <- SplineBasis %*% coefficients_model  # fitted values
  SplineBasis_difference <- (SplineBasis_plus_eps - SplineBasis) / eps   # first derivative spline basis
  SplineBasis_difference_2 <- (SplineBasis_plus_eps + SplineBasis_minus_eps - 2*SplineBasis)  / eps^2  # second derivative spline basis
  
  derivs <- SplineBasis_difference %*% coefficients_model  # fitted values first derivative
  derivs_second <- SplineBasis_difference_2 %*% coefficients_model  # fitted values second derivative
  
  # Compute the confidence bands
  if (CI){
    simultaneous1 <- SimultaneousIntervals(SplineBasis, coefficients_model, CovMatrix, plot_ = plot_, 
                                           title = "Fitted values", x = "time", y = "fitted values", x_axis = Time,...)
    simultaneous2 <- SimultaneousIntervals(SplineBasis_difference, coefficients_model, CovMatrix, plot_ = plot_, intercept_plot = 0, 
                                           title = "First derivative", x = "time", y = "first derivative", x_axis = Time, ...)
    simultaneous3 <- SimultaneousIntervals(SplineBasis_difference_2, coefficients_model, CovMatrix, plot_ = plot_, intercept_plot = 0, 
                                           title = "Second derivative", x = "time", y = "second derivative", x_axis = Time, ...)
    
    if (plot_){
      plot(plot_grid(simultaneous1$plot, simultaneous2$plot, simultaneous3$plot, nrow = 3, align = 'hv'))
    }
  }
  # curvature score (not used)
  if (curvature_score){
    ind1 <- c(2:length(Time))
    ind2 <- c(1:(length(Time)-1))
    diff <- Time[ind1] - Time[ind2]
    score <- sum(diff*(derivs_second[ind2]**2 + derivs_second[ind1]*derivs_second[ind2] + derivs_second[ind1]**2))/6
  }else{
    score <- NULL
  }
  
  return(list("fitted" = result_gam, "first_deriv" = derivs, "second_deriv" = derivs_second, "curvature" = score))
}


#' Find best lambda for regularization, based on results obtained with gam model
#' @param observed_data observed incidence curve
#' @param k_spline dimension of spline basis
#' @param prob_vector delay distribution
#' @param range_values range of values for lambda
#' @param which_regularization type of matrix for regularization (tri diagonal or diagonal)
#' @param smoothing time window for smoothing the convolved curve that is used for choosing the optimal parameter. If NULL, no smoothing is applied
#' @param diagnostic_plot TRUE to show the errors as function of the lambda
#' @param ... parameters to be passed to the function 'DeconvolutionGAM'
findBestLambda <- function(observed_data, k_spline, prob_vector, range_values = 10**(-1:-10), which_regularization = 'tridiag',
                           smoothing = NULL, diagnostic_plot = FALSE, ...){
  indices_init <- c(1:length(observed_data))
  # Apply a smoothing window on the observed data
  if (!is.null(smoothing)){
    data_inc_after <- rollapply(observed_data, width=smoothing, FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                by.column=TRUE, partial=TRUE, fill=NA, align="center")
  }
  errors <- numeric(length(range_values))
  B <- construct_B(length(observed_data), prob_vector)
  
  # iterate over the range of values for lambda
  pb <- txtProgressBar(min = 0, max = length(range_values), style = 3)
  for (j in c(1:length(range_values))){
    # regularization matrix
    if (which_regularization == 'tridiag'){
      delta <- construct_delta(k_spline + 4)
    }else{
      delta <- diag(1+numeric(k_spline + 4))
    }
    delta <- range_values[j] * delta
    
    # Deconvolution
    results_deconv <- DeconvolutionGAM(observed_data, prob_vector, k_spline, delta = delta, 
                                       plot_ = FALSE, ...)
    fitted <- results_deconv$fitted
    fitted[which(is.na(fitted))] <- 0
    # compute MSE between the convolved deconvolved data and the (smoothed) observed data
    errors[j] <- mean((B %*% fitted[indices_init] - data_inc_after[indices_init])**2, na.rm = TRUE)
    setTxtProgressBar(pb, j)
  }
  close(pb)
  
  # plot of the error as function of the lambda
  if (diagnostic_plot){
    plot_ <- ggplot(data = NULL, aes(x = range_values, y = errors)) + geom_point() + 
      geom_vline(xintercept = range_values[which(errors == min(errors))], color = "red") + 
      labs(title="Error as function of lambda", x="lambda", y="MSE") + 
      theme_for_the_plots + guides(color=guide_legend(override.aes=list(fill=NA))) + scale_x_log10()
    
    plot(plot_)
  }
  
  return(range_values[which(errors == min(errors))])
}

#' Find best dimension basis, based on results obtained with gam model
#' @param observed_data observed incidence curve
#' @param prob_vector delay distribution
#' @param range_values range of values for lambda
#' @param smoothing time window for smoothing the convolved curve that is used for choosing the optimal parameter. If NULL, no smoothing is applied
#' @param diagnostic_plot TRUE to show the errors as function of the dimension basis
#' @param ... parameters to be passed to the function 'DeconvolutionGAM'
findBestDimension <- function(observed_data, prob_vector, range_values = c(10:200), smoothing = NULL, diagnostic_plot = FALSE, ...){
  indices_init <- c(1:length(observed_data))
  # Apply a smoothing window on the observed data
  if (!is.null(smoothing)){
    data_inc_after <- rollapply(observed_data, width=smoothing, FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                by.column=TRUE, partial=TRUE, fill=NA, align="center")
  }
  errors <- numeric(length(range_values))
  B <- construct_B(length(observed_data), prob_vector)
  
  # iterate over the range of values for the basis dimension
  pb <- txtProgressBar(min = 0, max = length(range_values), style = 3)
  for (j in c(1:length(range_values))){
    # Deconvolution
    results_deconv <- DeconvolutionGAM(observed_data, prob_vector, range_values[j], delta = NULL, 
                                       plot_ = FALSE, ...)
    fitted <- results_deconv$fitted
    fitted[which(is.na(fitted))] <- 0
    # compute MSE between the convolved deconvolved data and the (smoothed) observed data
    errors[j] <- mean((B %*% fitted[indices_init] - data_inc_after[indices_init])**2, na.rm = TRUE)
    setTxtProgressBar(pb, j)
  }
  close(pb)
  
  # plot of the error as function of the basis dimension
  if (diagnostic_plot){
    plot_ <- ggplot(data = NULL, aes(x = range_values, y = errors)) + geom_point() + 
      geom_vline(xintercept = range_values[which(errors == min(errors))], color = "red") + 
      labs(title="Error as function of the dimension basis", x="dimension basis", y="MSE") + 
      theme_for_the_plots + guides(color=guide_legend(override.aes=list(fill=NA)))
    plot(plot_)
  }
  
  return(range_values[which(errors == min(errors))])
}



# ============================ Estimating R ==================
# -------------------------------- 1) GAM model --------------------------------
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
  return(list("Point_wise_CI_upper" = q_0.975, "Point_wise_CI_lower" = q_0.025, 
              "Simultaneous_upper"=q_0.975_crit, "Simultaneous_lower"=q_0.025_crit, "fits" = stackFits, "plot" = p, 
              "fitted" = result_gam))
}

# -------------------------------- 2) SIR model --------------------------------

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


# ============================ Construct dataset =======================

# Function to get the data for the corresponding country, including the daily number of infected, recovered and deaths
#' @param country country for which the data need to be collected
#' @param option 1 to get the data from bag.admin.ch (only for Switzerland), 2 for another source
#' @param plot_ TRUE to plot the collected data
ConstructDataset <- function(country, option, plot_ = F, DATA_already_provided = NULL){
  
  if (option == "OFSP"){
    option = 1
  }else if(option == "Our World in Data"){
    option = 2
  }else if (option == "HDX"){
    option = 3
  }else{
    option = 3
  }
  
  # -------------- Some functions --------------
  # Function to get the recovered cases for the corresponding country
  #' @param country country for which the data need to be collected
  getRecovered <- function(country, gotData = NULL){
    if (country == "United States"){
      country <- "US"
    }
    if (is.null(gotData)){
      data_recovered <- as.data.frame(read_csv(file = "https://data.humdata.org/hxlproxy/api/data-preview.csv?url=https%3A%2F%2Fraw.githubusercontent.com%2FCSSEGISandData%2FCOVID-19%2Fmaster%2Fcsse_covid_19_data%2Fcsse_covid_19_time_series%2Ftime_series_covid19_recovered_global.csv&filename=time_series_covid19_recovered_global.csv"))
    }else{
      data_recovered <- gotData
    }
    data_recovered_matrix <- data.matrix(data_recovered[5:length(data_recovered[1,])])
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
  BAGData <- function(gotData = NULL){
    if (is.null(gotData)){
      data <- read.xlsx(sep=",",startRow = 8, detectDates = TRUE,
                        "https://www.bag.admin.ch/dam/bag/de/dokumente/mt/k-und-i/aktuelle-ausbrueche-pandemien/2019-nCoV/covid-19-datengrundlage-lagebericht.xlsx.download.xlsx/200325_Datengrundlage_Grafiken_COVID-19-Bericht.xlsx")
      names(data) <- c("date","cases","casesCumul","hospitalized","hospitalizedCumul",
                       "deaths","deathsCumul")
    }else{
      data <- gotData
    }
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
    
    # d1 <- which(data_country$confirmed!=0)
    # d2 <- which(data_country$deaths!=0)
    # d3 <- which(data_country$recovered!=0)
    # date_first_not_0 <- intersect(d1, d2)[1]
    # return(list(data_country, date_first_not_0))
    return(data_country)
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
    geom_text(data=NULL, aes(x= times[length(indices_)] + (times[length(extension_results$DATA$time)] - times[length(indices_)])/2, y = 0, 
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


add_class <- function(x, class) {
  x$attribs <- append(x$attribs, list(class = class))
  x
}


# ------------ theme plots -------------
unique_countries <- c("Switzerland", "Germany", "Italy", "France", "China", "United States", "United Kingdom", "Belgium", "Spain", "Russia")
colors <- c("#00AFBB", "#E7B800", "#FC4E07")
dark_theme <- FALSE
if (dark_theme){
  font_color_plot <- '#1E1D1E'
  points_color <- "grey85"
  theme_for_the_plots <- theme(plot.background = element_rect(fill = font_color_plot, colour = "black"),
                               panel.background = element_rect(fill = "#1E1D1E", colour = "black"),
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



# ------------ Theme ------------
getTheme2 <- function(){
  
  # color gradient: top: #C73B59
  #                 middle: #D56252
  #                 bottom: #E59A51
  
  # white grey: #FBFBFB
  # white pure: #FFFFFF
  # black grey: #363636
  theme_purple_gradient <- shinyDashboardThemeDIY(
    
    ### general
    appFontFamily = "Courier New"
    ,appFontColor = "#363636"   
    ,primaryFontColor = "#363636"   # font color box title
    ,infoFontColor = "rgb(255,255,255)"
    ,successFontColor = "rgb(255,255,255)"
    ,warningFontColor = "rgb(255,255,255)"
    ,dangerFontColor = "rgb(255,255,255)"
    
    ,bodyBackColor = "#FBFBFB"    # background body
    
    # ,bodyBackColor = cssGradientThreeColors(
    #   direction = "down"
    #   ,colorStart = "rgb(49,56,107)"
    #   ,colorMiddle = "rgb(71,59,109)"
    #   ,colorEnd = "rgb(78,88,149)"
    #   ,colorStartPos = 0
    #   ,colorMiddlePos = 70
    #   ,colorEndPos = 100
    # )
    
    ### header
    ,logoBackColor = "#C73B59"  # coin haut gauche
    
    ,headerButtonBackColor = "#FBFBFB"
    ,headerButtonIconColor = "#FBFBFB"
    ,headerButtonBackColorHover = "#FBFBFB"
    ,headerButtonIconColorHover = "#FBFBFB"
    
    ,headerBackColor = "#FBFBFB"    # header background color
    ,headerBoxShadowColor = "#363636"
    ,headerBoxShadowSize = "0px 0px 0px"
    
    ### sidebar
    ,sidebarBackColor = "#C73B59"
    # ,sidebarBackColor = cssGradientThreeColors(
    #   direction = "down"
    #   ,colorStart = "#C73B59"
    #   ,colorMiddle = "#D56252"
    #   ,colorEnd = "#E59A51"
    #   ,colorStartPos = 0
    #   ,colorMiddlePos = 50
    #   ,colorEndPos = 100
    # )
    # color gradient: top: #C73B59
    #                 middle: #D56252
    #                 bottom: #E59A51
    
    ,sidebarShadowRadius = ""
    ,sidebarPadding = 10
    ,sidebarShadowColor = "0px 0px 0px"
    
    ,sidebarMenuBackColor = "rgba(235, 118, 101, 0)"
    # ,sidebarMenuBackColor = cssGradientThreeColors(
    #   direction = "right"
    #   ,colorStart = "rgb(48,103,157)"
    #   ,colorMiddle = "rgb(65,79,129)"
    #   ,colorEnd = "rgb(55,70,120)"
    #   ,colorStartPos = 0
    #   ,colorMiddlePos = 30
    #   ,colorEndPos = 100
    # )
    ,sidebarMenuPadding = 10
    ,sidebarMenuBorderRadius = 20
    
    ,sidebarUserTextColor = "rgb(128,177,221)"
    
    ,sidebarSearchBackColor = "rgb(40,70,115)"
    ,sidebarSearchIconColor = "rgb(50,115,145)"
    ,sidebarSearchBorderColor = "rgb(30,60,105)"
    
    ,sidebarTabTextColor = "rgb(255, 255, 255)"
    ,sidebarTabTextSize = 13
    ,sidebarTabBorderStyle = "none none solid none"
    ,sidebarTabBorderColor = "none"
    ,sidebarTabBorderWidth = 0
    
    # ,sidebarTabBackColorSelected = cssGradientThreeColors(
    #   direction = "right"
    #   ,colorStart = "rgb(56,137,189)"
    #   ,colorMiddle = "rgb(65,95,145)"
    #   ,colorEnd = "rgb(68,84,137)"
    #   ,colorStartPos = 0
    #   ,colorMiddlePos = 50
    #   ,colorEndPos = 100
    # )
    ,sidebarTabBackColorSelected = "rgb(237, 163, 153)"
    
    ,sidebarTabTextColorSelected = "rgb(255,255,255)"
    ,sidebarTabRadiusSelected = "30px"
    
    # ,sidebarTabBackColorHover = cssGradientThreeColors(
    #   direction = "right"
    #   ,colorStart = "rgb(56,137,189)"
    #   ,colorMiddle = "rgb(65,95,145)"
    #   ,colorEnd = "rgb(68,84,137)"
    #   ,colorStartPos = 0
    #   ,colorMiddlePos = 50
    #   ,colorEndPos = 100
    # )
    ,sidebarTabBackColorHover = "rgb(237, 163, 0)"  # "rgb(237, 163, 153)"
    
    ,sidebarTabTextColorHover = "rgb(255,255,255)"
    ,sidebarTabBorderStyleHover = "none"
    ,sidebarTabBorderColorHover = "none"
    ,sidebarTabBorderWidthHover = 0
    ,sidebarTabRadiusHover = "30px"
    
    ### boxes
    # ,boxBackColor = cssGradientThreeColors(
    #   direction = "right"
    #   ,colorStart = "rgb(70,75,125)"
    #   ,colorMiddle = "rgb(65,79,129)"
    #   ,colorEnd = "rgb(55,70,120)"
    #   ,colorStartPos = 0
    #   ,colorMiddlePos = 30
    #   ,colorEndPos = 100
    # )
    ,boxBackColor = "rgb(255,255,255)"
    
    ,boxBorderRadius = 15
    
    ,boxShadowSize = "0px 5px 10px"
    ,boxShadowColor = "#C16154"
    
    ,boxTitleSize = 16
    ,boxDefaultColor = "rgb(49,56,107)"
    ,boxPrimaryColor = "rgb(237, 163, 153)"  # change header box
    ,boxInfoColor = "rgb(20,100,160)"
    ,boxSuccessColor = "rgb(64,186,170)"
    ,boxWarningColor = "rgb(255,217,144)"
    ,boxDangerColor = "rgb(249,144,144)"
    
    ,tabBoxTabColor = "rgb(237, 163, 153)"
    ,tabBoxTabTextSize = 14
    ,tabBoxTabTextColor = "rgb(255,177,221)"
    
    ,tabBoxTabTextColorSelected = "rgb(0,255,255)"
    
    # ,tabBoxBackColor = cssGradientThreeColors(
    #   direction = "right"
    #   ,colorStart = "rgb(70,75,125)"
    #   ,colorMiddle = "rgb(65,79,129)"
    #   ,colorEnd = "rgb(55,70,120)"
    #   ,colorStartPos = 0
    #   ,colorMiddlePos = 30
    #   ,colorEndPos = 100
    # )
    ,tabBoxBackColor = "rgb(255,0,0)"
    ,tabBoxHighlightColor = "rgb(80,95,155)"
    ,tabBoxBorderRadius = 15
    
    ### inputs
    ,buttonBackColor = "rgb(72,190,229)"
    ,buttonTextColor = "rgb(40,63,106)"
    ,buttonBorderColor = "rgb(72,190,229)"
    ,buttonBorderRadius = 20
    
    ,buttonBackColorHover = "rgb(115,210,240)"
    ,buttonTextColorHover = "rgb(255,255,255)"
    ,buttonBorderColorHover = "rgb(115,210,240)"
    
    ,textboxBackColor = "rgb(237, 163, 153)"
    ,textboxBorderColor = "rgb(255,255,255)"
    ,textboxBorderRadius = 20
    ,textboxBackColorSelect = "rgb(237, 163, 153)"
    ,textboxBorderColorSelect = "rgb(237, 163, 153)"
    
    ### tables
    ,tableBackColor = "transparent"
    ,tableBorderColor = "rgb(80,95,155)"
    ,tableBorderTopSize = 1
    ,tableBorderRowSize = 1
    
  )
  return(theme_purple_gradient)
}


# ------------ Fluid page ------------
#                    Header --------
customLogo <- shinyDashboardLogoDIY(
  boldText = "COVID-19"
  ,mainText = "- Estimating R - "
  ,textSize = 16
  ,badgeText = NULL
  ,badgeTextColor = "white"
  ,badgeTextSize = 2
  ,badgeBackColor = rgb(0,0,0,0)
  ,badgeBorderRadius = 3
)
dbHeader <- dashboardHeaderPlus(title = customLogo, titleWidth = 350, enable_rightsidebar = TRUE
                                # ,left_menu = tagList(
                                #   actionButton(
                                #     inputId = "help",
                                #     label = "",
                                #     icon = icon("info"),
                                #     width = NULL
                                #   )
                                # )
                                )


#                    Sidebar -------
dbSidebar <- dashboardSidebar(
  introjsUI(),

  width = 350,
  sidebarMenu(
    chooseSliderSkin("Square", color = "#112446"),
    HTML('<center><img src="covid.png", height="100px", style="float:center"/></center>','<p style="color:black"></p>'),
    HTML("<br>"),
    menuItem("Welcome", tabName = "Description", icon = icon("home")) %>% add_class("Description"),
    menuItem("Estimating R", tabName = "EstimatingR", icon = icon("dashboard")) %>% add_class("EstimatingR"),
    menuItem("More about the method", tabName = "More", icon = icon("info-circle")) %>% add_class("More"),

    HTML("<br>")
    
    ,id = "sideBarContent"
  )
  ,uiOutput("contentSelectedTabSideBar")
)
#                    Body ----------

dbBody <- dashboardBody(
  
  useShinyjs(),
  getTheme2(),
  setShadow(class = "box"),
  setShadow(class = "dropdown-menu"),
  
  disconnectMessage(
    text = "Something went wrong! Try refreshing the page.",
    refresh = "Refresh",
    background = "rgba(225, 71, 54, 100)",
    colour = "#FFFFFF",
    refreshColour = "#337AB7",
    overlayColour = "rgba(225, 71, 54, 0.2)",
    overlayOpacity = 0.3,
    width = 450,
    top = "center",
    size = 24,
    css = "padding: 15px !important; box-shadow: none !important; background: rgba(199, 59, 89, 0.2);"
  ),
  
  # tags$head(
  #   tags$style(".direct-chat-contacts {background: rgb(255, 0, 255);}
  #              .content {
  #                   background-image: url(https://raw.githubusercontent.com/AntoineB0urret/EstimatingR/master/exampleIMG.png?token=AGGC2VRH7CWZBPKCLTXDEH3ABQ2QY);
  #                   object-fit: cover;
  #              }")
  #   # background-image: url(https://www.pngfind.com/pngs/m/461-4619276_hello-lettering-hd-png-download.png);
  # ),
  
  withSpinner(
  tabItems(
    tabItem(tabName = "EstimatingR",
            chooseSliderSkin("Square", color = "#112446"),
            
            tags$head(
              tags$style("#globalPlots {margin: 10px 10px 10px;}
                         ")
              ),
            
            fluidRow(id = "globalPlots",
              # h1(HTML("<center>Estimating the reproductive number of the COVID-19 epidemic \n </center>")),
              # HTML('<center><img src="covid.png", height="100px", style="float:center"/></center>','<p style="color:black"></p>'),
              # HTML("<br><br>"),
              
              boxPlus(
                tags$head(tags$style(".direct-chat-contacts {background: #ffffff;}
                                     ")),

                width = 12,
                title = "Deconvolution",
                closable = FALSE,
                status = "primary",
                solidHeader = FALSE,
                collapsible = TRUE,
                enable_sidebar = TRUE,
                sidebar_width = 35,
                sidebar_start_open = TRUE,
                id = "DeconvBox",

                sidebar_content = tagList(
                  fluidRow(id = "settings1",
                    # fluidRow(
                    #   column(6, numericInput("shape_deconv", "Shape:", min = 0, max = NA, value = 10.2, step = 0.1)),
                    #   column(6, numericInput("scale_deconv", "Scale:", min =  0, max = NA, value = 0.5, step = 0.1))
                    # ),
                    fluidRow(
                      column(6, numericInput("mean_deconv", "Mean:", min = 0, max = NA, value = 5.1, step = 0.1)),
                      column(6, numericInput("variance_deconv", "Variance:", min =  0, max = NA, value = 2.55, step = 0.05))
                    ),
                    # sliderInput("extension_length", "Extension:",
                    #             min = 1, max = 50,
                    #             value = 20),
                    conditionalPanel(
                      condition = "input.useBasisDim==1",
                      sliderInput("Basis_dimension", "Basis dimension :",
                                  min = 10, max = 200,
                                  value = 50)
                    ),
                    conditionalPanel(
                      condition = "input.useBasisDim==0",
                      numericInput("Frequency_dim", "Frequency :",
                                  min = 1, max = 200,
                                  value = 7)
                    ),
                    numericInput("Lambda_delta", "Regularization coefficient:",
                                 min = 1e-20, max = NA, value = 100, step = 10)
                  )
                ),
                withSpinner(plotOutput("DeconvolutionPlot",  width = "100%", height = "500px", hover = "plot_hover_deconv"), size = 1, color = "#E14736")
              )

              ,boxPlus(
                width = 12,
                title = "Reproductive number",
                closable = FALSE,
                status = "primary",
                solidHeader = FALSE,
                collapsible = TRUE,
                enable_sidebar = TRUE,
                sidebar_width = 35,
                sidebar_start_open = TRUE,
                id = "EstimBox",

                sidebar_content = tagList(
                  fluidRow(id = "settings2",
                           fluidRow(
                             column(6, numericInput("mean_generation", "Mean:", min = 0, max = NA, value = 5.2, step = 0.1)),
                             column(6, numericInput("var_generation", "Variance:", min =  0, max = NA, value = 2.95, step = 0.1))
                           ),
                    selectInput("family_estim", "Family:",
                                c("Poisson" = "poisson",
                                  "quasi-Poisson" = "quasipoisson",
                                  "Negative binomial" = "negbin"),
                                selected = "negbin"),
                    
                    conditionalPanel(
                      condition = "input.useBasisDim==1",
                      sliderInput("Basis_dimension_estim", "Basis dimension :",
                                  min = 10, max = 200,
                                  value = 70)
                    ),
                    conditionalPanel(
                      condition = "input.useBasisDim==0",
                      numericInput("Frequency_dim_estim", "Frequency :",
                                   min = 1, max = 200,
                                   value = 5)
                    )
                    # selectInput("weights", "Weights:", c("No weight", "1/var", "estimate/var", "sqrt(estimate)/var"), selected = "sqrt(estimate)/var")
                    # plotOutput("plotGenerationInterval")
                  )
                ),
                withSpinner(plotOutput("PlotRt", width = "100%", height = "600px", hover = "plot_hover_estim"), size = 1, color = "#E14736")
                # withSpinner(plotOutput("RtPlot", width = "100%", height = "300px", hover = "plot_hover_estim"), size = 1, color = colors[2]),
                # withSpinner(plotOutput("RtDerivPlot",  width = "100%", height = "300px", hover = "plot_hover_estim_deriv"), size = 1, color = colors[2])
              )
          

            )
            
    ),
    
    tabItem(tabName = "Description",
            fluidRow(id = "descriptionMethod",
                     h2("How to estimate the reproductive number ?"),
                     uiOutput('detailsMethod')
                     )
    ),
    tabItem(tabName = "More",
            includeHTML('www/EstimatingRtCOVID19.html')
    )
  )
  , size = 1, color = "#E14736")
)


#                    RightSideBar ---------
rightsidebar = rightSidebar(
  background = "dark",
  width = 250,
  rightSidebarTabContent(
    id  = "rightSideBarSettings",
    icon = NULL,
    active = TRUE,
    title = "Super settings"

    ,
    # checkboxInput(inputId = "useBasisDim", label = "Use basis dimension", value = FALSE)
    
    materialSwitch(inputId = "NightMode", label = "Night Mode", status = "danger", value = FALSE, right = TRUE)
    ,HTML("<br>")
    
    ,materialSwitch(inputId = "useBasisDim", label = "Use basis dimension", status = "danger", value = FALSE, right = TRUE)
    ,HTML("<br>")
    ,materialSwitch(inputId = "Show_deriv", label = "Show derivative", status = "danger", value = TRUE, right = TRUE)
    ,HTML("<br>")
    ,materialSwitch(inputId = "simple_graphs", label = "Simple graphs", status = "danger", value = TRUE, right = TRUE)
    ,HTML("<br>")
    # ,sliderInput("extension_length", "Extension:",
    #             min = 1, max = 50,
    #             value = 20)
    ,HTML("<br>")
    ,selectInput("weights", "Weights:", c("No weight", "1/var", "estimate/var", "sqrt(estimate)/var"), selected = "sqrt(estimate)/var")
    ,HTML("<br>")
    ,column(width = 12,
            actionButton(inputId = "EstimateDimAndRegularization",
                          label = "Choose best",
                          icon = icon("rocket"),
                          width = NULL,
                    )
            ,align = "center"
    )
    # ,checkboxInput(inputId = "Show_deriv", label = "Show derivative", value = TRUE)
  )
)
#                    ALL -----------
ui <- shinyUI(dashboardPagePlus(
  # includeCSS("www/style_day.css"),
  uiOutput("style"),
  header = dbHeader,
  sidebar = dbSidebar,
  body = dbBody,
  rightsidebar = rightsidebar
))

# ------------ Server ------------

server <- shinyServer(function(input, output, session) {
  
  # ----------------- TO BE REMOVED: test deconnection ------------- 
  observeEvent(input$disconnect, {
    session$close()
  })
  
  # observeEvent(input$bttn2, {
  #   js$collapse("settings")
  #   js$collapse("deconv")
  #   js$collapse("estim")
  #   js$collapse("estimDeriv")
  # })
  
  # ----- Methods for titles and description 
  # formulaText <- reactive({
  #   paste("COVID-19", "\n Selected country: ", input$location)
  # })
  # output$caption2 <- renderUI({
  #   h3(HTML(paste("Estimating the reproductive number of COVID-19", sep="<br/>")))
  # })
  # output$selected_country <- renderUI({
  #   h5(HTML(paste("Selected country: ", input$location), sep="<br/>"))
  # })
  # output$caption <- renderText({
  #   # formulaText()
  #   paste("COVID-19", "\n Selected country: ", input$location)
  # })
  
  # ============================ Styles and night/day modes ============================
  # Inport CSS style, dark and day mode
  output$style <- renderUI({
    if (!is.null(input$NightMode)) {
      if (input$NightMode) {
        includeCSS("www/style_dark.css")
      } else {
        includeCSS("www/style_day.css")
      }
    }
  })
  
  # Night/Day mode them for plots
  observeEvent(input$NightMode, {
    colors <- c("#00AFBB", "#E7B800", "#FC4E07")
    if (input$NightMode == 1){
      font_color_plot <- '#1E1D1E'
      points_color <- "grey85"
      PlotsData$theme <- theme(plot.background = element_rect(fill = font_color_plot, colour = "black"),
                               panel.background = element_rect(fill = font_color_plot, colour = "black"),
                               panel.grid.major = element_line(colour = "grey55"),
                               panel.grid.minor = element_line(colour = "grey40"),
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
      # theme_shiny <- shinytheme('slate')
    }else{
      font_color_plot <- "white"
      points_color <- "black"
      PlotsData$theme <- theme(plot.background = element_rect(fill = font_color_plot),
                               panel.background = element_rect(fill = "white", colour = "black"),
                               panel.grid.major = element_line(colour = "grey85"),
                               panel.grid.minor = element_line(colour = "grey90"),
                               legend.position = c(.95, .95),
                               legend.justification = c("right", "top"),
                               legend.box.just = "right",
                               legend.margin = margin(6, 6, 6, 6),
                               plot.margin = unit(c(1,1,1,1), "cm"),
                               legend.key = element_rect(colour = "transparent", fill = "transparent"))
      # theme_shiny <- shinytheme('paper')
    }
  })
  
  # ============================ Some auxiliary functions ============================
  is.integer0 <- function(x){ is.integer(x) && length(x) == 0L }
  
  removeNA <- function(x){
    x[which(is.na(x))] <- 0
    x[which(x < 0)] <- 0
    return(x)
  }
  
  #' Get all the results for a specific date and put them in a list
  getResults <- function(date_of_interest){
    dates_vector <- GetDATA_List$body$data$Date[GetDATA_List$body$index]
    index_ = which(dates_vector == date_of_interest)
    return(list("Deconv_fitted" = Deconvolution_List$body$deconvolution$fitted_values[Deconvolution_List$body$index_deconv][index_],
                "Deconv_simultaneous" = list("lower" = Deconvolution_List$body$confidence_bands$Simultaneous_lower[Deconvolution_List$body$index_deconv][index_],
                                             "upper" = Deconvolution_List$body$confidence_bands$Simultaneous_upper[Deconvolution_List$body$index_deconv][index_]),
                "Deconv_point_wise" = list("lower" = Deconvolution_List$body$confidence_bands$Point_wise_CI_lower[Deconvolution_List$body$index_deconv][index_],
                                           "upper" = Deconvolution_List$body$confidence_bands$Point_wise_CI_upper[Deconvolution_List$body$index_deconv][index_]),
                "R_fitted" = EstimatingR_List$body$modelFit$fitted.values[index_],
                "R_simultaneous" = list("lower" = EstimatingR_List$body$modelFit$derivPlots$simultaneousIntervals[[1]]$Simultaneous_lower[index_],
                                        "upper" = EstimatingR_List$body$modelFit$derivPlots$simultaneousIntervals[[1]]$Simultaneous_upper[index_]),
                "R_point_wise" = list("lower" = EstimatingR_List$body$modelFit$derivPlots$simultaneousIntervals[[1]]$Point_wise_CI_lower[index_],
                                      "upper" = EstimatingR_List$body$modelFit$derivPlots$simultaneousIntervals[[1]]$Point_wise_CI_upper[index_])
    ))
  }
  
  
  # ============================ side bar and header ============================
  # Hide sidebar
  shinyjs::runjs("document.getElementsByClassName('sidebar-toggle')[0].style.visibility = 'hidden';")
  observeEvent(input$hide,{
    shinyjs::runjs("document.getElementsByClassName('sidebar-toggle')[0].style.visibility = 'hidden';")
  })
  observeEvent(input$show,{
    shinyjs::runjs("document.getElementsByClassName('sidebar-toggle')[0].style.visibility = 'visible';")
  })
  
  # When the 'Estimating R' tab is selected (left side bar), add inputs for data and results below
  output$contentSelectedTabSideBar <- renderUI({
    if (input$sideBarContent == "EstimatingR"){
      if (!is.null(Deconvolution_List$body)){
        download_button <- downloadButton('downloadData', '')
      }else{
        download_button <- HTML("<br>")
      }
      fluidRow(
        HTML("<br>"),
        fluidRow(id = "sidebarsettings", inputId = "sidebarsettings",
                 HTML("<br>"),
                 HTML('<center>DATA</center>','<p style="color:black"></p>'),
                 tags$div(id='locationAndSource',
                          class='locationAndSource_class',
                          selectInput("location", "Location:",
                                      unique_countries, selected = "Switzerland", selectize = T),
                          selectInput("source", "Source:", c("OFSP", "Our World in Data", "HDX"), selected = "OFSP", selectize = T),
                          sliderInput("datesStartEnd", "Range:",
                                      min = as.Date("2020-01-01"), max = Sys.Date(), value = c(as.Date("2020-01-01"), Sys.Date()))
                 )
        )
        ,HTML("<br><br>")
        
        ,icon("arrow-alt-circle-down")
        ,HTML("<br><br>")
        
        ,actionButton(
          inputId = "bttn2",
          label = "Let's estimate !",
          icon = icon("rocket"),
          width = NULL
        )
        ,HTML("<br>")
        # ,actionButton('disconnect', 'Disconnect the app')
        # ,actionButton("help", "Press for instructions")
        
        ,actionButton(
          inputId = "help",
          label = "Need help !",
          icon = icon("info"),
          width = NULL
        )
        
        ,HTML("<br><br>")
        ,icon("arrow-alt-circle-down")
        
        ,HTML("<br><br>")
        
        ,fluidRow(id = "Results", inputId = "Results",
                  HTML("<br>"),
                  HTML('<center>Results</center>','<p style="color:black"></p>'),
                  sliderInput("resultInterestDate", "Date of interest:",
                              min = as.Date("2020-01-01"), max = Sys.Date(), value = Sys.Date() - 10)
                  ,column(width = 12,
                          # uiOutput("dynamic"),
                          tableOutput("dynamic"),
                          align = "center")
                  #,uiOutput("dynamic")
                  
                  ,HTML("<br>")
                  ,download_button
                  ,HTML("<br>")
        )
        ,HTML("<br><br>")
        ,align="center"
      )

    }
  })
  
  # Update source input when country is changed
  updateSource <- reactive({
    if (!is.null(input$location)){
      if (input$location == "Switzerland"){
        c("OFSP", "Our World in Data", "HDX")
      }else{
        c("Our World in Data", "HDX")
      }
    }
  })
  
  observe({
    updateSelectInput(session, "source", choices = updateSource())
  })
  
  # ============================ Reactive values and parameters ============================
  GetDATA_List <- reactiveValues("params_hash" = "", "changed" = TRUE, "plot"= NULL, "body" = NULL)
  Deconvolution_List <- reactiveValues("params_hash" = "", "changed" = TRUE, "plot"= NULL, "body" = NULL)
  EstimatingR_List <- reactiveValues("params_hash" = "", "changed" = TRUE, "plot"= NULL, "body" = NULL)
  # plotInfo <- reactiveValues("max_width_deconv" = NULL, "max_width_estim" = NULL, "max_width_deriv" = NULL)
  HoverInfo <- reactiveValues("Info" = NULL)
  PlotsData <- reactiveValues("theme" = NULL)
  DataAllCountries <- NULL
  SourceData <- NULL
  status <- reactiveVal()
  
  # Some Hyperparameters:
  threshold_beginning <- 10

  
  # ============================ Observe events ============================
  
  # Enable button when a parameter is changed
  observeEvent({list("location" = input$location,
                     "source" = input$source,
                     "datesStartEnd" = input$datesStartEnd,
                     "shape_deconv" = input$shape_deconv,
                     "scale_deconv" = input$scale_deconv,
                     "extension_length" = input$extension_length,
                     "Basis_dimension" = input$Basis_dimension,
                     "Frequency_dim" = input$Frequency_dim,
                     "Lambda_delta" = input$Lambda_delta,
                     "family_estim" = input$family_estim,
                     "Basis_dimension_estim" = input$Basis_dimension_estim,
                     "Frequency_dim_estim" = input$Frequency_dim_estim,
                     "weights" = input$weights,
                     "useBasisDim" = input$useBasisDim,
                     "mean_generation"=input$mean_generation, 
                     "var_generation"=input$var_generation,
                     "variance_deconv"= input$variance_deconv,
                     "mean_deconv" = input$mean_deconv)},
               {
                 status("Needs recalculation")
                 shinyjs::enable("bttn2")
                 shinyjs::enable("EstimateDimAndRegularization")})
  
  # Update panel -> data, deconv, estim
  observeEvent(input$bttn2, {
    suppressWarnings({
      # DATA
      obtainData()
      
      # Deconvolution
      getTrueIncidence()
      
      # Estimatnig Rt
      estimateRt()
    })
  })
  
  # ------------------------------------ Button estimating basis dimension and regularization -----------
  observeEvent(input$EstimateDimAndRegularization, {
    
    # DATA
    obtainData()
    
    shinyjs::disable("EstimateDimAndRegularization")
    
    if (!is.null(GetDATA_List$body)){
      print("OK, computing...")
      DATA <- GetDATA_List$body
      index <- DATA$index
      data_country <- DATA$data
      
      w_conv <- discrete_gamma(n = 20, normalize = TRUE, shape = (input$mean_deconv^2)/input$variance_deconv, scale = input$variance_deconv/input$mean_deconv)
      length_extension <- max(20, 2*input$mean_deconv)
      # w_conv <- discrete_gamma(n = 20, normalize = TRUE, shape = input$shape_deconv, scale = input$scale_deconv)
      
      # Extension before:
      extension_first_length <- 20
      data_cases_extension_start <- c(rep(1, extension_first_length), data_country$new_confirmed[index])
      
      smoothing_dim <- 7; smoothing_lambda <- 7;   # --> small order, small lambda   # Report: Figure 23a
      # smoothing_dim <- 1; smoothing_lambda <- 7;   # --> select high dimension basis, compensate by chosing large lambda for the penalization   # Report: Figure 23b
      # smoothing_dim <- 1; smoothing_lambda <- 1;   # high dimension basis, low lambda    # Report: Figure 23c
      
      
      # Now we search for the best basis dimension, as well as the best penalization based on the chosen dimension.
      runjs(paste0('$("#EstimateDimAndRegularization.btn").css("background","linear-gradient(to right, #333, #333 ', 0, '%, #eee ', 0, '%, #eee 100%)")'))
      
      observed_data <- data_cases_extension_start
      prob_vector <- w_conv
      range_values_basis = c(30:150)
      range_values_penalization = 10**(20:-10)
      total_length <- length(range_values_basis) + length(range_values_penalization)
      which_regularization <- 'tridiag'
      
      # Part 1: Find basis dimension:
      indices_init <- c(1:length(observed_data))
      if (!is.null(smoothing_dim)){
        data_inc_after <- rollapply(observed_data, width=smoothing_dim, FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                    by.column=TRUE, partial=TRUE, fill=NA, align="center")
      }
      errors <- numeric(length(range_values_basis))
      B <- construct_B(length(observed_data), prob_vector)
      
      suppressWarnings({
        for (j in c(1:length(range_values_basis))){
          # Deconvolution
          results_deconv <- DeconvolutionGAM(observed_data, prob_vector, range_values_basis[j], delta = NULL, 
                                             plot_ = FALSE, 
                                             extension = list('len' = length_extension, 'family' = nb(link = 'log'), 'keep_original' = TRUE, 'min_weights' = 0.1),
                                             fm = "cases ~ s(time, k = 50, bs = \'tp\', m = c(0, 0)) + s(day, k = 5)"
          )
          fitted <- results_deconv$fitted
          fitted[which(is.na(fitted))] <- 0
          # compute MSE between the convolved deconvolved data and the (smoothed) observed data
          errors[j] <- mean((B %*% fitted[indices_init] - data_inc_after[indices_init])**2, na.rm = TRUE)
          runjs(paste0('$("#EstimateDimAndRegularization.btn").css("background","linear-gradient(to right, #333, #333 ', j/total_length*100, '%, #eee ', j/total_length*100, '%, #eee 100%)")'))
        }
      })
      dimensionBasis <- range_values_basis[which(errors == min(errors))]
      
      
      # Part 2: regularization parameter
      
      indices_init <- c(1:length(observed_data))
      # Apply a smoothing window on the observed data
      if (!is.null(smoothing_lambda)){
        data_inc_after <- rollapply(observed_data, width=smoothing_lambda, FUN=function(x) mean(x, na.rm=TRUE), by=1,
                                    by.column=TRUE, partial=TRUE, fill=NA, align="center")
      }
      errors_reg <- numeric(length(range_values_penalization))
      B <- construct_B(length(observed_data), prob_vector)
      
      # iterate over the range of values for lambda
      suppressWarnings({
        for (j in c(1:length(range_values_penalization))){
          if (which_regularization == 'tridiag'){
            delta <- construct_delta(dimensionBasis + 4)
          }else{
            delta <- diag(1+numeric(dimensionBasis + 4))
          }
          delta <- range_values_penalization[j] * delta
          
          # Deconvolution
          results_deconv_reg <- DeconvolutionGAM(observed_data, prob_vector, dimensionBasis, delta = delta, 
                                                 plot_ = FALSE, 
                                                 extension = list('len' = length_extension, 'family' = nb(link = 'log'), 'keep_original' = TRUE, 'min_weights' = 0.1),
                                                 fm = "cases ~ s(time, k = 50, bs = \'tp\', m = c(0, 0)) + s(day, k = 5)",
                                                 family = quasipoisson(link = 'log'))
          fitted <- results_deconv_reg$fitted
          fitted[which(is.na(fitted))] <- 0
          # compute MSE between the convolved deconvolved data and the (smoothed) observed data
          errors_reg[j] <- mean((B %*% fitted[indices_init] - data_inc_after[indices_init])**2, na.rm = TRUE)
          runjs(paste0('$("#EstimateDimAndRegularization.btn").css("background","linear-gradient(to right, #333, #333 ', (j+length(range_values_basis))/total_length*100, '%, #eee ', (j+length(range_values_basis))/total_length*100, '%, #eee 100%)")'))
          
        }
      })
      lambda_delta <- range_values_penalization[which(errors_reg == min(errors_reg))]
      
      print(paste("Basis dimension:", dimensionBasis))
      print(paste("Regularization:", lambda_delta))
      runjs(paste0('$("#EstimateDimAndRegularization.btn").css("background","linear-gradient(to right, #333, #333 ', 100, '%, #eee ', 100, '%, #eee 100%)")'))
      
      updateNumericInput(session, "Lambda_delta", value = lambda_delta)
      updateSliderInput(session, "Basis_dimension", value = dimensionBasis)
      frequency_used <- round(length(observed_data)/dimensionBasis, digits = 2)
      updateNumericInput(session, "Frequency_dim", value = frequency_used)
      # if (input$useBasisDim == 0){  # use frequency
      #   frequency_used <- round(length(observed_data)/dimensionBasis, digits = 2)
      #   updateNumericInput(session, "Frequency_dim", value = frequency_used)
      # }else{
      #   updateSliderInput(session, "Basis_dimension", value = dimensionBasis)
      #   updateNumericInput(session, "Frequency_dim", value = frequency_used)
      # }
      
      # Deconvolution
      getTrueIncidence(basisDim = dimensionBasis, lambda = lambda_delta)
      
      # Estimatnig Rt
      estimateRt()
      
      
    }
  })
  
  # ============================ Results ============================
  output$dynamic <- renderTable({
    input$resultInterestDate
    hover <- input$plot_hover_deconv
    
    if(!is.null(EstimatingR_List$body)){
      if ((!is.null(hover) && (as.integer(hover$x) > 0))){
        date_of_interest <- as.Date(as.integer(hover$x))
        updateSliderInput(session, "resultInterestDate", value = date_of_interest)
      }else{
        date_of_interest <- input$resultInterestDate
      }
      results <- getResults(date_of_interest)
      if (!is.integer0(which(GetDATA_List$body$data$Date[GetDATA_List$body$index] == date_of_interest))){
        M <- matrix(c(ceiling(results$Deconv_fitted), 
                      paste0("[", ceiling(results$Deconv_point_wise$lower), ";", ceiling(results$Deconv_point_wise$upper), "]"),
                      paste0("[", ceiling(results$Deconv_simultaneous$lower), ";", ceiling(results$Deconv_simultaneous$upper), "]"),
                      
                      round(results$R_fitted, digits = 3),
                      paste0("[", round(results$R_point_wise$lower, digits = 3), ";", round(results$R_point_wise$upper, digits = 3), "]"),
                      paste0("[", round(results$R_simultaneous$lower, digits = 3), ";", round(results$R_simultaneous$upper, digits = 3), "]")
        ),nrow=3)
        rownames(M) <- c('Estimates','95% CI (p)', '95% CI (s)')
        colnames(M) <- c("I", "R")
        as.data.frame(M)
      }
    }
  }, rownames = TRUE, align = "ccc")
  
  # ============================ Titles and other =======================
  # Display the image
  output$renderLogo <- renderImage({
    outfile <- tempfile(fileext = '.png')
    
    # Generate the PNG
    png(outfile, width = 400, height = 300)
    hist(rnorm(input$obs), main = "Generated in renderImage()")
    dev.off()
    
    outfile <- "./www/covid.png"
    
    # Return a list containing the filename
    list(src = outfile,
         contentType = 'image/png',
         width = 70,
         height = 70,
         alt = "This is alternate text")
  }, deleteFile = FALSE)
  
  
  # ============================ Computation ============================
  # ------------------------------------ Fonctions ----------------------
  # Get data function
  obtainData <- function(){
    shinyjs::disable("bttn2")
    params <- list("location" = input$location,
                   "source" = input$source,
                   "datesStartEnd" = input$datesStartEnd)
    
    
    # Get database
    if ((is.null(DataAllCountries)) || (is.null(SourceData)) || (SourceData != params$source)){
      # Need to retreive the data 
      SourceData <<- params$source
      
      if (SourceData == "OFSP"){
        data <- read.xlsx(sep=",",startRow = 8, detectDates = TRUE,
                          "https://www.bag.admin.ch/dam/bag/de/dokumente/mt/k-und-i/aktuelle-ausbrueche-pandemien/2019-nCoV/covid-19-datengrundlage-lagebericht.xlsx.download.xlsx/200325_Datengrundlage_Grafiken_COVID-19-Bericht.xlsx")
        names(data) <- c("date","cases","casesCumul","hospitalized","hospitalizedCumul",
                         "deaths","deathsCumul")
        DataAllCountries <<- data
        
      }else if (SourceData == "Our World in Data"){
        data <- as.data.frame(read_csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv"))
        DataAllCountries <<- data
        
      }else if (SourceData == "HDX"){
        DataAllCountries <<- as.data.frame(read_csv(file = "https://data.humdata.org/hxlproxy/api/data-preview.csv?url=https%3A%2F%2Fraw.githubusercontent.com%2FCSSEGISandData%2FCOVID-19%2Fmaster%2Fcsse_covid_19_data%2Fcsse_covid_19_time_series%2Ftime_series_covid19_confirmed_global.csv&filename=time_series_covid19_confirmed_global.csv"))
        
      }else{
        # handle this situation....
      }
      
    }
    
    # Get data for specific country
    if (GetDATA_List$params_hash != paste(params$location, params$datesStartEnd[1], params$datesStartEnd[2], params$source)){
      print("(DATA) Getting data...")
      
      selected_country <- params$location
      # if (selected_country == "Switzerland"){ option <- 1 }else{ option <- 2 }
      # data <- ConstructDataset(selected_country, params$source, F)
      
      # Get data for specific country:
      if (SourceData == "OFSP"){
        data_cases <- DataAllCountries$cases
        data_dates <- DataAllCountries$date
        data_deaths <- DataAllCountries$deaths
        data_dates <- DataAllCountries$date
        
        inter_dates <- as.Date(intersect(data_dates, data_dates))
        data_deaths <- removeNA(data_deaths)
        data_cases <- removeNA(data_cases)
        
        # Put all data into a data frame
        data_country <- data.frame(matrix(ncol = 7, nrow = length(inter_dates)))
        colnames(data_country) <- c("Date", "confirmed", "deaths", "recovered", "new_confirmed", "new_deaths", "new_recovered")
        data_country$Date <- inter_dates
        data_country$new_confirmed <- data_cases[which(data_dates %in% inter_dates)]
        data_country$new_deaths <- data_deaths[which(data_dates %in% inter_dates)]
        
        data <- data_country
        
      }else if (SourceData == "Our World in Data"){
        data_c <- DataAllCountries[which(DataAllCountries[,"location"] == selected_country), ]
        
        # Some preprocessing, remove first zeros
        first_value = 1
        data_c$new_cases <- removeNA(data_c$new_cases)
        while(data_c$new_cases[first_value] < threshold_beginning){ first_value = first_value+1 }
        data_c <- data_c[first_value:length(data_c$date),]
        data_c$new_cases[which(data_c$new_cases < 1)] <- 1
        
        data_cases <- data_c$new_cases
        data_deaths <- data_c$new_deaths
        data_dates <- data_c$date
        
        inter_dates <- as.Date(data_dates)
        data_deaths <- removeNA(data_deaths)
        data_cases <- removeNA(data_cases)
        
        # Put all data into a data frame
        data_country <- data.frame(matrix(ncol = 7, nrow = length(inter_dates)))
        colnames(data_country) <- c("Date", "confirmed", "deaths", "recovered", "new_confirmed", "new_deaths", "new_recovered")
        data_country$Date <- inter_dates
        data_country$new_confirmed <- data_cases[which(data_dates %in% inter_dates)]
        data_country$new_deaths <- data_deaths[which(data_dates %in% inter_dates)]
        
        data <- data_country
        
      }else if (SourceData == "HDX"){
        
        data_confirmed_matrix <- data.matrix(DataAllCountries[5:length(DataAllCountries[1,])])
        
        l <- length(DataAllCountries[1, 5:length(DataAllCountries[1,])])
        
        data_country <- data.frame(matrix(ncol = 7, nrow = l))
        colnames(data_country) <- c("Date", "confirmed", "deaths", "recovered", "new_confirmed", "new_deaths", "new_recovered")
        data_country$Date <- as.Date(colnames(DataAllCountries)[5:length(DataAllCountries[1,])], "%m/%d/%y")
        # data_country$Date <- colnames(data_confirmed)[5:length(data_confirmed[1,])]
        
        if (length(which(DataAllCountries[,2] == selected_country)) <= 1){
          data_country$confirmed <- data_confirmed_matrix[which(DataAllCountries[,2] == selected_country),]
        }else{
          data_country$confirmed <- colSums(data_confirmed_matrix[which(DataAllCountries[,2] == selected_country),])
        }
        data_country$new_confirmed <- data_country$confirmed - c(0, data_country$confirmed[1:(l-1)])
        
        # Some preprocessing, remove first zeros
        first_value = 1
        data_country$new_confirmed <- removeNA(data_country$new_confirmed)
        while(data_country$new_confirmed[first_value] < threshold_beginning){ first_value = first_value+1 }
        data_country$new_confirmed[which(data_country$new_confirmed < 1)] <- 1
        
        data <- data_country[first_value:length(data_country$Date), ]
      }else{
        # handle this situation ...
      }
      
      data$new_confirmed <- removeNA(data$new_confirmed)
      
      start_index <- which(data$Date == params$datesStartEnd[1])
      if(is.integer0(start_index)){start_index <- 1}
      end_index <- which(data$Date == params$datesStartEnd[2])
      if(is.integer0(end_index)){end_index <- length(data$Date)}
      index <- intersect(c(1:(length(data$Date))), seq(start_index, end_index, by=1))
      
      GetDATA_List$body <- list("data" = data, "index" = index)
      GetDATA_List$params_hash <- paste(params$location, params$datesStartEnd[1], params$datesStartEnd[2], params$source)
      GetDATA_List$changed <- TRUE
    }else{
      GetDATA_List$changed <- FALSE
    }
  }
  
  # Estimate true incidence events function
  getTrueIncidence <- function(basisDim = NULL, lambda = NULL){
    shinyjs::disable("bttn2")
    
    if (is.null(basisDim)){
      basisDim <- input$Basis_dimension
    }
    if (is.null(lambda)){
      lambda <- input$Lambda_delta
    }
    params <- list("shape_deconv" = input$shape_deconv,
                   "scale_deconv" = input$scale_deconv,
                   "extension_length" = input$extension_length,
                   "Basis_dimension" = basisDim,
                   "Frequency_dim" = input$Frequency_dim,
                   "Lambda_delta" = lambda,
                   "useBasisDim" = input$useBasisDim,
                   "variance_deconv"= input$variance_deconv,
                   "mean_deconv" = input$mean_deconv)
    
    if ((Deconvolution_List$params_hash != paste(params$shape_deconv, params$scale_deconv, params$extension_length, basisDim, lambda, 
                                                 params$useBasisDim, params$Frequency_dim, params$variance_deconv, params$mean_deconv)) ||
        GetDATA_List$changed){
      
      print("(Deconvolution) Estimating new incidence events...")
      
      DATA <- GetDATA_List$body
      index <- DATA$index
      data_country <- DATA$data
      
      w_conv <- discrete_gamma(n = 20, normalize = TRUE, shape = (params$mean_deconv^2)/params$variance_deconv, scale = params$variance_deconv/params$mean_deconv)
      # w_conv <- discrete_gamma(n = 20, normalize = TRUE, shape = params$shape_deconv, scale = params$scale_deconv)
      
      # Extension before:
      extension_first_length <- 20
      data_cases_extension_start <- c(rep(1, extension_first_length), data_country$new_confirmed[index])
      
      if (params$useBasisDim == 0){  # use frequency
        basisDim_used <- ceiling(length(data_cases_extension_start)/params$Frequency_dim)
      }else{
        basisDim_used <- basisDim
      }
      
      length_extension <- max(20, 2*params$mean_deconv)
      
      delta <- construct_delta(basisDim_used + 4)
      delta <- lambda * delta
      results_IWLS_extended <- DeconvolutionGAM(data_cases_extension_start, w_conv, basisDim_used, plot_ = FALSE, delta = delta,
                                                get_derivatives = FALSE, family = quasipoisson(link = 'log'),
                                                extension = list('len' = length_extension, 'family' = nb(link = 'log'), 'keep_original' = TRUE, 'min_weights' = 0.1),
                                                fm = "cases ~ s(time, k = 100, bs = \'tp\', m = c(0, 0)) + s(day, k = 5)")
      ConfidenceBands_extended <- SimultaneousIntervals(results_IWLS_extended$SplineBasis, plot_ = TRUE,
                                                        coefficients_model = results_IWLS_extended$fit$coefficients,
                                                        CovMatrix = vcov(results_IWLS_extended$fit),
                                                        sample_splines = TRUE, ilink_fct = family(results_IWLS_extended$fit)$linkinv, N_keep = 100)
      
      Deconvolution_List$body <- list("deconvolution" = results_IWLS_extended,
                                      "confidence_bands" = ConfidenceBands_extended,
                                      "extension_original" = data_cases_extension_start,
                                      "parameters" = params,
                                      "data_original" = DATA,
                                      "index_deconv" = seq(extension_first_length+1, length(data_cases_extension_start), by = 1),
                                      "basisDim_used" = basisDim_used)
      Deconvolution_List$params_hash <- paste(params$shape_deconv, params$scale_deconv, params$extension_length, basisDim, lambda, 
                                              params$useBasisDim, params$Frequency_dim, params$variance_deconv, params$mean_deconv)
      Deconvolution_List$changed <- TRUE
      
      
      # Deconvolution plot
      deconv_res <- Deconvolution_List$body
      index <- deconv_res$data_original$index
      dates_vector <- GetDATA_List$body$data$Date[GetDATA_List$body$index]
      indices_used <- Deconvolution_List$body$index_deconv
      DATA <- GetDATA_List$body
      index <- DATA$index
      data_country <- DATA$data
      
      if (input$NightMode == 1){
        color_fitted <- "white"
      }else{
        color_fitted <- "black"
      }
      
      plot_deconvolution <- ggplot(data = NULL, aes(x = dates_vector, y = GetDATA_List$body$data$new_confirmed[GetDATA_List$body$index], colour = "Observed new cases")) + geom_line() +
        geom_ribbon(aes(ymin = Deconvolution_List$body$confidence_bands$Point_wise_CI_lower[indices_used], ymax = Deconvolution_List$body$confidence_bands$Point_wise_CI_upper[indices_used], fill = "Point-wise confidence bands", colour = NA), alpha = 0.3, fill = colors[3]) +
        geom_ribbon(aes(ymin = Deconvolution_List$body$confidence_bands$Simultaneous_lower[indices_used], ymax = Deconvolution_List$body$confidence_bands$Simultaneous_upper[indices_used], fill = "Simlutaneous confidence bands", colour = NA), alpha = 0.3, fill = colors[3]) +
        geom_line(aes(y = Deconvolution_List$body$deconvolution$fitted_values[indices_used], colour = 'Deconvolution')) +
        scale_colour_manual("",
                            breaks = c("incidence", 'Observed new cases', 'Deconvolution', "Point-wise confidence bands", "Simlutaneous confidence bands"),
                            values = c(color_fitted, color_fitted, "#f4e04d", colors[3], colors[3])) +
        scale_fill_manual("",
                          breaks = c("Point-wise confidence bands", "Simlutaneous confidence bands"),
                          values = c(colors[3], colors[3])) +
        labs(title="", x="Time", y="number of new cases") +
        theme_for_the_plots +
        coord_cartesian(ylim = c(0, 1.1 * max(data_country$new_confirmed[index])), xlim = c(dates_vector[1], dates_vector[length(dates_vector)])) +
        guides(color=guide_legend(override.aes=list(fill=NA))) +
        # theme(legend.position="bottom") +
        scale_x_date(date_breaks = "1 month", date_labels = "%b")    # or : date_labels = "%Y (%b)"
      
      # plotInfo$max_width_deconv <- grid::unit.pmax(ggplot_gtable(ggplot_build(plot_deconvolution))$widths[2:3])
      
      Deconvolution_List$plot <- plot_deconvolution
    }else{
      Deconvolution_List$changed <- FALSE
    }
  }
  
  # Estimate the reproductive number
  estimateRt <- function(){
    
    shinyjs::disable("bttn2")
    params <- list("family_estim" = input$family_estim,
                   "Basis_dimension_estim" = input$Basis_dimension_estim,
                   "Frequency_dim_estim" = input$Frequency_dim_estim,
                   "useBasisDim" = input$useBasisDim,
                   "weights" = input$weights,
                   "mean_generation" = input$mean_generation,
                   "var_generation" = input$var_generation)
    
    # getTrueIncidence()
    deconvolution_results <- Deconvolution_List$body
    
    if ((EstimatingR_List$params_hash != paste(params$family_estim, params$Basis_dimension_estim, params$weights,
                                               params$mean_generation, params$var_generation, params$useBasisDim, params$Frequency_dim_estim)) || Deconvolution_List$changed){
      
      print("(Reproductive number) Estimating R_t...")
      indices_used <- deconvolution_results$index_deconv
      Swiss_deconvolution <- deconvolution_results$deconvolution$fitted_values[indices_used]
      
      results_deconvolution <- deconvolution_results$deconvolution
      condidence_bands <- deconvolution_results$confidence_bands
      
      # W_gam <- discrete_normal()$prob_vector
      W_gam <- discrete_normal(30, mean_ = params$mean_generation + 1, se_ = sqrt(params$var_gen))$prob_vector
      
      offsetValues_IWLS <- convolve(c(0* 1:(length(W_gam)-1), Swiss_deconvolution), rev(W_gam), type = "filter")
      offsetValues_IWLS[which(offsetValues_IWLS <= 0)] <- 1
      
      if (params$weights == "No weight"){
        weights <- NULL
      }else if (params$weights == "1/var"){
        weights <- abs(condidence_bands$Point_wise_CI_upper[indices_used] - condidence_bands$Point_wise_CI_lower[indices_used])/4
        weights <- weights**2
        weights <- 1/weights
        weights <- weights/mean(weights)
      }else if (params$weights == "estimate/var"){
        weights <- abs(condidence_bands$Point_wise_CI_upper[indices_used] - condidence_bands$Point_wise_CI_lower[indices_used])/4
        weights <- weights**2
        weights <- weights /results_deconvolution$fitted_values[indices_used]
        weights <- 1/weights
        weights <- weights/mean(weights)
      }else if (params$weights == "sqrt(estimate)/var"){
        weights <- abs(condidence_bands$Point_wise_CI_upper[indices_used] - condidence_bands$Point_wise_CI_lower[indices_used])/4
        weights <- weights**2
        weights <- weights /sqrt(results_deconvolution$fitted_values[indices_used])
        weights <- 1/weights
        weights <- weights/mean(weights)
      }else{
        weights <- NULL
      }
      
      if (params$family_estim == "poisson"){
        which_model <- poisson(link = "log")
      }else if (params$family_estim == "quasipoisson"){
        which_model <- quasipoisson(link = "log")
      }else{
        which_model <- nb(link = "log")
      }
      
      if (params$useBasisDim == 0){  # use frequency
        basisDim_used <- ceiling(length(Swiss_deconvolution)/params$Frequency_dim_estim)
      }else{
        basisDim_used <- params$Basis_dimension_estim
      }
      
      formula <- paste('resp ~ s(x, k =', basisDim_used, ")")
      
      DAT_IWLS <- data.frame(resp = Swiss_deconvolution,
                             offs = offsetValues_IWLS,
                             x = c(1:length(Swiss_deconvolution)))
      fit_IWLS <- fitModel(DAT_IWLS, which_model = which_model, weights_ = weights, fm = formula, ignore_first_max = 10,
                           deriv_plots = c(0,1), detect_flat = TRUE,
                           axis_x = GetDATA_List$body$data$Date[GetDATA_List$body$index])
      
      EstimatingR_List$body <- list("modelFit" = fit_IWLS,
                                    "results_deconvolution" = deconvolution_results,
                                    "data" = deconvolution_results$data_original$data,
                                    "parameters" = params)
      EstimatingR_List$params_hash <- paste(params$family_estim, params$Basis_dimension_estim, params$weights, params$mean_generation, params$var_generation,
                                            params$useBasisDim, params$Frequency_dim_estim)
      EstimatingR_List$changed <- TRUE
      
      dates_vector <- GetDATA_List$body$data$Date[GetDATA_List$body$index]
      plot_estim_r <- EstimatingR_List$body$modelFit$derivPlots$list_plots[[1]] +
        # labs(title = "", y = latex2exp::TeX(paste("Estimates of", "$R_t$")), x = "time") +
        coord_cartesian(ylim = c(0, 4.2), xlim = c(dates_vector[1], dates_vector[length(dates_vector)])) + theme_for_the_plots
      plot_estim_r_deriv <- EstimatingR_List$body$modelFit$derivPlots$list_plots[[2]] +
        # labs(title = "", y = latex2exp::TeX(paste("Estimated first derivative of", "$R_t$")), x = "time") +
        coord_cartesian(ylim = c(-0.2, 0.2), xlim = c(dates_vector[1], dates_vector[length(dates_vector)])) + theme_for_the_plots
      
      EstimatingR_List$plot <- list("estim" = plot_estim_r, "deriv" = plot_estim_r_deriv)
      
      # plotInfo$max_width_estim <- grid::unit.pmax(ggplot_gtable(ggplot_build(plot_estim_r))$widths[2:3])
      # plotInfo$max_width_deriv <- grid::unit.pmax(ggplot_gtable(ggplot_build(plot_estim_r_deriv))$widths[2:3])
      
      # showNotification(latex2exp::TeX("Estimation of $R_t$ done."), type = 'message', duration = 2)
      # showNotification("Estimates of the reproductive number ready !", type = 'message', duration = 2)
      # shinyalert(title = "Estimates of the reproductive number ready !", type = "success")
      # sendSweetAlert(
      #   session = session,
      #   title = "",
      #   text = "Estimates of the reproductive number ready !",
      #   type = "success"
      # )
    }else{
      EstimatingR_List$changed <- FALSE
    }
    
  }
  
  # plot function for simple graphs
  customPlotEstimR <- function(timeAxis, fitted_values, lowerBP, upperBP, lowerBS, upperBS, th = 1,...){
    p <- ggplot(NULL, aes(x = timeAxis)) +
      geom_ribbon(aes(x = timeAxis, ymin = lowerBS, ymax = upperBS, fill = "95% simultaneous CI"), alpha = 1, linetype = 2) +
      geom_ribbon(aes(x = timeAxis, ymin = lowerBP, ymax = upperBP, fill = "95% point-wise CI"), alpha = 1, linetype = 1) +
      theme_for_the_plots +
      labs(...) +
      coord_cartesian(ylim = c(min(0, lowerBS[10:length(lowerBS)]), max(upperBS[10:length(upperBS)]))) + 
      scale_colour_manual("", 
                          breaks = c("Estimates", "95% point-wise CI", "95% simultaneous CI"),
                          values = c("#f4e04d", "#E14736", "#F2D5D1")) +
      scale_fill_manual("", 
                        breaks = c("Estimates", "95% point-wise CI", "95% simultaneous CI"),
                        values = c("#f4e04d", "#E14736", "#F2D5D1")) +
      geom_hline(yintercept = th, col = "red") +
      scale_x_date(date_breaks = "1 month", date_labels = "%b")
    
    
    if (inherits(timeAxis, "Date")){
      p <- p + scale_x_date(date_breaks = "1 month", date_labels = "%b")
    }
    
    p <- p + geom_line(aes(y = fitted_values, colour = "Estimates"), linetype = 1)   # #e09f3e
  }
  
  # ------------------------------------ Outputs ------------------------
  
  output$DeconvolutionPlot <- renderPlot({
    if (!is.null(Deconvolution_List$plot)){
      
      # Deconvolution plot
      deconv_res <- Deconvolution_List$body
      index <- deconv_res$data_original$index
      dates_vector <- GetDATA_List$body$data$Date[GetDATA_List$body$index]
      indices_used <- Deconvolution_List$body$index_deconv
      
      if (input$NightMode == 1){
        color_fitted <- "white"
      }else{
        color_fitted <- "black"
      }
      
      plot_deconvolution <- ggplot(data = NULL, aes(x = dates_vector, y = GetDATA_List$body$data$new_confirmed[GetDATA_List$body$index], colour = "Observed new cases")) + geom_line() +
        geom_ribbon(aes(ymin = Deconvolution_List$body$confidence_bands$Point_wise_CI_lower[indices_used], ymax = Deconvolution_List$body$confidence_bands$Point_wise_CI_upper[indices_used], fill = "Point-wise confidence bands", colour = NA), alpha = 0.3, fill = colors[3]) +
        geom_ribbon(aes(ymin = Deconvolution_List$body$confidence_bands$Simultaneous_lower[indices_used], ymax = Deconvolution_List$body$confidence_bands$Simultaneous_upper[indices_used], fill = "Simlutaneous confidence bands", colour = NA), alpha = 0.3, fill = colors[3]) +
        geom_line(aes(y = Deconvolution_List$body$deconvolution$fitted_values[indices_used], colour = 'Deconvolution')) +
        scale_colour_manual("",
                            breaks = c("incidence", 'Observed new cases', 'Deconvolution', "Point-wise confidence bands", "Simlutaneous confidence bands"),
                            values = c(color_fitted, color_fitted, "#f4e04d", colors[3], colors[3])) +
        scale_fill_manual("",
                          breaks = c("Point-wise confidence bands", "Simlutaneous confidence bands"),
                          values = c(colors[3], colors[3])) +
        labs(title="", x="Time", y="number of new cases") +
        theme_for_the_plots +
        coord_cartesian(ylim = c(0, 1.1 * max(GetDATA_List$body$data$new_confirmed[GetDATA_List$body$index])), xlim = c(dates_vector[1], dates_vector[length(dates_vector)])) +
        guides(color=guide_legend(override.aes=list(fill=NA))) +
        # theme(legend.position="bottom") +
        scale_x_date(date_breaks = "1 month", date_labels = "%b")    # or : date_labels = "%Y (%b)"
      
      
      # ggplotly(Deconvolution_List$plot)
      plot_deconvolution + 
        geom_vline(xintercept = input$datesStartEnd[1], color = "red") + 
        geom_vline(xintercept = input$datesStartEnd[2], color = "red") + PlotsData$theme + 
        theme(legend.position="bottom")
      
      # gtable_plot1 <- ggplot_gtable(ggplot_build(plt))
      # gtable_plot1$widths[2:3] <- grid::unit.pmax(plotInfo$max_width_deconv, plotInfo$max_width_estim, plotInfo$max_width_deriv)
      # plot(gtable_plot1)
    }else{
      ggplot() + theme_void() + theme(
        panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA)
      )
    }
  }, bg="transparent")
  
  output$PlotRt <- renderPlot({
    input$NightMode
    if (is.null(EstimatingR_List$plot)){
      ggplot() + theme_void() + theme(
        panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA)
      )
    }else{
      if (input$simple_graphs == 1){
        plot_deriv <- customPlotEstimR(GetDATA_List$body$data$Date[GetDATA_List$body$index], 
                                       EstimatingR_List$body$modelFit$derivPlots$simultaneousIntervals[[2]]$fitted,
                                       EstimatingR_List$body$modelFit$derivPlots$simultaneousIntervals[[2]]$Point_wise_CI_lower,
                                       EstimatingR_List$body$modelFit$derivPlots$simultaneousIntervals[[2]]$Point_wise_CI_upper,
                                       EstimatingR_List$body$modelFit$derivPlots$simultaneousIntervals[[2]]$Simultaneous_lower,
                                       EstimatingR_List$body$modelFit$derivPlots$simultaneousIntervals[[2]]$Simultaneous_upper,
                                       th = 0, title = "", y = latex2exp::TeX(paste("Derivatives of", "$R_t$")), x = "Time") + PlotsData$theme + 
          theme(legend.position="bottom")
        
        plot_estim <- customPlotEstimR(GetDATA_List$body$data$Date[GetDATA_List$body$index], 
                                       EstimatingR_List$body$modelFit$fitted.values,
                                       EstimatingR_List$body$modelFit$lower,
                                       EstimatingR_List$body$modelFit$upper,
                                       EstimatingR_List$body$modelFit$derivPlots$simultaneousIntervals[[1]]$Simultaneous_lower,
                                       EstimatingR_List$body$modelFit$derivPlots$simultaneousIntervals[[1]]$Simultaneous_upper,
                                       title = "", y = latex2exp::TeX(paste("Estimates of", "$R_t$")), x = "time") + PlotsData$theme + 
          theme(legend.position="bottom")
        
      }else{
        plot_deriv <- EstimatingR_List$plot$deriv + PlotsData$theme
        plot_estim <- EstimatingR_List$plot$estim + PlotsData$theme
      }
      
      if (input$Show_deriv == 1){
        plot_grid(plot_estim + theme(plot.margin = unit(c(0.5, 1, 0, 1), "cm")), 
                  plot_deriv + theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm")), 
                  nrow = 2, align = 'hv')
      }else{
        plot_grid(plot_estim + theme(plot.margin = unit(c(0.5, 1, 0, 1), "cm")), 
                  nrow = 1, align = 'hv')
      }
    }
  }, bg="transparent")
  
  output$plotGenerationInterval <- renderPlot({
    W_gam <- discrete_normal()$prob_vector
    ggplot(data = NULL, aes(x = c(1:length(W_gam)), y = W_gam)) + geom_bar(stat="identity", color="blue", fill="white") + theme_for_the_plots
  })

  # ============================ Help system ============================
  # Help system
  steps <- reactive(
    # data.frame(
    #   element=c(".sidebar-menu", ".main-header", ".sidebar-toggle", ".active", "#help", "#DeconvBox", "#EstimBox", "#bttn2", "#sidebarsettings"),
    #   intro=c(
    #     "This is a sidebar. Note that we access it with '.' instead of '#', because we track its class and not its id.",
    #     "This is a header.",
    #     "This is a button that allows to close and open the sidebar.",
    #     "This is the active element of the sidebar.",
    #     "This is a button that I added just to show the normal way to point to elements: with their id.",
    #     "Deconvolution box",
    #     "Estimating R box",
    #     "Launch simulation",
    #     "Modify data"
    #     
    #   ),
    #   position=c("right", "bottom", "bottom", "right", "top", "bottom", "top", "top", "right")
    # )
    data.frame(
      element=c(".Description",".EstimatingR", ".More", "#sidebarsettings", "#bttn2", "#DeconvBox", "#EstimBox", ".skin-blue .main-header .navbar .nav>li>a", "#Results"),
      intro=c(
        paste("Here you will find a small introduction about the project and how the reproductive number is estimated."),
        paste("In here you will be able to estimate the reproductive number and the true incidence curve."),
        paste("Here you can learn more about the methods that are used to provide the estimates. This includes a small R Markdown file to to see the code and the 
              results in action, and a small analysis of the results."),
        
        paste("These are the settings for the data. You can choose the country for which you want to estimate the reproductive number, as well as the source
              from which the data is retreived. <br/> 
              You can also select the time interval of interest."),
        
        paste("Press this button to launch the estimation. This will recover the reported incidence curve from the selected source and country, estimate 
              the true incidence curve and then estimate the reproductive number."),
        
        paste("Estimates of true incidence events are shown here. <br/> To access the settings for the deconvolution, press the", icon("cogs"), "icon.",
              "The incubation distribution is defined by the shape and scale parameters for the discretized gamma distribution. 
              Increasing the shape will increase the mean incubation period.",
              "The frequency defines the occurrence of the splines with respect to the period of time. For example, a frequency of 7 roughly 
              corresponds to a deconvolution process with a new parameter every 7 days.",
              "The regularization parameter incluences the smoothness of the curve. The bigger the regularization is, the smoother the 
              estimated incidence curve is."),
        
        paste("Estimates of the reproductive number are shown in this panel. <br/> To access the settings, press the", icon("cogs"), "icon.",
              "The distribution for the generation intervals is defined by the mean and the variance parameters.",
              "As with the deconvolution step, the frequency defines the occurrence of the splines with respect to the period of time.",
              "One can choose between a Poisson, quasi-Poisson and negative binomial model. The negative binomial model has proven 
              to give the best results and allow for over-dispersion."),
        
        paste("Press the", icon("cogs"), "icon to access the super settings. <br/> 
              You can directly specify the degree of the spline basis instead of the frequency, and by pressing the 'Choose best' button, 
              the method searches for the best basis dimension and regularization terms. This process may take a few minutes, so please do not 
              log out of the page.",
              "The simple graph option allows for more interpretable plots and removes sampled splines.",
              "Finally, you can specify how the uncertainty of the deconvolution step is used to estimate the reproductive number. If 'No weight' is chosen,
              the estimate of R does not take into account the uncertainty of the deconvolution."),
        
        paste("Here are shown the estimated incidence events (I) and reproductive number (R) for a specific date. The 95% point wise (95% CI (p)) and simultaneous (95% CI (s))
              can also be found in the table.")
      ),
      position=c("right","right","right","right", "right", "top", "top", "bottom", "right")
    )
  )
  
  observeEvent(input$help,{
    introjs(session,
            options = list(steps=steps(),
                           "nextLabel"="Next",
                           "prevLabel"="Previous",
                           "skipLabel"="Skip"
            )
            # events = list("oncomplete"=I('alert("Done")'))
    )
  }
  )
  
  # ============================ Description of the method ============================
  # Description of the method
  output$detailsMethod <- renderUI({
    withMathJax(
      HTML('The evolution of an epidemic can be tracked using its reproductive number, denoted \\(R_t\\), 
                        which represents the average number of secondary cases caused by the infection of an individual at time \\(t\\). 
                        In order for sanitary measures to effectively contain the spread of the epidemic, a good estimate of 
                        the reproduction number is required, as well as the uncertainties for the estimates. Monitoring \\(R_t\\) can 
                        help detect changes in the transmission of the disease. Estimates of \\(R_t\\) are mostly derived from reported 
                        quantities that are daily available, such as the incidence, death and recovered case curves. 
                        As a result, the symptom onset and reporting delay distributions must be taken into account in order 
                        to avoid inaccuracies and lagging errors. <br/> <br/>'),
      
      h4(HTML("Deconvolution <br/>")),
      HTML('We will use the number of daily new cases to estimate Rt, using generalized additive models with smoothing splines. 
                        In the following, we denote the reported curve (also called the observed cases) by \\((Y_t)_{t=1}^T\\) 
                        and the true incidence curve by \\((It)_{t=1}^T\\). These curves are assumed to be discretized in time, and are monitored and estimated on a daily basis. 
                        These curves are not necessarily equal, as infected individuals are usually not reported on the same day as the infection occurred. 
                        As such, the reported number of cases at some point in time might not be representative of the number of individuals that were infected, 
                        due to delay between the infection event and the symptom onset. Thus, we need to deconvolve the incidence curve before estimating \\(R_t\\), that is, 
                        estimate the true infection events from the observed ones, given a distribution for the incubation period. This step is carried out thanks to the use of a generalized additive model. <br/>
                        We assume that the observed event \\(Y_{t+1}\\) follows, conditionally on the incidence events \\(I_1, \\ldots, I_t\\), 
                        a Poisson distribution with a weighted sum of these incidence events. More formally, 
                        we assume that \\(Y_{t+1} | I_1, \\ldots, I_{t} \\sim \\text{Pois}(\\mu_{t+1})\\), with \\(\\mu_{t+1} = \\sum_{s=1}^K {w_s I_{t+1-s}}\\). 
                        To allow for overdispersion, we use a negative binomial distribution instead of the Poisson one. 
                        Also, instead of estimating directly the \\((It)_{t=1}^T\\) (which would result in one parameter for each day in our original model), 
                        we assume that these true incidence events are smooth function of time, and can be modelled with spline functions. <br/> <br/> '),
      
      h4(HTML("Estimating \\(R_t\\) <br/>")),
      HTML("Analogously to the deconvolution step, we estimate \\(R_t\\) with spline functions. 
                        In our model, we assume that the mean number of newly infected individuals at time \\(t\\) is a 
                        weighted average of the new infected individuals at previous time step, multiplied by some quantity \\(R_t\\) 
                        that varies across time. More formally, let \\(I_t\\) be the incidence curve for time step \\(t = 1, \\ldots, T\\), 
                        and let \\(\\{\\tilde{w}_j\\}_{j=1}^{K}\\) be the infectivity profile, which is different from the incubation period 
                        that we introduced earlier. We suppose that each \\(I_t\\) is Poisson distributed, with mean 
                        $$ \\mu_t = R_t \\sum_{j=1}^{K}\\tilde{w}_j I_{t-j} = \\exp \\left[ \\log(R_t) + \\log \\left (\\sum_{j=1}^{K}\\tilde{w}_j I_{t-j} \\right ) \\right ], $$
                        where the generation interval is assumed to be known. As such, the second term in the exponential (the offset) can 
                        be explicitly computed, and the only unknown is the log of the reproductive number. In order to estimate this, we suppose 
                        that \\(\\log(R_t)\\) is smoothed over the time, that is, $$ \\log(R_t) = \\sum_{i = 1}^{q_1}\\beta_j B_j(t), $$
                        where \\(B = \\{B_1,\\ldots,B_{q_1}\\}\\) is a spline basis of dimension \\(q_1\\), 
                        defined over the interval \\([0, T]\\). <br/> <br/><br/> <br/><br/> <br/>")
   
      # textFunctionTheme(dark_theme,'Let \\( y_1, \\ldots, y_n \\) be the observed events (for example the new number of cases per day). The Cori et. al estimation for \\( R_t \\) 
      #                   (the reproductive number) is as follows: $$ R_t = \\frac{y_t}{\\sum{w_s y_{t-s}}}. $$ However, it assumes that the observed events \\( y_t \\)
      #                   correspond to the true number of infections incident, which is not the case in practice. '),
      # textFunctionTheme(dark_theme,'Thus, the observations are transformed to get a new set of events \\( \\hat{x}_1, \\ldots, \\hat{x}_n \\), such that for each date 
      #                   \\( y_{it}, i = 1, \\ldots, y_t, \\), we have \\( \\hat{x}_{it} = y_{it} - d \\), where \\( d \\backsim DiscretizedGamma(scale, shape, n).\\)'),
      # textFunctionTheme(dark_theme,'Another way to estimate the reproductive number \\( R_t \\) is to assume that given \\( \\hat{x}_{1}, \\ldots, \\hat{x}_{t-1} \\), 
      #                   we have that \\( \\hat{x}_{t} \\) follows a Poisson distribution with mean \\( \\mu_{t} \\), with
      #                   $$ \\mu_{t} = R_t \\sum{w_s \\hat{x}_{t-s}} = \\exp( \\log(R_t) + \\log(\\sum{w_s \\hat{x}_{t-s}})) $$ where \\( w_s \\) are the generation 
      #                   intervals, and follow here a discretized gamma distribution as well.'),
      # textFunctionTheme(dark_theme,'Also, \\( \\log(R_t) \\) is assumed to be a smoothed function of the time. There is also the possibility of adding \'lockdown\' effects,
      #                   which are added in the linear predictor of the model with \\( I_{date >= lockdown}\\).'),
      # textFunctionTheme(dark_theme,'Using a Poisson family for the model fitting implies that \\( \\mathbb{E}(\\hat{x}_{t} | \\hat{x}_{t-1}, \\ldots) = \\mu_{t} \\), and 
      #                   \\( \\mathbb{V}(\\hat{x}_{t} | \\hat{x}_{t-1}, \\ldots) = \\mu_{t} \\), so the variance is a linear function of the mean.'),
      # textFunctionTheme(dark_theme,'One may relax this assumption by using a Quasi-poisson model, for which we have \\( \\mathbb{V}(\\hat{x}_{t} | \\hat{x}_{t-1}, \\ldots) = \\phi \\mu_{t} \\)
      #                   with $$ \\hat{\\phi} = \\frac{1}{n-k} \\sum_{i=1}^{n}{\\frac{(\\hat{x}_i - \\hat{\\mu}_i)^2}{\\hat{\\mu}_i}}  .$$'),
      # textFunctionTheme(dark_theme,'Another possibility is the use a negative binomial model, for which the density writes as
      #                   $$ f(y; k, \\mu) = \\frac{\\Gamma(y+k)}{\\Gamma(k)y!} \\left( \\frac{k}{\\mu + k} \\right)^k \\left(1-\\frac{k}{\\mu + k} \\right)^y $$
      #                   and \\( \\mathbb{E}(\\hat{x}_{t} | \\hat{x}_{t-1}, \\ldots) = \\mu_{t} \\), \\( \\mathbb{V}(\\hat{x}_{t} | \\hat{x}_{t-1}, \\ldots) = \\mu_{t}(\\mu_t / k +1) \\), 
      #                   which implies that the variance is now a quadratic function of the mean.')
      )
  })
  
  
  
  # ============================ Download data ======================
  gatherDataAndResults <- function(){
    dates_vector <- GetDATA_List$body$data$Date[GetDATA_List$body$index]
    Deconv_fitted <- Deconvolution_List$body$deconvolution$fitted_values[Deconvolution_List$body$index_deconv]
    Deconv_simultaneous <- list("lower" = Deconvolution_List$body$confidence_bands$Simultaneous_lower[Deconvolution_List$body$index_deconv],
                                 "upper" = Deconvolution_List$body$confidence_bands$Simultaneous_upper[Deconvolution_List$body$index_deconv])
    Deconv_point_wise <- list("lower" = Deconvolution_List$body$confidence_bands$Point_wise_CI_lower[Deconvolution_List$body$index_deconv],
                               "upper" = Deconvolution_List$body$confidence_bands$Point_wise_CI_upper[Deconvolution_List$body$index_deconv])
    R_fitted <- EstimatingR_List$body$modelFit$fitted.values
    R_simultaneous <- list("lower" = EstimatingR_List$body$modelFit$derivPlots$simultaneousIntervals[[1]]$Simultaneous_lower,
                            "upper" = EstimatingR_List$body$modelFit$derivPlots$simultaneousIntervals[[1]]$Simultaneous_upper)
    R_point_wise <- list("lower" = EstimatingR_List$body$modelFit$derivPlots$simultaneousIntervals[[1]]$Point_wise_CI_lower,
                          "upper" = EstimatingR_List$body$modelFit$derivPlots$simultaneousIntervals[[1]]$Point_wise_CI_upper)
    results <- data.frame("dates" = dates_vector,
                          "Estimate incidence" = Deconv_fitted,
                          "Estimate incidence 95% point wise lower" = Deconv_point_wise$lower,
                          "Estimate incidence 95% point wise upper" = Deconv_point_wise$upper,
                          "Estimate incidence 95% simultaneous lower" = Deconv_simultaneous$lower,
                          "Estimate incidence 95% simultaneous upper" = Deconv_simultaneous$upper,
                          "Estimate incidence" = R_fitted,
                          "Estimate incidence 95% point wise lower" = R_point_wise$lower,
                          "Estimate incidence 95% point wise upper" = R_point_wise$upper,
                          "Estimate incidence 95% simultaneous lower" = R_simultaneous$lower,
                          "Estimate incidence 95% simultaneous upper" = R_simultaneous$upper)
    return(results)
  }
  
  output$downloadData <- downloadHandler(
    filename = function() { 
      paste("results-", Sys.Date(), ".csv", sep="")
    }
    ,content = function(file) {
      if (!is.null(Deconvolution_List$body)){
        write.csv(gatherDataAndResults(), file)
      }
    }
  )
  
})


# ------------ App -------------
shinyApp(ui = ui, server = server)




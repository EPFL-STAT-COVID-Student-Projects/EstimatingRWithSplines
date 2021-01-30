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
# -------------------------- Auxiliary functions --------------------------------

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



# -------------------------- 1) Roll back ---------------------------------------

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


# Example of use
if(FALSE){
  source("/Users/Antoine/Desktop/Projet Semestre 2020/CODE/simulation.R")
  
  sim_list <- simulation()
  df <- sim_list$sim_df
  
  data_init <- df$incidence[2:(length(df$incidence)-1)]
  
  w_conv<- discrete_gamma(n = 100, normalize = TRUE, shape = 10, scale = 0.5)
  
  data_convolved <- convolveData(w_conv, data_init)
  data_country_roll_back <- rollBack(data_convolved, w_conv)
  data_country_roll_back <- data_country_roll_back[(100+1):length(data_country_roll_back)]
  
  results <- rollBackCI(data_convolved, w_conv, 100)
  
  plot(data_init, type = 'l')
  lines(data_convolved, col = "blue")
  lines(results$mean, col = "red")
  lines(results$q2.5, col = "red")
  lines(results$q97.5, col = "red")
}

# -------------------------- 2) Optimization ------------------------------------

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

# Example of use
if(FALSE){
  source("/Users/Antoine/Desktop/Projet Semestre 2020/CODE/simulation.R")
  
  sim_list <- simulation()
  df <- sim_list$sim_df
  data_init <- df$incidence[2:(length(df$incidence)-1)]
  w_conv<- discrete_gamma(n = 100, normalize = TRUE, shape = 10, scale = 0.5)
  data_convolved <- convolveData(w_conv, data_init)
  
  results <- optimPoissonDeconvolution(data_convolved, w_conv, smooth_observations = TRUE, maxIter = 100, method = "CG", data_start = NULL)
  
  # print(diag(solve(-results$result$hessian)))
  # results$result
  # 
  # hess <- log_lik_hessian(results$fitted_values, data_convolved, w_conv)
  # print(diag(solve(hess[1:(dim(hess)[1]-1), 1:(dim(hess)[1]-1)])))
  
  plot(data_init, type = 'l')
  lines(data_convolved, col = "blue")
  lines(results$fitted_values, col = "red")
  lines(rollapply(results$fitted_values, width=3, FUN=function(x) mean(x, na.rm=TRUE), by=1, 
                  by.column=TRUE, partial=TRUE, fill=NA, align="center"), col = "green")
  
}


# -------------------------- 3) IWLS --------------------------

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



# --------------------------
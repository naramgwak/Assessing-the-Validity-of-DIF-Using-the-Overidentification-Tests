set.seed(1234)

library(boot)
library(dplyr)
library(ggplot2)

sim_dat <- function(a = rep(1, 10), # discrimination
                    b = c(-2.5, -2, -1.5, -1, -.5, 0, .5 , 1  , 1.5, 2), # difficulty
                    uTheta = 1, gTheta = 1, zTheta = 1, 
                    uY = 1, gY = 0, zY = 0,
                    n = 100){
  U <- rbinom(n, 1, .5)  
  G <- rbinom(n, 1, .5)   
  Z <- rbinom(n, 1, .5)   # a valid iv -> there is no item bias
  
  theta <- rnorm(n, 0, 1)   # true ability
  
  # U, G, Z  -> theta
  if (uTheta == 1){
    theta <- ifelse(U == 1, theta + 1, theta)
  }
  
  if (gTheta == 1){
    theta <- ifelse(G == 1, theta + 1, theta)
  }
  
  if (zTheta == 1){
    theta <- ifelse(Z == 1, theta + 1, theta)
  }
  
  # p_ij
  pij <- matrix(NA, nrow = n, ncol = length(a))
  for (j in 1:length(a)){
    for (i in 1:n){
      pij[i, j] <- (exp(a[j] * (theta[i] - b[j])) / (1 + exp(a[j] * (theta[i] - b[j]))))
    }
  }
  
  # U, G -> Y (by changing item difficulty of #10; more easier)
  if (uY == 1){
    pij[, 10] <- ifelse(U == 1, (exp(a[10] * (theta - (b[10] - 1))) / (1 + exp(a[10] * (theta - (b[10] - 1))))), pij[, 10])
  }
  
  if (gY != 0){  # if gY is not 0, adjust the difficulty based on gY value
    pij[, 10] <- ifelse(G == 1, (exp(a[10] * (theta - (b[10] - gY))) / (1 + exp(a[10] * (theta - (b[10] - gY))))), pij[, 10])
  } 
  
  Y <- matrix(NA, nrow = n, ncol = length(a))
  for (j in 1:length(a)){
    for (i in 1:n){
      Y[i, j] <- rbinom(1, 1, pij[i, j])
    }
  }
  score <- rowSums(Y)
  
  data <- data.frame(U = U, G = G, Z = Z, theta = theta, score = score)
  Y_df <- as.data.frame(Y)
  colnames(Y_df) <- paste0("Y", 1:ncol(Y_df))
  data <- cbind(data, Y_df)
  return(data)
}

sim_mod <- function(data, out){
  Y_col <- data[, out + 5]
  U <- data$U
  G <- data$G
  Z <- data$Z
  theta <- data$theta
  
  m1_g <- lm(theta ~ G, data = data)
  m1_z <- lm(theta ~ Z, data = data)
  
  summary(m1_g)
  
  m2_g <- lm(Y_col ~ fitted(m1_g))
  m2_z <- lm(Y_col ~ fitted(m1_z))
  
  coef_g <- coef(m2_g)['fitted(m1_g)']
  coef_z <- coef(m2_z)['fitted(m1_z)']
  
  difm <- coef_g - coef_z
  
  return(difm)
}

boot_sim <- function(data, out, n.boot) {
  boot_func <- function(data, indices) {
    d <- data[indices, ]
    mod_res <- sim_mod(d, out)
    return(mod_res)
  }
  
  boot_com <- boot(data, boot_func, R = n.boot)
  boot_ci <- boot.ci(boot_com, type = "perc")$percent[4:5]
  zero_in_ci <- boot_ci[1] < 0 & boot_ci[2] > 0
  return(list(boot_com = boot_com, boot_ci = boot_ci, zero = zero_in_ci))
}

calculate_power <- function(n_reps, n.boot, sample_sizes, gY_values) {
  power_results <- data.frame()
  
  for (n in sample_sizes) {
    for (gY in gY_values) {
      detection_count <- 0
      
      pb <- txtProgressBar(min = 0, max = n_reps, style = 3)
      
      for (i in 1:n_reps) {
        data <- sim_dat(gTheta = 1, gY = gY, uTheta = 1, uY = 1, n = n)
        boot_res <- boot_sim(data, 10, n.boot = n.boot)
        
        if (!boot_res$zero) {
          detection_count <- detection_count + 1
        }
        
        setTxtProgressBar(pb, i)
      }
      close(pb)
      
      power <- detection_count / n_reps
      power_results <- rbind(power_results, data.frame(sample_size = n, gY = gY, power = power))
    }
  }
  
  return(power_results)
}

sample_sizes <- c(300, 500, 1000)
gY_values <- seq(0, 2, by = 0.5)

n_reps <- 100
n.boot <- 1000

power_results_df <- calculate_power(n_reps, n.boot, sample_sizes, gY_values)

print(power_results_df)

power_results_long <- power_results_df %>%
  pivot_longer(cols = starts_with("power"), names_to = "power", values_to = "value") %>%
  mutate(gY = as.numeric(gY))

ggplot(power_results_long, aes(x = sample_size, y = value, color = as.factor(gY), group = gY)) +
  geom_line() +
  geom_point() +
  labs(x = "Sample Size",
       y = "Power",
       color = "gY") +
  theme_minimal() +
  theme(legend.position = "right")

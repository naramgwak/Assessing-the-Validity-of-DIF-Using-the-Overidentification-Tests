set.seed(1234)

library(boot)
library(writexl)

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
  
  if (gY == 1){
    pij[, 10] <- ifelse(G == 1, (exp(a[10] * (theta - (b[10] - 1))) / (1 + exp(a[10] * (theta - (b[10] - 1))))), pij[, 10])
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
  G1 <- data$G1
  G2 <- data$G2
  theta <- data$theta
  
  m1_g1 <- lm(theta ~ G1, data = data)
  m1_g2 <- lm(theta ~ G2, data = data)
  
  summary(m1_g1)
  
  m2_g1 <- lm(Y_col ~ fitted(m1_g1))
  m2_g2 <- lm(Y_col ~ fitted(m1_g2))
  
  coef_g1 <- coef(m2_g1)['fitted(m1_g1)']
  coef_g2 <- coef(m2_g2)['fitted(m1_g2)']
  
  difm <- coef_g1 - coef_g2
  
  return(difm)
}

boot_sim <- function(data, out, n.boot) {
  boot_func <- function(data, indices) {
    d <- data[indices, ]
    mod_res <- sim_mod(d, out)
    return(mod_res)
  }
  
  boot_com <- boot(data, boot_func, R = n.boot)
  boot_ci <- boot.ci(boot_com, type = "perc")$percent[4:5] # bca, NULL
  zero_in_ci <- boot_ci[1] < 0 & boot_ci[2] > 0
  return(list(boot_com = boot_com, boot_ci = boot_ci, zero = zero_in_ci))
}

n.boot <- 1000
n_reps <- c(100, 500, 1000)

sDIF_results_df <- data.frame() # g1 X g2 X (no DIF), both are ivs
teDIF_results_df <- data.frame() # g1 O g2 X (g1 DIF), g2 is iv
tbDIF_results_df <- data.frame() # g1 O g2 O (g1 & g2 DIF), both are not ivs

for (n_rep in n_reps) {
  pb <- txtProgressBar(min = 0, max = n_rep, style = 3)
  
  sDIF_mod_res <- numeric(n_rep)
  sDIF_ci_lower <- numeric(n_rep)
  sDIF_ci_upper <- numeric(n_rep)
  sDIF_zero <- numeric(n_rep)
  
  teDIF_mod_res <- numeric(n_rep)
  teDIF_ci_lower <- numeric(n_rep)
  teDIF_ci_upper <- numeric(n_rep)
  teDIF_zero <- numeric(n_rep)
  
  tbDIF_mod_res <- numeric(n_rep)
  tbDIF_ci_lower <- numeric(n_rep)
  tbDIF_ci_upper <- numeric(n_rep)
  tbDIF_zero <- numeric(n_rep)
  
  for (i in 1:n_rep) {
    # sDIF
    data <- sim_dat(g1Theta = 1, g1Y = 0, g2Theta = 1, g2Y = 0, uTheta = 1, uY = 1, n = 1000)
    mod_res <- sim_mod(data, 10)
    boot_res <- boot_sim(data, 10, n.boot = n.boot)
    
    sDIF_mod_res[i] <- mod_res
    sDIF_ci_lower[i] <- boot_res$boot_ci[1]
    sDIF_ci_upper[i] <- boot_res$boot_ci[2]
    sDIF_zero[i] <- boot_res$zero
    
    # teDIF
    data <- sim_dat(g1Theta = 1, g1Y = 1, g2Theta = 1, g2Y = 0, uTheta = 1, uY = 1, n = 1000)
    mod_res <- sim_mod(data, 10)
    boot_res <- boot_sim(data, 10, n.boot = n.boot)
    
    teDIF_mod_res[i] <- mod_res
    teDIF_ci_lower[i] <- boot_res$boot_ci[1]
    teDIF_ci_upper[i] <- boot_res$boot_ci[2]
    teDIF_zero[i] <- boot_res$zero
    
    # tbDIF
    data <- sim_dat(g1Theta = 1, g1Y = 1, g2Theta = 1, g2Y = 1, uTheta = 1, uY = 1, n = 1000)
    mod_res <- sim_mod(data, 10)
    boot_res <- boot_sim(data, 10, n.boot = n.boot)
    
    tbDIF_mod_res[i] <- mod_res
    tbDIF_ci_lower[i] <- boot_res$boot_ci[1]
    tbDIF_ci_upper[i] <- boot_res$boot_ci[2]
    tbDIF_zero[i] <- boot_res$zero
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  sDIF_results_df <- rbind(sDIF_results_df, data.frame(
    n_rep = n_rep,
    mod_res_mean = mean(sDIF_mod_res),
    ci_lower_mean = mean(sDIF_ci_lower),
    ci_upper_mean = mean(sDIF_ci_upper),
    zero_mean = mean(sDIF_zero)
  ))
  
  teDIF_results_df <- rbind(teDIF_results_df, data.frame(
    n_rep = n_rep,
    mod_res_mean = mean(teDIF_mod_res),
    ci_lower_mean = mean(teDIF_ci_lower),
    ci_upper_mean = mean(teDIF_ci_upper),
    zero_mean = mean(teDIF_zero)
  ))
  
  tbDIF_results_df <- rbind(tbDIF_results_df, data.frame(
    n_rep = n_rep,
    mod_res_mean = mean(tbDIF_mod_res),
    ci_lower_mean = mean(tbDIF_ci_lower),
    ci_upper_mean = mean(tbDIF_ci_upper),
    zero_mean = mean(tbDIF_zero)
  ))
}

write_xlsx(list(sDIF_results = sDIF_results_df, teDIF_results = teDIF_results_df, tbDIF_results = tbDIF_results_df), "scenario2.xlsx")
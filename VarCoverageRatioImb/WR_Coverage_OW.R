library(doMC)
library(MASS)
library(nnet)
#library(VGAM)
library(ggplot2)
library(parallel)
library(foreach)
library(doParallel)
#numCores <- detectCores()

#############################################################
################# Variance Coverage of OW ###################
#############################################################

args<-commandArgs(trailingOnly = TRUE)
k<-as.integer(args[1])
if (is.na(k)) k <- 1
paste("Scenario:",k)

numCores<-as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK",8))
if (is.na(numCores)) numCores<-8
registerDoMC(cores=numCores)

# define scenarios

scenarios <- read.table("WR_Params.txt", header=TRUE, sep="")
scenarios <- subset(scenarios, scenario == k)


scenario <- k

sds <- c(scenarios$sd_x1, scenarios$sd_x2, scenarios$sd_x3)
bern_param <- c(scenarios$param_x4, scenarios$param_x5, scenarios$param_x6)


trt_eff1 <- scenarios$trt_eff1

bi_trt <- c(1,-1,1,-1,1,-1)*scenarios$bi_trt
bi_ctrl <- c(1,-1,1,-1,1,-1)*scenarios$bi_ctrl

b01 <- 1
b02 <- 0.05
step_size <- scenarios$step_size #(sample size increase step)
re_count <- 3 #sample size count
sim_num <- 1000 #simulation iteration
# re_count <- scenarios$re_count #sample size count
# sim_num <- scenarios$sim_num #simulation iteration


start_time <- Sys.time()

##################################Simulation##################################


#vector of outcomes as a factor with ordered levels
outcomes_3lvl <- factor(c("first", "second", "third"), 
                        levels = c("first", "second", "third"), 
                        ordered = TRUE)

#Order: first < second < third, the larger the better

inv_logit <- function(logit) exp(logit)/(1 + exp(logit))

WP_trt_OW_coverage <- numeric(re_count)
WP_ctrl_OW_coverage <- numeric(re_count)
WR_OW_coverage <- numeric(re_count)

WP_trt_OW_var_ratio <- numeric(re_count)
WP_ctrl_OW_var_ratio <- numeric(re_count)
WR_OW_var_ratio <- numeric(re_count)

WP_trt_OW_theoretical_var_mean <- numeric(re_count)
WP_trt_OW_empirical_var <- numeric(re_count)

WP_ctrl_OW_theoretical_var_mean <- numeric(re_count)
WP_ctrl_OW_empirical_var <- numeric(re_count)

WR_OW_theoretical_var_mean <- numeric(re_count)
WR_OW_empirical_var <- numeric(re_count)

WP_trt_sim_true <- 0.2959808
WP_ctrl_sim_true <- 0.2048627
WR_sim_true <- 1.45124



var_ratio = function(uR,uS,sigR2,sigS2,cov_RS){
  WR_appro_var = (uR^2/uS^2)*(sigR2/uR^2-2*cov_RS/(uR*uS)+sigS2/uS^2)
  return(WR_appro_var)
}




sample_size_list <- c(200,300,400)




for(size_count in 1:re_count){
  
WP_trt_sim_OW <- numeric(sim_num)
WP_ctrl_sim_OW <- numeric(sim_num)
WR_sim_OW <- numeric(sim_num)

theory_var_OW_trt <- numeric(sim_num)
theory_var_OW_ctrl <- numeric(sim_num)
cov_trt_ctrl <- numeric(sim_num)
theory_var_OW_WR <- numeric(sim_num)

for (count_temp in 1:sim_num){

  set.seed(count_temp)
  n_count <- sample_size_list[size_count]

# covariates
x1 <- rnorm(n_count, mean = 1, sd = sds[1])
x2 <- rnorm(n_count, mean = 0.9, sd = sds[2])
x3 <- rnorm(n_count, mean = 0.8, sd = sds[3])

x4 <- rbinom(n_count,1,bern_param[1])
x5 <- rbinom(n_count,1,bern_param[2])
x6 <- rbinom(n_count,1,bern_param[3])

df_cov <- data.frame(x1, x2, x3, x4, x5, x6)
treatment_assignment <- rbinom(n_count, 1, 0.7)

trt_cov <- df_cov[treatment_assignment == 1, ]
ctrl_cov <- df_cov[treatment_assignment == 0, ]

trt_cov_quad <- trt_cov^2
ctrl_cov_quad <- ctrl_cov^2

bi_trt_quad <- 2*bi_trt
bi_ctrl_quad <- 2*bi_ctrl


combined_trt_quad <- cbind(trt_cov, trt_cov_quad)
combined_ctrl_quad <- cbind(ctrl_cov, ctrl_cov_quad)

logodds1_trt <- b01 + as.matrix(combined_trt_quad) %*% c(bi_trt, bi_trt_quad) + trt_eff1
logodds2_trt <- b02 + as.matrix(combined_trt_quad) %*% c(bi_trt, bi_trt_quad) + trt_eff1

logodds1_ctrl <- b01 + as.matrix(combined_ctrl_quad) %*% c(bi_ctrl, bi_ctrl_quad)
logodds2_ctrl <- b02 + as.matrix(combined_ctrl_quad) %*% c(bi_ctrl, bi_ctrl_quad)


## Probability Trt

prob_2to3_trt <- inv_logit(logodds1_trt)
prob_3_trt <- inv_logit(logodds2_trt)
prob_1_trt <- 1 - prob_2to3_trt
prob_2_trt <- prob_2to3_trt - prob_3_trt



## Probability Ctrl

prob_2to3_ctrl <- inv_logit(logodds1_ctrl)
prob_3_ctrl <- inv_logit(logodds2_ctrl)
prob_1_ctrl <- 1 - prob_2to3_ctrl
prob_2_ctrl <- prob_2to3_ctrl - prob_3_ctrl

#generate random outcomes
outcomes_trt <- c()
for (i in 1:nrow(trt_cov)) {
  outcomes_trt[i] <- sample(
    outcomes_3lvl, 
    size = 1,
    prob = c(prob_1_trt[i], prob_2_trt[i], prob_3_trt[i])
  )
}



outcomes_ctrl <- c()
for (i in 1:nrow(ctrl_cov)) {
  outcomes_ctrl[i] <- sample(
    outcomes_3lvl, 
    size = 1,
    prob = c(prob_1_ctrl[i], prob_2_ctrl[i], prob_3_ctrl[i])
  )
}







##########################OW########################################

df_trt <- data.frame(outcomes_trt)
df_ctrl <- data.frame(outcomes_ctrl)

colnames(df_trt) <- "outcomes_comb"
colnames(df_ctrl) <- "outcomes_comb"

df_comb <- rbind(df_trt, df_ctrl)
df_comb$treatment <- c(rep(1,nrow(trt_cov)), rep(0,nrow(ctrl_cov)))

df_comb$x1 <- c(x1[treatment_assignment == 1], x1[treatment_assignment == 0])
df_comb$x2 <- c(x2[treatment_assignment == 1], x2[treatment_assignment == 0])
df_comb$x3 <- c(x3[treatment_assignment == 1], x3[treatment_assignment == 0])
df_comb$x4 <- c(x4[treatment_assignment == 1], x4[treatment_assignment == 0])
df_comb$x5 <- c(x5[treatment_assignment == 1], x5[treatment_assignment == 0])
df_comb$x6 <- c(x6[treatment_assignment == 1], x6[treatment_assignment == 0])

# Calculate propensity scores
PropScore <- glm(treatment~x1 + x2 + x3 + x4 + x5 + x6, 
                 data = df_comb, family=binomial)
pi_func <- fitted(PropScore)
#pi_func <- rep(0.5, nrow(df_comb))



compare_rows_OW <- function(i) {
  g_hfunc_OW1_vec <- numeric(n_count)
  g_hfunc_OW2_vec <- numeric(n_count)
  
  for (j in 1:n_count) {
    if (df_comb[i,1] > df_comb[j,1]) {
      g_hfunc_OW1_vec[j] = (df_comb[i,"treatment"] * (1-df_comb[j,"treatment"]) * (1-pi_func[i]) * pi_func[j])
      g_hfunc_OW2_vec[j] = 0
    } else if (df_comb[i,1] < df_comb[j,1]) {
      g_hfunc_OW1_vec[j] = 0
      g_hfunc_OW2_vec[j] = (df_comb[i,"treatment"] * (1-df_comb[j,"treatment"]) * (1-pi_func[i]) * pi_func[j])
    }
  }
  
  return(list(OW1 = g_hfunc_OW1_vec, OW2 = g_hfunc_OW2_vec))
}

# Parallelize the outer loop
results_OW <- mclapply(1:n_count, compare_rows_OW, mc.cores = numCores)

# Aggregate results
g_hfunc_OW1 <- do.call(rbind, lapply(results_OW, function(x) x$OW1))
g_hfunc_OW2 <- do.call(rbind, lapply(results_OW, function(x) x$OW2))

# Calculate tau values
tau1_OW <- sum(as.vector(g_hfunc_OW1), na.rm = TRUE) /
  (sum(df_comb[,"treatment"]*(1-pi_func), na.rm = TRUE) * sum((1-df_comb[,"treatment"])*pi_func, na.rm = TRUE))
tau2_OW <- sum(as.vector(g_hfunc_OW2), na.rm = TRUE) /
  (sum(df_comb[,"treatment"]*(1-pi_func), na.rm = TRUE) * sum((1-df_comb[,"treatment"])*pi_func, na.rm = TRUE))

# Update WD_sim_OW and WR_sim_OW
WP_trt_sim_OW[count_temp] <- tau1_OW
WP_ctrl_sim_OW[count_temp] <- tau2_OW
WR_sim_OW[count_temp] <- tau1_OW / tau2_OW





###########OW Theoretical Variance###########


###g1_OW_trt###

n <- nrow(df_comb)
adjusted_sum_i_neq_j <- sum(sapply(1:n, function(j) {sum(df_comb$treatment[-j]*(1-pi_func[-j])) * 
    (1 - df_comb$treatment[j])*pi_func[j]}))
g1_constant_denominator_OW <- (1/(n*(n-1))) * adjusted_sum_i_neq_j

g1_OW_trt_fun_parallel <- function(i, outcomes, treatment, pi_func) {
  pairwise_comparisons_trt <- (1/2) * (
    (((1 - pi_func[i]) * pi_func) * treatment[i] * (1 - treatment) * as.numeric(outcomes[i] > outcomes)) +
      (((1 - pi_func) * pi_func[i]) * treatment * (1 - treatment[i]) * as.numeric(outcomes > outcomes[i]))
  )
  
  #pairwise_comparisons_trt[i] <- 0  # Avoid self-comparison
  sum_g1_OW_trt <- sum(pairwise_comparisons_trt) / g1_constant_denominator_OW
  return((1/(n-1)) * sum_g1_OW_trt)
}

results_g1_OW_trt <- mclapply(1:n, g1_OW_trt_fun_parallel, df_comb$outcomes_comb, 
                              df_comb$treatment, pi_func, mc.cores = numCores - 1)


g1_OW_trt <- unlist(results_g1_OW_trt)




###g1_OW_ctrl###

g1_OW_ctrl_fun_parallel <- function(i, outcomes, treatment, pi_func) {
  pairwise_comparisons_ctrl <- (1/2) * (
    (((1 - pi_func[i]) * pi_func) * treatment[i] * (1 - treatment) * as.numeric(outcomes[i] < outcomes)) +
      (((1 - pi_func) * pi_func[i]) * treatment * (1 - treatment[i]) * as.numeric(outcomes < outcomes[i]))
  )
  
  #pairwise_comparisons_ctrl[i] <- 0  # Avoid self-comparison
  sum_g1_OW_ctrl <- sum(pairwise_comparisons_ctrl) / g1_constant_denominator_OW
  return((1/(n-1)) * sum_g1_OW_ctrl)
}

results_g1_OW_ctrl <- mclapply(1:n, g1_OW_ctrl_fun_parallel, df_comb$outcomes_comb, 
                               df_comb$treatment, pi_func, mc.cores = numCores - 1)


g1_OW_ctrl <- unlist(results_g1_OW_ctrl)







###OW A Treatment### (Dim: 6*1)


dBeta_A_trt_OW <- function(i, j, pi_func, df_comb) {
  n <- nrow(df_comb)
  pi_i = pi_func[i]
  pi_j = pi_func[j]
  Z_i = df_comb[i, "treatment"]
  Z_j = df_comb[j, "treatment"]
  Y_i = df_comb[i, "outcomes_comb"]
  Y_j = df_comb[j, "outcomes_comb"]
  X_i = df_comb[i, c("x1", "x2", "x3", "x4", "x5", "x6")]
  X_j = df_comb[j, c("x1", "x2", "x3", "x4", "x5", "x6")]
  
  comparison_i_j = as.numeric(Y_i > Y_j)  # Convert logical TRUE/FALSE to 1/0
  comparison_j_i = as.numeric(Y_j > Y_i)  # Convert logical TRUE/FALSE to 1/0
  
  A = 0.5 * (
    ((1-pi_j) * X_j - pi_i * X_i) * (1-pi_i) * pi_j * Z_i * (1 - Z_j) * comparison_i_j  +
      ((1-pi_i) * X_i - pi_j * X_j) * (1-pi_j) * pi_i * Z_j * (1 - Z_i) * comparison_j_i
  )
  #A = sum(A)
  
  return(A)
  
}




###OW A Control### (Dim: 6*1)


dBeta_A_ctrl_OW <- function(i, j, pi_func, df_comb) {
  n <- nrow(df_comb)
  pi_i = pi_func[i]
  pi_j = pi_func[j]
  Z_i = df_comb[i, "treatment"]
  Z_j = df_comb[j, "treatment"]
  Y_i = df_comb[i, "outcomes_comb"]
  Y_j = df_comb[j, "outcomes_comb"]
  X_i = df_comb[i, c("x1", "x2", "x3", "x4", "x5", "x6")]
  X_j = df_comb[j, c("x1", "x2", "x3", "x4", "x5", "x6")]
  
  comparison_i_j = as.numeric(Y_i < Y_j)  # Convert logical TRUE/FALSE to 1/0
  comparison_j_i = as.numeric(Y_j < Y_i)  # Convert logical TRUE/FALSE to 1/0
  
  A = 0.5 * (
    ((1-pi_j) * X_j - pi_i * X_i) * (1-pi_i) * pi_j * Z_i * (1 - Z_j) * comparison_i_j  +
      ((1-pi_i) * X_i - pi_j * X_j) * (1-pi_j) * pi_i * Z_j * (1 - Z_i) * comparison_j_i
  )
  #A = sum(A)
  
  return(A)
  
}









### B ### (Dim: 6*1)


dBeta_B_OW <- function(pi_func, df_comb) {
  n <- nrow(df_comb)
  B <- matrix(0, nrow = 1, ncol = 6) # Initialize B as a zero matrix with 1 row and 6 columns
  
  # Pre-compute the values that are not dependent on i
  pi = pi_func
  Z = df_comb[, "treatment"]
  X = as.matrix(df_comb[, c("x1", "x2", "x3", "x4", "x5", "x6")])
  
  # Compute sum_j1 and sum_j2 outside the loop
  sum_j1 <- sum((1 - Z) * pi)
  sum_j2 <- colSums((1 - Z) * pi * (1 - pi) * X)
  
  for (i in 1:n) {
    pi_i = pi[i]
    Z_i = Z[i]
    X_i = X[i, ]
    
    term1 <- (-Z_i * pi_i * (1-pi_i)) * X_i * sum_j1
    term2 <- (Z_i * (1 - pi_i)) * sum_j2
    
    B <- B + term1 + term2
  }
  
  B <- B / (n * (n - 1))
  return(B)
}
dBeta_B_OW_val <- dBeta_B_OW(pi_func, df_comb)




### C ### (Dim: 6*1)

dBeta_C_OW <- function(pi_func, df_comb) {
  n <- nrow(df_comb)
  pi <- pi_func
  Z <- df_comb[, "treatment"]
  
  sum_i <- sum(Z * (1 - pi))
  sum_j <- sum((1 - Z) * pi)
  
  C <- (1 / (n * (n - 1))) * sum_i * sum_j
  
  return(C)
}

dBeta_C_OW_val <- dBeta_C_OW(pi_func, df_comb)


###dBeta_gg_OW Treatment### (Dim: 6*1)

dBeta_gg_OW_trt <- function(i, j, pi_func, df_comb) {
  temp_gg_OW_trt <- (1/2)*
    ((((1-pi_func[i]) * pi_func[j]) * 
        df_comb[i,"treatment"]*(1-df_comb[j,"treatment"])*
        as.numeric(df_comb[i,"outcomes_comb"] > df_comb[j,"outcomes_comb"]))
     + 
       (((1-pi_func[j]) * pi_func[i]) * 
          df_comb[j,"treatment"]*(1-df_comb[i,"treatment"])*
          as.numeric(df_comb[j,"outcomes_comb"] > df_comb[i,"outcomes_comb"]))
    ) /
    (
      (1/(nrow(df_comb)*(nrow(df_comb)-1)))*
        sum(df_comb[,"treatment"] * (1-pi_func))*
        sum((1-df_comb[,"treatment"]) * pi_func)
    )
  
  dBeta_gg_OW_trt_output <- dBeta_A_trt_OW(i, j, pi_func, df_comb)/dBeta_C_OW_val - 
    temp_gg_OW_trt * dBeta_B_OW_val/dBeta_C_OW_val
  return(dBeta_gg_OW_trt_output)
}



###dBeta_gg_OW Control### (Dim: 6*1)

dBeta_gg_OW_ctrl <- function(i, j, pi_func, df_comb) {
  temp_gg_OW_ctrl <- (1/2)*
    ((((1-pi_func[i]) * pi_func[j]) * 
        df_comb[i,"treatment"]*(1-df_comb[j,"treatment"])*
        as.numeric(df_comb[i,"outcomes_comb"] < df_comb[j,"outcomes_comb"]))
     + 
       (((1-pi_func[j]) * pi_func[i]) * 
          df_comb[j,"treatment"]*(1-df_comb[i,"treatment"])*
          as.numeric(df_comb[j,"outcomes_comb"] < df_comb[i,"outcomes_comb"]))
    ) /
    (
      (1/(nrow(df_comb)*(nrow(df_comb)-1)))*
        sum(df_comb[,"treatment"] * (1-pi_func))*
        sum((1-df_comb[,"treatment"]) * pi_func)
    )
  
  dBeta_gg_OW_ctrl_output <- dBeta_A_ctrl_OW(i, j, pi_func, df_comb)/dBeta_C_OW_val - 
    temp_gg_OW_ctrl * dBeta_B_OW_val/dBeta_C_OW_val
  return(dBeta_gg_OW_ctrl_output)
}



### An Treatment ### (Dim: 6*1)

calculate_An_trt_OW_parallel <- function(pi_func, df_comb) {
  n <- nrow(df_comb)
  comb_factor <- 2 / (n * (n - 1))
  
  # Define a function for the computation that needs to be parallelized
  compute_An_trt_OW <- function(i) {
    An_i <- matrix(0, nrow = 6, ncol = 1)  # Initialize a matrix to hold the results for this iteration
    for (j in (i + 1):n) {
      An_i <- An_i + dBeta_gg_OW_trt(i, j, pi_func, df_comb)
    }
    return(An_i)
  }
  
  results <- mclapply(1:(n - 1), compute_An_trt_OW, mc.cores = numCores - 1)
  
  An_matrix <- do.call(rbind, results)
  An_vector <- colSums(An_matrix) * comb_factor
  return(An_vector)
}

An_trt_OW_vector <- calculate_An_trt_OW_parallel(pi_func, df_comb)


### An Control ### (Dim: 6*1)

calculate_An_ctrl_OW_parallel <- function(pi_func, df_comb) {
  n <- nrow(df_comb)
  comb_factor <- 2 / (n * (n - 1))
  
  # Define a function for the computation that needs to be parallelized
  compute_An_ctrl_OW <- function(i) {
    An_i <- matrix(0, nrow = 6, ncol = 1)  # Initialize a matrix to hold the results for this iteration
    for (j in (i + 1):n) {
      An_i <- An_i + dBeta_gg_OW_ctrl(i, j, pi_func, df_comb)
    }
    return(An_i)
  }
  
  results <- mclapply(1:(n - 1), compute_An_ctrl_OW, mc.cores = numCores - 1)
  
  An_matrix <- do.call(rbind, results)
  An_vector <- colSums(An_matrix) * comb_factor
  return(An_vector)
}

An_ctrl_OW_vector <- calculate_An_ctrl_OW_parallel(pi_func, df_comb)




### E_beta0beta0 ### (Dim: 6*6)


Ebeta0 <- function() {
  n <- nrow(df_comb)
  X <- as.matrix(df_comb[, c("x1", "x2", "x3", "x4", "x5", "x6")])
  ps <- pi_func
  Ebeta = crossprod(sqrt(ps*(1-ps)) * X) / n
  return(Ebeta)
  
}

### l_beta ### (Dim: 6*1)


l_beta_fun <- function(i) {
  X_i <- as.numeric(df_comb[i, c("x1", "x2", "x3", "x4", "x5", "x6")])
  Z_i <- df_comb[i, "treatment"]
  pi_i <- pi_func[i]
  l_beta <- (solve(Ebeta0())) %*% X_i * (Z_i - pi_i)
  return(l_beta)
}



OW_Variance_trt <- function() {
  n <- nrow(df_comb)
  var_sum <- 0
  An_t <- t(An_trt_OW_vector)
  for (i in 1:n) {
    var_sum <- var_sum + (2 * (g1_OW_trt[i] - WP_trt_sim_OW[count_temp]) +
                            An_t %*% l_beta_fun(i))^2

  }
  
  OW_Var_trt <- var_sum / n
  return(OW_Var_trt)
}


OW_Variance_ctrl <- function() {
  n <- nrow(df_comb)
  var_sum <- 0
  An_t <- t(An_ctrl_OW_vector)
  for (i in 1:n) {
    var_sum <- var_sum + (2 * (g1_OW_ctrl[i] - WP_ctrl_sim_OW[count_temp]) +
                            An_t %*% l_beta_fun(i))^2

  }
  
  OW_Var_ctrl <- var_sum / n
  return(OW_Var_ctrl)
}

### OW Covariance 
OW_cov <- function(){
  cov_sum <- 0
  An_t_trt <- t(An_trt_OW_vector)
  An_t_ctrl <- t(An_ctrl_OW_vector)
  for (i in 1:n) {
    l_beta_fun_i = l_beta_fun(i)
    temp_trt = (2 * (g1_OW_trt[i] - WP_trt_sim_OW[count_temp]) + An_t_trt %*% l_beta_fun_i)
    temp_ctrl = (2 * (g1_OW_ctrl[i] - WP_ctrl_sim_OW[count_temp]) + An_t_ctrl %*% l_beta_fun_i)
    cov_sum = cov_sum + temp_trt*temp_ctrl
  }
  OW_cov_results <- cov_sum / n
  return(OW_cov_results)
}

Var_trt_OW_calculated <- OW_Variance_trt()/nrow(df_comb)
Var_ctrl_OW_calculated <- OW_Variance_ctrl()/nrow(df_comb)
cov_OW_calculated <- OW_cov()/nrow(df_comb)







#################OW Estimator Output#########################


###Theoretical Variance###

theory_var_OW_trt[count_temp] <- Var_trt_OW_calculated
theory_var_OW_ctrl[count_temp] <- Var_ctrl_OW_calculated
cov_trt_ctrl[count_temp] <- cov_OW_calculated
theory_var_OW_WR[count_temp] <- var_ratio(WP_trt_sim_OW[count_temp],WP_ctrl_sim_OW[count_temp],
                                          theory_var_OW_trt[count_temp],theory_var_OW_ctrl[count_temp],
                                          cov_trt_ctrl[count_temp])




}#open on line




theory_sd_OW_WR <- sqrt(theory_var_OW_WR)


WR_OW_coverage[size_count] <- mean((WR_sim_true > (WR_sim_OW - qnorm(0.975)*theory_sd_OW_WR)) & 
                                     (WR_sim_true < (WR_sim_OW + qnorm(0.975)*theory_sd_OW_WR)))
WR_OW_var_ratio[size_count] <- mean(theory_var_OW_WR)/var(WR_sim_OW)



WR_OW_theoretical_var_mean[size_count] <- mean(theory_var_OW_WR)
WR_OW_empirical_var[size_count] <- var(WR_sim_OW)




}


end_time <- Sys.time()
run_time <- end_time - start_time
run_time
#20 mins per 2 sim_num



######Coverage Output######




VarCoverage_OW_df <- as.data.frame(cbind(WR_OW_var_ratio,
                                         WR_OW_coverage,
                                         WR_OW_theoretical_var_mean,
                                         WR_OW_empirical_var))

colnames(VarCoverage_OW_df) <- c("Var_Ratio_WR",
                                 "Coverage_WR",
                                 "theoretical_var_mean_WR",
                                 "empirical_var_WR")

write.table(VarCoverage_OW_df,
            file=paste("results/OW_VarCoverage.txt",sep=""), 
            sep="\t", row.names=F)


















library(doMC)
library(MASS)
library(nnet)
#library(VGAM)
library(ggplot2)
library(parallel)
library(foreach)
library(doParallel)
#numCores <- detectCores()

######################################################
#########Theoretical Variance of Win Ratio############
######################################################

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
re_count <- scenarios$re_count #sample size count
#sim_num <- scenarios$sim_num #simulation iteration


start_time <- Sys.time()

##################################Simulation##################################


#vector of outcomes as a factor with ordered levels
outcomes_3lvl <- factor(c("first", "second", "third"), 
                        levels = c("first", "second", "third"), 
                        ordered = TRUE)

#Order: first < second < third, the larger the better

inv_logit <- function(logit) exp(logit)/(1 + exp(logit))


theory_var_Ori_trt <- numeric(re_count)
theory_var_Ori_ctrl <- numeric(re_count)
theory_var_IPW_trt <- numeric(re_count)
theory_var_IPW_ctrl <- numeric(re_count)
theory_var_OW_trt <- numeric(re_count)
theory_var_OW_ctrl <- numeric(re_count)
theory_var_AIPW_trt <- numeric(re_count)
theory_var_AIPW_ctrl <- numeric(re_count)
theory_var_AOW_trt <- numeric(re_count)
theory_var_AOW_ctrl <- numeric(re_count)


trueWD <- numeric(re_count)
trueWR <- numeric(re_count)








for (count_temp in 1:re_count){
  set.seed(count_temp)
n_count <- count_temp*step_size


# covariates
x1_trt <- rnorm(n_count, mean = 1, sd = sds[1])
x1_ctrl <- rnorm(n_count, mean = 1, sd = sds[1])
x2_trt <- rnorm(n_count, mean = 0.9, sd = sds[2])
x2_ctrl <- rnorm(n_count, mean = 0.9, sd = sds[2])
x3_trt <- rnorm(n_count, mean = 0.8, sd = sds[3])
x3_ctrl <- rnorm(n_count, mean = 0.8, sd = sds[3])

x4_trt <- rbinom(n_count,1,bern_param[1])
x4_ctrl <- rbinom(n_count,1,bern_param[1])
x5_trt <- rbinom(n_count,1,bern_param[2])
x5_ctrl <- rbinom(n_count,1,bern_param[2])
x6_trt <- rbinom(n_count,1,bern_param[3])
x6_ctrl <- rbinom(n_count,1,bern_param[3])

trt_cov <- data.frame(x1_trt, x2_trt, x3_trt, x4_trt, x5_trt, x6_trt)
ctrl_cov <- data.frame(x1_ctrl, x2_ctrl, x3_ctrl, x4_ctrl, x5_ctrl, x6_ctrl)

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
for (i in 1:n_count) {
  outcomes_trt[i] <- sample(
    outcomes_3lvl, 
    size = 1,
    prob = c(prob_1_trt[i], prob_2_trt[i], prob_3_trt[i])
  )
}



outcomes_ctrl <- c()
for (i in 1:n_count) {
  outcomes_ctrl[i] <- sample(
    outcomes_3lvl, 
    size = 1,
    prob = c(prob_1_ctrl[i], prob_2_ctrl[i], prob_3_ctrl[i])
  )
}







############################True Win Ratio##############################

# ctrl but assigned in trt group
logodds1_trt_assigned <- b01 + as.matrix(combined_ctrl_quad) %*% c(bi_trt, bi_trt_quad) + trt_eff1
logodds2_trt_assigned <- b02 + as.matrix(combined_ctrl_quad) %*% c(bi_trt, bi_trt_quad) + trt_eff1

# trt but assigned in ctrl group
logodds1_ctrl_assigned <- b01 + as.matrix(combined_trt_quad) %*% c(bi_ctrl, bi_ctrl_quad)
logodds2_ctrl_assigned <- b02 + as.matrix(combined_trt_quad) %*% c(bi_ctrl, bi_ctrl_quad)

## Probability

prob_2to3_trt_assigned <- inv_logit(logodds1_trt_assigned)
prob_3_trt_assigned <- inv_logit(logodds2_trt_assigned)
prob_1_trt_assigned <- 1 - prob_2to3_trt_assigned
prob_2_trt_assigned <- prob_2to3_trt_assigned - prob_3_trt_assigned

## Probability

prob_2to3_ctrl_assigned <- inv_logit(logodds1_ctrl_assigned)
prob_3_ctrl_assigned <- inv_logit(logodds2_ctrl_assigned)
prob_1_ctrl_assigned <- 1 - prob_2to3_ctrl_assigned
prob_2_ctrl_assigned <- prob_2to3_ctrl_assigned - prob_3_ctrl_assigned

## outcomes
outcomes_trt_assigned <- c()
for (i in 1:n_count) {
  outcomes_trt_assigned[i] <- sample(
    outcomes_3lvl, 
    size = 1,
    prob = c(prob_1_trt_assigned[i], prob_2_trt_assigned[i], prob_3_trt_assigned[i])
  )
}

outcomes_ctrl_assigned <- c()
for (i in 1:n_count) {
  outcomes_ctrl_assigned[i] <- sample(
    outcomes_3lvl, 
    size = 1,
    prob = c(prob_1_ctrl_assigned[i], prob_2_ctrl_assigned[i], prob_3_ctrl_assigned[i])
  )
}



df_trt_true <- data.frame(c(outcomes_trt, outcomes_trt_assigned))
df_ctrl_true <- data.frame(c(outcomes_ctrl_assigned, outcomes_ctrl))

# Initialize matrices
trtwin_true <- matrix(NA, nrow = nrow(df_trt_true), ncol = nrow(df_ctrl_true))
ctrlwin_true <- matrix(NA, nrow = nrow(df_trt_true), ncol = nrow(df_ctrl_true))

# Define the comparison function for true values
compare_rows_true <- function(i) {
  trt_row_true <- df_trt_true[i, 1]
  trt_vec_true <- numeric(nrow(df_ctrl_true))
  ctrl_vec_true <- numeric(nrow(df_ctrl_true))
  
  for (j in 1:nrow(df_ctrl_true)) {
    if (trt_row_true > df_ctrl_true[j, 1]) {
      trt_vec_true[j] = 1
      ctrl_vec_true[j] = 0
    } else if (trt_row_true < df_ctrl_true[j, 1]) {
      trt_vec_true[j] = 0
      ctrl_vec_true[j] = 1
    }
  }
  
  return(list(trt = trt_vec_true, ctrl = ctrl_vec_true))
}

# Parallelize the outer loop for true values
results_true <- mclapply(1:nrow(df_trt_true), compare_rows_true, mc.cores = numCores)

# Aggregate results for true values
for (i in 1:length(results_true)) {
  trtwin_true[i, ] <- results_true[[i]]$trt
  ctrlwin_true[i, ] <- results_true[[i]]$ctrl
}

# Calculate win proportions for true values
trt_winpr_true <- sum(as.vector(trtwin_true), na.rm = TRUE) / (nrow(df_trt_true) * nrow(df_ctrl_true))
ctrl_winpr_true <- sum(as.vector(ctrlwin_true), na.rm = TRUE) / (nrow(df_trt_true) * nrow(df_ctrl_true))

WD_sim_trueWR <- trt_winpr_true - ctrl_winpr_true
WR_sim_trueWR <- trt_winpr_true / ctrl_winpr_true






##############Original Unadjusted Estimator##########################

df_trt <- data.frame(outcomes_trt)
df_ctrl <- data.frame(outcomes_ctrl)

# Initialize matrices
trtwin <- matrix(NA, nrow = nrow(df_trt), ncol = nrow(df_ctrl))
ctrlwin <- matrix(NA, nrow = nrow(df_trt), ncol = nrow(df_ctrl))

# Define the comparison function
compare_rows_unadj <- function(i) {
  trt_row <- df_trt[i, 1]
  trt_vec <- numeric(nrow(df_ctrl))
  ctrl_vec <- numeric(nrow(df_ctrl))
  
  for (j in 1:nrow(df_ctrl)) {
    if (trt_row > df_ctrl[j, 1]) {
      trt_vec[j] = 1
      ctrl_vec[j] = 0
    } else if (trt_row < df_ctrl[j, 1]) {
      trt_vec[j] = 0
      ctrl_vec[j] = 1
    }
  }
  
  return(list(trt = trt_vec, ctrl = ctrl_vec))
}

# Parallelize the outer loop
results <- mclapply(1:nrow(df_trt), compare_rows_unadj, mc.cores = numCores)

# Aggregate results
for (i in 1:length(results)) {
  trtwin[i, ] <- results[[i]]$trt
  ctrlwin[i, ] <- results[[i]]$ctrl
}

# Calculate win proportions
trt_winpr <- sum(as.vector(trtwin), na.rm = TRUE) / (nrow(df_trt) * nrow(df_ctrl))
ctrl_winpr <- sum(as.vector(ctrlwin), na.rm = TRUE) / (nrow(df_trt) * nrow(df_ctrl))

# Update WD_sim and WR_sim
WP_trt_sim_ori <- trt_winpr
WP_ctrl_sim_ori <- ctrl_winpr
WD_sim <- trt_winpr - ctrl_winpr
WR_sim <- trt_winpr / ctrl_winpr









###########Unadjusted Theoretical Variance###########

colnames(df_trt) <- "outcomes_comb"
colnames(df_ctrl) <- "outcomes_comb"

df_comb <- rbind(df_trt, df_ctrl)
df_comb$treatment <- c(rep(1,n_count), rep(0,n_count))

sum_treatment <- sum(df_comb[,"treatment"])
sum_control <- sum(1 - df_comb[,"treatment"])

# Function to parallelize for the treatment group
p_ijik_trt_parallel <- function(i, df_comb) {
  n <- nrow(df_comb)
  temp_inner <- numeric(n)
  
  for (j in 1:n) {
    wij_unadjusted <- (df_comb[i,"treatment"] * (1 - df_comb[j,"treatment"])) /
      ((1/n) * sum_treatment * (1/n) * sum_control)
    temp_inner[j] <- wij_unadjusted * as.numeric(df_comb[i,"outcomes_comb"] > df_comb[j,"outcomes_comb"])
  }
  
  return((sum(temp_inner))^2)
}

# Parallel computation of the term1 for the treatment group
term1_trt_list <- mclapply(1:nrow(df_comb), p_ijik_trt_parallel, df_comb, mc.cores = numCores - 1)
term1_trt <- sum(unlist(term1_trt_list))/(nrow(df_comb)^2 * (nrow(df_comb) - 1))
term2_trt <- (1/(nrow(df_comb)-1)) * WP_trt_sim_ori

p_ijik_trt_results <- term1_trt - term2_trt


# Function to parallelize for the control group
p_ijik_ctrl_parallel <- function(i, df_comb) {
  n <- nrow(df_comb)
  temp_inner <- numeric(n)
  
  for (j in 1:n) {
    wij_unadjusted <- (df_comb[i,"treatment"] * (1 - df_comb[j,"treatment"])) /
      ((1/n) * sum_treatment * (1/n) * sum_control)
    temp_inner[j] <- wij_unadjusted * as.numeric(df_comb[i,"outcomes_comb"] < df_comb[j,"outcomes_comb"])
  }
  
  return((sum(temp_inner))^2)
}

# Parallel computation of the term1 for the control group
term1_ctrl_list <- mclapply(1:nrow(df_comb), p_ijik_ctrl_parallel, df_comb, mc.cores = numCores - 1)
term1_ctrl <- sum(unlist(term1_ctrl_list))/(nrow(df_comb)^2 * (nrow(df_comb) - 1))
term2_ctrl <- (1/(nrow(df_comb)-1)) * WP_ctrl_sim_ori

p_ijik_ctrl_results <- term1_ctrl - term2_ctrl


Var_trt_ori_calculated <- 1/(nrow(df_comb)) * (
  1/(nrow(df_trt)/nrow(df_comb)) * p_ijik_trt_results + 
    1/(1-nrow(df_trt)/nrow(df_comb)) * p_ijik_trt_results - 
    2*(WP_trt_sim_ori)^2
)

Var_ctrl_ori_calculated <- 1/(nrow(df_comb)) * (
  1/(nrow(df_ctrl)/nrow(df_comb)) * p_ijik_ctrl_results + 
    1/(1-nrow(df_ctrl)/nrow(df_comb)) * p_ijik_ctrl_results - 
    2*(WP_ctrl_sim_ori)^2
)









##########################IPW########################################




df_comb$x1 <- c(x1_trt, x1_ctrl)
df_comb$x2 <- c(x2_trt, x2_ctrl)
df_comb$x3 <- c(x3_trt, x3_ctrl)
df_comb$x4 <- c(x4_trt, x4_ctrl)
df_comb$x5 <- c(x5_trt, x5_ctrl)
df_comb$x6 <- c(x6_trt, x6_ctrl)

# Calculate propensity scores
PropScore <- glm(treatment~x1 + x2 + x3 + x4 + x5 + x6, 
                 data = df_comb, family=binomial)
pi_func <- fitted(PropScore)




# Define the IPW comparison function
compare_rows_IPW <- function(i) {
  g_hfunc_IPW1_vec <- numeric(2*n_count)
  g_hfunc_IPW2_vec <- numeric(2*n_count)
  
  for (j in 1:(2*n_count)) {
    if (df_comb[i,1] > df_comb[j,1]) {
      g_hfunc_IPW1_vec[j] = (df_comb[i,"treatment"] * (1-df_comb[j,"treatment"]) * 1 / (pi_func[i] * (1-pi_func[j])))
      g_hfunc_IPW2_vec[j] = 0
    } else if (df_comb[i,1] < df_comb[j,1]) {
      g_hfunc_IPW1_vec[j] = 0
      g_hfunc_IPW2_vec[j] = (df_comb[i,"treatment"] * (1-df_comb[j,"treatment"]) * 1 / (pi_func[i] * (1-pi_func[j])))
    }
  }
  
  return(list(IPW1 = g_hfunc_IPW1_vec, IPW2 = g_hfunc_IPW2_vec))
}

# Parallelize the outer loop
results_IPW <- mclapply(1:(2*n_count), compare_rows_IPW, mc.cores = numCores)

# Aggregate results
g_hfunc_IPW1 <- do.call(rbind, lapply(results_IPW, function(x) x$IPW1))
g_hfunc_IPW2 <- do.call(rbind, lapply(results_IPW, function(x) x$IPW2))

# Calculate tau values
tau1_IPW <- sum(as.vector(g_hfunc_IPW1), na.rm = TRUE) /
  (sum(df_comb[,"treatment"]/pi_func, na.rm = TRUE) * sum((1-df_comb[,"treatment"])/(1-pi_func)))
tau2_IPW <- sum(as.vector(g_hfunc_IPW2), na.rm = TRUE) /
  (sum(df_comb[,"treatment"]/pi_func, na.rm = TRUE) * sum((1-df_comb[,"treatment"])/(1-pi_func)))

# Update WD_sim_IPW and WR_sim_IPW
WP_trt_sim_IPW <- tau1_IPW
WP_ctrl_sim_IPW <- tau2_IPW
WD_sim_IPW <- tau1_IPW - tau2_IPW
WR_sim_IPW <- tau1_IPW / tau2_IPW




###########IPW Theoretical Variance###########




###g1_IPW_trt###

# Pre-calculate constants
sum_trt_over_pi <- sum(df_comb[,"treatment"]/pi_func)
sum_ctrl_over_one_minus_pi <- sum((1-df_comb[,"treatment"])/(1-pi_func))
n <- nrow(df_comb)
constant_denominator <- (1/(n*(n-1))) * sum_trt_over_pi * sum_ctrl_over_one_minus_pi

g1_IPW_trt_fun <- function(i, outcomes, treatment, pi_func, constant_denominator) {
  # Vectorized operation to calculate pairwise comparisons
  pairwise_comparisons <- (1/2) * (
    ((treatment[i] * (1 - treatment) * as.numeric(outcomes[i] > outcomes)) / 
      (pi_func[i] * (1 - pi_func))) + 
      (((1-treatment[i]) * treatment * as.numeric(outcomes[i] < outcomes)) / 
          (pi_func * (1 - pi_func[i])))
  )
  
  # Avoid comparison of the subject with itself
  pairwise_comparisons[i] <- 0
  
  # Sum and adjust according to the formula
  sum_g1_IPW_trt <- sum(pairwise_comparisons) / constant_denominator
  return((1/(n-1)) * sum_g1_IPW_trt)
}

# Parallelize the computation for g1_IPW_trt
results_g1_IPW_trt <- mclapply(1:n, g1_IPW_trt_fun, df_comb$outcomes_comb, 
                               df_comb$treatment, pi_func, constant_denominator, 
                               mc.cores = numCores)
# Aggregate results
g1_IPW_trt <- unlist(results_g1_IPW_trt)




###g1_IPW_ctrl###

g1_IPW_ctrl_fun <- function(i, outcomes, treatment, pi_func, constant_denominator) {
  # Vectorized operation to calculate pairwise comparisons
  pairwise_comparisons <- (1/2) * (
    ((treatment[i] * (1 - treatment) * as.numeric(outcomes[i] < outcomes)) / 
       (pi_func[i] * (1 - pi_func))) + 
      (((1-treatment[i]) * treatment * as.numeric(outcomes[i] > outcomes)) / 
         (pi_func * (1 - pi_func[i])))
  )
  
  # Avoid comparison of the subject with itself
  pairwise_comparisons[i] <- 0
  
  # Sum and adjust according to the formula
  sum_g1_IPW_ctrl <- sum(pairwise_comparisons) / constant_denominator
  return((1/(n-1)) * sum_g1_IPW_ctrl)
}

# Parallelize the computation for g1_IPW_trt
results_g1_IPW_ctrl <- mclapply(1:n, g1_IPW_ctrl_fun, df_comb$outcomes_comb, 
                               df_comb$treatment, pi_func, constant_denominator, 
                               mc.cores = numCores)

# Aggregate results
g1_IPW_ctrl <- unlist(results_g1_IPW_ctrl)




###A Treatment### (Dim: 6*1)




dBeta_A_trt_IPW <- function(i, j) {
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
    (pi_j * X_j - (1 - pi_i) * X_i) * Z_i * (1 - Z_j) * comparison_i_j / (pi_i * (1 - pi_j)) +
      (pi_i * X_i - (1 - pi_j) * X_j) * Z_j * (1 - Z_i) * comparison_j_i / (pi_j * (1 - pi_i))
  )
  #A = sum(A)
  
  return(A)
  
}


###A Control### (Dim: 6*1)


dBeta_A_ctrl_IPW <- function(i, j) {
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
    (pi_j * X_j - (1 - pi_i) * X_i) * Z_i * (1 - Z_j) * comparison_i_j / (pi_i * (1 - pi_j)) +
      (pi_i * X_i - (1 - pi_j) * X_j) * Z_j * (1 - Z_i) * comparison_j_i / (pi_j * (1 - pi_i))
  )
  #A = sum(A)
  
  return(A)
  
}




### B ### (Dim: 6*1)


dBeta_B <- function() {
  n <- nrow(df_comb)
  B <- matrix(0, nrow = 1, ncol = 6) # Initialize B as a zero matrix with 1 row and 6 columns
  
  # Pre-compute the values that are not dependent on i
  pi = pi_func
  Z = df_comb[, "treatment"]
  X = as.matrix(df_comb[, c("x1", "x2", "x3", "x4", "x5", "x6")])
  
  # Compute sum_j1 and sum_j2 outside the loop
  sum_j1 <- sum((1 - Z) / (1 - pi))
  sum_j2 <- colSums((1 - Z) * pi * X / (1 - pi))
  
  for (i in 1:n) {
    pi_i = pi[i]
    Z_i = Z[i]
    X_i = X[i, ]
    
    term1 <- (-Z_i * (1 - pi_i) / pi_i) * X_i * sum_j1
    term2 <- (Z_i / pi_i) * sum_j2
    
    B <- B + term1 + term2
  }
  
  B <- B / (n * (n - 1))
  return(B)
}





### C ### (Dim: 1*1)

dBeta_C <- function() {
  n <- nrow(df_comb)
  pi <- pi_func
  Z <- df_comb[, "treatment"]
  
  sum_i <- sum(Z / pi)
  sum_j <- sum((1 - Z) / (1 - pi))
  
  C <- (1 / (n * (n - 1))) * sum_i * sum_j
  
  return(C)
}



###dBeta_gg_IPW Treatment### (Dim: 6*1)


dBeta_gg_IPW_trt <- function(i,j) {
  temp_gg_IPW_trt <- (1/2)*
    (((df_comb[i,"treatment"]*(1-df_comb[j,"treatment"])*
         as.numeric(df_comb[i,"outcomes_comb"] > df_comb[j,"outcomes_comb"]))/
        (pi_func[i] * (1-pi_func[j]))) + 
       ((df_comb[j,"treatment"]*(1-df_comb[i,"treatment"])*
           as.numeric(df_comb[j,"outcomes_comb"] > df_comb[i,"outcomes_comb"]))/
          (pi_func[j] * (1-pi_func[i])))) /
    (
      (1/(nrow(df_comb)*(nrow(df_comb)-1)))*
        sum(df_comb[,"treatment"]/pi_func)*
        sum((1-df_comb[,"treatment"])/(1-pi_func))
    )
  
  dBeta_gg_IPW_trt_output <- dBeta_A_trt_IPW(i,j)/dBeta_C() - temp_gg_IPW_trt * dBeta_B()/dBeta_C()
  return(dBeta_gg_IPW_trt_output)
}



###dBeta_gg_IPW Control### (Dim: 6*1)


dBeta_gg_IPW_ctrl <- function(i,j) {
  temp_gg_IPW_ctrl <- (1/2)*
    (((df_comb[i,"treatment"]*(1-df_comb[j,"treatment"])*
         as.numeric(df_comb[i,"outcomes_comb"] < df_comb[j,"outcomes_comb"]))/
        (pi_func[i] * (1-pi_func[j]))) + 
       ((df_comb[j,"treatment"]*(1-df_comb[i,"treatment"])*
           as.numeric(df_comb[j,"outcomes_comb"] < df_comb[i,"outcomes_comb"]))/
          (pi_func[j] * (1-pi_func[i])))) /
    (
      (1/(nrow(df_comb)*(nrow(df_comb)-1)))*
        sum(df_comb[,"treatment"]/pi_func)*
        sum((1-df_comb[,"treatment"])/(1-pi_func))
    )
  
  dBeta_gg_IPW_ctrl_output <- dBeta_A_ctrl_IPW(i,j)/dBeta_C() - temp_gg_IPW_ctrl * dBeta_B()/dBeta_C()
  return(dBeta_gg_IPW_ctrl_output)
}





### An Treatment ### (Dim: 6*1)

calculate_An_trt_IPW_parallel <- function() {
  n <- nrow(df_comb)
  comb_factor <- 2 / (n * (n - 1))

  compute_An_trt <- function(i) {
    An_i <- matrix(0, nrow = 6, ncol = 1)
    for (j in (i + 1):n) {
      An_i <- An_i + dBeta_gg_IPW_trt(i, j)
    }
    return(An_i)
  }
  results <- mclapply(1:(n - 1), compute_An_trt, mc.cores = numCores - 1)
  
  An_matrix <- do.call(rbind, results)
  An_vector <- colSums(An_matrix) * comb_factor
  return(An_vector)
}

An_trt_IPW <- calculate_An_trt_IPW_parallel()




### An Control ### (Dim: 6*1)

calculate_An_ctrl_IPW_parallel <- function() {
  n <- nrow(df_comb)
  comb_factor <- 2 / (n * (n - 1))
  
  compute_An_ctrl <- function(i) {
    An_i <- matrix(0, nrow = 6, ncol = 1)
    for (j in (i + 1):n) {
      An_i <- An_i + dBeta_gg_IPW_ctrl(i, j)
    }
    return(An_i)
  }
  results <- mclapply(1:(n - 1), compute_An_ctrl, mc.cores = numCores - 1)
  
  An_matrix <- do.call(rbind, results)
  An_vector <- colSums(An_matrix) * comb_factor
  return(An_vector)
}

An_ctrl_IPW <- calculate_An_ctrl_IPW_parallel()




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



IPW_Variance_trt <- function() {
  n <- nrow(df_comb)
  var_sum <- 0
  An_t <- t(An_trt_IPW)
  for (i in 1:n) {
    var_sum <- var_sum + (2 * (g1_IPW_trt[i] - WP_trt_sim_IPW) + An_t %*% l_beta_fun(i))^2
  }
  
  IPW_Var_trt <- var_sum / n
  return(IPW_Var_trt)
}




IPW_Variance_ctrl <- function() {
  n <- nrow(df_comb)
  var_sum <- 0
  An_t <- t(An_ctrl_IPW)
  for (i in 1:n) {
    var_sum <- var_sum + (2 * (g1_IPW_ctrl[i] - WP_ctrl_sim_IPW) + An_t %*% l_beta_fun(i))^2
  }
  
  IPW_Var_ctrl <- var_sum / n
  return(IPW_Var_ctrl)
}



Var_trt_IPW_calculated <- IPW_Variance_trt()/nrow(df_comb)
Var_ctrl_IPW_calculated <- IPW_Variance_ctrl()/nrow(df_comb)









#################Unadjusted Estimator Output#########################


###Theoretical Variance###

theory_var_Ori_trt[count_temp] <- Var_trt_ori_calculated
theory_var_Ori_ctrl[count_temp] <- Var_ctrl_ori_calculated



#################IPW Estimator Output#########################


###Theoretical Variance###

theory_var_IPW_trt[count_temp] <- Var_trt_IPW_calculated
theory_var_IPW_ctrl[count_temp] <- Var_ctrl_IPW_calculated



}#open on line 85



end_time <- Sys.time()
run_time <- end_time - start_time




run_time
#51.23171 mins from 30-300 by 30
#36.90234 secs from 10-30 by 10
#14.40525 secs from 10-30 by 10






######Variance Output######




Ori_Var_df <- as.data.frame(cbind(theory_var_Ori_trt,
                                  theory_var_Ori_ctrl))

colnames(Ori_Var_df) <- c("TheoreticalVariance_trt",
                          "TheoreticalVariance_ctrl")

IPW_Var_df <- as.data.frame(cbind(theory_var_IPW_trt,
                                  theory_var_IPW_ctrl))

colnames(IPW_Var_df) <- c("IPW_TheoreticalVariance_trt",
                          "IPW_TheoreticalVariance_ctrl")




write.table(Ori_Var_df,
            file=paste("results/Unadjusted_theoreticalVar.txt",sep=""), 
            sep="\t", row.names=F)


write.table(IPW_Var_df,
            file=paste("results/IPW_theoreticalVar.txt",sep=""), 
            sep="\t", row.names=F)








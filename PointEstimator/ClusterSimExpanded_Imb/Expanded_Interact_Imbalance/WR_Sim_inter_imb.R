library(doMC)
library(MASS)
library(nnet)
#library(VGAM)
library(ggplot2)
library(parallel)
library(foreach)
library(doParallel)
#numCores <- detectCores()

#################################################################
######## Expanded with interaction terms (Imbalance) ############
#################################################################

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
sim_num <- scenarios$sim_num #simulation iteration


start_time <- Sys.time()

##################################Simulation##################################


#vector of outcomes as a factor with ordered levels
outcomes_3lvl <- factor(c("first", "second", "third"), 
                    levels = c("first", "second", "third"), 
                    ordered = TRUE)

#Order: first < second < third, the larger the better

inv_logit <- function(logit) exp(logit)/(1 + exp(logit))









empirical_var_Ori_WD <- numeric(re_count)
empirical_var_IPW_WD <- numeric(re_count)
empirical_var_OW_WD <- numeric(re_count)
empirical_var_DR_WD <- numeric(re_count)
empirical_var_DROW_WD <- numeric(re_count)

empirical_var_Ori_WR <- numeric(re_count)
empirical_var_IPW_WR <- numeric(re_count)
empirical_var_OW_WR <- numeric(re_count)
empirical_var_DR_WR <- numeric(re_count)
empirical_var_DROW_WR <- numeric(re_count)

bias_Ori_WD <- numeric(re_count)
bias_IPW_WD <- numeric(re_count)
bias_OW_WD <- numeric(re_count)
bias_DR_WD <- numeric(re_count)
bias_DROW_WD <- numeric(re_count)

bias_Ori_WR <- numeric(re_count)
bias_IPW_WR <- numeric(re_count)
bias_OW_WR <- numeric(re_count)
bias_DR_WR <- numeric(re_count)
bias_DROW_WR <- numeric(re_count)

pct_bias_Ori_WD <- numeric(re_count)
pct_bias_IPW_WD <- numeric(re_count)
pct_bias_OW_WD <- numeric(re_count)
pct_bias_DR_WD <- numeric(re_count)
pct_bias_DROW_WD <- numeric(re_count)

pct_bias_Ori_WR <- numeric(re_count)
pct_bias_IPW_WR <- numeric(re_count)
pct_bias_OW_WR <- numeric(re_count)
pct_bias_DR_WR <- numeric(re_count)
pct_bias_DROW_WR <- numeric(re_count)

mse_Ori_WD <- numeric(re_count)
mse_IPW_WD <- numeric(re_count)
mse_OW_WD <- numeric(re_count)
mse_DR_WD <- numeric(re_count)
mse_DROW_WD <- numeric(re_count)

mse_Ori_WR <- numeric(re_count)
mse_IPW_WR <- numeric(re_count)
mse_OW_WR <- numeric(re_count)
mse_DR_WR <- numeric(re_count)
mse_DROW_WR <- numeric(re_count)

trueWD <- 0.2888708 - 0.2000238
trueWR <- 1.474961





for (count_temp in 1:re_count){

n_count <- 2*count_temp*step_size


WD_sim <- numeric(sim_num)
WD_sim_IPW <- numeric(sim_num)
WD_sim_OW <- numeric(sim_num)
WD_sim_DR <- numeric(sim_num)
WD_sim_DROW <- numeric(sim_num)

WR_sim <- numeric(sim_num)
WR_sim_IPW <- numeric(sim_num)
WR_sim_OW <- numeric(sim_num)
WR_sim_DR <- numeric(sim_num)
WR_sim_DROW <- numeric(sim_num)




for (sim_count in 1:sim_num){
set.seed(sim_count)

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

trt_cov_interactions <- data.frame(matrix(nrow = nrow(trt_cov), ncol = 0))
ctrl_cov_interactions <- data.frame(matrix(nrow = nrow(ctrl_cov), ncol = 0))

for (i in 1:(ncol(trt_cov)-1)) {
  for (j in (i+1):ncol(trt_cov)) {
    trt_cov_interactions[, paste0(names(trt_cov)[i], "_x_", names(trt_cov)[j])] <- trt_cov[,i] * trt_cov[,j]
    ctrl_cov_interactions[, paste0(names(ctrl_cov)[i], "_x_", names(ctrl_cov)[j])] <- ctrl_cov[,i] * ctrl_cov[,j]
  }
}

bi_trt_interactions <- numeric()
bi_ctrl_interactions <- numeric()



for (i in 1:(length(bi_trt)-1)) {
  for (j in (i+1):length(bi_trt)) {
    if(sign(bi_trt[i]) == sign(bi_trt[j])){
      bi_trt_interactions <- c(bi_trt_interactions, bi_trt[i] + bi_trt[j])
      bi_ctrl_interactions <- c(bi_ctrl_interactions, bi_ctrl[i] + bi_ctrl[j])
    }
    else {
      bi_trt_interactions <- c(bi_trt_interactions, 0.25*(sample(c(-1, 1), 1))*bi_trt[i] * bi_trt[j])
      bi_ctrl_interactions <- c(bi_ctrl_interactions, 0.25*(sample(c(-1, 1), 1))*bi_ctrl[i] * bi_ctrl[j])
    }
  }
}



combined_trt_interact <- cbind(trt_cov, trt_cov_interactions)
combined_ctrl_interact <- cbind(ctrl_cov, ctrl_cov_interactions)

logodds1_trt <- b01 + as.matrix(combined_trt_interact) %*% c(bi_trt, bi_trt_interactions) + trt_eff1
logodds2_trt <- b02 + as.matrix(combined_trt_interact) %*% c(bi_trt, bi_trt_interactions) + trt_eff1

logodds1_ctrl <- b01 + as.matrix(combined_ctrl_interact) %*% c(bi_ctrl, bi_ctrl_interactions)
logodds2_ctrl <- b02 + as.matrix(combined_ctrl_interact) %*% c(bi_ctrl, bi_ctrl_interactions)


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
WD_sim[sim_count] <- trt_winpr - ctrl_winpr
WR_sim[sim_count] <- trt_winpr / ctrl_winpr











##########################IPW########################################




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

# Define the IPW comparison function
compare_rows_IPW <- function(i) {
  g_hfunc_IPW1_vec <- numeric(n_count)
  g_hfunc_IPW2_vec <- numeric(n_count)
  
  for (j in 1:n_count) {
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
results_IPW <- mclapply(1:n_count, compare_rows_IPW, mc.cores = numCores)

# Aggregate results
g_hfunc_IPW1 <- do.call(rbind, lapply(results_IPW, function(x) x$IPW1))
g_hfunc_IPW2 <- do.call(rbind, lapply(results_IPW, function(x) x$IPW2))

# Calculate tau values
tau1_IPW <- sum(as.vector(g_hfunc_IPW1), na.rm = TRUE) /
  (sum(df_comb[,"treatment"]/pi_func, na.rm = TRUE) * sum((1-df_comb[,"treatment"])/(1-pi_func)))
tau2_IPW <- sum(as.vector(g_hfunc_IPW2), na.rm = TRUE) /
  (sum(df_comb[,"treatment"]/pi_func, na.rm = TRUE) * sum((1-df_comb[,"treatment"])/(1-pi_func)))

# Update WD_sim_IPW and WR_sim_IPW
WD_sim_IPW[sim_count] <- tau1_IPW - tau2_IPW
WR_sim_IPW[sim_count] <- tau1_IPW / tau2_IPW









######################Overlap Weighting##############################




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
WD_sim_OW[sim_count] <- tau1_OW - tau2_OW
WR_sim_OW[sim_count] <- tau1_OW / tau2_OW







###############Doubly Robust (IPW) Estimator################################


###DR mu-function

if((!(1 %in% df_comb[df_comb$treatment==1,]$outcomes_comb) |
    !(2 %in% df_comb[df_comb$treatment==1,]$outcomes_comb) |
    !(3 %in% df_comb[df_comb$treatment==1,]$outcomes_comb))|
   (!(1 %in% df_comb[df_comb$treatment==0,]$outcomes_comb) |
    !(2 %in% df_comb[df_comb$treatment==0,]$outcomes_comb) |
    !(3 %in% df_comb[df_comb$treatment==0,]$outcomes_comb))){
  next
}

else {
  #dr_trt <- multinom(factor(outcomes_comb) ~ (x1 + x2 + x3 + x4 + x5 + x6)^2,
  #                   data = df_comb[df_comb$treatment == 1,], trace = FALSE)
  lp_start <- rep(0, 21)
  th_start <- c(-1, 1)
  start_values <- c(lp_start, th_start)
  dr_trt <- polr(factor(outcomes_comb) ~ (x1 + x2 + x3 + x4 + x5 + x6)^2,
                 data = df_comb[df_comb$treatment == 1,],
                 start = start_values)
  cond_prob_trt <- predict(dr_trt, newdata = df_comb, type = "probs")
  
  #dr_ctrl <- multinom(factor(outcomes_comb) ~ (x1 + x2 + x3 + x4 + x5 + x6)^2,
  #                    data = df_comb[df_comb$treatment == 0,], trace = FALSE)
  dr_ctrl <- polr(factor(outcomes_comb) ~ (x1 + x2 + x3 + x4 + x5 + x6)^2,
                  data = df_comb[df_comb$treatment == 0,],
                  start = start_values)
  cond_prob_ctrl <- predict(dr_ctrl, newdata = df_comb, type = "probs")
}



cl <- makeCluster(numCores - 1)  # Use all cores except one
registerDoParallel(cl)

# Define the parallelized function
compare_rows_DR <- function(i, df_comb, pi_func, cond_prob_trt, cond_prob_ctrl) {
  n_count <- nrow(df_comb)
  
  DR_num1_vec <- numeric(n_count)
  DR_denom1_vec <- numeric(n_count)
  DR_num2_vec <- numeric(n_count)
  DR_denom2_vec <- numeric(n_count)
  DR_mu1_vec <- numeric(n_count)
  DR_mu2_vec <- numeric(n_count)
  DR_mu_denom1_vec <- numeric(n_count)
  DR_mu_denom2_vec <- numeric(n_count)
  
  for(j in 1:n_count){
    if(i != j){
      DR_num1_vec[j] = (df_comb[i,'treatment'] * (1 - df_comb[j,'treatment'])) * (1/(pi_func[i] * (1 - pi_func[j]))) *
        ((df_comb[i,1] > df_comb[j,1]) - 
           (cond_prob_trt[i, 2] * cond_prob_ctrl[j, 1] + 
              cond_prob_trt[i, 3] * cond_prob_ctrl[j, 1] +
              cond_prob_trt[i, 3] * cond_prob_ctrl[j, 2]))
      
      DR_denom1_vec[j] = (df_comb[i,'treatment'] * (1 - df_comb[j,'treatment'])) * (1/(pi_func[i] * (1 - pi_func[j])))
      
      DR_num2_vec[j] = (df_comb[i,'treatment'] * (1 - df_comb[j,'treatment'])) * (1/(pi_func[i] * (1 - pi_func[j]))) *
        ((df_comb[i,1] < df_comb[j,1]) - 
           (cond_prob_trt[i, 1] * cond_prob_ctrl[j, 2] + 
              cond_prob_trt[i, 1] * cond_prob_ctrl[j, 3] +
              cond_prob_trt[i, 2] * cond_prob_ctrl[j, 3]))
      
      DR_denom2_vec[j] = (df_comb[i,'treatment'] * (1 - df_comb[j,'treatment'])) * (1/(pi_func[i] * (1 - pi_func[j])))
      
      DR_mu1_vec[j] = 1 * (cond_prob_trt[i, 2] * cond_prob_ctrl[j, 1] + 
                             cond_prob_trt[i, 3] * cond_prob_ctrl[j, 1] + 
                             cond_prob_trt[i, 3] * cond_prob_ctrl[j, 2])
      
      DR_mu2_vec[j] = 1 * (cond_prob_trt[i, 1] * cond_prob_ctrl[j, 2] + 
                             cond_prob_trt[i, 1] * cond_prob_ctrl[j, 3] + 
                             cond_prob_trt[i, 2] * cond_prob_ctrl[j, 3])
      
      DR_mu_denom1_vec[j] <- 1
      DR_mu_denom2_vec[j] <- 1
    }
  }
  
  return(list(DR_num1 = DR_num1_vec, DR_denom1 = DR_denom1_vec, 
              DR_num2 = DR_num2_vec, DR_denom2 = DR_denom2_vec,
              DR_mu1 = DR_mu1_vec, DR_mu_denom1 = DR_mu_denom1_vec,
              DR_mu2 = DR_mu2_vec, DR_mu_denom2 = DR_mu_denom2_vec))
  
}


#Use the function to compute the matrices
# results_DR <- lapply(1:n_count, compare_rows_DR,
#                      df_comb=df_comb, pi_func=pi_func,
#                      cond_prob_trt=cond_prob_trt, cond_prob_ctrl=cond_prob_ctrl)

results_DR <- foreach(i=1:n_count) %dopar% {
  compare_rows_DR(i, df_comb, pi_func, cond_prob_trt, cond_prob_ctrl)
}


# Combine the results
tau_DR_num1 <- do.call(rbind, lapply(results_DR, function(x) x$DR_num1))
tau_DR_denom1 <- do.call(rbind, lapply(results_DR, function(x) x$DR_denom1))
tau_DR_num2 <- do.call(rbind, lapply(results_DR, function(x) x$DR_num2))
tau_DR_denom2 <- do.call(rbind, lapply(results_DR, function(x) x$DR_denom2))
tau_DR_mu1 <- do.call(rbind, lapply(results_DR, function(x) x$DR_mu1))
tau_DR_mu_denom1 <- do.call(rbind, lapply(results_DR, function(x) x$DR_mu_denom1))
tau_DR_mu2 <- do.call(rbind, lapply(results_DR, function(x) x$DR_mu2))
tau_DR_mu_denom2 <- do.call(rbind, lapply(results_DR, function(x) x$DR_mu_denom2))



tau1_DR <- (sum(as.vector(tau_DR_num1), na.rm = T) / sum(as.vector(tau_DR_denom1), na.rm = T)) +
  (sum(as.vector(tau_DR_mu1), na.rm = T)/sum(as.vector(tau_DR_mu_denom1), na.rm = T))

tau2_DR <- (sum(as.vector(tau_DR_num2), na.rm = T) / sum(as.vector(tau_DR_denom2), na.rm = T)) +
  (sum(as.vector(tau_DR_mu2), na.rm = T)/sum(as.vector(tau_DR_mu_denom2), na.rm = T))



WD_sim_DR[sim_count] <- tau1_DR - tau2_DR
WR_sim_DR[sim_count] <- tau1_DR/tau2_DR



stopCluster(cl)




###############Doubly Robust (OW) Estimator################################




cl <- makeCluster(numCores - 1)  # Use all cores except one
registerDoParallel(cl)



compare_rows_DROW <- function(i, df_comb, pi_func, cond_prob_trt, cond_prob_ctrl) {
  n_count <- nrow(df_comb)
  
  DROW_num1_vec <- numeric(n_count)
  DROW_denom1_vec <- numeric(n_count)
  DROW_num2_vec <- numeric(n_count)
  DROW_denom2_vec <- numeric(n_count)
  DROW_mu1_vec <- numeric(n_count)
  DROW_mu2_vec <- numeric(n_count)
  DROW_mu_denom1_vec <- numeric(n_count)
  DROW_mu_denom2_vec <- numeric(n_count)
  
  for(j in 1:n_count){
    if(i != j){
      DROW_num1_vec[j] = (df_comb[i,'treatment'] * (1 - df_comb[j,'treatment'])) * ((1-pi_func[i]) * pi_func[j]) *
        ((df_comb[i,1] > df_comb[j,1]) - 
           (cond_prob_trt[i, 2] * cond_prob_ctrl[j, 1] + 
              cond_prob_trt[i, 3] * cond_prob_ctrl[j, 1] +
              cond_prob_trt[i, 3] * cond_prob_ctrl[j, 2]))
      
      DROW_denom1_vec[j] = (df_comb[i,'treatment'] * (1 - df_comb[j,'treatment'])) * ((1-pi_func[i]) * pi_func[j])
      
      DROW_num2_vec[j] = (df_comb[i,'treatment'] * (1 - df_comb[j,'treatment'])) * ((1-pi_func[i]) * pi_func[j]) *
        ((df_comb[i,1] < df_comb[j,1]) - 
           (cond_prob_trt[i, 1] * cond_prob_ctrl[j, 2] + 
              cond_prob_trt[i, 1] * cond_prob_ctrl[j, 3] +
              cond_prob_trt[i, 2] * cond_prob_ctrl[j, 3]))
      
      DROW_denom2_vec[j] = (df_comb[i,'treatment'] * (1 - df_comb[j,'treatment'])) * ((1-pi_func[i]) * pi_func[j])
      
      DROW_mu1_vec[j] = (pi_func[i] * pi_func[j] * (1-pi_func[i]) *(1-pi_func[j])) * 
        (cond_prob_trt[i, 2] * cond_prob_ctrl[j, 1] + 
           cond_prob_trt[i, 3] * cond_prob_ctrl[j, 1] + 
           cond_prob_trt[i, 3] * cond_prob_ctrl[j, 2])
      
      DROW_mu2_vec[j] = (pi_func[i] * pi_func[j] * (1-pi_func[i]) *(1-pi_func[j])) * 
        (cond_prob_trt[i, 1] * cond_prob_ctrl[j, 2] + 
           cond_prob_trt[i, 1] * cond_prob_ctrl[j, 3] + 
           cond_prob_trt[i, 2] * cond_prob_ctrl[j, 3])
      
      DROW_mu_denom1_vec[j] = pi_func[i] * pi_func[j] * (1-pi_func[i]) *(1-pi_func[j])
      DROW_mu_denom2_vec[j] = pi_func[i] * pi_func[j] * (1-pi_func[i]) *(1-pi_func[j])
    }
  }
  
  return(list(DROW_num1 = DROW_num1_vec, DROW_denom1 = DROW_denom1_vec, 
              DROW_num2 = DROW_num2_vec, DROW_denom2 = DROW_denom2_vec, 
              DROW_mu1 = DROW_mu1_vec, DROW_mu_denom1 = DROW_mu_denom1_vec,
              DROW_mu2 = DROW_mu2_vec, DROW_mu_denom2 = DROW_mu_denom2_vec))
}





# Use the function to compute the matrices
# results_DROW <- lapply(1:n_count, compare_rows_DROW, 
#                        df_comb=df_comb, pi_func=pi_func, 
#                        cond_prob_trt=cond_prob_trt, cond_prob_ctrl=cond_prob_ctrl)

results_DROW <- foreach(i=1:n_count) %dopar% {
  compare_rows_DROW(i, df_comb, pi_func, cond_prob_trt, cond_prob_ctrl)
}


# Combine the results
tau_DROW_num1 <- do.call(rbind, lapply(results_DROW, function(x) x$DROW_num1))
tau_DROW_denom1 <- do.call(rbind, lapply(results_DROW, function(x) x$DROW_denom1))
tau_DROW_num2 <- do.call(rbind, lapply(results_DROW, function(x) x$DROW_num2))
tau_DROW_denom2 <- do.call(rbind, lapply(results_DROW, function(x) x$DROW_denom2))
tau_DROW_mu1 <- do.call(rbind, lapply(results_DROW, function(x) x$DROW_mu1))
tau_DROW_mu_denom1 <- do.call(rbind, lapply(results_DROW, function(x) x$DROW_mu_denom1))
tau_DROW_mu2 <- do.call(rbind, lapply(results_DROW, function(x) x$DROW_mu2))
tau_DROW_mu_denom2 <- do.call(rbind, lapply(results_DROW, function(x) x$DROW_mu_denom2))

tau1_DROW <- (sum(as.vector(tau_DROW_num1), na.rm = T) / sum(as.vector(tau_DROW_denom1), na.rm = T)) +
  (sum(as.vector(tau_DROW_mu1), na.rm = T)/sum(as.vector(tau_DROW_mu_denom1), na.rm = T))

tau2_DROW <- (sum(as.vector(tau_DROW_num2), na.rm = T) / sum(as.vector(tau_DROW_denom2), na.rm = T)) +
  (sum(as.vector(tau_DROW_mu2), na.rm = T)/sum(as.vector(tau_DROW_mu_denom2), na.rm = T))





WD_sim_DROW[sim_count] <- tau1_DROW - tau2_DROW
WR_sim_DROW[sim_count] <- tau1_DROW/tau2_DROW



stopCluster(cl)




}#open on line 145









#################Unadjusted Estimator Output#########################


###Difference###
WD_sim_noninf <- WD_sim[WD_sim!=Inf &
                        !is.na(WD_sim)]

empirical_var_Ori_WD[count_temp] <- var(WD_sim_noninf)
bias_Ori_WD[count_temp] <- mean(WD_sim) - trueWD
pct_bias_Ori_WD[count_temp] <- (bias_Ori_WD[count_temp]/trueWD)*100
mse_Ori_WD[count_temp] <- bias_Ori_WD[count_temp]^2 + empirical_var_Ori_WD[count_temp]




###Ratio###
WR_sim_noninf <- WR_sim[WR_sim!=Inf &
                        WR_sim!=0 &
                        !is.na(WR_sim)]

empirical_var_Ori_WR[count_temp] <- var(WR_sim_noninf)
bias_Ori_WR[count_temp] <- mean(WR_sim) - trueWR
pct_bias_Ori_WR[count_temp] <- (bias_Ori_WR[count_temp]/trueWR)*100
mse_Ori_WR[count_temp] <- bias_Ori_WR[count_temp]^2 + empirical_var_Ori_WR[count_temp]






###############IPW Estimator Output###########################

###Difference###
WD_sim_noninf_IPW <- WD_sim_IPW[WD_sim_IPW!=Inf &
                                !is.na(WD_sim_IPW)]

empirical_var_IPW_WD[count_temp] <- var(WD_sim_noninf_IPW)
bias_IPW_WD[count_temp] <- mean(WD_sim_IPW) - trueWD
pct_bias_IPW_WD[count_temp] <- (bias_IPW_WD[count_temp]/trueWD)*100
mse_IPW_WD[count_temp] <- bias_IPW_WD[count_temp]^2 + empirical_var_IPW_WD[count_temp]



###Ratio###
WR_sim_noninf_IPW <- WR_sim_IPW[WR_sim_IPW!=Inf &
                                WR_sim_IPW!=0 &
                                !is.na(WR_sim_IPW)]

empirical_var_IPW_WR[count_temp] <- var(WR_sim_noninf_IPW)
bias_IPW_WR[count_temp] <- mean(WR_sim_IPW) - trueWR
pct_bias_IPW_WR[count_temp] <- (bias_IPW_WR[count_temp]/trueWR)*100
mse_IPW_WR[count_temp] <- bias_IPW_WR[count_temp]^2 + empirical_var_IPW_WR[count_temp]








###############OW Estimator Output###########################


###Difference###
WD_sim_noninf_OW <- WD_sim_OW[WD_sim_OW!=Inf &
                              !is.na(WD_sim_OW)]

empirical_var_OW_WD[count_temp] <- var(WD_sim_noninf_OW)
bias_OW_WD[count_temp] <- mean(WD_sim_OW) - trueWD
pct_bias_OW_WD[count_temp] <- (bias_OW_WD[count_temp]/trueWD)*100
mse_OW_WD[count_temp] <- bias_OW_WD[count_temp]^2 + empirical_var_OW_WD[count_temp]




###Ratio###
WR_sim_noninf_OW <- WR_sim_OW[WR_sim_OW!=Inf &
                              WR_sim_OW!=0 &
                              !is.na(WR_sim_OW)]

empirical_var_OW_WR[count_temp] <- var(WR_sim_noninf_OW)
bias_OW_WR[count_temp] <- mean(WR_sim_OW) - trueWR
pct_bias_OW_WR[count_temp] <- (bias_OW_WR[count_temp]/trueWR)*100
mse_OW_WR[count_temp] <- bias_OW_WR[count_temp]^2 + empirical_var_OW_WR[count_temp]






###############DR(IPW) Estimator Output###########################


###Difference###
WD_sim_noninf_DR <- WD_sim_DR[WD_sim_DR!=Inf &
                              WD_sim_DR!=0 &
                              !is.na(WD_sim_DR)]

empirical_var_DR_WD[count_temp] <- var(WD_sim_noninf_DR)
bias_DR_WD[count_temp] <- mean(WD_sim_DR) - trueWD
pct_bias_DR_WD[count_temp] <- (bias_DR_WD[count_temp]/trueWD)*100
mse_DR_WD[count_temp] <- bias_DR_WD[count_temp]^2 + empirical_var_DR_WD[count_temp]



###Ratio###
WR_sim_noninf_DR <- WR_sim_DR[WR_sim_DR!=Inf &
                              WR_sim_DR!=0 &
                              !is.na(WR_sim_DR)]

empirical_var_DR_WR[count_temp] <- var(WR_sim_noninf_DR)
bias_DR_WR[count_temp] <- mean(WR_sim_DR) - trueWR
pct_bias_DR_WR[count_temp] <- (bias_DR_WR[count_temp]/trueWR)*100
mse_DR_WR[count_temp] <- bias_DR_WR[count_temp]^2 + empirical_var_DR_WR[count_temp]






###############DR(OW) Estimator Output###########################


###Difference###
WD_sim_noninf_DROW <- WD_sim_DROW[WD_sim_DROW!=Inf &
                                  WD_sim_DROW!=0 &
                                  !is.na(WD_sim_DROW)]

empirical_var_DROW_WD[count_temp] <- var(WD_sim_noninf_DROW)
bias_DROW_WD[count_temp] <- mean(WD_sim_DROW) - trueWD
pct_bias_DROW_WD[count_temp] <- (bias_DROW_WD[count_temp]/trueWD)*100
mse_DROW_WD[count_temp] <- bias_DROW_WD[count_temp]^2 + empirical_var_DROW_WD[count_temp]



###Ratio###
WR_sim_noninf_DROW <- WR_sim_DROW[WR_sim_DROW!=Inf &
                                  WR_sim_DROW!=0 &
                                  !is.na(WR_sim_DROW)]

empirical_var_DROW_WR[count_temp] <- var(WR_sim_noninf_DROW)
bias_DROW_WR[count_temp] <- mean(WR_sim_DROW) - trueWR
pct_bias_DROW_WR[count_temp] <- (bias_DROW_WR[count_temp]/trueWR)*100
mse_DROW_WR[count_temp] <- bias_DROW_WR[count_temp]^2 + empirical_var_DROW_WR[count_temp]






}#open on line 125




relative_eff_IPW_WD <- empirical_var_Ori_WD/empirical_var_IPW_WD
relative_eff_IPW_WR <- empirical_var_Ori_WR/empirical_var_IPW_WR

relative_eff_OW_WD <- empirical_var_Ori_WD/empirical_var_OW_WD
relative_eff_OW_WR <- empirical_var_Ori_WR/empirical_var_OW_WR

relative_eff_DR_WD <- empirical_var_Ori_WD/empirical_var_DR_WD
relative_eff_DR_WR <- empirical_var_Ori_WR/empirical_var_DR_WR

relative_eff_DROW_WD <- empirical_var_Ori_WD/empirical_var_DROW_WD
relative_eff_DROW_WR <- empirical_var_Ori_WR/empirical_var_DROW_WR


re_df_WD <- as.data.frame(cbind(relative_eff_IPW_WD,
                            relative_eff_OW_WD,
                            relative_eff_DR_WD,
                            relative_eff_DROW_WD))


re_df_WR <- as.data.frame(cbind(relative_eff_IPW_WR,
                            relative_eff_OW_WR,
                            relative_eff_DR_WR,
                            relative_eff_DROW_WR))


Ori_df <- as.data.frame(cbind(c(empirical_var_Ori_WD, empirical_var_Ori_WR),
                          c(pct_bias_Ori_WD, pct_bias_Ori_WR),
                          c(mse_Ori_WD, mse_Ori_WR)))


IPW_df <- as.data.frame(cbind(c(empirical_var_IPW_WD, empirical_var_IPW_WR),
                          c(pct_bias_IPW_WD, pct_bias_IPW_WR),
                          c(mse_IPW_WD, mse_IPW_WR)))


OW_df <- as.data.frame(cbind(c(empirical_var_OW_WD, empirical_var_OW_WR),
                         c(pct_bias_OW_WD, pct_bias_OW_WR),
                         c(mse_OW_WD, mse_OW_WR)))


DR_df <- as.data.frame(cbind(c(empirical_var_DR_WD, empirical_var_DR_WR),
                         c(pct_bias_DR_WD, pct_bias_DR_WR),
                         c(mse_DR_WD, mse_DR_WR)))


DROW_df <- as.data.frame(cbind(c(empirical_var_DROW_WD, empirical_var_DROW_WR),
                           c(pct_bias_DROW_WD, pct_bias_DROW_WR),
                           c(mse_DROW_WD, mse_DROW_WR)))


end_time <- Sys.time()
run_time <- end_time - start_time
run_time#3.068355 per sim_num


write.table(re_df_WD,
        file=paste("results/sim_re_WD_newscenario",scenario,".txt",sep=""), 
        sep="\t", row.names=F)


write.table(re_df_WR,
        file=paste("results/sim_re_WR_newscenario",scenario,".txt",sep=""), 
        sep="\t", row.names=F)


write.table(Ori_df,
        file=paste("results/Unadjusted_newscenario",scenario,".txt",sep=""), 
        sep="\t", row.names=F)


write.table(IPW_df,
        file=paste("results/IPW_newscenario",scenario,".txt",sep=""), 
        sep="\t", row.names=F)


write.table(OW_df,
        file=paste("results/OW_newscenario",scenario,".txt",sep=""), 
        sep="\t", row.names=F)


write.table(DR_df,
        file=paste("results/DR_newscenario",scenario,".txt",sep=""), 
        sep="\t", row.names=F)


write.table(DROW_df,
        file=paste("results/DROW_newscenario",scenario,".txt",sep=""), 
        sep="\t", row.names=F)








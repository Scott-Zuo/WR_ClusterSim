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
################# Variance Coverage of AOW ##################
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

WP_trt_AOW_coverage_classIII <- numeric(re_count)
WP_ctrl_AOW_coverage_classIII <- numeric(re_count)
WR_AOW_coverage <- numeric(re_count)

WP_trt_AOW_var_ratio_classIII <- numeric(re_count)
WP_ctrl_AOW_var_ratio_classIII <- numeric(re_count)
WR_AOW_var_ratio <- numeric(re_count)

WP_trt_AOW_theoretical_var_mean_classIII <- numeric(re_count)
WP_trt_AOW_empirical_var <- numeric(re_count)

WP_ctrl_AOW_theoretical_var_mean_classIII <- numeric(re_count)
WP_ctrl_AOW_empirical_var <- numeric(re_count)

WR_AOW_theoretical_var_mean <- numeric(re_count)
WR_AOW_empirical_var <- numeric(re_count)

WP_trt_sim_true <- 0.2953701
WP_ctrl_sim_true <- 0.2050298
WR_sim_true <- 1.447066






var_ratio = function(uR,uS,sigR2,sigS2,cov_RS){
  WR_appro_var = (uR^2/uS^2)*(sigR2/uR^2-2*cov_RS/(uR*uS)+sigS2/uS^2)
  return(WR_appro_var)
}

sample_size_list <- c(200,300,400)




for(size_count in 1:re_count){

WP_trt_sim_AOW <- numeric(sim_num)
WP_ctrl_sim_AOW <- numeric(sim_num)
WR_sim_AOW <- numeric(sim_num)

theory_var_AOW_trt_classIII <- numeric(sim_num)
theory_var_AOW_ctrl_classIII <- numeric(sim_num)
cov_trt_ctrl <- numeric(sim_num)
theory_var_AOW_WR <- numeric(sim_num)

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
treatment_assignment <- rbinom(n_count, 1, 0.5)

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





##########################AOW########################################

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


### mu-function

if((!(1 %in% df_comb[df_comb$treatment==1,]$outcomes_comb) |
    !(2 %in% df_comb[df_comb$treatment==1,]$outcomes_comb) |
    !(3 %in% df_comb[df_comb$treatment==1,]$outcomes_comb))|
   (!(1 %in% df_comb[df_comb$treatment==0,]$outcomes_comb) |
    !(2 %in% df_comb[df_comb$treatment==0,]$outcomes_comb) |
    !(3 %in% df_comb[df_comb$treatment==0,]$outcomes_comb))){
  next
}


else {
  lp_start <- rep(0, 6)
  th_start <- c(-1, 1)
  start_values <- c(lp_start, th_start)
  dr_trt <- polr(factor(outcomes_comb) ~ x1 + x2 + x3 + x4 + x5 + x6, 
                 data = df_comb[df_comb$treatment == 1,],
                 start = start_values)
  cond_prob_trt <- predict(dr_trt, newdata = df_comb, 
                           type = "probs")
  
  #dr_ctrl <- multinom(outcomes_comb ~ x1 + x2 + x3 + x4 + x5 + x6 + 
  #                      I(x1^2) + I(x2^2) + I(x3^2) + I(x4^2) + I(x5^2) + I(x6^2),
  #                    data = df_comb[df_comb$treatment == 0,], trace = FALSE)
  dr_ctrl <- polr(factor(outcomes_comb) ~ x1 + x2 + x3 + x4 + x5 + x6, 
                  data = df_comb[df_comb$treatment == 0,],
                  start = start_values)
  
  cond_prob_ctrl <- predict(dr_ctrl, newdata = df_comb, 
                            type = "probs")
}



cl <- makeCluster(numCores - 1)  # Use all cores except one
registerDoParallel(cl)



compare_rows_AOW <- function(i, df_comb, pi_func, cond_prob_trt, cond_prob_ctrl) {
  n_count <- nrow(df_comb)
  
  AOW_num1_vec <- numeric(n_count)
  AOW_denom1_vec <- numeric(n_count)
  AOW_num2_vec <- numeric(n_count)
  AOW_denom2_vec <- numeric(n_count)
  AOW_mu1_vec <- numeric(n_count)
  AOW_mu2_vec <- numeric(n_count)
  AOW_mu_denom1_vec <- numeric(n_count)
  AOW_mu_denom2_vec <- numeric(n_count)
  
  for(j in 1:(n_count)){
    if(i != j){
      AOW_num1_vec[j] = (df_comb[i,'treatment'] * (1 - df_comb[j,'treatment'])) * ((1-pi_func[i]) * pi_func[j]) *
        ((df_comb[i,1] > df_comb[j,1]) - 
           (cond_prob_trt[i, 2] * cond_prob_ctrl[j, 1] + 
              cond_prob_trt[i, 3] * cond_prob_ctrl[j, 1] +
              cond_prob_trt[i, 3] * cond_prob_ctrl[j, 2]))
      
      AOW_denom1_vec[j] = (df_comb[i,'treatment'] * (1 - df_comb[j,'treatment'])) * ((1-pi_func[i]) * pi_func[j])
      
      AOW_num2_vec[j] = (df_comb[i,'treatment'] * (1 - df_comb[j,'treatment'])) * ((1-pi_func[i]) * pi_func[j]) *
        ((df_comb[i,1] < df_comb[j,1]) - 
           (cond_prob_trt[i, 1] * cond_prob_ctrl[j, 2] + 
              cond_prob_trt[i, 1] * cond_prob_ctrl[j, 3] +
              cond_prob_trt[i, 2] * cond_prob_ctrl[j, 3]))
      
      AOW_denom2_vec[j] = (df_comb[i,'treatment'] * (1 - df_comb[j,'treatment'])) * ((1-pi_func[i]) * pi_func[j])
      
      AOW_mu1_vec[j] = (pi_func[i] * pi_func[j] * (1-pi_func[i]) *(1-pi_func[j])) * 
        (cond_prob_trt[i, 2] * cond_prob_ctrl[j, 1] + 
           cond_prob_trt[i, 3] * cond_prob_ctrl[j, 1] + 
           cond_prob_trt[i, 3] * cond_prob_ctrl[j, 2])
      
      AOW_mu2_vec[j] = (pi_func[i] * pi_func[j] * (1-pi_func[i]) *(1-pi_func[j])) * 
        (cond_prob_trt[i, 1] * cond_prob_ctrl[j, 2] + 
           cond_prob_trt[i, 1] * cond_prob_ctrl[j, 3] + 
           cond_prob_trt[i, 2] * cond_prob_ctrl[j, 3])
      
      AOW_mu_denom1_vec[j] = pi_func[i] * pi_func[j] * (1-pi_func[i]) *(1-pi_func[j])
      AOW_mu_denom2_vec[j] = pi_func[i] * pi_func[j] * (1-pi_func[i]) *(1-pi_func[j])
    }
  }
  
  return(list(AOW_num1 = AOW_num1_vec, AOW_denom1 = AOW_denom1_vec, 
              AOW_num2 = AOW_num2_vec, AOW_denom2 = AOW_denom2_vec, 
              AOW_mu1 = AOW_mu1_vec, AOW_mu_denom1 = AOW_mu_denom1_vec,
              AOW_mu2 = AOW_mu2_vec, AOW_mu_denom2 = AOW_mu_denom2_vec))
}





# Use the function to compute the matrices
# results_AOW <- lapply(1:(2*n_count), compare_rows_AOW, 
#                        df_comb=df_comb, pi_func=pi_func, 
#                        cond_prob_trt=cond_prob_trt, cond_prob_ctrl=cond_prob_ctrl)

results_AOW <- foreach(i=1:n_count) %dopar% {
  compare_rows_AOW(i, df_comb, pi_func, cond_prob_trt, cond_prob_ctrl)
}


# Combine the results
tau_AOW_num1 <- do.call(rbind, lapply(results_AOW, function(x) x$AOW_num1))
tau_AOW_denom1 <- do.call(rbind, lapply(results_AOW, function(x) x$AOW_denom1))
tau_AOW_num2 <- do.call(rbind, lapply(results_AOW, function(x) x$AOW_num2))
tau_AOW_denom2 <- do.call(rbind, lapply(results_AOW, function(x) x$AOW_denom2))
tau_AOW_mu1 <- do.call(rbind, lapply(results_AOW, function(x) x$AOW_mu1))
tau_AOW_mu_denom1 <- do.call(rbind, lapply(results_AOW, function(x) x$AOW_mu_denom1))
tau_AOW_mu2 <- do.call(rbind, lapply(results_AOW, function(x) x$AOW_mu2))
tau_AOW_mu_denom2 <- do.call(rbind, lapply(results_AOW, function(x) x$AOW_mu_denom2))

tau1_AOW <- (sum(as.vector(tau_AOW_num1), na.rm = T) / sum(as.vector(tau_AOW_denom1), na.rm = T)) +
  (sum(as.vector(tau_AOW_mu1), na.rm = T)/sum(as.vector(tau_AOW_mu_denom1), na.rm = T))

tau2_AOW <- (sum(as.vector(tau_AOW_num2), na.rm = T) / sum(as.vector(tau_AOW_denom2), na.rm = T)) +
  (sum(as.vector(tau_AOW_mu2), na.rm = T)/sum(as.vector(tau_AOW_mu_denom2), na.rm = T))


# Update WD_sim_AOW and WR_sim_AOW
WP_trt_sim_AOW[count_temp] <- tau1_AOW
WP_ctrl_sim_AOW[count_temp] <- tau2_AOW
WR_sim_AOW[count_temp] <- tau1_AOW / tau2_AOW




stopCluster(cl)




########### AOW Theoretical Variance###########




mu_func_trt <- function(a, b, cond_prob_trt, cond_prob_ctrl) {
mu_a_b <- (
  cond_prob_trt[a, 2] * cond_prob_ctrl[b, 1] +
    cond_prob_trt[a, 3] * (cond_prob_ctrl[b, 1] + cond_prob_ctrl[b, 2])
)

  return(mu_a_b)
}


# mu_func_trt <- function(a, b, cond_prob_trt, cond_prob_ctrl) {
#   ca <- if (df_comb$treatment[a] == 1) cond_prob_trt[a, ] else cond_prob_ctrl[a, ]
#   cb <- if (df_comb$treatment[b] == 1) cond_prob_trt[b, ] else cond_prob_ctrl[b, ]
#   mu_a_b <- (
#     (1 - ca[1]) * cb[1] + ca[3] * cb[2]
#   )
#   return(mu_a_b)
# }


mu_func_ctrl <- function(a, b, cond_prob_trt, cond_prob_ctrl) {
  mu_a_b <- (
    cond_prob_trt[a, 2] * cond_prob_ctrl[b, 3]+
      cond_prob_trt[a, 1] * (cond_prob_ctrl[b, 2] + cond_prob_ctrl[b, 3])
  )
  
  return(mu_a_b)
}


# mu_func_ctrl <- function(a, b, cond_prob_trt, cond_prob_ctrl) {
#   ca <- if (df_comb$treatment[a] == 1) cond_prob_trt[a, ] else cond_prob_ctrl[a, ]
#   cb <- if (df_comb$treatment[b] == 1) cond_prob_trt[b, ] else cond_prob_ctrl[b, ]
#   mu_a_b <- (
#     ca[1] * (cb[2] + cb[3]) + ca[2] * cb[3]
#   )
#   return(mu_a_b)
# }







###g1_AOW_trt###

n <- nrow(df_comb)
g1_AOW_trt_fun <- function(i, outcomes, treatment, pi_func) {
  n <- nrow(df_comb)
  pairwise_comparisons_AOW_term1 <- (1/2) * (
    ((treatment[i] * (1 - treatment) * (1 - pi_func[i]) * pi_func) * 
       (as.numeric(outcomes[i] > outcomes) - sapply(1:n, function(b) mu_func_trt(i, b, cond_prob_trt, cond_prob_ctrl)))) +
      (((1 - treatment[i]) * treatment * pi_func[i] * (1 - pi_func)) * 
         (as.numeric(outcomes > outcomes[i]) - sapply(1:n, function(a) mu_func_trt(a, i, cond_prob_trt, cond_prob_ctrl))))
  )
  pairwise_comparisons_AOW_term2 <- (1/2) * pi_func[i] * (1 - pi_func[i]) *
    pi_func * (1 - pi_func) * (sapply(1:n, function(b) mu_func_trt(i, b, cond_prob_trt, cond_prob_ctrl)) +
                                 sapply(1:n, function(a) mu_func_trt(a, i, cond_prob_trt, cond_prob_ctrl)))
  #pairwise_comparisons_AOW_term1[i] <- 0
  #pairwise_comparisons_AOW_term2[i] <- 0# condition i!=j
  
  # Sum and adjust according to the formula
  adjusted_denom1_term <- sum(sapply(1:n, function(j) {sum(treatment[-j]*(1-pi_func[-j])) * 
      (1 - treatment[j]) * pi_func[j]}))
  g1_denominator1_AOW <- (1/(n*(n-1))) * adjusted_denom1_term
  
  adjusted_denom2_term <- sum(sapply(1:n, function(j) {sum(pi_func[-j]*(1-pi_func[-j])) * 
      (pi_func[j]) * (1-pi_func[j])}))
  g1_denominator2_AOW <- (1/(n*(n-1))) * adjusted_denom2_term
  
  sum_g1_AOW_trt <- (sum(pairwise_comparisons_AOW_term1) / g1_denominator1_AOW) +
    (sum(pairwise_comparisons_AOW_term2) / g1_denominator2_AOW)
  return((1/(n-1)) * sum_g1_AOW_trt)
  
}


# Parallelize the computation for g1_AOW_trt
results_g1_AOW_trt <- mclapply(1:n, g1_AOW_trt_fun, df_comb$outcomes_comb, 
                               df_comb$treatment, pi_func, mc.cores = numCores - 1)

g1_AOW_trt <- unlist(results_g1_AOW_trt)


g1_AOW_ctrl_fun <- function(i, outcomes, treatment, pi_func) {
  n <- nrow(df_comb)
  pairwise_comparisons_AOW_term1 <- (1/2) * (
    ((treatment[i] * (1 - treatment) * (1 - pi_func[i]) * pi_func) * 
       (as.numeric(outcomes[i] < outcomes) - sapply(1:n, function(b) mu_func_ctrl(i, b, cond_prob_trt, cond_prob_ctrl)))) +
      (((1 - treatment[i]) * treatment * pi_func[i] * (1 - pi_func)) * 
         (as.numeric(outcomes < outcomes[i]) - sapply(1:n, function(a) mu_func_ctrl(a, i, cond_prob_trt, cond_prob_ctrl))))
  )
  pairwise_comparisons_AOW_term2 <- (1/2) * pi_func[i] * (1 - pi_func[i]) *
    pi_func * (1 - pi_func) * (sapply(1:n, function(b) mu_func_ctrl(i, b, cond_prob_trt, cond_prob_ctrl)) +
                                 sapply(1:n, function(a) mu_func_ctrl(a, i, cond_prob_trt, cond_prob_ctrl)))
  #pairwise_comparisons_AOW_term1[i] <- 0
  #pairwise_comparisons_AOW_term2[i] <- 0# condition i!=j
  
  # Sum and adjust according to the formula
  adjusted_denom1_term <- sum(sapply(1:n, function(j) {sum(treatment[-j]*(1-pi_func[-j])) * 
      (1 - treatment[j]) * pi_func[j]}))
  g1_denominator1_AOW <- (1/(n*(n-1))) * adjusted_denom1_term
  
  adjusted_denom2_term <- sum(sapply(1:n, function(j) {sum(pi_func[-j]*(1-pi_func[-j])) * 
      (pi_func[j]) * (1-pi_func[j])}))
  g1_denominator2_AOW <- (1/(n*(n-1))) * adjusted_denom2_term
  
  sum_g1_AOW_trt <- (sum(pairwise_comparisons_AOW_term1) / g1_denominator1_AOW) +
    (sum(pairwise_comparisons_AOW_term2) / g1_denominator2_AOW)
  return((1/(n-1)) * sum_g1_AOW_trt)
  
}


results_g1_AOW_ctrl <- mclapply(1:n, g1_AOW_ctrl_fun, df_comb$outcomes_comb, 
                                df_comb$treatment, pi_func, mc.cores = numCores - 1)

g1_AOW_ctrl <- unlist(results_g1_AOW_ctrl)






### B ###


### B1_trt & B1_ctrl ### (Dim: 6*1)

B1_AOW_trt <- function(df, pi_func, cond_prob_trt, cond_prob_ctrl) {
  
  covariate_names <- c("x1", "x2", "x3", "x4", "x5", "x6")
  n <- nrow(df)
  B1_values <- mclapply(1:n, function(j) {
    Xj <- as.matrix(df[j, covariate_names])
    numerator_vec <- numeric(length(covariate_names))
    denominator <- 0
    # Loop through all control indices and calculate contributions to numerator and denominator
    for (i in 1:n) {
      if (i != j) {
        Xi <- as.matrix(df[i, covariate_names])
        e_Xi <- pi_func[i]
        e_Xj <- pi_func[j]
        mu_ij <- mu_func_trt(i, j, cond_prob_trt, cond_prob_ctrl)
        term <- (df$treatment[i] * (1 - df$treatment[j])) * (1 - e_Xi) * e_Xj * 
          ((1-e_Xj) * Xj - e_Xi * Xi) * (as.numeric(df$outcomes[i] > df$outcomes[j]) - mu_ij)
        numerator_vec <- numerator_vec + term
        denominator <- denominator + df$treatment[i] * (1-e_Xi)
      }
    }
    denominator <- denominator * (1 - df$treatment[j]) * pi_func[j]
    return(list(numerator = numerator_vec, denominator = denominator))
  }, mc.cores = numCores - 1)
  
  total_numerator <- Reduce("+", lapply(B1_values, `[[`, "numerator"))
  denominator <- sum(sapply(B1_values, `[[`, "denominator"))
  B1_vector <- total_numerator / denominator
  
  return(B1_vector)
}




B1_AOW_ctrl <- function(df, pi_func, cond_prob_trt, cond_prob_ctrl) {
  
  covariate_names <- c("x1", "x2", "x3", "x4", "x5", "x6")
  n <- nrow(df)
  B1_values <- mclapply(1:n, function(j) {
    Xj <- as.matrix(df[j, covariate_names])
    numerator_vec <- numeric(length(covariate_names))
    denominator <- 0
    # Loop through all control indices and calculate contributions to numerator and denominator
    for (i in 1:n) {
      if (i != j) {
        Xi <- as.matrix(df[i, covariate_names])
        e_Xi <- pi_func[i]
        e_Xj <- pi_func[j]
        mu_ij <- mu_func_ctrl(i, j, cond_prob_trt, cond_prob_ctrl)
        term <- (df$treatment[i] * (1 - df$treatment[j])) * (1 - e_Xi) * e_Xj * 
          ((1-e_Xj) * Xj - e_Xi * Xi) * (as.numeric(df$outcomes[i] < df$outcomes[j]) - mu_ij)
        numerator_vec <- numerator_vec + term
        denominator <- denominator + df$treatment[i] * (1-e_Xi)
      }
    }
    denominator <- denominator * (1 - df$treatment[j]) * pi_func[j]
    return(list(numerator = numerator_vec, denominator = denominator))
  }, mc.cores = numCores - 1)
  
  total_numerator <- Reduce("+", lapply(B1_values, `[[`, "numerator"))
  denominator <- sum(sapply(B1_values, `[[`, "denominator"))
  B1_vector <- total_numerator / denominator
  
  return(B1_vector)
}


### B2_trt & B2_ctrl ### (Dim: 1*1)

B2_AOW_trt <- function(df, pi_func, cond_prob_trt, cond_prob_ctrl) {
  n <- nrow(df)
  B2_values <- mclapply(1:n, function(j) {
    numerator_vec <- 0
    denominator <- 0
    # Loop through all control indices and calculate contributions to numerator and denominator
    for (i in 1:n) {
      if (i != j) {
        e_Xi <- pi_func[i]
        e_Xj <- pi_func[j]
        mu_ij <- mu_func_trt(i, j, cond_prob_trt, cond_prob_ctrl)
        term <- (df$treatment[i] * (1 - df$treatment[j])) * (1-e_Xi) * e_Xj *
          (as.numeric(df$outcomes[i] > df$outcomes[j]) - mu_ij)
        numerator_vec <- numerator_vec + term
        denominator <- denominator + df$treatment[i] * (1-e_Xi)
      }
    }
    denominator <- denominator * (1 - df$treatment[j]) * pi_func[j]
    return(list(numerator = numerator_vec, denominator = denominator))
  }, mc.cores = numCores - 1)
  
  total_numerator <- Reduce("+", lapply(B2_values, `[[`, "numerator"))
  denominator <- sum(sapply(B2_values, `[[`, "denominator"))
  B2_vector <- total_numerator / denominator
  return(B2_vector)
}




B2_AOW_ctrl <- function(df, pi_func, cond_prob_trt, cond_prob_ctrl) {
  n <- nrow(df)
  B2_values <- mclapply(1:n, function(j) {
    numerator_vec <- 0
    denominator <- 0
    # Loop through all control indices and calculate contributions to numerator and denominator
    for (i in 1:n) {
      if (i != j) {
        e_Xi <- pi_func[i]
        e_Xj <- pi_func[j]
        mu_ij <- mu_func_ctrl(i, j, cond_prob_trt, cond_prob_ctrl)
        term <- (df$treatment[i] * (1 - df$treatment[j])) * (1-e_Xi) * e_Xj *
          (as.numeric(df$outcomes[i] < df$outcomes[j]) - mu_ij)
        numerator_vec <- numerator_vec + term
        denominator <- denominator + df$treatment[i] * (1-e_Xi)
      }
    }
    denominator <- denominator * (1 - df$treatment[j]) * pi_func[j]
    return(list(numerator = numerator_vec, denominator = denominator))
  }, mc.cores = numCores - 1)
  
  total_numerator <- Reduce("+", lapply(B2_values, `[[`, "numerator"))
  denominator <- sum(sapply(B2_values, `[[`, "denominator"))
  B2_vector <- total_numerator / denominator
  return(B2_vector)
}




### B3 ### (Dim: 6*1)

B3_AOW_fun <- function(df, pi_func) {
  n <- nrow(df)
  covariate_names <- c("x1", "x2", "x3", "x4", "x5", "x6")
  B3_term1 <- mclapply(1:n, function(j) {
    Xj <- as.matrix(df[j, covariate_names])
    term1_vec <- numeric(length(covariate_names))
    for (i in 1:n) {
      Xi <- as.matrix(df[i, covariate_names])
      term1_vec <- term1_vec + (-df$treatment[i]) * pi_func[i] * (1-pi_func[i]) * Xi
    }
    term1_vec <- term1_vec * (1-df$treatment[j]) * pi_func[j]
    return(list(term1 = term1_vec))
  }, mc.cores = numCores - 1)
  
  B3_term2 <- mclapply(1:n, function(j) {
    Xj <- as.matrix(df[j, covariate_names])
    term2_vec <- numeric(length(covariate_names))
    for (i in 1:n) {
      term2_vec <- term2_vec + (df$treatment[i] * (1-pi_func[i]))
    }
    term2_vec <- term2_vec * (1 - df$treatment[j]) * pi_func[j] * (1 - pi_func[j]) * Xj
    return(list(term2 = term2_vec))
  }, mc.cores = numCores - 1)
  
  total_term1 <- Reduce("+", lapply(B3_term1, `[[`, "term1"))
  total_term2 <- Reduce("+", lapply(B3_term2, `[[`, "term2"))
  return(total_term1 + total_term2)  
  
}


### B4 ### (Dim: 1*1)

B4_AOW_fun <- function(df, pi_func) {
  n <- nrow(df)
  B4_term <- mclapply(1:n, function(j) {
    term_val <- 0
    for (i in 1:n) {
      if (i != j) {
        term_val <- term_val + df$treatment[i] * (1-pi_func[i])
      }
    }
    term_val <- term_val * (1 - df$treatment[j]) * pi_func[j]
    return(list(B4_term_val = term_val))
  }, mc.cores = numCores - 1)
  
  B4_AOW <- Reduce("+", lapply(B4_term, `[[`, "B4_term_val"))
  return(B4_AOW)  
}


### B5 ### (Dim: 6*1)

B5_AOW_trt <- function(df, pi_func, cond_prob_trt, cond_prob_ctrl) {
  covariate_names <- c("x1", "x2", "x3", "x4", "x5", "x6")
  n <- nrow(df)
  B5_values <- mclapply(1:n, function(j) {
    numerator_vec <- numeric(length(covariate_names))
    denominator_val <- 0
    # Loop through all control indices and calculate contributions to numerator and denominator
    for (i in 1:n) {
      if (i != j) {
        Xi <- as.matrix(df[i, covariate_names])
        e_Xi <- pi_func[i]
        e_Xj <- pi_func[j]
        h_ij <- e_Xi * (1-e_Xi) * e_Xj * (1-e_Xj)
        mu_ij <- mu_func_trt(i, j, cond_prob_trt, cond_prob_ctrl)
        term <- h_ij * ((1 - 2 * e_Xi) * Xi) * mu_ij
        numerator_vec <- numerator_vec + term
        denominator_val <- denominator_val + h_ij
      }
    }
    return(list(numerator = numerator_vec, denominator = denominator_val))
  }, mc.cores = numCores - 1)
  
  total_numerator <- Reduce("+", lapply(B5_values, `[[`, "numerator"))
  denominator <- sum(sapply(B5_values, `[[`, "denominator"))
  B5_vector <- total_numerator / denominator
  
  return(B5_vector)
  
}

B5_AOW_ctrl <- function(df, pi_func, cond_prob_trt, cond_prob_ctrl) {
  covariate_names <- c("x1", "x2", "x3", "x4", "x5", "x6")
  n <- nrow(df)
  B5_values <- mclapply(1:n, function(j) {
    numerator_vec <- numeric(length(covariate_names))
    denominator_val <- 0
    # Loop through all control indices and calculate contributions to numerator and denominator
    for (i in 1:n) {
      if (i != j) {
        Xi <- as.matrix(df[i, covariate_names])
        e_Xi <- pi_func[i]
        e_Xj <- pi_func[j]
        h_ij <- e_Xi * (1-e_Xi) * e_Xj * (1-e_Xj)
        mu_ij <- mu_func_ctrl(i, j, cond_prob_trt, cond_prob_ctrl)
        term <- h_ij * ((1 - 2 * e_Xi) * Xi) * mu_ij
        numerator_vec <- numerator_vec + term
        denominator_val <- denominator_val + h_ij
      }
    }
    return(list(numerator = numerator_vec, denominator = denominator_val))
  }, mc.cores = numCores - 1)
  
  total_numerator <- Reduce("+", lapply(B5_values, `[[`, "numerator"))
  denominator <- sum(sapply(B5_values, `[[`, "denominator"))
  B5_vector <- total_numerator / denominator
  
  return(B5_vector)
  
}




### B6 ### (Dim: 1*1)

B6_AOW_trt <- function(df, pi_func, cond_prob_trt, cond_prob_ctrl) {
  n <- nrow(df)
  B6_values <- mclapply(1:n, function(j) {
    numerator_val <- 0
    denominator_val <- 0
    # Loop through all control indices and calculate contributions to numerator and denominator
    for (i in 1:n) {
      if (i != j) {
        e_Xi <- pi_func[i]
        e_Xj <- pi_func[j]
        h_ij <- e_Xi * (1-e_Xi) * e_Xj * (1-e_Xj)
        mu_ij <- mu_func_trt(i, j, cond_prob_trt, cond_prob_ctrl)
        term <- h_ij * mu_ij
        numerator_val <- numerator_val + term
        denominator_val <- denominator_val + h_ij
      }
    }
    return(list(numerator = numerator_val, denominator = denominator_val))
  }, mc.cores = numCores - 1)
  
  total_numerator <- Reduce("+", lapply(B6_values, `[[`, "numerator"))
  denominator <- sum(sapply(B6_values, `[[`, "denominator"))
  B6_vector <- total_numerator / denominator
  
  return(B6_vector)
}




B6_AOW_ctrl <- function(df, pi_func, cond_prob_trt, cond_prob_ctrl) {
  n <- nrow(df)
  B6_values <- mclapply(1:n, function(j) {
    numerator_val <- 0
    denominator_val <- 0
    # Loop through all control indices and calculate contributions to numerator and denominator
    for (i in 1:n) {
      if (i != j) {
        e_Xi <- pi_func[i]
        e_Xj <- pi_func[j]
        h_ij <- e_Xi * (1-e_Xi) * e_Xj * (1-e_Xj)
        mu_ij <- mu_func_ctrl(i, j, cond_prob_trt, cond_prob_ctrl)
        term <- h_ij * mu_ij
        numerator_val <- numerator_val + term
        denominator_val <- denominator_val + h_ij
      }
    }
    return(list(numerator = numerator_val, denominator = denominator_val))
  }, mc.cores = numCores - 1)
  
  total_numerator <- Reduce("+", lapply(B6_values, `[[`, "numerator"))
  denominator <- sum(sapply(B6_values, `[[`, "denominator"))
  B6_vector <- total_numerator / denominator
  
  return(B6_vector)
}




### B7 ### (Dim: 6*1)

B7_AOW_fun <- function(df, pi_func) {
  n <- nrow(df)
  covariate_names <- c("x1", "x2", "x3", "x4", "x5", "x6")
  B7_term1 <- mclapply(1:n, function(j) {
    term1_vec <- numeric(length(covariate_names))
    for (i in 1:n) {
      if (i != j) {
        Xi <- as.matrix(df[i, covariate_names])
        term1_vec <- term1_vec + pi_func[i] * (1-pi_func[i]) * ((1-2*pi_func[i]) * Xi)
      }
    }
    term1_vec <- term1_vec * pi_func[j] * (1 - pi_func[j])
    return(list(term1 = term1_vec))
  }, mc.cores = numCores - 1)
  
  B7_term2 <- mclapply(1:n, function(j) {
    Xj <- as.matrix(df[j, covariate_names])
    term2_vec <- numeric(length(covariate_names))
    for (i in 1:n) {
      if (i != j) {
        term2_vec <- term2_vec + (pi_func[i] * (1-pi_func[i]))
      }
    }
    term2_vec <- term2_vec * pi_func[j] * (1 - pi_func[j]) * (1 - 2*pi_func[j]) * Xj
    return(list(term2 = term2_vec))
  }, mc.cores = numCores - 1)
  
  total_term1 <- Reduce("+", lapply(B7_term1, `[[`, "term1"))
  total_term2 <- Reduce("+", lapply(B7_term2, `[[`, "term2"))
  
  return(total_term1 + total_term2)  
  
  
}


### B8 ### (Dim: 1*1)

B8_AOW_fun <- function(df, pi_func) {
  n <- nrow(df)
  B8_term <- mclapply(1:n, function(j) {
    term_val <- 0
    for (i in 1:n) {
      if (i != j) {
        term_val <- term_val + pi_func[i] * (1-pi_func[i])
      }
    }
    term_val <- term_val * pi_func[j] * (1-pi_func[j])
    return(list(B8_term_val = term_val))
  }, mc.cores = numCores - 1)
  
  B8_AOW <- Reduce("+", lapply(B8_term, `[[`, "B8_term_val"))
  return(B8_AOW)  
}




### Bn ### (Dim: 6*1)

Bn_AOW_trt_fun <- function(){
  Bn_trt <- B1_AOW_trt(df_comb, pi_func, cond_prob_trt, cond_prob_ctrl) -
    B3_AOW_fun(df_comb, pi_func) * 
    (B2_AOW_trt(df_comb, pi_func, cond_prob_trt, cond_prob_ctrl) / 
       B4_AOW_fun(df_comb, pi_func)) + 
    B5_AOW_trt(df_comb, pi_func, cond_prob_trt, cond_prob_ctrl) - 
    B7_AOW_fun(df_comb, pi_func) * 
    (B6_AOW_trt(df_comb, pi_func, cond_prob_trt, cond_prob_ctrl) /
       B8_AOW_fun(df_comb, pi_func))
  return(t(Bn_trt))
}


Bn_AOW_ctrl_fun <- function(){
  Bn_ctrl <- B1_AOW_ctrl(df_comb, pi_func, cond_prob_trt, cond_prob_ctrl) -
    B3_AOW_fun(df_comb, pi_func) * 
    (B2_AOW_ctrl(df_comb, pi_func, cond_prob_trt, cond_prob_ctrl) / 
       B4_AOW_fun(df_comb, pi_func)) + 
    B5_AOW_ctrl(df_comb, pi_func, cond_prob_trt, cond_prob_ctrl) - 
    B7_AOW_fun(df_comb, pi_func) * 
    (B6_AOW_ctrl(df_comb, pi_func, cond_prob_trt, cond_prob_ctrl) /
       B8_AOW_fun(df_comb, pi_func))
  return(t(Bn_ctrl))
}


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






### E_theta_theta ### (Dim: 8*8)

thresholds_trt <- dr_trt$zeta
thresholds_ctrl <- dr_ctrl$zeta

E_thetatheta_fun <- function(df, cond_prob_trt, cond_prob_ctrl, thresholds_trt, thresholds_ctrl) {
  n <- nrow(df)
  F_matrices <- lapply(1:n, function(i) {
    Xi <- t(as.matrix(df[i, c("x1", "x2", "x3", "x4", "x5", "x6")]))
    
    thresholds <- if(df$treatment[i] == 1) thresholds_trt else thresholds_ctrl
    cond_prob <- if(df$treatment[i] == 1) cond_prob_trt else cond_prob_ctrl
    
    # Calculate delta values based on the outcome for individual i
    delta <- as.numeric(df[i, ]$outcomes_comb == 1:3)
    
    # Calculate the c values for individual i
    c1 <- cond_prob[i, 1]
    c2 <- cond_prob[i, 1] + cond_prob[i, 2]
    
    
    # Calculate the F matrix components for individual i
    F11 <- -(delta[1] + delta[2]) * (c1 - c1^2) - delta[2] * exp(thresholds[2] - thresholds[1]) *
      ((exp(thresholds[2] - thresholds[1]) - 1)^(-2))
    F12 <- delta[2] * exp(thresholds[2] - thresholds[1]) * ((exp(thresholds[2] - thresholds[1]) - 1)^(-2))
    F21 <- F12
    F13 <- (delta[1] + delta[2]) * (c1 - c1^2) * t(Xi)
    F31 <- t(F13)
    F22 <- -(delta[2] + delta[3]) * (c2 - c2^2) - delta[2] * exp(thresholds[1] - thresholds[2]) * 
      (1 - exp(thresholds[1] - thresholds[2]))^(-2)
    F23 <- (delta[2] + delta[3]) * (c2 - c2^2) * t(Xi)
    F32 <- t(F23)
    F33 <- -((delta[1] + delta[2]) * (c1 - c1^2) + (delta[2] + delta[3]) * (c2 - c2^2)) * (Xi %*% t(Xi))
    
    # Construct the F matrix for individual i
    top_left <- rbind(cbind(F11, F12), cbind(F21, F22))
    top_2row <- cbind(top_left, rbind(F13, F23))
    bottom_row <- cbind(cbind(F31, F32), F33)
    F_matrix <- as.matrix(rbind(top_2row, bottom_row))
    
    return(F_matrix)
  })
  
  # Compute the average of the F matrices across all individuals
  E_thetatheta <- (-1/n) * Reduce("+", F_matrices)
  return(E_thetatheta)
}



E_thetatheta <- E_thetatheta_fun(df_comb, cond_prob_trt, cond_prob_ctrl, thresholds_trt, thresholds_ctrl)




### mu(Xi, Xj) ### (Dim: 8*1)

mu_tilde_AOW_fun <- function(i, j, df, cond_prob_trt, cond_prob_ctrl) {
  # Extract Xi and Xj as row vectors
  Xi <- t(as.matrix(df[i, c("x1", "x2", "x3", "x4", "x5", "x6")]))
  Xj <- t(as.matrix(df[j, c("x1", "x2", "x3", "x4", "x5", "x6")]))
  
  # Determine which set of conditional probabilities to use based on treatment status
  # cond_prob_i <- if (df$treatment[i] == 1) cond_prob_trt else cond_prob_ctrl
  # cond_prob_j <- if (df$treatment[j] == 1) cond_prob_trt else cond_prob_ctrl
  
  # Calculate the c values for individuals i and j
  ci_1 <- cond_prob_trt[i, 1]
  ci_2 <- cond_prob_trt[i, 1] + cond_prob_trt[i, 2]
  cj_1 <- cond_prob_ctrl[j, 1]
  cj_2 <- cond_prob_ctrl[j, 1] + cond_prob_ctrl[j, 2]
  
  # Construct the mu_tilde vector according to the formula
  mu_tilde_vector_1 <- -ci_1 * (1 - ci_1) * cj_1 + cj_1 * (1 - cj_1) * (ci_2 - ci_1)
  mu_tilde_vector_2 <- -ci_2 * (1 - ci_2) * (cj_2 - cj_1) + (1 - ci_2) * cj_2 * (1 - cj_2)
  mu_tilde_vector_3 <- (cj_1 * (ci_1 - ci_1^2) + (ci_2 - ci_2^2) * (cj_2 - cj_1)) * Xi - 
    ((ci_2 - ci_1) * (cj_1 - cj_1^2) + (1 - ci_2) * (cj_2 - cj_2^2)) * Xj
  mu_tilde_vector <- rbind(mu_tilde_vector_1, mu_tilde_vector_2,
                           mu_tilde_vector_3)
  
  return(mu_tilde_vector)
}



### Cn ### (Dim: 8*1)

Cn_AOW_fun <- function(df, pi_func, cond_prob_trt, cond_prob_ctrl) {
  n <- nrow(df)
  # Calculate the first term of Cn
  Cn_first_term <- mclapply(1:n, function(j) {
    numerator_vec <- 0
    denominator <- 0
    # Loop through all control indices and calculate contributions to numerator and denominator
    for (i in 1:n) {
      if (i != j) {
        e_Xi <- pi_func[i]
        e_Xj <- pi_func[j]
        mu_tilde_ij <- mu_tilde_AOW_fun(i, j, df, cond_prob_trt, cond_prob_ctrl)
        term <- (df$treatment[i] * (1 - df$treatment[j])) * (1-e_Xi) * e_Xj * (-mu_tilde_ij)
        numerator_vec <- numerator_vec + term
        denominator <- denominator + df$treatment[i] * (1-e_Xi)
      }
    }
    denominator <- denominator * (1 - df$treatment[j]) * pi_func[j]
    return(list(numerator = numerator_vec, denominator = denominator))
  }, mc.cores = numCores - 1)
  
  term1_numerator <- Reduce("+", lapply(Cn_first_term, `[[`, "numerator"))
  term1_denominator <- sum(sapply(Cn_first_term, `[[`, "denominator"))
  first_term <- term1_numerator / term1_denominator
  
  
  # Second term of Cn
  Cn_second_term <- mclapply(1:n, function(j) {
    numerator_vec <- 0
    denominator <- 0
    # Loop through all control indices and calculate contributions to numerator and denominator
    for (i in 1:n) {
      if (i != j) {
        e_Xi <- pi_func[i]
        e_Xj <- pi_func[j]
        h_ij <- e_Xi * (1-e_Xi) * e_Xj * (1-e_Xj)
        mu_tilde_ij <- mu_tilde_AOW_fun(i, j, df, cond_prob_trt, cond_prob_ctrl)
        term <- h_ij * (mu_tilde_ij)
        numerator_vec <- numerator_vec + term
        denominator <- denominator + e_Xi * (1-e_Xi)
      }
    }
    denominator <- denominator * pi_func[j] * (1-pi_func[j])
    return(list(numerator = numerator_vec, denominator = denominator))
  }, mc.cores = numCores - 1)
  
  term2_numerator <- Reduce("+", lapply(Cn_second_term, `[[`, "numerator"))
  term2_denominator <- sum(sapply(Cn_second_term, `[[`, "denominator"))
  second_term <- term2_numerator / term2_denominator
  
  Cn_value <- first_term + second_term
  return(Cn_value)
  
}



# Use the function to calculate Cn
Cn_value <- Cn_AOW_fun(df_comb, pi_func, cond_prob_trt, cond_prob_ctrl)




### l_tilde_theta ### (Dim: 8*1)


l_theta_AOW_fun <- function(i, df_comb, cond_prob_trt, cond_prob_ctrl, 
                             thresholds_trt, thresholds_ctrl, E_thetatheta) {
  E_thetatheta_inv <- ginv(E_thetatheta)
  
  Xi <- t(as.matrix(df_comb[i, c("x1", "x2", "x3", "x4", "x5", "x6")]))
  Yi <- df_comb$outcomes[i]
  
  cond_prob <- if (df_comb$treatment[i] == 1) cond_prob_trt else cond_prob_ctrl
  
  c1 <- cond_prob[i, 1]
  c2 <- c1 + cond_prob[i, 2]
  delta <- as.numeric(df_comb[i, ]$outcomes_comb == 1:3)
  thresholds <- if (df_comb$treatment[i] == 1) thresholds_trt else thresholds_ctrl
  
  # Construct the l_tilde vector according to the formula
  l_tilde_vector_1 <- delta[1] - (delta[1] + delta[2]) * c1 - delta[2] * (exp(thresholds[2] - thresholds[1])-1)^(-1)
  l_tilde_vector_2 <- -(delta[2] + delta[3]) * c2 + delta[2] * (1 - exp(thresholds[1] - thresholds[2]))^(-1)
  l_tilde_vector_3 <- (-(delta[1] + delta[2]) + (delta[1] + delta[2]) * c1 + (delta[2] + delta[3]) * c2) * Xi
  l_tilde_vector <- rbind(l_tilde_vector_1, l_tilde_vector_2, l_tilde_vector_3)
  l_tilde_value <- E_thetatheta_inv %*% l_tilde_vector
  
  return(l_tilde_value)
}




AOW_classIII_Variance_trt <- function() {
  n <- nrow(df_comb)
  var_sum <- 0
  Bn_t <- t(Bn_AOW_trt_fun())
  Cn_t <- t(Cn_value)
  for (i in 1:n) {
    l_theta <- l_theta_AOW_fun(i, df_comb, cond_prob_trt, cond_prob_ctrl, 
                               thresholds_trt, thresholds_ctrl, E_thetatheta)
    var_sum <- var_sum + (2 * (g1_AOW_trt[i] - WP_trt_sim_AOW[count_temp]) +
                            Bn_t %*% l_beta_fun(i) +
                            Cn_t %*% l_theta)^2
  }
  
  AOW_Var_trt <- var_sum / n
  return(AOW_Var_trt)
}




AOW_classIII_Variance_ctrl <- function() {
  n <- nrow(df_comb)
  var_sum <- 0
  Bn_t <- t(Bn_AOW_ctrl_fun())
  Cn_t <- t(Cn_value)
  for (i in 1:n) {
    l_theta <- l_theta_AOW_fun(i, df_comb, cond_prob_trt, cond_prob_ctrl, 
                               thresholds_trt, thresholds_ctrl, E_thetatheta)
    var_sum <- var_sum + (2 * (g1_AOW_ctrl[i] - WP_ctrl_sim_AOW[count_temp]) +
                            Bn_t %*% l_beta_fun(i) +
                            Cn_t %*% l_theta)^2
  }
  
  AOW_Var_ctrl <- var_sum / n
  return(AOW_Var_ctrl)
}



### AOW Covariance 
AOW_classIII_cov <- function() {
  n <- nrow(df_comb)
  cov_sum <- 0
  Bn_t_trt <- t(Bn_AOW_trt_fun())
  Bn_t_ctrl <- t(Bn_AOW_ctrl_fun())
  Cn_t <- t(Cn_value)
  for (i in 1:n) {
    l_beta_fun_i = l_beta_fun(i)
    l_theta <- l_theta_AOW_fun(i, df_comb, cond_prob_trt, cond_prob_ctrl, 
                                thresholds_trt, thresholds_ctrl, E_thetatheta)
    temp_trt = (2 * (g1_AOW_trt[i] - WP_trt_sim_AOW[count_temp]) + 
                  Bn_t_trt %*% l_beta_fun_i + 
                  Cn_t %*% l_theta)
    temp_ctrl = (2 * (g1_AOW_ctrl[i] - WP_ctrl_sim_AOW[count_temp]) + 
                   Bn_t_ctrl %*% l_beta_fun_i +
                   Cn_t %*% l_theta)
    cov_sum = cov_sum + temp_trt*temp_ctrl
  }
  AOW_cov_results <- cov_sum / n
  return(AOW_cov_results)
}


Var_trt_AOW_classIII <- AOW_classIII_Variance_trt()/nrow(df_comb)
Var_ctrl_AOW_classIII <- AOW_classIII_Variance_ctrl()/nrow(df_comb)
cov_AOW_calculated <- AOW_classIII_cov()/nrow(df_comb)




#################AOW Estimator Output#########################


###Theoretical Variance###

theory_var_AOW_trt_classIII[count_temp] <- Var_trt_AOW_classIII
theory_var_AOW_ctrl_classIII[count_temp] <- Var_ctrl_AOW_classIII
cov_trt_ctrl[count_temp] <- cov_AOW_calculated
theory_var_AOW_WR[count_temp] <- var_ratio(WP_trt_sim_AOW[count_temp], WP_ctrl_sim_AOW[count_temp],
                                           theory_var_AOW_trt_classIII[count_temp],
                                           theory_var_AOW_ctrl_classIII[count_temp], 
                                           cov_trt_ctrl[count_temp])



}#open on line




theory_sd_AOW_WR <- sqrt(theory_var_AOW_WR)


WR_AOW_coverage[size_count] <- mean((WR_sim_true > (WR_sim_AOW - qnorm(0.975)*theory_sd_AOW_WR)) & 
                                       (WR_sim_true < (WR_sim_AOW + qnorm(0.975)*theory_sd_AOW_WR)))
WR_AOW_var_ratio[size_count] <- mean(theory_var_AOW_WR)/var(WR_sim_AOW)

WR_AOW_theoretical_var_mean[size_count] <- mean(theory_var_AOW_WR)
WR_AOW_empirical_var[size_count] <- var(WR_sim_AOW)

}


end_time <- Sys.time()
run_time <- end_time - start_time
run_time
#7.26111 mins per 2 sim_num



######Coverage Output######




VarCoverage_AOW_df_classIII <- as.data.frame(cbind(WR_AOW_var_ratio,
                                                   WR_AOW_coverage,
                                                   WR_AOW_theoretical_var_mean,
                                                   WR_AOW_empirical_var))

colnames(VarCoverage_AOW_df_classIII) <- c("Var_Ratio_WR",
                                           "Coverage_WR",
                                           "theoretical_var_mean_WR",
                                           "empirical_var_WR")

write.table(VarCoverage_AOW_df_classIII,
            file=paste("results/AOW_VarCoverage_classIII.txt",sep=""), 
            sep="\t", row.names=F)













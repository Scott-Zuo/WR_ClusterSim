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
######### Variance Coverage of Unadjusted ############
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

WP_trt_Ori_coverage <- numeric(re_count)
WP_ctrl_Ori_coverage <- numeric(re_count)
WR_Ori_coverage <- numeric(re_count)

WP_trt_Ori_var_ratio <- numeric(re_count)
WP_ctrl_Ori_var_ratio <- numeric(re_count)
WR_Ori_var_ratio <- numeric(re_count)

WP_trt_Ori_theoretical_var_mean <- numeric(re_count)
WP_trt_Ori_empirical_var <- numeric(re_count)

WP_ctrl_Ori_theoretical_var_mean <- numeric(re_count)
WP_ctrl_Ori_empirical_var <- numeric(re_count)

WR_Ori_theoretical_var_mean <- numeric(re_count)
WR_Ori_empirical_var <- numeric(re_count)

WP_trt_sim_true <- 0.2953701
WP_ctrl_sim_true <- 0.2050298
WR_sim_true <- 1.447066



var_ratio = function(uR,uS,sigR2,sigS2,cov_RS){
  WR_appro_var = (uR^2/uS^2)*(sigR2/uR^2-2*cov_RS/(uR*uS)+sigS2/uS^2)
  return(WR_appro_var)
}




sample_size_list <- c(200,300,400)




for(size_count in 1:re_count){

WP_trt_sim_Ori <- numeric(sim_num)
WP_ctrl_sim_Ori <- numeric(sim_num)
WR_sim_Ori <- numeric(sim_num)

theory_var_Ori_trt <- numeric(sim_num)
theory_var_Ori_ctrl <- numeric(sim_num)
cov_trt_ctrl <- numeric(sim_num)
theory_var_WR <- numeric(sim_num)


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
WP_trt_sim_Ori[count_temp] <- trt_winpr
WP_ctrl_sim_Ori[count_temp] <- ctrl_winpr
WR_sim_Ori[count_temp] <- trt_winpr / ctrl_winpr

colnames(df_trt) <- "outcomes_comb"
colnames(df_ctrl) <- "outcomes_comb"

df_comb <- rbind(df_trt, df_ctrl)
df_comb$treatment <- c(rep(1,nrow(df_trt)), rep(0,nrow(df_ctrl)))









###########Unadjusted Theoretical Variance###########

n1 <- nrow(df_trt)
n2 <- nrow(df_ctrl)

sum_treatment <- sum(df_comb[,"treatment"])
sum_control <- sum(1 - df_comb[,"treatment"])

# Function to parallelize for the treatment group ijik
p_ijik_trt_parallel <- function(i) {

  temp_inner <- numeric(n2)
  
  for (j in 1:n2) {
    wij_unadjusted <- 1
    temp_inner[j] <- wij_unadjusted * as.numeric(df_trt[i,"outcomes_comb"] > df_ctrl[j,"outcomes_comb"])
  }
  
  return((sum(temp_inner))^2)
}

# Parallel computation of the term1 for the treatment group
p_ijik_term1_trt_list <- mclapply(1:n1, p_ijik_trt_parallel, mc.cores = numCores - 1)
p_ijik_term1_trt <- sum(unlist(p_ijik_term1_trt_list))/(n1 * n2 * (n2 - 1))
p_ijik_term2_trt <- (1/(n2-1)) * WP_trt_sim_Ori[count_temp]

p_ijik_trt_results <- p_ijik_term1_trt - p_ijik_term2_trt



# Function to parallelize for the treatment group ijkj
p_ijkj_trt_parallel <- function(j) {

  temp_inner <- numeric(n1)
  
  for (i in 1:n1) {
    wij_unadjusted <- 1
    temp_inner[i] <- wij_unadjusted * as.numeric(df_trt[i,"outcomes_comb"] > df_ctrl[j,"outcomes_comb"])
  }
  
  return((sum(temp_inner))^2)
}



# Parallel computation of the term1 for the treatment group
p_ijkj_term1_trt_list <- mclapply(1:n2, p_ijkj_trt_parallel, mc.cores = numCores - 1)
p_ijkj_term1_trt <- sum(unlist(p_ijkj_term1_trt_list))/(n2 * n1 * (n1 - 1))
p_ijkj_term2_trt <- (1/(n1-1)) * WP_trt_sim_Ori[count_temp]

p_ijkj_trt_results <- p_ijkj_term1_trt - p_ijkj_term2_trt




# Function to parallelize for the control group ijik
p_ijik_ctrl_parallel <- function(i) {

  temp_inner <- numeric(n2)
  
  for (j in 1:n2) {
    wij_unadjusted <- 1
    temp_inner[j] <- wij_unadjusted * as.numeric(df_trt[i,"outcomes_comb"] < df_ctrl[j,"outcomes_comb"])
  }
  
  return((sum(temp_inner))^2)
}


# Parallel computation of the term1 for the control group
p_ijik_term1_ctrl_list <- mclapply(1:n1, p_ijik_ctrl_parallel, mc.cores = numCores - 1)
p_ijik_term1_ctrl <- sum(unlist(p_ijik_term1_ctrl_list))/(n1 * n2 * (n2 - 1))
p_ijik_term2_ctrl <- (1/(n2-1)) * WP_ctrl_sim_Ori[count_temp]


p_ijik_ctrl_results <- p_ijik_term1_ctrl - p_ijik_term2_ctrl




# Function to parallelize for the control group ijkj
p_ijkj_ctrl_parallel <- function(j) {

  temp_inner <- numeric(n1)
  
  for (i in 1:n1) {
    wij_unadjusted <- 1
    temp_inner[i] <- wij_unadjusted * as.numeric(df_trt[i,"outcomes_comb"] < df_ctrl[j,"outcomes_comb"])
  }
  
  return((sum(temp_inner))^2)
}


# Parallel computation of the term1 for the control group
p_ijkj_term1_ctrl_list <- mclapply(1:n2, p_ijkj_ctrl_parallel, mc.cores = numCores - 1)
p_ijkj_term1_ctrl <- sum(unlist(p_ijkj_term1_ctrl_list))/(n2 * n1 * (n1 - 1))
p_ijkj_term2_ctrl <- (1/(n1-1)) * WP_ctrl_sim_Ori[count_temp]


p_ijkj_ctrl_results <- p_ijkj_term1_ctrl - p_ijkj_term2_ctrl




Var_trt_Ori_calculated <- (1/(n1+n2)) * (((n1+n2)/n1) * p_ijik_trt_results -
                                         ((n1+n2)/n1) * (WP_trt_sim_Ori[count_temp])^2 +
                                         ((n1+n2)/n2) * p_ijkj_trt_results -
                                         ((n1+n2)/n2) * (WP_trt_sim_Ori[count_temp])^2)

Var_ctrl_Ori_calculated <- (1/(n1+n2)) * (((n1+n2)/n1) * p_ijik_ctrl_results -
                                           ((n1+n2)/n1) * (WP_ctrl_sim_Ori[count_temp])^2 +
                                           ((n1+n2)/n2) * p_ijkj_ctrl_results -
                                           ((n1+n2)/n2) * (WP_ctrl_sim_Ori[count_temp])^2)






########### Computing Covariance ###########

p_ijik_cov <- function(i) {
  temp_inner1 <- temp_inner2 <- numeric(n2)
  for (j in 1:n2) {
    wij_unadjusted <- 1
    temp_inner1[j] <- wij_unadjusted * as.numeric(df_trt[i,"outcomes_comb"] > df_ctrl[j,"outcomes_comb"])
    temp_inner2[j] <- wij_unadjusted * as.numeric(df_trt[i,"outcomes_comb"] < df_ctrl[j,"outcomes_comb"])
  }
  return((sum(temp_inner1))*(sum(temp_inner2)))
}
term1_cov_list1 <- apply(matrix(1:(n1),nrow=1),2, p_ijik_cov)
p_ijik_cov_results <- sum(unlist(term1_cov_list1))/(n1 * n2 * (n2 - 1))

p_ijkj_cov <- function(j) {
  temp_inner1 <- temp_inner2 <- numeric(n1)
  for (i in 1:n1) {
    wij_unadjusted <- 1
    temp_inner1[i] <- wij_unadjusted * as.numeric(df_trt[i,"outcomes_comb"] > df_ctrl[j,"outcomes_comb"])
    temp_inner2[i] <- wij_unadjusted * as.numeric(df_trt[i,"outcomes_comb"] < df_ctrl[j,"outcomes_comb"])
  }
  return((sum(temp_inner1))*(sum(temp_inner2)))
}
term2_cov_list1 <- apply(matrix(1:(n2),nrow=1),2, p_ijkj_cov)
p_ijkj_cov_results <- sum(unlist(term2_cov_list1))/(n2 * n1 * (n1 - 1))
cov_Ori_calculated <- (1/(n1+n2)) * (((n1+n2)/n1) * p_ijik_cov_results -
                                       ((n1+n2)/n1) * (WP_trt_sim_Ori[count_temp])*(WP_ctrl_sim_Ori[count_temp]) +
                                       ((n1+n2)/n2) * p_ijkj_cov_results -
                                       ((n1+n2)/n2) * (WP_trt_sim_Ori[count_temp])*(WP_ctrl_sim_Ori[count_temp]))

theory_var_Ori_trt[count_temp] <- Var_trt_Ori_calculated
theory_var_Ori_ctrl[count_temp] <- Var_ctrl_Ori_calculated
cov_trt_ctrl[count_temp] <- cov_Ori_calculated
theory_var_WR[count_temp] <- var_ratio(WP_trt_sim_Ori[count_temp],WP_ctrl_sim_Ori[count_temp],
                                       theory_var_Ori_trt[count_temp],theory_var_Ori_ctrl[count_temp],
                                       cov_trt_ctrl[count_temp])




}#open on line
  
  
  

theory_sd_WR <- sqrt(theory_var_WR)


WR_Ori_coverage[size_count] <- mean((WR_sim_true > (WR_sim_Ori - qnorm(0.975)*theory_sd_WR)) & 
                                      (WR_sim_true < (WR_sim_Ori + qnorm(0.975)*theory_sd_WR)))

WR_Ori_var_ratio[size_count] <- mean(theory_var_WR)/var(WR_sim_Ori)
WR_Ori_theoretical_var_mean[size_count] <- mean(theory_var_WR)
WR_Ori_empirical_var[size_count] <- var(WR_sim_Ori)

}


end_time <- Sys.time()
run_time <- end_time - start_time
run_time
#32.42666 secs per 5 sim_num













######Coverage Output######




VarCoverage_Ori_df <- as.data.frame(cbind(WR_Ori_var_ratio,
                                          WR_Ori_coverage,
                                          WR_Ori_theoretical_var_mean,
                                          WR_Ori_empirical_var))

colnames(VarCoverage_Ori_df) <- c("Var_Ratio_WR",
                                  "Coverage_WR",
                                  "theoretical_var_mean_WR",
                                  "empirical_var_WR")

write.table(VarCoverage_Ori_df,
            file=paste("results/Unadjusted_VarCoverage.txt",sep=""), 
            sep="\t", row.names=F)


















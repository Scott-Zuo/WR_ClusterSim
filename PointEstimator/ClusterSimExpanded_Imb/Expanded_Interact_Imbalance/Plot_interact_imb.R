library(ggplot2)


######################################################
#####Plot with the output from cluster simulation#####
######################################################

###Treatment/Control: 70/30###

step_size <- 60
re_count <- 10
#setwd("/Users/sizhezuo/Desktop/YCASResearch/WR_Cluster/Expanded_Interact_Imbalance/results/results_0226")

Ori_df <- read.table("Unadjusted_newscenario1.txt", header=TRUE, sep="")
IPW_df <- read.table("IPW_newscenario1.txt", header=TRUE, sep="")
OW_df <- read.table("OW_newscenario1.txt", header=TRUE, sep="")
DR_df <- read.table("DR_newscenario1.txt", header=TRUE, sep="")
DROW_df <- read.table("DROW_newscenario1.txt", header=TRUE, sep="")

colnames(Ori_df) <- c("variance", "bias", "mse")
colnames(IPW_df) <- c("variance", "bias", "mse")
colnames(OW_df) <- c("variance", "bias", "mse")
colnames(DR_df) <- c("variance", "bias", "mse")
colnames(DROW_df) <- c("variance", "bias", "mse")


sim_df_WD <- rbind(Ori_df[1:re_count,], 
                   IPW_df[1:re_count,], 
                   OW_df[1:re_count,],
                   DR_df[1:re_count,],
                   DROW_df[1:re_count,])

sim_df_WR <- rbind(Ori_df[(re_count+1):(2*re_count),], 
                   IPW_df[(re_count+1):(2*re_count),], 
                   OW_df[(re_count+1):(2*re_count),],
                   DR_df[(re_count+1):(2*re_count),],
                   DROW_df[(re_count+1):(2*re_count),])


##############Win Difference#########

sim_df_WD$estimator <- c(rep("unadjusted",re_count), rep("IPW",re_count),
                         rep("OW",re_count), rep("AIPW",re_count),
                         rep("AOW",re_count))
sim_df_WD$estimator <- as.factor(sim_df_WD$estimator)
sim_df_WD$samplesize <- c(c(1:re_count)*step_size, c(1:re_count)*step_size, 
                          c(1:re_count)*step_size, c(1:re_count)*step_size,
                          c(1:re_count)*step_size)

sim_df_WD <- sim_df_WD[is.finite(rowSums(cbind(sim_df_WD$variance, 
                                               sim_df_WD$bias,
                                               sim_df_WD$mse))), ]


#################Win Ratio###########
sim_df_WR$estimator <- c(rep("unadjusted",re_count), rep("IPW",re_count),
                         rep("OW",re_count), rep("AIPW",re_count),
                         rep("AOW",re_count))
sim_df_WR$estimator <- as.factor(sim_df_WR$estimator)
sim_df_WR$samplesize <- c(c(1:re_count)*step_size, c(1:re_count)*step_size, 
                          c(1:re_count)*step_size, c(1:re_count)*step_size,
                          c(1:re_count)*step_size)

sim_df_WR <- sim_df_WR[is.finite(rowSums(cbind(sim_df_WR$variance, 
                                               sim_df_WR$bias,
                                               sim_df_WR$mse))), ]








############################## Variance Plot ################################




###WD###
var_plot_WD <- ggplot(sim_df_WD, aes(x=samplesize, y = variance,
                                     color = estimator)) + 
  geom_point() + 
  geom_line(aes(group = estimator)) + 
  scale_color_manual(values = c("unadjusted" = "black", "IPW" = "gold",
                                "OW" = "orange", "AIPW" = "red", "AOW" = "brown")) +
  ggtitle("Variance comparison from Simulation of Win Difference") +
  xlab("Total Sample Size") + ylab("Variance") +
  scale_x_continuous(breaks=seq(0,step_size*re_count,step_size)) + 
  theme(
    plot.title = element_text(color="red", size=12, face="bold.italic"),
    axis.title.x = element_text(size=12),
    axis.text.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    legend.text = element_text(size=12)
  )



###WR###
var_plot_WR <- ggplot(sim_df_WR, aes(x=samplesize, y = variance,
                                     color = estimator)) + 
  geom_point() + 
  geom_line(aes(group = estimator)) + 
  scale_color_manual(values = c("unadjusted" = "black", "IPW" = "gold",
                                "OW" = "orange", "AIPW" = "red", "AOW" = "brown")) +
  ggtitle("Variance comparison from Simulation of Win Ratio") +
  xlab("Total Sample Size") + ylab("Variance") +
  scale_x_continuous(breaks=seq(0,step_size*re_count,step_size)) + 
  theme(
    plot.title = element_text(color="red", size=12, face="bold.italic"),
    axis.title.x = element_text(size=12),
    axis.text.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    legend.text = element_text(size=12)
  )



png(filename = "fig_interactImb_Var_WD.png", units = "in",
    width = 6, height = 4, res = 300)
var_plot_WD
dev.off()

png(filename = "fig_interactImb_Var_WR.png", units = "in",
    width = 6, height = 4, res = 300)
var_plot_WR
dev.off()






############################### MSE Plot #################################




mse_plot_WD <- ggplot(sim_df_WD, aes(x=samplesize, y = mse,
                                     color = estimator)) + 
  geom_point() + 
  geom_line(aes(group = estimator)) + 
  scale_color_manual(values = c("unadjusted" = "black", "IPW" = "gold",
                                "OW" = "orange", "AIPW" = "red", "AOW" = "brown"))+
  ggtitle("MSE comparison from Simulation for Win Difference") +
  xlab("Total Sample Size") + ylab("MSE") +
  scale_x_continuous(breaks=seq(0,step_size*re_count,step_size)) + 
  theme(
    plot.title = element_text(color="red", size=12, face="bold.italic"),
    axis.title.x = element_text(size=12),
    axis.text.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    legend.text = element_text(size=12)
  )




mse_plot_WR <- ggplot(sim_df_WR, aes(x=samplesize, y = mse,
                                     color = estimator)) + 
  geom_point() + 
  geom_line(aes(group = estimator)) + 
  scale_color_manual(values = c("unadjusted" = "black", "IPW" = "gold",
                                "OW" = "orange", "AIPW" = "red", "AOW" = "brown"))+
  ggtitle("MSE comparison from Simulation for Win Ratio") +
  xlab("Total Sample Size") + ylab("MSE") +
  scale_x_continuous(breaks=seq(0,step_size*re_count,step_size)) + 
  theme(
    plot.title = element_text(color="red", size=12, face="bold.italic"),
    axis.title.x = element_text(size=12),
    axis.text.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    legend.text = element_text(size=12)
  )




png(filename = "fig_interactImb_MSE_WD.png", units = "in",
    width = 6, height = 4, res = 300)
mse_plot_WD
dev.off()

png(filename = "fig_interactImb_MSE_WR.png", units = "in",
    width = 6, height = 4, res = 300)
mse_plot_WR
dev.off()






############################### Bias Plot #################################




bias_plot_WD <- ggplot(sim_df_WD, aes(x = samplesize, y = abs(bias),
                                      color = estimator)) + 
  geom_point() + 
  geom_line(aes(group = estimator)) +
  scale_color_manual(values = c("unadjusted" = "black", "IPW" = "gold",
                                "OW" = "orange", "AIPW" = "red", "AOW" = "brown")) +
  ggtitle("Bias from Simulation for Win Difference") +
  xlab("Total Sample Size") + ylab("Percentage of Bias") +
  scale_x_continuous(breaks=seq(0,step_size*re_count,step_size)) +
  theme(
    plot.title = element_text(color="red", size=12, face="bold.italic"),
    axis.title.x = element_text(size=12),
    axis.text.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    legend.text = element_text(size=12)
  )




bias_plot_WR <- ggplot(sim_df_WR, aes(x = samplesize, y = abs(bias),
                                      color = estimator)) + 
  geom_point() + 
  geom_line(aes(group = estimator)) +
  scale_color_manual(values = c("unadjusted" = "black", "IPW" = "gold",
                                "OW" = "orange", "AIPW" = "red", "AOW" = "brown")) +
  ggtitle("Bias from Simulation for Win Ratio") +
  xlab("Total Sample Size") + ylab("Percentage of Bias") +
  scale_x_continuous(breaks=seq(0,step_size*re_count,step_size)) +
  theme(
    plot.title = element_text(color="red", size=12, face="bold.italic"),
    axis.title.x = element_text(size=12),
    axis.text.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    legend.text = element_text(size=12)
  )




png(filename = "fig_interactImb_Bias_WD.png", units = "in",
    width = 6, height = 4, res = 300)
bias_plot_WD
dev.off()

png(filename = "fig_interactImb_Bias_WR.png", units = "in",
    width = 6, height = 4, res = 300)
bias_plot_WR
dev.off()






###################### Relative Efficiency Plot ########################


re_df_WD <- read.table("sim_re_WD_newscenario1.txt", header=TRUE, sep="")
re_df_WR <- read.table("sim_re_WR_newscenario1.txt", header=TRUE, sep="")





re_df_WD <- data.frame(c(re_df_WD$relative_eff_IPW_WD, re_df_WD$relative_eff_OW_WD,
                         re_df_WD$relative_eff_DR_WD, re_df_WD$relative_eff_DROW_WD))
colnames(re_df_WD) <- "relative_eff"
re_df_WD$samplesize <- c(c(1:re_count)*step_size, c(1:re_count)*step_size,
                         c(1:re_count)*step_size, c(1:re_count)*step_size)
re_df_WD$estimator <- c(rep("IPW",re_count), rep("OW",re_count),
                        rep("AIPW",re_count), rep("AOW",re_count))
re_df_WD$estimator <- as.factor(re_df_WD$estimator)


re_df_WR <- data.frame(c(re_df_WR$relative_eff_IPW_WR, re_df_WR$relative_eff_OW_WR,
                         re_df_WR$relative_eff_DR_WR, re_df_WR$relative_eff_DROW_WR))
colnames(re_df_WR) <- "relative_eff"
re_df_WR$samplesize <- c(c(1:re_count)*step_size, c(1:re_count)*step_size,
                         c(1:re_count)*step_size, c(1:re_count)*step_size)
re_df_WR$estimator <- c(rep("IPW",re_count), rep("OW",re_count),
                        rep("AIPW",re_count), rep("AOW",re_count))
re_df_WR$estimator <- as.factor(re_df_WR$estimator)







re_plot_WD <- ggplot(re_df_WD, aes(x=samplesize, y = relative_eff,
                                   color = estimator)) + 
  geom_point() + 
  geom_line(aes(group = estimator)) +
  scale_color_manual(values = c("IPW" = "gold", "OW" = "orange",
                                "AIPW" = "red", "AOW" = "brown")) +
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  ggtitle("Relative Efficiency from Simulation for Win Difference") +
  xlab("Total Sample Size") + ylab("Relative Efficiency") +
  scale_x_continuous(breaks=seq(0,step_size*re_count,step_size)) +
  coord_cartesian(ylim = c(0.8, 1.8)) +
  theme(
    plot.title = element_text(color="red", size=12, face="bold.italic"),
    axis.title.x = element_text(size=12),
    axis.text.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    legend.text = element_text(size=12)
  )


re_plot_WR <- ggplot(re_df_WR, aes(x=samplesize, y = relative_eff,
                                   color = estimator)) + 
  geom_point() + 
  geom_line(aes(group = estimator)) +
  scale_color_manual(values = c("IPW" = "gold", "OW" = "orange",
                                "AIPW" = "red", "AOW" = "brown")) +
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  ggtitle("Relative Efficiency from Simulation for Win Ratio") +
  xlab("Total Sample Size") + ylab("Relative Efficiency") +
  scale_x_continuous(breaks=seq(0,step_size*re_count,step_size)) +
  coord_cartesian(ylim = c(0.8, 1.8)) +
  theme(
    plot.title = element_text(color="red", size=12, face="bold.italic"),
    axis.title.x = element_text(size=12),
    axis.text.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    legend.text = element_text(size=12)
  )


png(filename = "fig_interactImb_RE_WD.png", units = "in",
    width = 6, height = 4, res = 300)
re_plot_WD
dev.off()

png(filename = "fig_interactImb_RE_WR.png", units = "in",
    width = 6, height = 4, res = 300)
re_plot_WR
dev.off()





ggsave("fig_interactImb_RE_WD.eps", plot = re_plot_WD, device = "eps", width = 6, height = 4)
ggsave("fig_interactImb_RE_WR.eps", plot = re_plot_WR, device = "eps", width = 6, height = 4)















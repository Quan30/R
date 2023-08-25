

### Load the estimated parameters as data tables
library(data.table)
# Set the working directory
mydir <- "C:/Users/JiaqiLu/Documents/PhD@Helmholtz/proj0.5_ProbGen/code/parameter_estimations"
setwd(mydir)
# List all .csv files
csv_files <- list.files(pattern = "*.csv")
# Read each .csv file into a list of data frames using fread
list_of_dataframes <- lapply(csv_files, function(file) fread(file, header = TRUE))
# Remove the .csv extension from the file names to use as dataframe names
names(list_of_dataframes) <- sub("\\.csv$", "", csv_files)
# Assign each dataframe to the global environment with its name
list2env(list_of_dataframes, envir = .GlobalEnv)
# n_genes & n_perturbations
n_pert <- dim(perturb_disp_lfc_mu)[1]
n_gene <- dim(perturb_disp_lfc_mu)[2]
print(paste("We observe", n_gene, "genes, where", n_pert, "perturbations are made."))



########################################
# plot log_var_mean, size_factor_mean, batch_effect, covariate_effects, nb_log_mean_ctrl, nb_log_disp_ctrl

par(mfrow=c(3,2))
# blue
barplot(log_var_mean, las=2, col=rgb(0.25, 0.41, 0.88), border=FALSE, axisnames=FALSE,
        main="log_var_mean", xlab = "Gene", ylab = "Value")
abline(h=mean(log_var_mean), lty=2, col="black")
mtext(side=2, at=mean(log_var_mean), text=sprintf("%.2f", mean(log_var_mean)), las=2, line=0.2, cex=0.8, col="black")
# blue
barplot(size_factor_mean, las=2, col=rgb(0.25, 0.41, 0.88), border=FALSE, axisnames=FALSE,
        main="size_factor_mean", xlab = "Gene", ylab = "Value")
abline(h=size_factor_mean, lty=2, col="black")
mtext(side=2, at=size_factor_mean, text=sprintf("%.2f", size_factor_mean), las=2, line=0.2, cex=0.8, col="black")
# blue
barplot(batch_effect, las=2, col=rgb(0.25, 0.41, 0.88), border=FALSE, axisnames=FALSE,
        main="batch_effect", xlab = "Gene", ylab = "Value")
abline(h=mean(batch_effect), lty=2, col="black")
mtext(side=2, at=mean(batch_effect), text=sprintf("%.2f", mean(batch_effect)), las=2, line=0.2, cex=0.8, col="black")
# blue
barplot(covariate_effects, las=2, col=rgb(0.25, 0.41, 0.88), border=FALSE, axisnames=FALSE,
        main="covariate_effects", xlab = "Gene", ylab = "Value")
abline(h=mean(covariate_effects), lty=2, col="black")
mtext(side=2, at=mean(covariate_effects), text=sprintf("%.2f", mean(covariate_effects)), las=2, line=0.2, cex=0.8, col="black")
# blue
barplot(nb_log_mean_ctrl, las=2, col=rgb(0.25, 0.41, 0.88), border=FALSE, axisnames=FALSE,
        main="nb_log_mean_ctrl", xlab = "Gene", ylab = "Value")
abline(h=mean(nb_log_mean_ctrl), lty=2, col="black")
mtext(side=2, at=mean(nb_log_mean_ctrl), text=sprintf("%.2f", mean(nb_log_mean_ctrl)), las=2, line=0.2, cex=0.8, col="black")
# green
barplot(nb_log_disp_ctrl, las=2, col=rgb(0.13, 0.55, 0.13), border=FALSE, axisnames=FALSE,
        main="nb_log_disp_ctrl", xlab = "Gene", ylab = "Value")
abline(h=mean(nb_log_disp_ctrl), lty=2, col="black")
mtext(side=2, at=mean(nb_log_disp_ctrl), text=sprintf("%.2f", mean(nb_log_disp_ctrl)), las=2, line=0.2, cex=0.8, col="black")

# Reset plotting parameters to default
par(mfrow=c(1,1))





#####################
# Visualize Posterior Distribution of each Parameter
## - Draw the approximated posterior density (refer to guide about what distirbution is assumed for the posterior)
## - Mark the 95% credible interval (all posterior expect log_pooling&multiplicative_noise are assumed to be normal, all are symmetric)
## - Mark the mean of the posterior & the MAP estimation


## How different is the MAP estimation from the mean of the posterior distribution
# Correlation matrix
corr_perturb_mean_lfc <- cor(perturb_disp_lfc_mu, perturb_disp_lfc_mu_map)
View(corr_perturb_mean_lfc)
# take perturb_mean_lfc as example
perturb_mean_lfc_mu1 <- as.numeric(perturb_mean_lfc_mu[1,])
perturb_mean_lfc_mu_map1 <- as.numeric(perturb_mean_lfc_mu_map[1,])
plot(perturb_mean_lfc_mu1, perturb_mean_lfc_mu_map1, 
     xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2),
     xlab = "perturb_mean_lfc_mu", ylab = "perturb_mean_lfc_mu_map", main = "Compare Perturb Mean LFC")
# Add a line for better visualization
abline(a=0, b=1, col="red", lty=2)
# RESULT: returns an almost horizontal line instead of a line with slope 1. soooooo weird!

# let us try again with cont_cov_effect
cont_cov_effect_mu1 <- as.numeric(cont_cov_effect_mu[1,])
cont_cov_effect_mu_map1 <- as.numeric(cont_cov_effect_mu_map[1,])
plot(cont_cov_effect_mu1, cont_cov_effect_mu_map1, 
     xlim = c(-1, 1), ylim = c(-1, 1),
     xlab = "perturb_mean_lfc_mu", ylab = "perturb_mean_lfc_mu_map", main = "Compare Perturb Mean LFC")
# Add a line for better visualization
abline(a=0, b=1, col="red", lty=2)
# RESULT: looks much better than perturb_mean_lfc, even though still not on the slope-1-line, at least already similar
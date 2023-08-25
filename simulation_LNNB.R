## Our plan in general
# 1. based on the MAP that we estimated, sample from the corresponding LNNB distribution
# 2. adjust n_cells, compare the power within each methode (MKmisc, scPower, powersimR, muscat)

################
## sample from LNNB

## Toy example: sample 100 data from LNNB
n <- 1000        # number of samples
norm_mean <- 0
norm_sd <- 1
r <- 5            # total_count: r=5
p <- 0.2          # set 1-p as the porbability
l <- log(p/(1-p)) # logit argument of the LNNB.


# nb samples
samples_nb <- rnbinom(n = n, size = r, prob = 1-p)

nb_mean_theo <- r * exp(l)
nb_var_theo <- nb_mean_theo + (nb_mean_theo)^2 / r
nb_mean_emp <- mean(samples_nb)
nb_var_emp <- var(samples_nb)

# print results
print(paste("The theoretical mean is", round(nb_mean_theo,4), ". The theoretical variance is", round(nb_var_theo,4), "."))
print(paste("The empirical mean is", round(nb_mean_emp,4), ". The empirical variance is", round(nb_var_emp,4), "."))


# Sample n observations from a normal distribution with mean 0 and standard deviation 1
samples_normal <- rnorm(n = n, mean = norm_mean, sd = norm_sd)

# Sample n observations from a negative binomial distribution with size 5 and probability 0.5
samples_lnnb <- c()
for (e in samples_normal) {            # here e = epsilon is the normal distributed 'noise'
  l1 <- l + e
  p1 <- 1 / (1 + exp(-l1))
  #p1 <- p*exp(e) / (1 - p + p*exp(e))  # the formula compute from l1 = l + epsilon
  new_sample <- rnbinom(n = 1, size = r, prob = 1-p1)
  samples_lnnb <- c(samples_lnnb, new_sample)
}


# test whether it is the distribution we expected by computing empirical with theoretical mean&var 
# theoretical mean & variance
lnnb_mean_theo <- exp(l + log(r) + 1/2 * norm_sd^2)
lnnb_var_theo <- lnnb_mean_theo + (exp(norm_sd^2) * (1+1/r) - 1) * (lnnb_mean_theo^2)
# empirical mean & variance
lnnb_mean_emp <- mean(samples_lnnb)
lnnb_var_emp <- var(samples_lnnb)

# print results
print(paste("The theoretical mean is", round(lnnb_mean_theo,4), ". The theoretical variance is", round(lnnb_var_theo,4), "."))
print(paste("The empirical mean is", round(lnnb_mean_emp,4), ". The empirical variance is", round(lnnb_var_emp,4), "."))


## write it as a function as in pyro
rlnnb <- function(n, logits, total_count, multiplicative_noise_scale){
  # total_count = r
  # logits = l
  # multiplicative_noise_scale = sigma = norm_sd
  
  # generate normal distributed ultiplicative noise
  samples_normal <- rnorm(n = n, mean = 0, sd = multiplicative_noise_scale)
  
  logits1 <- logits + samples_normal             
  p1 <- 1 / (1 + exp(-logits1))
  samples_lnnb <- rnbinom(n = n, size = total_count, prob = 1-p1)
  
  return(samples_lnnb)
}

# test 
lnnb_samples1 <- rlnnb(n = 1000, total_count = 5, logits = log(0.2/0.8), multiplicative_noise_scale = 1)
mean(lnnb_samples1)




#########################
## sample from LNNB with parameter set to our MAP estimations
# Goal: input perturbation indicator | number_of_cells to simulate | batch_numer(or the mean?)


### Load the estimated parameters as data tables
library(data.table)
# Set the working directory
mydir <- "~/work/proj0.5_ProbGen/code/parameter_estimations"
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


### set all parameters

## parameters for nb_log_mean
log_var_mean <- mapply(function(mu, sigma) rnorm(1, mu, sigma), array(log_var_mean_mu)[[1]], array(log_var_mean_sigma)[[1]])  # array, length = 5859, estimated | (input default)

size_factor_mean <- mean(array(size_factor[[1]]))  # single number, length = 1, mean of given size factor over all cells | (input default)

batch_effect <- mapply(function(mu, sigma) rnorm(1, mu, sigma), rapply(batch_effect_mu, mean), rapply(batch_effect_sigma, mean)) # array, length = 5859, mean over estimated batch_effect of 6 batches | (input default)

cont_cov_effect <- mapply(function(mu, sigma) rnorm(1, mu, sigma), array(unlist(cont_cov_effect_mu)), array(unlist(cont_cov_effect_sigma)))
cont_cov_effect <- as.data.frame(matrix(cont_cov_effect, ncol = n_gene))  # dataframe, dim = 3*5859, estimated
covariate_effects <- colMeans(as.matrix(cont_covariates) %*% as.matrix(cont_cov_effect))  # array, length = 5859, mean of covariate effects over all cells | (input default)

nb_log_mean_ctrl <- log_var_mean + size_factor_mean + batch_effect + covariate_effects  # array, length = 5859

element_mean_lfc <- mapply(function(mu, sigma) rnorm(1, mu, sigma), array(unlist(element_mean_lfc_mu)), array(unlist(element_mean_lfc_sigma)))
element_mean_lfc <- as.data.frame(matrix(element_mean_lfc, ncol = n_gene))  # dataframe, dim = 1243*5859, estimated elementwise lfc

interm_mean <- perturb_mean_lfc_mu + as.matrix(guide_by_element) %*% as.matrix(element_mean_lfc)  # a intermediate mean of siye 1243 * 5859
perturb_mean_lfc <- mapply(function(mu, sigma) rnorm(1, mu, sigma), array(unlist(interm_mean)), array(unlist(perturb_mean_lfc_sigma)))
perturb_mean_lfc <- as.data.frame(matrix((perturb_mean_lfc), ncol = n_gene))   # dataframe, dim = 1243*5859

perturbations <- rep(0,n_pert)  # array, length = 1243, define by myself as input | (input default)
#perturbations <- rep(1,n_pert) 

nb_log_mean <- nb_log_mean_ctrl + t(as.matrix(perturbations)) %*% as.matrix(perturb_mean_lfc)  # matrix, dim = 1*5859, an input to the LNNB

## parameters for nb_log_disp
log_var_disp <- mapply(function(mu, sigma) rnorm(1, mu, sigma), array(log_var_disp_mu)[[1]], array(log_var_disp_sigma)[[1]])  # array, length = 5859, estimated | (input default)
#nb_log_disp_ctrl <- matrix(rep(log_var_disp, n_pert), nrow=n_pert)
nb_log_disp_ctrl <- log_var_disp   # array, length = 5859

element_disp_lfc <- mapply(function(mu, sigma) rnorm(1, mu, sigma), array(unlist(element_disp_lfc_mu)), array(unlist(element_disp_lfc_sigma)))
element_disp_lfc <- as.data.frame(matrix(element_disp_lfc, ncol = n_gene))  # dataframe, dim = 1243*5859, estimated elementwise lfc

interm_mean_disp <- perturb_disp_lfc_mu + as.matrix(guide_by_element) %*% as.matrix(element_disp_lfc)
perturb_disp_lfc <- mapply(function(mu, sigma) rnorm(1, mu, sigma), array(unlist(interm_mean_disp)), array(unlist(perturb_disp_lfc_sigma)))
perturb_disp_lfc <- as.data.frame(matrix((perturb_disp_lfc), ncol = n_gene))   # dataframe, dim = 1243*5859

nb_log_dispersion <- nb_log_disp_ctrl + t(as.matrix(perturbations)) %*% as.matrix(perturb_disp_lfc)  # matrix, dim = 1*5859, an input to the LNNB

## multiplicative noise
multiplicative_noise <- as.matrix(multiplicative_noise_mu)  # data.table, dim = 5859*1, the estimated sigma
multiplicative_noise <- matrix(multiplicative_noise, ncol = n_gene)   # matrix, dim = 1*5859, reshaped


### simulate from LNNB

## sample from multivariate LNNB (refer to 1-dim in the toy example)

# Number of samples to simulate
n_samples <- 1000

rmvlnnb <- function(n = 100, 
                    logits = nb_log_mean - nb_log_dispersion - multiplicative_noise^2/2, 
                    total_count = exp(nb_log_dispersion), 
                    multiplicative_noise_scale = multiplicative_noise){
  # inputs: need to be a matrix, need dim = 1*n_vars.
  # output: a n*n_vars matrix. n is the number of samples to simulate. n_vars is the dim of variates
  
  n_vars <- dim(logits)[2]  # determine dimension of variates
  samples <- matrix(rep(0, n*n_vars), nrow = n)  # create a n*n_vars zero matrix to save the sampling results
  
  samples <- mapply(function(logits, total_count, multiplicative_noise_scale) rlnnb(n=n, logits, total_count, multiplicative_noise_scale), 
                    logits[1,], total_count[1,], multiplicative_noise_scale[1,])
  
  return(samples)
}

# test
# set params as our estimation, remember that here we set all perturbation to be zero
logits <- nb_log_mean - nb_log_dispersion - multiplicative_noise^2/2
total_count <- exp(nb_log_dispersion)

lnnb_samples2 <- rmvlnnb(n = n_samples, logits = logits, total_count = total_count, multiplicative_noise_scale = multiplicative_noise)
View(lnnb_samples2)



## compare distribution with or without perturbation
# 1. simulate 1000 samples each from mvlnnb, when set perturbation to be all 0 or all 1 (cannot set all 1, cause problem).
# 2. compute empirical mean(sum) of each gene's expression level (count)
# 3. draw the histogram in same plot
# 4. plot other entries that is not affected by perturbation: log_var_mean, log_var_disp, size_factor_mean, batch_effect, covariate_effects

# 1~3 seem to not make sense
# 4. plot log_var_mean, size_factor_mean, batch_effect, covariate_effects, nb_log_mean_ctrl, nb_log_disp_ctrl

par(mfrow=c(3,2))
# blue
barplot(log_var_mean, las=2, col=rgb(0.25, 0.41, 0.88), border=NA, axisnames=FALSE,
        main="log_var_mean", xlab = "Gene", ylab = "Value")
abline(h=mean(log_var_mean), lty=2, col="black")
mtext(side=2, at=mean(log_var_mean), text=sprintf("%.2f", mean(log_var_mean)), las=2, line=0.2, cex=0.8, col="black")
# blue
barplot(size_factor_mean, las=2, col=rgb(0.25, 0.41, 0.88), border=NA, axisnames=FALSE,
        main="size_factor_mean", xlab = "Gene", ylab = "Value")
abline(h=size_factor_mean, lty=2, col="black")
mtext(side=2, at=size_factor_mean, text=sprintf("%.2f", size_factor_mean), las=2, line=0.2, cex=0.8, col="black")
# blue
barplot(batch_effect, las=2, col=rgb(0.25, 0.41, 0.88), border=NA, axisnames=FALSE,
        main="batch_effect", xlab = "Gene", ylab = "Value")
abline(h=mean(batch_effect), lty=2, col="black")
mtext(side=2, at=mean(batch_effect), text=sprintf("%.2f", mean(batch_effect)), las=2, line=0.2, cex=0.8, col="black")
# blue
barplot(covariate_effects, las=2, col=rgb(0.25, 0.41, 0.88), border=NA, axisnames=FALSE,
        main="covariate_effects", xlab = "Gene", ylab = "Value")
abline(h=mean(covariate_effects), lty=2, col="black")
mtext(side=2, at=mean(covariate_effects), text=sprintf("%.2f", mean(covariate_effects)), las=2, line=0.2, cex=0.8, col="black")
# blue
barplot(nb_log_mean_ctrl, las=2, col=rgb(0.25, 0.41, 0.88), border=NA, axisnames=FALSE,
        main="nb_log_mean_ctrl", xlab = "Gene", ylab = "Value")
abline(h=mean(nb_log_mean_ctrl), lty=2, col="black")
mtext(side=2, at=mean(nb_log_mean_ctrl), text=sprintf("%.2f", mean(nb_log_mean_ctrl)), las=2, line=0.2, cex=0.8, col="black")
# green
barplot(nb_log_disp_ctrl, las=2, col=rgb(0.13, 0.55, 0.13), border=NA, axisnames=FALSE,
        main="nb_log_disp_ctrl", xlab = "Gene", ylab = "Value")
abline(h=mean(nb_log_disp_ctrl), lty=2, col="black")
mtext(side=2, at=mean(nb_log_disp_ctrl), text=sprintf("%.2f", mean(nb_log_disp_ctrl)), las=2, line=0.2, cex=0.8, col="black")

# Reset plotting parameters to default
par(mfrow=c(1,1))




###################################
# create a new function to sample from mvlnnb (rmvlnnb_diy)
# with input: 
# - nb_log_mean_ctrl
# - perturbations
# - perturb_mean_lfc
# - nb_log_disp_ctrl
# - perturb_disp_lfc 
###################################
nb_log_mean <- nb_log_mean_ctrl + t(as.matrix(perturbations)) %*% as.matrix(perturb_mean_lfc)
nb_log_dispersion <- nb_log_disp_ctrl + t(as.matrix(perturbations)) %*% as.matrix(perturb_disp_lfc)

rmvlnnb_diy <- function(n = 100, 
                    nb_log_mean_ctrl,  # array, length = 5895
                    perturbations,     # array, length = 5859
                    perturb_mean_lfc,  # data.frame, dim = 1243 * 5859
                    nb_log_disp_ctrl,  # array, length = 5895
                    perturb_disp_lfc,  # data.frame, dim = 1243 * 5859
                    multiplicative_noise_scale = multiplicative_noise # matrix, dim = 1 * 5859
                    ){ 
  # output: a n*n_vars matrix. n is the number of samples to simulate. n_vars is the dim of variates
  
  nb_log_mean <- nb_log_mean_ctrl + t(as.matrix(perturbations)) %*% as.matrix(perturb_mean_lfc)  # matrix, dim = 1*5859, an input to the LNNB
  nb_log_dispersion <- nb_log_disp_ctrl + t(as.matrix(perturbations)) %*% as.matrix(perturb_disp_lfc)  # matrix, dim = 1*5859, an input to the LNNB
  
  logits <- nb_log_mean - nb_log_dispersion - multiplicative_noise^2/2
  total_count <- exp(nb_log_dispersion)
  
  samples <- rmvlnnb(n = n, logits = logits, total_count = total_count, multiplicative_noise_scale = multiplicative_noise_scale)
  
  return(samples)
}

# test
n_samples <- 1000
lnnb_samples3 <- rmvlnnb_diy(n = n_samples,
                             nb_log_mean_ctrl = nb_log_mean_ctrl, perturbations = rep(0, n_pert), perturb_mean_lfc = perturb_mean_lfc,
                             nb_log_disp_ctrl = nb_log_mean_ctrl, perturb_disp_lfc = perturb_disp_lfc, 
                             multiplicative_noise_scale = multiplicative_noise)

# a function to return the mean, variance, dispersion (later maybe other properties) of LNNB, same input as rmvlnnb_diy()
mvlnnb_diy_prop <- function(nb_log_mean_ctrl,  # array, length = 5895
                            perturbations,     # array, length = 5859
                            perturb_mean_lfc,  # data.frame, dim = 1243 * 5859
                            nb_log_disp_ctrl,  # array, length = 5895
                            perturb_disp_lfc,  # data.frame, dim = 1243 * 5859
                            multiplicative_noise_scale = multiplicative_noise # matrix, dim = 1 * 5859
                            ){
  # empty list to save results
  results <- list()
  
  nb_log_mean <- nb_log_mean_ctrl + t(as.matrix(perturbations)) %*% as.matrix(perturb_mean_lfc)
  results$mean <- exp(nb_log_mean)
  
  nb_log_dispersion <- nb_log_disp_ctrl + t(as.matrix(perturbations)) %*% as.matrix(perturb_disp_lfc)  # matrix, dim = 1*5859, an input to the LNNB
  results$disp <- exp(nb_log_dispersion)
  
  results$var <- exp(nb_log_mean) + (exp(multiplicative_noise_scale^2) * (1 + 1/exp(nb_log_dispersion)) - 1) * (exp(nb_log_mean))^2
  
  return(results)
}

# test
lnnb_samples3_prop <- mvlnnb_diy_prop(nb_log_mean_ctrl = nb_log_mean_ctrl, perturbations = rep(0, n_pert), perturb_mean_lfc = perturb_mean_lfc,
                                      nb_log_disp_ctrl = nb_log_disp_ctrl, perturb_disp_lfc = perturb_disp_lfc, 
                                      multiplicative_noise_scale = multiplicative_noise)
View(lnnb_samples3_prop$mean)
View(lnnb_samples3_prop$disp)
View(lnnb_samples3_prop$var)

## plot each gene's expression level's mean (&disp &var?)
# haha$mean is vector of length 5859, make it a bar plot and draw the dotted line as the mean of the means
# purpose: get a resonable range of pert to add up




### empirical power
# input: n_samples - number of simulations to make, lfc - perturbation effect size (for now just a constant added to nb_log_mean)
## methodology
# 0. simulate n_samples with perturbations all 0 (as control group)
# 1. randomly select 50 genes out of 5859 genes  --------------------------------------------------------------------------------------------------------- why 50(?)
# 2. add perturbation to the selected 50 gene, resimulate again (need to be able to change perturbation size)
# 3. compare the original simulation and perturbed simulation. do t-test to see which gene has significantly different expression. compute the type 2 error.
# 4. write nested loop. vary sample size & pert size. Draw 3 plots: sample_size vs power (fix pert_size); pert_size vs power (fix sample_size); heatmap vary both

## 0. simulate n_samples with perturbations all 0 (as control group)
n_samples <- 1000
perturbations <- rep(0, n_pert)
samples_orig <- rmvlnnb_diy(n = n_samples,
                            nb_log_mean_ctrl = nb_log_mean_ctrl, perturbations = perturbations, perturb_mean_lfc = perturb_mean_lfc,
                            nb_log_disp_ctrl = nb_log_disp_ctrl, perturb_disp_lfc = perturb_disp_lfc, 
                            multiplicative_noise_scale = multiplicative_noise)

## 1. randomly select 50 genes out of 5859 genes
set.seed(123)
n_pert_genes <- 50  # number of randomly chosen genes to perturb
pert_idx <- sample(1:5859, n_pert_genes, replace = FALSE)  # the indices of genes to add perturbation

## 2. add perturbation to the selected 50 gene, resimulate again
pert_size <- 0.1   # perturbation directly add to nb_log_mean (or nb_log_mean_ctrl, they do the same), leave nb_log_disp the same ------------------------ plus or minus(?) give both a try!

nb_log_mean_ctrl_add_pert <- copy(nb_log_mean_ctrl)  # make a copy to avoid changing the orig simulation
nb_log_mean_ctrl_add_pert[pert_idx] <- nb_log_mean_ctrl_add_pert[pert_idx] + pert_size

samples_pert <- rmvlnnb_diy(n = n_samples,
                            nb_log_mean_ctrl = nb_log_mean_ctrl_add_pert, perturbations = perturbations, perturb_mean_lfc = perturb_mean_lfc,
                            nb_log_disp_ctrl = nb_log_disp_ctrl, perturb_disp_lfc = perturb_disp_lfc, 
                            multiplicative_noise_scale = multiplicative_noise)

## 3. compare for each gene whether it is differently expressed, using t-test and look at the p-values
p_values <- sapply(1:n_gene, function(i) {
  t.test(samples_orig[,i], samples_pert[,i])$p.value
})

# get the indices of p-values < 0.05
sign_idx <- which(p_values <=0.05)

# compute type 1 error and type 2 error, compute
common_pert_idx <- sign_idx[sign_idx %in% pert_idx]

# power and specificity
false_detect_idx <- sign_idx[!sign_idx %in% pert_idx] # not perturbed, but detected

# empirical power from t-test
emp_power <- length(common_pert_idx) / length(pert_idx)  # (TPR) = 1 - type_1_error = number of correctly detected perts / number of actural perts
emp_spec <- 1 - (length(false_detect_idx)) / (n_gene - length(pert_idx)) # (TNG) = 1 - type_1_error = 1 - number of falsely rejected / number of unperturbed

print(paste("The power is", emp_power, "with type I error", round(1 - emp_spec, 4)))


## 4. write loop
library(tictoc)
#toy
sample_size_range <- c(100, 500)
pert_size_range <- c(0.1, 1)
#sample_size_range <- c(100, 500, 1000, 2000, 5000, 10000)
#pert_size_range <- c(0.01, 0.05, 0.1, 0.2, 0.5, 1, 2, 5)

# empty dataframe to same results
powers_df <- data.frame(sample_size = c(),
                       lfc_size = c(),
                       power = c())

tic()
for (n_samples in sample_size_range) {
  for (pert_size in pert_size_range) {
    # sample without perturbation
    perturbations <- rep(0, n_pert)
    samples_orig <- rmvlnnb_diy(n = n_samples,
                                nb_log_mean_ctrl = nb_log_mean_ctrl, perturbations = perturbations, perturb_mean_lfc = perturb_mean_lfc,
                                nb_log_disp_ctrl = nb_log_disp_ctrl, perturb_disp_lfc = perturb_disp_lfc, 
                                multiplicative_noise_scale = multiplicative_noise)
    
    # sample with perturbation
    nb_log_mean_ctrl_add_pert <- copy(nb_log_mean_ctrl)  # make a copy to avoid changing the orig simulation
    nb_log_mean_ctrl_add_pert[pert_idx] <- nb_log_mean_ctrl_add_pert[pert_idx] + pert_size
    
    samples_pert <- rmvlnnb_diy(n = n_samples,
                                nb_log_mean_ctrl = nb_log_mean_ctrl_add_pert, perturbations = perturbations, perturb_mean_lfc = perturb_mean_lfc,
                                nb_log_disp_ctrl = nb_log_disp_ctrl, perturb_disp_lfc = perturb_disp_lfc, 
                                multiplicative_noise_scale = multiplicative_noise)
    
    # compute p_value with t-test
    p_values <- sapply(1:n_gene, function(i) {
      t.test(samples_orig[,i], samples_pert[,i])$p.value
    })
    #get the indices of p-values < 0.05
    sign_idx <- which(p_values <=0.05)
    #compute type 1 error and type 2 error, compute
    common_pert_idx <- sign_idx[sign_idx %in% pert_idx]  # detected, really perturbed
    #empirical power from t-test
    emp_power <- length(common_pert_idx) / length(pert_idx)  # (TPR) = 1 - type_1_error = number of correctly detected perts / number of actural perts
    #put into result table
    powers_df <- rbind(powers_df, data.frame(sample_size = n_samples, lfc_size = pert_size, power = emp_power))
    
  }
}
toc()  # time similar with running on cluster or running locally
print(powers_df)


#toy
#sample_size_range <- c(100, 500)
#pert_size_range <- c(0.1, 1)
sample_size_range <- c(100, 500, 1000, 2000, 5000, 10000)
pert_size_range <- c(0.01, 0.05, 0.1, 0.2, 0.5, 1, 2)

# need Monte Carlo estimation and compute the mean of the powers
mc_powers_df <- data.frame(sample_size = rep(0, length(sample_size_range) * length(pert_size_range)),
                           lfc_size = rep(0, length(sample_size_range) * length(pert_size_range)),
                           power = rep(0, length(sample_size_range) * length(pert_size_range))
                           )

# add up all power estimations, then compute mean by divide by n_mc
n_mc <- 1000  # number of monte carlo estimations
tic()
for (k in 1:n_mc) {
  # empty dataframe to same results
  powers_df <- data.frame(sample_size = c(),
                          lfc_size = c(),
                          power = c())
  for (n_samples in sample_size_range) {
    for (pert_size in pert_size_range) {
      # sample without perturbation
      perturbations <- rep(0, n_pert)
      samples_orig <- rmvlnnb_diy(n = n_samples,
                                  nb_log_mean_ctrl = nb_log_mean_ctrl, perturbations = perturbations, perturb_mean_lfc = perturb_mean_lfc,
                                  nb_log_disp_ctrl = nb_log_disp_ctrl, perturb_disp_lfc = perturb_disp_lfc, 
                                  multiplicative_noise_scale = multiplicative_noise)
      
      # sample with perturbation
      nb_log_mean_ctrl_add_pert <- copy(nb_log_mean_ctrl)  # make a copy to avoid changing the orig simulation
      nb_log_mean_ctrl_add_pert[pert_idx] <- nb_log_mean_ctrl_add_pert[pert_idx] + pert_size
      
      samples_pert <- rmvlnnb_diy(n = n_samples,
                                  nb_log_mean_ctrl = nb_log_mean_ctrl_add_pert, perturbations = perturbations, perturb_mean_lfc = perturb_mean_lfc,
                                  nb_log_disp_ctrl = nb_log_disp_ctrl, perturb_disp_lfc = perturb_disp_lfc, 
                                  multiplicative_noise_scale = multiplicative_noise)
      
      # compute p_value with t-test
      p_values <- sapply(1:n_gene, function(i) {
        t.test(samples_orig[,i], samples_pert[,i])$p.value
      })
      #get the indices of p-values < 0.05
      sign_idx <- which(p_values <=0.05)
      #compute type 1 error and type 2 error, compute
      common_pert_idx <- sign_idx[sign_idx %in% pert_idx]  # detected, really perturbed
      #empirical power from t-test
      emp_power <- length(common_pert_idx) / length(pert_idx)  # (TPR) = 1 - type_1_error = number of correctly detected perts / number of actural perts
      #put into result table
      powers_df <- rbind(powers_df, data.frame(sample_size = n_samples, lfc_size = pert_size, power = emp_power))
    }
  }
  mc_powers_df["power"] <- mc_powers_df["power"] + powers_df["power"]
}
toc()
# finalize the output table
mc_powers_df["sample_size"] <- powers_df["sample_size"]
mc_powers_df["lfc_size"] <- powers_df["lfc_size"]
mc_powers_df["power"] <- mc_powers_df["power"] / n_mc
View(mc_powers_df)

write.csv(mc_powers_df, "powers.csv", row.names = FALSE)



# plot the lfc against gene expression level
# plot the effect size against gene expression level
# i.e. x-axis: effect size c(0.01, 0.1, 0.2, 0.5, 1, ...) | y-axis: orig gene expression level



a <- matrix(rep(0,6), nrow = 3)
a_df <- data.frame(a)
names(a_df) <- c("A", "B")
a_df

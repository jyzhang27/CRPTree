## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning=FALSE, message=FALSE)
library(CRPTree)
library(ggplot2)
set.seed(1234)


## -----------------------------------------------------------------------------
N <- 30
B <- 10
alpha <- 2
tree_ex <- rcrp_tree(N, B, alpha)


## -----------------------------------------------------------------------------
plot_crp_tree(tree_ex)
plot_crp_tree(tree_ex, TRUE)

## -----------------------------------------------------------------------------
pcrp_tree(tree_ex, 1)
pcrp_tree(tree_ex, 2)
pcrp_tree(tree_ex, 5)


## -----------------------------------------------------------------------------
pnull_tree(tree_ex)

## -----------------------------------------------------------------------------
count_same_attachments(tree_ex)

## -----------------------------------------------------------------------------
s_0 <- compute_S_mean(tree_ex, 500)
s_0

## -----------------------------------------------------------------------------
S_obs <- simulate_S_observed(tree_ex, 50)
S_mean <- mean(S_obs)

S_tree_null <- simulate_S_tree_null(tree_ex, 50)
pval_tree <- mean(rbind(outer(S_tree_null, S_obs, '>='), rep(TRUE, length(S_obs))))

ggplot(as.data.frame(S_tree_null), aes(S_tree_null)) + geom_bar(color='black', fill='gray') + theme_light()+ xlab('S') + ggtitle('Approximate null distribution of S conditional on T_N')
  

## -----------------------------------------------------------------------------
S_mean_null <- simulate_S_mean_tree(tree_ex, 50, 49)
pval_avg <- (1+ sum(S_mean_null >= S_mean))/(1+length(S_mean_null))

ggplot(as.data.frame(S_mean_null), aes(S_mean_null)) + 
  geom_histogram(color='black', fill='gray', bins=20) + theme_light()+
  scale_x_continuous(limits=c(floor(min(S_mean_null)), ceiling(max(c(max(S_mean_null), S_mean))))) +
  geom_vline(xintercept=S_mean, color='red') + xlab('S_0') + 
  ggtitle('Approximate null distribution of S_0, Red line is observed value')

## -----------------------------------------------------------------------------
one_tree_all_methods(tree_ex, 50, 49, 49)

## -----------------------------------------------------------------------------
one_tree_all_methods(tree_ex, 50, 49, 49, alt_methods = TRUE)

## -----------------------------------------------------------------------------
alphas <- c(2,5, 10, 20)
crp_tree_power(tree_ex, alphas)


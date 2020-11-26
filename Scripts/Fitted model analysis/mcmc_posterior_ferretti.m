% Import posterior parameter distributions for the Ferretti model obtained
% using data augmentation MCMC, and calculate the posterior distribution of
% the proportion of presymptomatic transmissions.

clear all; close all; clc;

addpath('../../Data')
addpath('../../Results')
addpath('../../Functions/Ferretti')

% Load known model parameters
load('../../Data/input_data.mat','inc_shape','inc_scale')

% Load output of MCMC fitting procedure
load('../../Results/param_fit_mcmc_ferretti.mat','theta_mat','ll_vec','acceptance_rate')

m_inc = inc_shape*inc_scale;

% Calculate posterior distributions of individual model parameters,
% assuming uniform priors for mu, sigma and alpha
mu_post = theta_mat(:,1);
sigma_post = exp(theta_mat(:,2));
alpha_post = exp(theta_mat(:,3));

% Calculate posterior distribution for the proportion of presymptomatic
% transmissions
prob_presymp_post = ((1+exp(mu_post./sigma_post)).^(-alpha_post)-(1+exp((m_inc+mu_post)./sigma_post)).^(-alpha_post))./(1-(1+exp((m_inc+mu_post)./sigma_post)).^(-alpha_post));

% Find parameter values giving the best fit to the data
[ll_best,posn_best] = max(ll_vec);
params_best = get_params_ferretti(theta_mat(posn_best,:),m_inc);
prob_presymp_best = prob_presymp_post(posn_best);

% Remove burn-in and thin posteriors
no_steps = length(mu_post);
indices_kept = (500001):100:no_steps;

mu_post = mu_post(indices_kept);
sigma_post = sigma_post(indices_kept);
alpha_post = alpha_post(indices_kept);
prob_presymp_post = prob_presymp_post(indices_kept);

% Median and 95% confidence bounds for the proportion of presymptomatic
% transmissions
quantile(prob_presymp_post,[0.025,0.5,0.975])

% Save results
save('../../Results/mcmc_posterior_ferretti','ll_best','params_best','prob_presymp_best','prob_presymp_post')

rmpath('../../Data')
rmpath('../../Results')
rmpath('../../Functions/Ferretti')
% Import posterior parameter distributions for the independent transmission
% and symptoms model obtained using data augmentation MCMC, and calculate
% the posterior distribution of the proportion of presymptomatic
% transmissions.

clear all; close all; clc;

addpath('../../Data')
addpath('../../Results')
addpath('../../Functions/Indep')

% Load known model parameters
load('../../Data/input_data.mat','inc_shape','inc_scale')

% Load output of MCMC fitting procedure
load('../../Results/param_fit_mcmc_indep.mat','theta_mat','ll_vec','acceptance_rate')

% Calculate posterior distributions of individual model parameters,
% assuming uniform priors
gen_mean_post = exp(theta_mat(:,1));
gen_var_post = exp(theta_mat(:,2));
gen_shape_post = (gen_mean_post.^2)./gen_var_post;
gen_scale_post = gen_var_post./gen_mean_post;

% Find parameter values giving the best fit to the data
[ll_best,posn_best] = max(ll_vec);
params_best = [gen_shape_post(posn_best),gen_scale_post(posn_best)];
prob_presymp_best = get_presymp_trans_probs_indep(params_best(1),params_best(2),inc_shape,inc_scale);

% Remove burn-in and thin posteriors
no_steps = length(gen_mean_post);
indices_kept = (500001):100:no_steps;

gen_mean_post = gen_mean_post(indices_kept);
gen_var_post = gen_var_post(indices_kept);
gen_shape_post = gen_shape_post(indices_kept);
gen_scale_post = gen_scale_post(indices_kept);

% Calculate posterior distribution for the proportion of presymptomatic
% transmissions
prob_presymp_post = get_presymp_trans_probs_indep(gen_shape_post,gen_scale_post,inc_shape,inc_scale);

% Median and 95% confidence bounds for the proportion of presymptomatic
% transmissions
quantile(prob_presymp_post,[0.025,0.5,0.975])

% Save results
save('../../Results/mcmc_posterior_indep','ll_best','params_best','prob_presymp_best','prob_presymp_post')

rmpath('../../Data')
rmpath('../../Results')
rmpath('../../Functions/Indep')
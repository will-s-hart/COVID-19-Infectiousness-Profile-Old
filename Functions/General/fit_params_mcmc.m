function [theta_mat,ll_vec,acceptance_rate] = fit_params_mcmc(ll_form,data_struct_observed,no_steps,theta_prop_sd,theta_init,t_s_nochange_prob)
    
    % Use data augmentation MCMC to fit the parameters, theta, of the
    % model of infectiousness under consideration, to the observed
    % transmission pair data contained in data_struct.
    
    % Number of fitted parameters
    no_params_fitted = length(theta_init);
    
    % Create matrices containing each possible symptom onset time for each
    % source and recipient
    t_s1_mat = cell2mat((data_struct_observed.t_s1_data)')';
    t_s2_mat = cell2mat((data_struct_observed.t_s2_data)')';
    
    % Number of transmission pairs
    no_pairs = size(t_s1_mat,1);
    
    % Length of grid of possible onset times
    t_s_grid_length = size(t_s1_mat,2);

    % Initialise vector of parameters, theta
    theta = theta_init;
    
    % Initialise vectors t_s1 and t_s2 containing onset times for each
    % source and recipient
    t_s1_inds = randi(t_s_grid_length,no_pairs,1);
    t_s2_inds = randi(t_s_grid_length,no_pairs,1);
    t_s1 = t_s1_mat(sub2ind([no_pairs,t_s_grid_length],(1:no_pairs)',t_s1_inds));
    t_s2 = t_s2_mat(sub2ind([no_pairs,t_s_grid_length],(1:no_pairs)',t_s2_inds));
    
    % Structure array containing initial augmented data
    data_struct_augmented_init = rmfield(data_struct_observed,{'t_s1_data','t_s2_data'});
    data_struct_augmented_init.t_s1_data = t_s1; data_struct_augmented_init.t_s2_data = t_s2;
    
    % Calculate initial likelihood
    ll = ll_form(theta,data_struct_augmented_init);

    % Create vectors to hold output of fitting procedure
    theta_mat = zeros(no_steps,no_params_fitted); %parameters at each step
    acceptance_vec = zeros(no_steps,1); %to record which steps are accepted
    ll_vec = zeros(no_steps,1); %likelihood at each step
    
    % Structure array to populate with augmented data at each step in the
    % chain
    data_struct_augmented_prop = data_struct_augmented_init;
    
    % Break steps into 100 groups to record progress
    step_no_mat = reshape(1:no_steps,no_steps/100,100);
    
    % Run chain
    for j = 1:size(step_no_mat,2)
        for i = 1:size(step_no_mat,1)

            step_no = step_no_mat(i,j);

            if mod(step_no,2)>0.5
                
                % For odd step numbers, update model parameters using
                % independent normal proposal distributions
                
                theta_prop = theta + randn(1,no_params_fitted).*theta_prop_sd;
                
                t_s1_inds_prop = t_s1_inds;
                t_s2_inds_prop = t_s2_inds;
                
                t_s1_prop = t_s1;
                t_s2_prop = t_s2;
            else
                
                % For even step numbers, update the symptom onset times
                % (each onset time index stays the same with probabilities
                % t_s_nochange_prob, and else jumps up or down by one with
                % equal probabilities)
                
                theta_prop = theta;
                
                t_s1_inds_prop = round(t_s1_inds-1/(2*t_s_nochange_prob)+rand(no_pairs,1)/t_s_nochange_prob);
                t_s1_inds_prop = max(min(t_s1_inds_prop,t_s_grid_length),1); %reset proposed values outside the grid to remain at endpoints
                t_s2_inds_prop = round(t_s2_inds-1/(2*t_s_nochange_prob)+rand(no_pairs,1)/t_s_nochange_prob);
                t_s2_inds_prop = max(min(t_s2_inds_prop,t_s_grid_length),1); %reset proposed values outside the grid to remain at endpoints

                t_s1_prop = t_s1_mat(sub2ind([no_pairs,t_s_grid_length],(1:no_pairs)',t_s1_inds_prop));
                t_s2_prop = t_s2_mat(sub2ind([no_pairs,t_s_grid_length],(1:no_pairs)',t_s2_inds_prop));
            end
            
            % Populate the structure array data_struct_augmented_prop with
            % the proposed onset times
            data_struct_augmented_prop.t_s1_data = t_s1_prop;
            data_struct_augmented_prop.t_s2_data = t_s2_prop;
            
            % Calculate the likelihood using the proposed parameters and
            % onset times
            ll_prop = ll_form(theta_prop,data_struct_augmented_prop);
            
            % Ratio between proposed and previous likelihoods
            a = exp(ll_prop-ll);
            
            % Accept the proposed parameters and onset times with
            % probability a
            if rand < a
                
                % Update the current parameters and onset times
                theta = theta_prop;
                t_s1 = t_s1_prop;
                t_s2 = t_s2_prop;
                t_s1_inds = t_s1_inds_prop;
                t_s2_inds = t_s2_inds_prop;

                ll = ll_prop;

                acceptance_vec(step_no) = 1;
            end
            
            % Record parameter values and likelihood
            theta_mat(step_no,:) = theta;
            ll_vec(step_no) = ll;
        end
        
        % Display progress when an integer percentage of steps has been
        % completed
        fprintf('%d%% complete\n',100*step_no/no_steps);
    end
    
    % Calculate acceptance rates: overall, and when either the parameters
    % or the onset times were updated.
    acceptance_rate_overall = mean(acceptance_vec)
    acceptance_rate_theta = mean(acceptance_vec(1:2:end))
    acceptance_rate_data = mean(acceptance_vec(2:2:end))

    acceptance_rate.overall = acceptance_rate_overall;
    acceptance_rate.theta = acceptance_rate_theta;
    acceptance_rate.data = acceptance_rate_data;
end
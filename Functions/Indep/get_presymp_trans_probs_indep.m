function p_vec = get_presymp_trans_probs_indep(gen_shape_vec,gen_scale_vec,inc_shape,inc_scale)

    % Numerically calculate the proportion of presymptomatic transmissions for
    % the independent transmission and symptom model, for each set of
    % generation time parameters given by the entries in gen_shape_vec and
    % gen_scale_vec, by integrating the generation time distribution weighted
    % by the proportion of hosts who have developed symptoms.

    % Grid of times since infection to integrate over
    t_max = 50;
    dt = 0.1;
    t_vec = ((dt/2):dt:(t_max-dt/2))';

    % Grids of pairs of values of the shape/scale parameter of the generation
    % time distribution, and of the time since infection
    [gen_shape_grid,t_grid] = ndgrid(gen_shape_vec,t_vec);
    [gen_scale_grid,~] = ndgrid(gen_scale_vec,t_vec);

    % Evaluate the density of the generation time and the cumulative
    % distribution of the incubation period, at each grid point
    F_inc_vec = gamcdf(t_vec,inc_shape,inc_scale);
    f_gen_grid = gampdf(t_grid,gen_shape_grid,gen_scale_grid);

    % Numerically integrate over the generation time distribution, weighted by
    % the proportion of hosts who have developed symptoms, to calculate the
    % proportion of presymptomatic transmissions for each parameter set.
    p_vec = 1-f_gen_grid*F_inc_vec*dt;

end
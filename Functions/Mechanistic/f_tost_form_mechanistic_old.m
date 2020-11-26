function f_tost = f_tost_form_mechanistic_old(t_tost,params)

    % TOST distribution for our mechanistic approach with parameters given
    % by params.

    gamma = params(1); mu = params(2);
    k_inc = params(3); k_E = params(4); k_I = params(5);
    alpha = params(6);
    k_P = k_inc-k_E;

    m_P = k_P/(k_inc*gamma);
    m_I = 1/mu;

    q_P = (alpha*m_P)/((alpha*m_P)+m_I);
    q_I = 1-q_P;

    fm = (q_P/m_P)*(1-gamcdf(-t_tost,k_P,1/(k_inc*gamma)));
    fp = (q_I/m_I)*(1-gamcdf(t_tost,k_I,1/(k_I*mu)));

    indm = (t_tost<0);
    ind0 = (t_tost==0);
    indp = logical(1-(indm+ind0));

    f_tost = zeros(size(t_tost));
    f_tost(indm) = fm(indm);
    f_tost(indp) = fp(indp);
    f_tost(ind0) = fp(ind0);

end
% Calculate both the proportion of contacts identified for different
% contact elicitation windows, and the reduction in transmission through
% isolation for different times between exposure and isolation.

% This script requires the Chebfun package to run (freely available at
% https://www.chebfun.org/download/).

clear all; close all; clc;
splitting on

addpath('../../Results')

load('../../Results/gen_tost_serial_indep.mat','f_gen_indep','f_tost_indep','f_serial_indep')
load('../../Results/gen_tost_serial_constinf.mat','f_gen_constinf','f_tost_constinf','f_serial_constinf')
load('../../Results/gen_tost_serial_varinf.mat','f_gen_varinf','f_tost_varinf','f_serial_varinf')
load('../../Results/gen_tost_serial_ferretti.mat','f_gen_ferretti','f_tost_ferretti','f_serial_ferretti')

% Proportion of infectious presymptomatic contacts found when
% contacts of a symptomatic host are traced

F_tost_indep = cumsum(f_tost_indep);
F_tost_constinf = cumsum(f_tost_constinf);
F_tost_varinf = cumsum(f_tost_varinf);
F_tost_ferretti = cumsum(f_tost_ferretti);

t_traced_from = chebfun('t',[0,7]);
contacts_traced_indep = 1-F_tost_indep(-t_traced_from);
contacts_traced_constinf = 1-F_tost_constinf(-t_traced_from);
contacts_traced_varinf = 1-F_tost_varinf(-t_traced_from);
contacts_traced_ferretti = 1-F_tost_ferretti(-t_traced_from);

% Transmissions prevented when contacts isolated at different
% times since exposure

F_gen_indep = cumsum(f_gen_indep);
F_gen_constinf = cumsum(f_gen_constinf);
F_gen_varinf = cumsum(f_gen_varinf);
F_gen_ferretti = cumsum(f_gen_ferretti);

t_isolated = chebfun('t',[0,7]);
transmissions_prevented_indep = 1-F_gen_indep(t_isolated);
transmissions_prevented_constinf = 1-F_gen_constinf(t_isolated);
transmissions_prevented_varinf = 1-F_gen_varinf(t_isolated);
transmissions_prevented_ferretti = 1-F_gen_ferretti(t_isolated);

save('../../Results/contact_tracing','contacts_traced_indep','contacts_traced_constinf','contacts_traced_varinf','contacts_traced_ferretti','transmissions_prevented_indep','transmissions_prevented_constinf','transmissions_prevented_varinf','transmissions_prevented_ferretti')

rmpath('../../Results')
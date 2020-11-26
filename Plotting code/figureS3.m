% Produce the panels in Fig. S3 of our manuscript.

% This script requires the Chebfun package to run (freely available at
% https://www.chebfun.org/download/), and requires the export_fig package
% (freely available at https://github.com/altmany/export_fig) to export the
% figure panels to pdf files.

clear all; close all; clc;

addpath('../Results')

load('../Results/contact_tracing','contacts_traced_indep','contacts_traced_constinf','contacts_traced_varinf','contacts_traced_ferretti','transmissions_prevented_indep','transmissions_prevented_constinf','transmissions_prevented_varinf','transmissions_prevented_ferretti')

% Colorscheme

c1 = [0, 0.4470, 0.7410]; c2 = [0.8500, 0.3250, 0.0980];
c3 = [0.9290, 0.6940, 0.1250]; c4 = [0.4940, 0.1840, 0.5560];

% Plot figures

for k = 1:6
figsetup(k)
end

% Time to go back when contact tracing, for different values of the contact
% tracing effectiveness

tracing_eff_vec = [0.8,0.6,0.4];

for k = 1:3

    tracing_eff = tracing_eff_vec(k);   

    figure(k); hold on;
    plot(100*tracing_eff*contacts_traced_constinf,'color',c2,'linewidth',3)
    plot(100*tracing_eff*contacts_traced_ferretti,'color',c3,'linewidth',3)
    plot(100*tracing_eff*contacts_traced_indep,'color',c4,'linewidth',3)
    plot(100*tracing_eff*contacts_traced_varinf,'-.','color',c1,'linewidth',3)
    xticks(0:7);
    ylim([0,80])
    xlabel('Time contacts traced before onset (days)')
    ylabel('Infectious contacts found (%)')
end

% Time in infection isolated, for different values of the isolation
% effectiveness

isolation_eff_vec = [0.8,0.6,0.4];

for k = 1:3

    isolation_eff = isolation_eff_vec(k);   

    figure(k+3); hold on;
    plot(100*isolation_eff*transmissions_prevented_constinf,'color',c2,'linewidth',3)
    plot(100*isolation_eff*transmissions_prevented_ferretti,'color',c3,'linewidth',3)
    plot(100*isolation_eff*transmissions_prevented_indep,'color',c4,'linewidth',3)
    plot(100*isolation_eff*transmissions_prevented_varinf,'-.','color',c1,'linewidth',3)
    xticks(0:7);
    ylim([0,80])
    xlabel('Time since infection isolated (days)')
    ylabel('Onward transmissions prevented (%)')
end

for k = 1:6
figsetup(k)
end

figure(1); export_fig Figures/figS3a.pdf -nocrop -transparent
figure(2); export_fig Figures/figS3b.pdf -nocrop -transparent
figure(3); export_fig Figures/figS3c.pdf -nocrop -transparent
figure(4); export_fig Figures/figS3d.pdf -nocrop -transparent
figure(5); export_fig Figures/figS3e.pdf -nocrop -transparent
figure(6); export_fig Figures/figS3f.pdf -nocrop -transparent

rmpath('../Results')
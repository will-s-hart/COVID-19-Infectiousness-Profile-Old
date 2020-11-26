% Produce the panels in Fig. 4 of our manuscript.

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

for k = 1:2
figsetup(k)
end

% Time to go back when contact tracing

figure(1); hold on;
plot(100*contacts_traced_constinf,'color',c2,'linewidth',3)
plot(100*contacts_traced_ferretti,'color',c3,'linewidth',3)
plot(100*contacts_traced_indep,'color',c4,'linewidth',3)
plot(100*contacts_traced_varinf,'-.','color',c1,'linewidth',3)
xticks(0:7);
ylim([20,100])
xlabel('Time contacts traced before onset (days)')
ylabel('Infectious contacts found (%)')

% Time in infection contacts isolated (assuming not asymptomatic)

figure(2); hold on;
plot(100*transmissions_prevented_constinf,'color',c2,'linewidth',3)
plot(100*transmissions_prevented_ferretti,'color',c3,'linewidth',3)
plot(100*transmissions_prevented_indep,'color',c4,'linewidth',3)
plot(100*transmissions_prevented_varinf,'-.','color',c1,'linewidth',3)
xticks(0:7);
ylim([20,100])
xlabel('Time since infection isolated (days)')
ylabel('Onward transmissions prevented (%)')

for k = 1:2
figsetup(k)
end

figure(1); export_fig Figures/fig4a.pdf -nocrop -transparent
figure(2); export_fig Figures/fig4b.pdf -nocrop -transparent

rmpath('../Results')
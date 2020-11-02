close all
clear all
clc

load('param_new_optim_31.mat')
load('cost_new_optim_31.mat')
ideale_ottim_23 = [0.6541    0.4364    0.5585    0.3132    0.9000    0.2417];
ideale_ottim_24 = [0.6541    0.1939    0.5585    0.3132    0.9000    0.2417];
ideale_ottim_25 = [0.6541    0.1939    0.5585    0.3132    0.7500    0.2417]; % anche per 26,27,28,29,30,32, 33
ideale_ottim_31 = [0.5739    0.0701    0.5585    0.3132    0.7500    0.2417];

names = ["K_a_d_a_p", "K2", "A2", "K1", "A1", "I_e"];
legend_values = [];
median_parameters = [];
for i=1:length(parametri)
    parameters = parametri{1,i};
    cost = cost_function{1,i};
    legend_values = [legend_values; string(cost(end))];
    for j=1:size(parameters,2)
        line_width = 2;
        median_parameters(i,j) = parameters(end,j);
        x=1:length(parameters);
        y=parameters(:,j);
        figure(1)
        subplot(2,3,j)
        plot(x,y, 'LineWidth',line_width)   
        axis([0 inf 0 1])
        title(names(j))
        hold on
        figure(2)
        plot(cost)
        hold on
    end
end
% figure(1)
% for j=1:size(parameters,2)
%     x = 1:100;
%     y = ideale_ottim_25(j)*ones(length(x),1);
%     subplot(2,3,j)
%     hold on
%     plot(x,y,'k--','Linewidth',3)
%     legend(legend_values)
%     xlabel('# optimization step')
%     ylabel('normalized value')
% end
figure(1)
for j=1:size(parameters,2)
    subplot(2,3,j)
    legend(legend_values)
    xlabel('# optimization step')
    ylabel('normalized value')
end
figure(2)
xlabel('# optimization step')
ylabel('Cost function value')
title('Cost function evolution during optimization')

median_parameters = median(median_parameters,1);

Cm = 104;
tau_m = -45.76;
m_IF = 0;
t_ref = 1.65;
% low = K_adap, K2, A2, K1, A1, I_e

low = [Cm/(tau_m^2)+0.000001,-1/tau_m+0.00001,0.001,3/((1/m_IF)*1000),0.001,-150];
up_2 = 10*low(2);
up = [((up_2-1/tau_m)^2)*Cm/4-0.000001,up_2,100,3/t_ref,600,-50];

% low = [Cm/(tau_m^2)+0.000001,-1/tau_m+0.025,0.001,3/((1/m_IF)*1000),0.001,-150];
% up_2 = 5*low(2);
% up = [((up_2-1/tau_m)^2)*Cm/4-0.000001,up_2,100,3/t_ref,600,-50];

median_value = median_parameters.*(up-low)+low
    
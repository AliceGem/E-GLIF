% Check solution


% Model definition - 3D linear ODE system
% y = Vm(t); Iadap(t); Idep2(t)          % State variables


clear
close all

time_window = [0 60]
initial_conditions = [-65; 1453; -0.03]
[t,y] = ode45(@eglif,time_window,initial_conditions);

plot(t,y(:,1),'-',t,y(:,2),'-',t,y(:,3),'-', 'Linewidth',2)
title('E-GLIF Solution with ODE45');
xlabel('Time t');
ylabel('Solution y');
legend('V_m','I_a_d_a_p', 'I_d_e_p','location','northeastoutside')

function dydt = eglif(t,y)
    tau_m = 15.58
    E_L = -65
    Ie = 0
    Ist = 200
    Cm = 189.79
    k_adap = 0.623
    k2 = 0.064184852
    k1 = -0.254
    dydt = [(1/tau_m)*y(1) - E_L/tau_m + Ie/Cm + Ist/Cm + y(3)/Cm - y(2)/Cm;...
            k_adap*y(1) - k2*y(2) - k_adap*E_L;...
            -k1*y(3)];
end


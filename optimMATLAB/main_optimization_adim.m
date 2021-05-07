% E-GLIF neuron model optimization
% LIF neuron + adaptive current for modelling subthreshold-related
% mechanisms and Spike-Frequency Adaptation (SFA) + spike-triggered current modelling
% Na channel-related fast mechanisms; therefore, a simple linear model (useful for
% large scale simulations) and biologically plausible.
% References: [Geminiani et al., Front Neuroinform, 2018]; [Geminiani et al.,  Front Comput Neurosci, 2019]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: the parameter k2 should be set to -1/tau_m, for neurons with self-sustained
% oscillations of the membrane potential (thus, oscillatory not damped solution of the model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

% Neurons to optimize.
% In this example:
% MLI = CA1 Pyramidal Cell,
% DCN = Deep Cerebellar Nuclei neurons (glutamatergic/GAD-negative large neurons)
neu_name = 'CA1-PC';
%i = 1;      % i-th neuron selected - 1 is CA1 PC, 2 is DCN in this case

% Step 1: set passive membrane parameters (e.g. Cm, tau_m, etc) from neurophysiology experiments
% for each of the neurons to be optimized. Reference values should be taken from papers
% (for consistency with reference stimulation protocol) or neuroelectro.org
% (if not available in papers).
% Each parameter is saved in an array of values for each neuron
% CA1 PC data from Helene MM paper, control young mice for protocol stim;
% passive properties from neuroelectro.org
Cm = 90.0;        % [pF]
tau_m = -15.0;   % [ms] should be the opposite of the given value - to fix
t_ref = 2.15;        % [ms]    ---> to check!
E_L = -65;      % [mV]
Vth = -48;      % [mV]          % to check if we want to stay below threshold
Vr = -55;       % = E_L-10 [mV] ---> could it be correct (Vth - AHP amplitude) for hippocampus?


% Step 2: set the target input-output (Istim-firing frequency) relationship from literature stimulation protocols
% For target frequencies, mean and Standard Deviation (SD) values are considered, in order to fit a distribution

% Intrinsic firing frequency (/autorhythm/spontaneous firing):
m_IF = 0;        % Mean intrinsic frequency
sd_IF = 0;     % Standard Deviation of intrinsic frequency (!SE is reported in some studies!)

% Depolarization phases:
% * input current Istim = [3xnn] in [pA], so we use 3 values of input current for 3 depolarazion phases/steps,
% for the nn neurons considered for optimization
% * mean target frequency during depolarization m_Fdep = [3xnn] in [Hz]
% * SD of target frequency during depolarization sd_Fdep = [3xnn] in [Hz]
Istim = [200.0; 400.0; 600.0];             % Annalisa sims                          %[100.0; 200.0; 300.0], ...     % Helene paper Table 1
         
%TODO: add reference spike times for hippocampus

m_Fdep = {[[22.2; 25; 37]...               % Frequency of first ISI - mean; time*tau_m = ISI from graphs Annalisa
          [50.0; 80.0; 110.0]], [[33.3; 45.4; 166.7]...               % Frequency of second ISI - mean
          [50.0; 80.0; 110.0]], [[13.3; 35.7; 100]...               % Frequency of third ISI - mean
          [50.0; 80.0; 110.0]],[[11.1; 19.2; 62.5]...               % Frequency of fourth ISI - mean
          [50.0; 80.0; 110.0]] };

sd_Fdep = {[[4.8; 6.6; 6.2]...
           [2.0; 5.0; 15.0]], [[4.8; 6.6; 6.2]...
           [2.0; 5.0; 15.0]], [[4.8; 6.6; 6.2]...
           [2.0; 5.0; 15.0]], [[4.8; 6.6; 6.2]...
           [2.0; 5.0; 15.0]]};

% Target frequency at the end of a depolarization step should take into account SFA,
% using the parameter SFA_gain = ratio between initial and steady-state firing rate [nnx3]
% (it should be set to 1 if no SFA is present)
SFA_gain = {[2.0, 2.4, 3.0; ...
            1.2, 1.2, 1.2], [2.0, 2.4, 3.0; ...
            1.2, 1.2, 1.2], [2.0, 2.4, 3.0; ...
            1.2, 1.2, 1.2], [2.0, 2.4, 3.0; ...
            1.2, 1.2, 1.2]};

% Hyperpolarization phase:
% * input current Iinh [pA]
% * minimum Vm value during hyperpolarization Vinh_min [mV]
% * steady-state Vm value during hyperpolarization Vinh_ss [mV]
Iinh = -100.0; % values to check
Vinh_min = -95;
Vinh_ss = -90;

% Following hyperpolarization, a rebound burst is present in some neuron types
% Burst frequency is equal to intrinsic frequency if no rebound burst is present (e.g. for DCN_gly)
m_Fburst = [m_IF(1), m_IF(2)*2];           % [Hz]
sd_Fburst = [sd_IF(2), sd_IF(2)];


% Step 3: deriving a spike times distribution of 'ne' samples to fit during optimization

ne = 10;      % 'ne' target values for each stimulation step

% Target spiking times during spontaneous firing (Istim = 0)
T_tonic(i,:) = (1./(sd_IF(i).*randn(1,ne) + m_IF(i)))*1000;

% Avoid negative values in the distribution
while ~isempty(find(T_tonic<0))
    T_tonic(i,find(T_tonic<0))=(1./(sd_IF(i).*randn(1,length(find(T_tonic<0))) + m_IF(i)))*1000;
end

% Target spiking times during depolarization phases (Istim > 0)
for sp = 1:length(m_Fdep)
T_dep1(i,:) = (1./(sd_Fdep(1,i).*randn(1,ne) + m_Fdep(1,i)))*1000;
T_dep2(i,:) = (1./(sd_Fdep(2,i).*randn(1,ne) + m_Fdep(2,i)))*1000;
T_dep3(i,:) = (1./(sd_Fdep(3,i).*randn(1,ne) + m_Fdep(3,i)))*1000;
end
Tdep = {T_dep1, T_dep2, T_dep3};

% Target spiking times following hyperpolarization (Istim < 0) - to fit rebound bursting if present
Tburst = (1./(sd_Fburst(i).*randn(1,ne) + m_Fburst(i)))*1000;
Tlb = 5.*randn(1,ne) + 1000*(1/mean(m_IF(i)));

% Target spiking times in the afterhyperpolarization (AHP) when returning to instrinsic firing after depolarization
% In the cerebellum, taken into account only for optimization tests on the Golgi cell
Tahp = [5.*randn(ne,1) + 80, 5.*randn(ne,1) + 100, 5.*randn(ne,1) + 120];


%% Model definition - 3D linear ODE system
syms Vm(t) Iadap(t) Idep(t)          % State variables
syms alfa beta delta IaA0 IdA0...            % Parameters to be optimized
     Ist


ode1 = diff(Vm) == alfa-beta*(Iadap-Idep)+delta*(1+Vm); (-1/tau_m(i))*Vm + Ie/Cm(i) + Ist/Cm(i) + Idep/Cm(i) - Iadap/Cm(i) + E_L(i)/tau_m(i);
ode2 = diff(Iadap) == 1-Iadap + Vm;
ode3 = diff(Idep) == -beta*Idep;


I = Ie+Ist;


% Analytic solution
Psi = sqrt(); %beta^2 + (beta-1)*delta;

Vm = @(t,alfa,beta,delta,IaA0,IdA0) (1/2).*(beta+(-1).*delta).^(-1).*(beta.^2+((-1)+beta).*delta).^(-1).*(4.*beta+( ...
  -1).*(1+delta).^2).^(-1).*Psi.*(2.*exp(1).^(((-1).*t+t0).*beta).*IdA0.*( ...
  (-1)+beta).*beta.*(beta+(-1).*delta).*Psi+(-2).*(alpha+(-1).*beta+delta).*(beta.^2+(( ...
  -1)+beta).*delta).*Psi+exp(1).^((1/2).*(t+(-1).*t0).*((-1)+delta+(-1).*Psi)) ...
  .*(IdA0.*beta.*(beta+(-1).*delta).*((-1)+(-1).*delta+beta.*(3+delta+(-1).*Psi)+Psi) ...
  +(-1).*(beta.^2+(-1).*delta+beta.*delta).*(alpha.*(1+(-2).*beta+delta+(-1).*Psi)+(beta+ ...
  (-1).*delta).*((-1)+2.*IaA0.*beta+(-1).*delta+Psi+V0.*((-1)+(-1).*delta+Psi)))) ...
  +exp(1).^((1/2).*(t+(-1).*t0).*((-1)+delta+Psi)).*((-1).*IdA0.*beta.*( ...
  beta+(-1).*delta).*((-1)+(-1).*delta+(-1).*Psi+beta.*(3+delta+Psi))+(beta.^2+(-1).* ...
  delta+beta.*delta).*(alpha.*(1+(-2).*beta+delta+Psi)+(beta+(-1).*delta).*((-1)+2.*IaA0.* ...
  beta+(-1).*delta+(-1).*Psi+(-1).*V0.*(1+delta+Psi)))));

Iadap = @(t,alfa,beta,delta,IaA0,IdA0) (1/2).*exp(1).^(t0+(-1).*t0.*delta+(-1/2).*t.*((-1)+2.*beta+delta+Psi)).*( ...
  beta+(-1).*delta).^(-1).*(beta.^2+((-1)+beta).*delta).^(-1).*(4.*beta+(-1).*(1+ ...
  delta).^2).^(-1).*(2.*exp(1).^(t0.*((-1)+beta+delta)+(1/2).*t.*((-1)+delta+ ...
  Psi)).*IdA0.*beta.*(beta+(-1).*delta).*(4.*beta+(-1).*(1+delta).^2)+(-2).*exp( ...
  1).^(t0.*((-1)+delta)+(1/2).*t.*((-1)+2.*beta+delta+Psi)).*alpha.*(beta.^2+((-1) ...
  +beta).*delta).*((-4).*beta+(1+delta).^2)+exp(1).^((1/2).*t0.*((-1)+delta+(-1) ...
  .*Psi)+t.*((-1)+beta+delta+Psi)).*((-1).*IdA0.*beta.*(beta+(-1).*delta).*((-1).* ...
  (1+delta).^2+((-1)+delta).*Psi+2.*beta.*(2+Psi))+(beta.^2+((-1)+beta).*delta).*( ...
  alpha.*(1+(-4).*beta+delta.*(2+delta+(-1).*Psi)+Psi)+(beta+(-1).*delta).*(4.*IaA0.* ...
  beta+(-2).*(1+V0).*Psi+IaA0.*(1+delta).*((-1)+(-1).*delta+Psi))))+exp(1).^( ...
  t.*((-1)+beta+delta)+(1/2).*t0.*((-1)+delta+Psi)).*(IdA0.*beta.*(beta+(-1).*delta) ...
  .*((1+delta).^2+2.*beta.*((-2)+Psi)+((-1)+delta).*Psi)+(beta.^2+((-1)+beta).*delta) ...
  .*(alpha.*((-4).*beta+(1+delta).^2+((-1)+delta).*Psi)+(beta+(-1).*delta).*(4.* ...
  IaA0.*beta+2.*(1+V0).*Psi+(-1).*IaA0.*(1+delta).*(1+delta+Psi)))));

Idep = @(t,alfa,beta,delta,IaA0,IdA0) exp(1).^(((-1).*t+t0).*beta).*IdA0;
%% Optimization
% Optimization uses a multi-objective strategy, minimizing an error that takes into account multiple features at
% the same time. So the found solution is a compromise between all the features, while also aiming at fulfilling the constraints
delta = (-1/tau_m(i)-k2)^2-4*(-k2/tau_m(i)+k_adap/Cm(i));     % [1/ms^2]
param3_low = 3/((1/(m_IF(i)+3*sd_IF(i)))*1000); % param3 = k1
param3_high = 3/t_ref(i);

% Global variables for saving optimization info
global par cf w con error_all
w = 1;          % The weight to consider whether Vm reaches the threshold or not
% low = K_adap, K2, A2, K1, A1, I_e
low = [Cm(i)/(tau_m(i)^2)+0.000001,-1/tau_m(i)+0.00001,0.001,3/((1/m_IF(i))*1000),0.001,-100.0];   %k_adap (sicuro >Cm/tau_m^2),k2 (sicuro > 1/tau_m),A2,k1,A1, Ie
up_2 = 10*low(2);
up = [((up_2-1/tau_m(i))^2)*Cm(i)/4-0.000001,up_2,500,3/t_ref(i),500,100];

% low = [Cm(i)/(tau_m(i)^2)+0.000001,-1/tau_m(i)+0.025,0.001,3/((1/m_IF(i))*1000),0.001,-150];   %k_adap (sicuro >Cm/tau_m^2),k2 (sicuro > 1/tau_m),A2,k1,A1, Ie
% up_2 = 5*low(2);
% up = [((up_2-1/tau_m(i))^2)*Cm(i)/4-0.000001,up_2,100,3/t_ref(i),600,-50];

% Linear inequality constraints: A2<A1; kadap>(Cm/tau_m)*k2 -> in normalized
% ranges!!
% Att: the first constraints should be modified if A1 and A2 are not in the
% same ranges (should be converted in the original range as done for Kadap
% and K2)
% A = [0 0 1 0 -1 0;(low(1)-up(1)) (-Cm(i)/tau_m(i))*(up(2)-low(2)) 0 0 0 0];
% b = [0;low(1)+(Cm(i)/tau_m(i))*low(2)-0.000001];
% A = [0 0 (up(3)-low(3)) 0 -(up(5)-low(5)) 0;(low(1)-up(1)) (-Cm(i)/tau_m(i))*(up(2)-low(2)) 0 0 0 0];
% b = [low(5)-low(3);low(1)+(Cm(i)/tau_m(i))*low(2)-0.000001];
A = [0 0 1 0 -1 0; (low(1)-up(1)) (-Cm(i)/tau_m(i))*(up(2)-low(2)) 0 0 0 0];
b = [0; low(1)+(Cm(i)/tau_m(i))*low(2)-0.000001];

% Linear equality constraints
Aeq = [];
beq = [];

% Lower and upper bounds of normalized parameters
lb = zeros(6,1);
ub = ones(6,1);

for nopt = 1:10
    par = [];
    cf = [];
    con = [];
    error_all = [];

    start_param = rand(6,1)'
    % Check that the constraints are satisfied at initial point
    con1 = (confun_eglif(start_param,low,up,i,Iinh,Cm,tau_m,E_L,Vth, Vinh_ss,t_ref,L2, L3, Sp, T_tonic, Istim,V1,V2,Vss, T_dep3, SFA_gain, Iahp1, Iahp2, Iahp3))
    while con1(1)>0 || con1(2)>0 %|| con1(3)>0  || con1(4)>0
        start_param = rand(6,1)'
        con1 = (confun_eglif(start_param,low,up,i,Iinh,Cm,tau_m,E_L,Vth, Vinh_ss,t_ref,L2, L3, Sp, T_tonic, Istim,V1,V2,Vss, T_dep3, SFA_gain, Iahp1, Iahp2, Iahp3))
    end

    equs = [equ1;equ2;equ3;equ4;equ5;equ6;equ7;equ8];

% Optimization algorithm options
    options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter-detailed','TolX',1e-3,'TolCon',1e-3,...
        'TolFun',1e-2,'ObjectiveLimit',0.1,'MaxFunEvals',100,'MaxIter',100);         %,'ScaleProblem','obj-and-constr');          % Test also different algorithms!!!!
    
% % Separating linear and non linear constraints - error on area Vm
    [param_all,fval_all] = ...
    fmincon(@(param)objfun_eglif(param,low,up,equs,[1.2 1/2 1/4 1/6 1.1 1 1.1 1.1 1 1 1 1 1],Istim,i,T_tonic, Tdep, Tburst, Tlb, t_ref, SFA_gain,tau_m, Cm,E_L,Vth,Vr,Sp,L2,L3,Iinh(i),Tahp,Vinh_ss),start_param,...
        A,b,Aeq,beq,lb,ub,@(param)confun_eglif(param,low,up,i,Iinh,Cm,tau_m,E_L,Vth, Vinh_ss,t_ref,L2, L3, Sp, T_tonic, Istim,V1,V2,Vss, T_dep3, SFA_gain, Iahp1, Iahp2, Iahp3),options);

    parametri{nopt} = par;
    
    cost_function{nopt} = cf;
    par_init(nopt,:) = start_param;
    constraints{nopt} = con;
    err_all{nopt} = error_all;
    nopt
    
    parameters = parametri{1,nopt};
    param_all = parameters(end,:);
    disp('valori parametri riportati ai range originali')
    optimal_parameters = param_all.*(up-low)+low

end


% Saving optimization data
save param_CA1PC.mat parametri
save cost_CA1PC.mat cost_function
save init_CA1PC.mat par_init
save constr_CA1PC.mat constraints
save err_all_CA1PC.mat err_all


%% Plot optim
figure(10)
hold on
nopt=10
for i=1:nopt
    
    subplot(1,2,1)
    subplot(1,2,2)
    plot(cost_function{1,i})
    legend
    
end
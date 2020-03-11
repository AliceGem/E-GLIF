% E-GLIF neuron model optimization
% LIF neuron + adaptive current for modelling subthreshold-related
% mechanisms and SFA + spike-triggered currents modelling io
% channel-related mechanisms; therefore, a simple linear model (useful for
% large scale simulations) and biologically plausible.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% k2 fissato a -1/tau_m; ottimizzati solo gli altri parametri %%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 clear
% clc

neu_names = { 'MLI', 'DCN', 'IO'};

% Neuron properties from neurophysiology: checked values from papers (for
% coherence on protocol) + neuroelectro.org
Cm = [14.6 142.0 197.0];
tau_m = [-9.125 -33.0 -20.0];   % [ms] should be the opposite of the given value - to fix
t_ref = [1.59 1.5 12.0];
E_L = [-68.0 -45.0 -55.0];    % [mV]
Vth = [-53.0 -36.0 -50.0];    % [mV]
Vr = [-78.0 -55.0 -65.0];    % EL-10;
nn = length(Cm);
fI_PC = 0.1;        % Masoli, 2015 - Fig.2

m_IF = [8.5 20.0 0.0];        % Mean intrinsic frequency - Lachamp, 2009
sd_IF = [2.7 0 3.0 0.0];     % SE o Standard deviation intrinsic frequency - Lachamp, 2009 (reports SE over n = 6)
            

    
% Protocolli di stimolazione dep: Istim = [3xnn] in [pA]; 3 fasi di stimolazione,
% nn neuroni
Istim = [[12.0; 24.0; 36.0]...     % MLI - Stellate example from Galliano et al., 2013 - Fig.S3
        [0.0; 0.0; 0.0]]; % IO - tonic bursting - da capire il protocollo di stimolazione
                                        
m_Fdep = [[30.0; 60.0; 90.0]...              %[[30.0; 70.0; 90.0]...         % f-I slope da Galliano  
          [45.0; 80.0; 115.0]...       % [Uusisaari et al., 2007]
          [0.0; 0.0; 0.0]];

sd_Fdep = 1*[[1.0; 5.0; 10.0]...   % La sd non è disponibile nei dati dal model
           [5.0; 7.0; 10.0]...
           [0.0; 0.0; 0.0]];
 

% Stimolazione hyp     
Iinh = [-24.0];       
Vinh_min = [-125];      
Vinh_ss = [-115];

SFA_gain = [1 1 1; ];

i = 1;      % i-th neuron


% Model - 3D linear ODE system
syms Vm(t) I1(t) I2(t)          % Variabili di stato
syms Ie k_adap k1 k2 A1 A2 ...            % Parametri da ottimizzare
     Ist


ode1(i) = diff(Vm) == (-1/tau_m(i))*Vm + Ie/Cm(i) + Ist/Cm(i) + I1/Cm(i) - I2/Cm(i) + E_L(i)/tau_m(i);
ode2(i) = diff(I1) == -k1*I1;
ode3(i) = diff(I2) == k_adap*Vm - k2*I2 - k_adap*E_L(i);



%%
I = Ie+Ist;  	% Ie resta da ottimizzare (param(6)) perch� le correnti sono attive anche nella fase tonic, in assenza di input esterno. 

% Eigenvalues
l1 = -k1;
D = (1/tau_m(i)+k2)^2-4*(k2/tau_m(i)+k_adap/Cm(i));   % Discriminante
l2 = 0.5*(-(1/tau_m(i)+k2)+sqrt(D));
l3 = 0.5*(-(1/tau_m(i)+k2)-sqrt(D));

% Eigenvectors
csi = (k2-k1)*tau_m(i)/((1-k1*tau_m(i))*(k2-k1)*Cm(i)+k_adap*tau_m(i));
csi2 = k_adap*tau_m(i)/((1-k1*tau_m(i))*(k2-k1)*Cm(i)+k_adap*tau_m(i));



x1 = [csi; csi2; 1];      % Associated to l1

L2 = Cm(i)*(-1/tau_m(i)-l2);
x2 = [1; L2; 0];            % Associated to l2

L3 = Cm(i)*(-1/tau_m(i)-l3);
x3 = [1; L3; 0];            % Associated to l3

% Specific solution
Sp = [(tau_m(i)*k2*I + E_L(i)*(tau_m(i)*k_adap+Cm(i)*k2))/(tau_m(i)*k_adap + Cm(i)*k2);...      % Vm
      I*k_adap*tau_m(i)/(tau_m(i)*k_adap + Cm(i)*k2);...                                        % I2
      0];                                                                                       % I1

%% Optimization

% clc
% clear param_all

ne = 10;

T_tonic(i,:) = (1./(sd_IF(i).*randn(1,ne) + m_IF(i)))*1000;

% Avoid negative values in the distribution
while ~isempty(find(T_tonic<0))
    T_tonic(i,find(T_tonic<0))=(1./(sd_IF(i).*randn(1,length(find(T_tonic<0))) + m_IF(i)))*1000;
end

T_dep1(i,:) = (1./(sd_Fdep(1,i).*randn(1,ne) + m_Fdep(1,i)))*1000;
T_dep2(i,:) = (1./(sd_Fdep(2,i).*randn(1,ne) + m_Fdep(2,i)))*1000;
T_dep3(i,:) = (1./(sd_Fdep(3,i).*randn(1,ne) + m_Fdep(3,i)))*1000;
 
Tdep = {T_dep1,T_dep2,T_dep3};

m_Fburst = [m_IF(i)];           % [Hz]     % No rebound burst
sd_Fburst = [sd_IF(i)];        %[1];
Tburst = (1./(sd_Fburst(i).*randn(1,ne) + m_Fburst(i)))*1000;
Tlb = 5.*randn(1,ne) + 1000*(1/mean(m_IF(i)));  

Tahp = [5.*randn(ne,1) + 80, 5.*randn(ne,1) + 100, 5.*randn(ne,1) + 120];  %[80 100 120];
%la soluzione seguente ipotizza che i parametri non siano tutti contemporaneamente
% nulli e i vincoli tengono conto di: soluzione IMMAGINARIA e STABILE


% c1, c2, c3 are vectors containing the constants' values during the
% different phases of optimization: latency, 1st spk, SS spk (tonic and exc1,2,3), lb, 1st burst ISI 

% Phase 1: time to first spike (latency)
syms c1
c1(1) = 0;
c2(1) = ((E_L(i)-Sp(1))*L3+Sp(2))/(L3-L2);
c3(1) = ((Sp(1)-E_L(i))*L2-Sp(2))/(L3-L2);

V1 = c1(1)*x1(1)*exp(l1*t)+c2(1)*x2(1)*exp(l2*t)+c3(1)*x3(1)*exp(l3*t)+Sp(1);
I2_lat = c1(1)*x1(2)*exp(l1*t)+c2(1)*x2(2)*exp(l2*t)+c3(1)*x3(2)*exp(l3*t)+Sp(2);

% Stessa soluzione di Hertag perch� le c.i. portano Idep ad essere nulla
% fino al primo spike

equ1 = V1 - Vth(i);

% Phase 2: time of the second spike (onset)
syms t1
dV1 = ((c1(1))*(x1(1))*(l1))*exp(l1*t1)+((c2(1))*(x2(1))*(l2))*exp(l2*t1)+((c3(1))*(x3(1))*(l3))*exp(l3*t1);

Vr_ = Vr(i)-Sp(1)-A1*csi;
beta1 = dV1-(Vr(i)-Vth(i))/tau_m(i)-A2/Cm(i)+A1/Cm(i);

c1(2) = A1;
c2(2) = (l3*Vr_-beta1+A1*l1*csi)/(l3-l2);
c3(2) = (beta1-A1*l1*csi-l2*Vr_)/(l3-l2);

V2 = ((c1(2))*(x1(1)))*exp(l1*t)+((c2(2))*(x2(1)))*exp(l2*t)+((c3(2))*(x3(1)))*exp(l3*t)+Sp(1);
I2_1spk = c1(2)*x1(2)*exp(l1*t)+c2(2)*x2(2)*exp(l2*t)+c3(2)*x3(2)*exp(l3*t)+Sp(2);

equ2 = V2-Vth(i);

% Phase 3: steady-state spiking time
syms tss

mi = x1(2)*exp(l1*tss);
eta = x2(2)*exp(l2*tss);
teta = x3(2)*exp(l3*tss);

c1(3) = A1;
c2(3) = (Vr_*L3-Vr_*teta-A1*mi-A2+A1*csi2)/(L3-L2+eta-teta);
c3(3) = (A1*mi+Vr_*(eta-L2)+A2-A1*csi2)/(L3-L2+eta-teta);

Vss = c1(3)*x1(1)*exp(l1*t)+c2(3)*x2(1)*exp(l2*t)+c3(3)*x3(1)*exp(l3*t)+Sp(1);
I2_ss = c1(3)*x1(2)*exp(l1*t)+c2(3)*x2(2)*exp(l2*t)+c3(3)*x3(2)*exp(l3*t)+Sp(2);

equ3 = Vss-Vth(i);


% Phase 4A: latency of rebound burst
Vhyp_ss = subs(Sp(1),Ist,Iinh(i));
Ihyp_ss = subs(Sp(2),Ist,Iinh(i));

phi = -Vhyp_ss/tau_m(i)-Ihyp_ss/Cm(i)+E_L(i)/tau_m(i)+Ie/Cm(i);
c1(4) = 0;
c2(4) = ((Vhyp_ss-Sp(1))*l3-phi)/(l3-l2);
c3(4) = (-(Vhyp_ss-Sp(1))*l2+phi)/(l3-l2);
Vlb = c1(4)*x1(1)*exp(l1*t)+c2(4)*x2(1)*exp(l2*t)+c3(4)*x3(1)*exp(l3*t)+Sp(1);
Ihyp_lb = c1(4)*x1(2)*exp(l1*t)+c2(4)*x2(2)*exp(l2*t)+c3(4)*x3(2)*exp(l3*t)+Sp(2);

equ4 = Vlb-Vth(i);

% Phase 4: rebound burst
syms tlb
% tlb = 5;           % IPOTESI!!! Meglio metterlo nell'ottimizzazione
dVlb = ((c1(4))*(x1(1))*(l1))*exp(l1*tlb)+((c2(4))*(x2(1))*(l2))*exp(l2*tlb)+((c3(4))*(x3(1))*(l3))*exp(l3*tlb);
beta2 = dVlb-(Vr(i)-Vth(i))/tau_m(i)-A2/Cm(i)+A1/Cm(i); % Con A2 sottratto dal valore raggiunto al tempo dello spike
            
c1(5) = A1;
c2(5) = (l3*Vr_-beta2+A1*l1*csi)/(l3-l2);
c3(5) = (beta2-A1*l1*csi-l2*Vr_)/(l3-l2);

Vb = c1(5)*x1(1)*exp(l1*t)+c2(5)*x2(1)*exp(l2*t)+c3(5)*x3(1)*exp(l3*t)+Sp(1);
equ5 = Vb-Vth(i);

% Phase 5 AHP1
% Vdep1_ss = subs(Sp(1),Ist,Istim(1,i));
% Idep1_ss = subs(Sp(2),Ist,Istim(1,i));
    %Vdep1_ss = subs(Sp(1),Ist,0);
Vdep1_ss = E_L(i);              % Punto di partenza: un valore di potenziale tra Vr e Vth perch� il neurone sta sparando
Idep1_ss = subs(Sp(2),Ist,Istim(1,i));              % TO BE CHECKED

phi = -Vdep1_ss/tau_m(i)-Idep1_ss/Cm(i)+E_L(i)/tau_m(i)+Ie/Cm(i);
c1(6) = 0;
c2(6) = ((Vdep1_ss-Sp(1))*l3-phi)/(l3-l2);
c3(6) = (-(Vdep1_ss-Sp(1))*l2+phi)/(l3-l2);
Vahp1 = c1(6)*x1(1)*exp(l1*t)+c2(6)*x2(1)*exp(l2*t)+c3(6)*x3(1)*exp(l3*t)+Sp(1);
Iahp1 = c1(6)*x1(2)*exp(l1*t)+c2(6)*x2(2)*exp(l2*t)+c3(6)*x3(2)*exp(l3*t)+Sp(2);

equ6 = Vahp1-Vth(i);

% Phase 6 AHP2
% Vdep2_ss = subs(Sp(1),Ist,Istim(2,i));
% Idep2_ss = subs(Sp(2),Ist,Istim(2,i));
% Vdep2_ss = subs(Sp(1),Ist,0);
Vdep2_ss = E_L(i);
Idep2_ss = subs(Sp(2),Ist,Istim(2,i));

phi = -Vdep2_ss/tau_m(i)-Idep2_ss/Cm(i)+E_L(i)/tau_m(i)+Ie/Cm(i);
c1(7) = 0;
c2(7) = ((Vdep2_ss-Sp(1))*l3-phi)/(l3-l2);
c3(7) = (-(Vdep2_ss-Sp(1))*l2+phi)/(l3-l2);
Vahp2 = c1(7)*x1(1)*exp(l1*t)+c2(7)*x2(1)*exp(l2*t)+c3(7)*x3(1)*exp(l3*t)+Sp(1);
Iahp2 = c1(7)*x1(2)*exp(l1*t)+c2(7)*x2(2)*exp(l2*t)+c3(7)*x3(2)*exp(l3*t)+Sp(2);

equ7 = Vahp2-Vth(i);

% Phase 7 AHP3
% Vdep3_ss = subs(Sp(1),Ist,Istim(3,i));
% Idep3_ss = subs(Sp(2),Ist,Istim(3,i));
% Vdep3_ss = subs(Sp(1),Ist,0);
Vdep3_ss = E_L(i);
Idep3_ss = subs(Sp(2),Ist,Istim(3,i));

phi = -Vdep3_ss/tau_m(i)-Idep3_ss/Cm(i)+E_L(i)/tau_m(i)+Ie/Cm(i);
c1(8) = 0;
c2(8) = ((Vdep3_ss-Sp(1))*l3-phi)/(l3-l2);
c3(8) = (-(Vdep3_ss-Sp(1))*l2+phi)/(l3-l2);
Vahp3 = c1(8)*x1(1)*exp(l1*t)+c2(8)*x2(1)*exp(l2*t)+c3(8)*x3(1)*exp(l3*t)+Sp(1);
Iahp3 = c1(8)*x1(2)*exp(l1*t)+c2(8)*x2(2)*exp(l2*t)+c3(8)*x3(2)*exp(l3*t)+Sp(2);

equ8 = Vahp3-Vth(i);
%% Plot delle 3 soluzioni per check
% close all
% param_all = [0.01, 0.03, 0.35*Cm(i), 0.05, 1.5*Cm(i)];
% T_tonic(i,1)=8;
% Tonic
Vm_tonic_lat = subs(V1,{Ist,k_adap,k2,A2,k1,A1,Ie},{0,param_all(1),param_all(2),param_all(3),param_all(4),param_all(5),param_all(6)})
Vm_tonic_1spk = subs(V2,{Ist,k_adap,k2,A2,k1,A1,Ie,t1},{0,param_all(1),param_all(2),param_all(3),param_all(4),param_all(5),param_all(6),T_tonic(i,1)});
Vm_tonic_ss = subs(Vss,{Ist,k_adap,k2,A2,k1,A1,Ie,tss},{0,param_all(1),param_all(2),param_all(3),param_all(4),param_all(5),param_all(6),T_tonic(i,1)});
I2_tonic_lat = subs(I2_lat,{Ist,k_adap,k2,A2,k1,A1,Ie},{0,param_all(1),param_all(2),param_all(3),param_all(4),param_all(5),param_all(6)});
I2_tonic_1spk = subs(I2_1spk,{Ist,k_adap,k2,A2,k1,A1,Ie,t1},{0,param_all(1),param_all(2),param_all(3),param_all(4),param_all(5),param_all(6),T_tonic(i,1)});
I2_tonic_ss = subs(I2_ss,{Ist,k_adap,k2,A2,k1,A1,Ie,tss},{0,param_all(1),param_all(2),param_all(3),param_all(4),param_all(5),param_all(6),T_tonic(i,1)});


figure
Xax = [0 1000]
hton = ezplot(Vm_tonic_lat,Xax);
set(hton,'Color','r','LineStyle','-.')
hold on
line(Xax,[Vth(i) Vth(i)],'LineStyle','--');
h1 = ezplot((Vm_tonic_1spk),Xax);
set(h1,'Color','b','LineStyle','-.')
hss = ezplot(Vm_tonic_ss,Xax);
set(hss,'Color','m','LineStyle','-.')
ylim([-80 120])
title(['param = ', num2str(param_all)])
%%
figure
Xax = [0 1000]
hton = ezplot(I2_tonic_lat,Xax);
set(hton,'Color','r','LineStyle','-.')
hold on
line(Xax,[Vth(i) Vth(i)],'LineStyle','--');
h1 = ezplot((I2_tonic_1spk),Xax);
set(h1,'Color','b','LineStyle','-.')
hss = ezplot(I2_tonic_ss,Xax);
set(hss,'Color','m','LineStyle','-.')
ylim([-80 120])
title(['Iadap - param = ', num2str(param_all)])
%%
V_1(1) = subs(V1,{Ist,k_adap,k2,A2,k1,A1,Ie},{Istim(1,i),param_all(1),param_all(2),param_all(3),param_all(4),param_all(5),param_all(6)});
V_1(2) = subs(V1,{Ist,k_adap,k2,A2,k1,A1,Ie},{Istim(2,i),param_all(1),param_all(2),param_all(3),param_all(4),param_all(5),param_all(6)});
V_1(3)= subs(V1,{Ist,k_adap,k2,A2,k1,A1,Ie},{Istim(3,i),param_all(1),param_all(2),param_all(3),param_all(4),param_all(5),param_all(6)});
V_1(4)= subs(V1,{Ist,k_adap,k2,A2,k1,A1,Ie},{-Istim(1,i),param_all(1),param_all(2),param_all(3),param_all(4),param_all(5),param_all(6)});

% I2_dep_lat(1) = subs(I2_lat,{Ist,k_adap,k2,A2,k1,A1,Ie},{Istim(1,i),param_all(1),-1/tau_m(i),param_all(2),param_all(3),param_all(4),param_all(5)});
% I2_dep_1spk(1) = subs(I2_1spk,{Ist,k_adap,k2,A2,k1,A1,Ie,t1},{Istim(1,i),param_all(1),-1/tau_m(i),param_all(2),param_all(3),param_all(4),param_all(5),T_tonic(i,1)});
% I2_dep_ss(1) = subs(I2_ss,{Ist,k_adap,k2,A2,k1,A1,Ie,tss},{Istim(1,i),param_all(1),-1/tau_m(i),param_all(2),param_all(3),param_all(4),param_all(5),T_tonic(i,1)});

figure;

Xax = [0 1000]
hlat = ezplot(V_1(1),Xax);
set(hlat,'Color','r')
hold on
line(Xax,[Vth(i) Vth(i)],'LineStyle','--');
h = ezplot(V_1(2),Xax);
set(h,'Color','r','LineStyle','--')
h = ezplot(V_1(3),Xax);
set(h,'Color','r','LineStyle',':')
h = ezplot(V_1(4),Xax);
set(h,'Color','r','LineStyle','-.')
ylim([-80 120])

V_2(1) = subs(V2,{Ist,k_adap,k2,A2,k1,A1,t1,Ie},{Istim(1,i),param_all(1),param_all(2),param_all(3),param_all(4),param_all(5),T_dep1(i,1),param_all(6)});
V_2(2) = subs(V2,{Ist,k_adap,k2,A2,k1,A1,t1,Ie},{Istim(2,i),param_all(1),param_all(2),param_all(3),param_all(4),param_all(5),T_dep2(i,1),param_all(6)});
V_2(3)= subs(V2,{Ist,k_adap,k2,A2,k1,A1,t1,Ie},{Istim(3,i),param_all(1),param_all(2),param_all(3),param_all(4),param_all(5),T_dep3(i,1),param_all(6)});
V_2(4)= subs(V2,{Ist,k_adap,k2,A2,k1,A1,t1,Ie},{-Istim(1,i),param_all(1),param_all(2),param_all(3),param_all(4),param_all(5),125,param_all(6)});


%figure;
hon1 = ezplot(V_2(1),Xax);
set(hon1,'Color','b')
hold on
h = ezplot(V_2(2),Xax);
set(h,'Color','b','LineStyle','--')
h = ezplot(V_2(3),Xax);
set(h,'Color','b','LineStyle',':')
h = ezplot(V_2(4),Xax);
set(h,'Color','b','LineStyle','-.')


V_3(1) = subs(Vss,{Ist,k_adap,k2,A2,k1,A1,tss,Ie},{Istim(1,i),param_all(1),param_all(2),param_all(3),param_all(4),param_all(5),T_dep1(i,1)*1.1,param_all(6)});
V_3(2) = subs(Vss,{Ist,k_adap,k2,A2,k1,A1,tss,Ie},{Istim(2,i),param_all(1),param_all(2),param_all(3),param_all(4),param_all(5),T_dep2(i,1)*1.1,param_all(6)});
V_3(3)= subs(Vss,{Ist,k_adap,k2,A2,k1,A1,tss,Ie},{Istim(3,i),param_all(1),param_all(2),param_all(3),param_all(4),param_all(5),T_dep3(i,1)*1.1,param_all(6)});
V_3(4) = subs(Vss,{Ist,k_adap,k2,A2,k1,A1,tss,Ie},{-Istim(1,i),param_all(1),param_all(2),param_all(3),param_all(4),param_all(5),T_dep1(i,1)*1.1,param_all(6)});

hss = ezplot(V_3(1),Xax);
set(hss,'Color','m')
hold on
h = ezplot(V_3(2),Xax);
set(h,'Color','m','LineStyle','--')
h = ezplot(V_3(3),Xax);
set(h,'Color','m','LineStyle',':')
h = ezplot(V_3(4),Xax);
set(h,'Color','m','LineStyle','-.')

% legend([hlat,hon1,hss],{'lat','on','ss'},'Location', 'northeastoutside')

V_4 = subs(Vlb,{Ist,k_adap,k2,A2,k1,A1,Ie},{0,param_all(1),param_all(2),param_all(3),param_all(4),param_all(5),param_all(6)});

h_lb = ezplot(V_4,Xax);
set(h_lb,'Color','g')
hold on

ylim([-80 120])

V_5 = subs(Vb,{Ist,k_adap,k2,A2,k1,A1,tlb,Ie},{0,param_all(1),param_all(2),param_all(3),param_all(4),param_all(5),Tlb(1),param_all(6)});
h_b = ezplot(V_5,Xax);
set(h_b,'Color','c')
ylim([-120 120])

V_6 = subs(Vahp1,{Ist,k_adap,k2,A2,k1,A1,Ie},{0,param_all(1),param_all(2),param_all(3),param_all(4),param_all(5),param_all(6)});
h_ahp1 = ezplot(V_6,Xax);
set(h_ahp1,'Color','k')
ylim([-120 120])

V_7 = subs(Vahp2,{Ist,k_adap,k2,A2,k1,A1,Ie},{0,param_all(1),param_all(2),param_all(3),param_all(4),param_all(5),param_all(6)});
h_ahp2 = ezplot(V_7,Xax);
set(h_ahp2,'Color','k')
ylim([-120 120])

V_8 = subs(Vahp3,{Ist,k_adap,k2,A2,k1,A1,Ie},{0,param_all(1),param_all(2),param_all(3),param_all(4),param_all(5),param_all(6)});
h_ahp3 = ezplot(V_8,Xax);
set(h_ahp3,'Color','k')
ylim([-200 200])
title(['param = ',num2str(param_all)])

figure;
ezplot(subs(Ihyp_lb,{Ist,k_adap,k2,A2,k1,A1},{0,param_all(1),param_all(2),param_all(3),param_all(4),param_all(5)}),Xax);

ylim([-120 120])

figure;
ezplot(subs(I2_lat,{Ist,k_adap,k2,A2,k1,A1,Ie},{Istim(3,i),param_all(1),param_all(2),param_all(3),param_all(4),param_all(5)}),Xax);
hold on
ezplot(subs(I2_1spk,{Ist,k_adap,k2,A2,k1,A1,t1,Ie},{Istim(3,i),param_all(1),param_all(2),param_all(3),param_all(4),param_all(5),T_tonic(i,1),param_all(6)}),Xax);
ezplot(subs(I2_1spk,{Ist,k_adap,k2,A2,k1,A1,tss,Ie},{Istim(3,i),param_all(1),param_all(2),param_all(3),param_all(4),param_all(5),T_tonic(i,1),param_all(6)}),Xax);



%%
% Att! Quello che otteniamo è il meglio per quei constraints!
delta = (-1/tau_m(i)-k2)^2-4*(-k2/tau_m(i)+k_adap/Cm(i));     % [1/ms^2]
param3_low = 3/((1/(m_IF(i)+3*sd_IF(i)))*1000);       
param3_high = 3/t_ref(i);

global par cf w con error_all
w = 1;          % The weight to consider whether Vm reaches the threshold or not

low = [Cm(i)/(tau_m(i)^2)+0.000001,-1/tau_m(i)+0.000001,0.0001,3/((1/m_IF(i))*1000),0.0001,0.0001];   %k_adap (sicuro >Cm/tau_m^2),k2 (sicuro > 1/tau_m),A2,k1,A1, Ie
up_2 = 10*low(2);
up = [((up_2-1/tau_m(i))^2)*Cm(i)/4-0.000001,up_2,10.0,3/t_ref(i),10.0,10.0];

A = [0 0 1 0 -1 0;(low(1)-up(1)) (-Cm(i)/tau_m(i))*(up(2)-low(2)) 0 0 0 0];       
        % Linear constraints: A2<A1; kadap>(Cm/tau_m)*k2 -> in normalized
        % ranges!! Att: il primo vincolo va modificato se A1 e A2 non hanno
        % gli stessi ranges!
b = [0;low(1)+(Cm(i)/tau_m(i))*low(2)-0.000001];

Aeq = [];
beq = [];

lb = zeros(6,1);
ub = ones(6,1);

for nopt = 1:10
    par = [];
    cf = [];
    con = [];
    error_all = [];
    
    start_param = rand(6,1)'
    %Check that the constraints are satisfied at initial point
    con1 = (confun_eglif(start_param,low,up,i,Iinh,Cm,tau_m,E_L,Vth, Vinh_ss,t_ref,L2, L3, Sp, T_tonic, Istim,V1,V2,Vss, T_dep3, SFA_gain, Iahp1, Iahp2, Iahp3))
    while con1(1)>0 || con1(2)>0 || con1(3)>0  || con1(4)>0
        start_param = rand(6,1)'
        con1 = (confun_eglif(start_param,low,up,i,Iinh,Cm,tau_m,E_L,Vth, Vinh_ss,t_ref,L2, L3, Sp, T_tonic, Istim,V1,V2,Vss, T_dep3, SFA_gain, Iahp1, Iahp2, Iahp3))
    end
    
    equs = [equ1;equ2;equ3;equ4;equ5;equ6;equ7;equ8];

    options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter-detailed','TolX',1e-3,'TolCon',1e-3,...
        'TolFun',1e-3,'ObjectiveLimit',0.1,'MaxFunEvals',200,'MaxIter',200);         %,'ScaleProblem','obj-and-constr');          % Test also different algorithms!!!!


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
end

save param.mat parametri
save cost.mat cost_function
save init_par.mat par_init
save constr.mat constraints
save err_all.mat err_all
%%
clear parametri cost_function par_init constraints


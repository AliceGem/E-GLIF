function [ sq_err ] = objfun_eglif( param, low, up, equs, weights, Istim, i, Tton, T, Tb, Tlb, tref, gainSFA, tau_m, Cm, E_L, Vth, Vr, Sp, L2, L3, Iinh,T_ahp, Vhyp_ss,n_optim )

% objfun_eglif compute the sum of squared errors over the depolarizing phases for
% latency, onset and steady-state firing period:
% error_'phase'N = error during 'phase' N (e.g. depol1)
% and rebound phase.
% For some neurons also AHP could be considered

global par cf error_all


syms k_adap k2 Ist A2 k1 A1 Ie

par = [par;param];
neu_ind = i;
fac = 1.2;    % To consider the interval including the desired spk time! Therefore we can not search the roots just up to the desired spk time
step = 0.1;   % Factor to compute the step in the t_vec of points where we search for zero crossing

delta = (-1/tau_m(neu_ind)-norma(param(2),low(2),up(2)))^2-4*(norma(param(2),low(2),up(2))/tau_m(neu_ind)+norma(param(1),low(1),up(1))/Cm(neu_ind));     % [1/ms^2]


c1 = 0;
c2 = subs((((E_L(i)-Sp(1))*L3+Sp(2))/(L3-L2)),{Ist,k_adap,k2,A2,k1,A1,Ie},{0,norma(param(1),low(1),up(1)),norma(param(2),low(2),up(2)),norma(param(3),low(3),up(3)),...
    norma(param(4),low(4),up(4)),norma(param(5),low(5),up(5)),norma(param(6),low(6),up(6))});
c3 = subs((((Sp(1)-E_L(i))*L2-Sp(2))/(L3-L2)),{Ist,k_adap,k2,A2,k1,A1,Ie},{0,norma(param(1),low(1),up(1)),norma(param(2),low(2),up(2)),norma(param(3),low(3),up(3)),...
    norma(param(4),low(4),up(4)),norma(param(5),low(5),up(5)),norma(param(6),low(6),up(6))});
% Eigenvalues
l1 = -norma(param(4),low(4),up(4));
D = delta;   % Discriminante
l2 = 0.5*(-(1/tau_m(neu_ind)+norma(param(2),low(2),up(2)))+sqrt(D));
l3 = 0.5*(-(1/tau_m(neu_ind)+norma(param(2),low(2),up(2)))-sqrt(D));
% Eigenvectors
csi = (norma(param(2),low(2),up(2))-norma(param(4),low(4),up(4)))*tau_m(neu_ind)/((1-norma(param(4),low(4),up(4))*tau_m(neu_ind))*(norma(param(2),low(2),up(2))-norma(param(4),low(4),up(4)))*Cm(neu_ind)+...
    norma(param(1),low(1),up(1))*tau_m(neu_ind));
csi2 = norma(param(1),low(1),up(1))*tau_m(neu_ind)/((1-norma(param(4),low(4),up(4))*tau_m(neu_ind))*(norma(param(2),low(2),up(2))-norma(param(4),low(4),up(4)))*Cm(neu_ind)+...
    norma(param(1),low(1),up(1))*tau_m(neu_ind));

v1 = [csi; csi2; 1];      % Associated to l1
x1 = v1(1);
L2 = Cm(neu_ind)*(-1/tau_m(neu_ind)-l2);
v2 = [1; L2; 0];            % Associated to l2
x2 = v2(1);
L3 = Cm(neu_ind)*(-1/tau_m(neu_ind)-l3);
v3 = [1; L3; 0];            % Associated to l3
x3 = v3(1);
Sp_1 = double(subs(Sp(1),{Ist,k_adap,k2,A2,k1,A1,Ie},{0,norma(param(1),low(1),up(1)),norma(param(2),low(2),up(2)),norma(param(3),low(3),up(3)),...
    norma(param(4),low(4),up(4)),norma(param(5),low(5),up(5)),norma(param(6),low(6),up(6))}));
        
% Depolarization steps
for dp = 1:3  % Number of depolarizing phases
    clear Sp_1
    Sp_1 = double(subs(Sp(1),{Ist,k_adap,k2,A2,k1,A1,Ie},{Istim(dp,i),norma(param(1),low(1),up(1)),norma(param(2),low(2),up(2)),norma(param(3),low(3),up(3)),...
        norma(param(4),low(4),up(4)),norma(param(5),low(5),up(5)),norma(param(6),low(6),up(6))}));

    for ind = 1:length(T{1}(i,:))

        % Phase 1
        sol1 = @(t,c1,c2,c3,x1,x2,x3,l1,l2,l3,Sp_1) Vth(i)-vpa(real(((c1*x1*exp(l1*t)+c2*x2*exp(l2*t)+c3*x3*exp(l3*t)+Sp_1))),2);

        c1 = 0;
        c2 = subs((((E_L(i)-Sp(1))*L3+Sp(2))/(L3-L2)),{Ist,k_adap,k2,A2,k1,A1,Ie},{Istim(dp,i),norma(param(1),low(1),up(1)),norma(param(2),low(2),up(2)),...
            norma(param(3),low(3),up(3)),...
            norma(param(4),low(4),up(4)),norma(param(5),low(5),up(5)),norma(param(6),low(6),up(6))});
        c3 = subs((((Sp(1)-E_L(i))*L2-Sp(2))/(L3-L2)),{Ist,k_adap,k2,A2,k1,A1,Ie},{Istim(dp,i),norma(param(1),low(1),up(1)),norma(param(2),low(2),up(2)),norma(param(3),low(3),up(3)),...
            norma(param(4),low(4),up(4)),norma(param(5),low(5),up(5)),norma(param(6),low(6),up(6))});

        Sp_1 = double(subs(Sp(1),{Ist,k_adap,k2,A2,k1,A1,Ie},{Istim(dp,i),norma(param(1),low(1),up(1)),norma(param(2),low(2),up(2)),norma(param(3),low(3),up(3)),...
            norma(param(4),low(4),up(4)),norma(param(5),low(5),up(5)),norma(param(6),low(6),up(6))}));


        area_act_dep(ind,1,dp) = integral(@(t) 0.5*(double(sol1(t,c1,c2,c3,x1,x2,x3,l1,l2,l3,Sp_1))).*(sign(double(sol1(t,c1,c2,c3,x1,x2,x3,l1,l2,l3,Sp_1)))+1),0.0,T{dp}(i,ind));

        error_dep(ind,1,dp) = abs(area_act_dep(ind,1,dp)-0.5*T{dp}(i,ind)*(Vth(i)-E_L(i)));
        %error_dep(ind,1,dp) = (abs(area_act_dep(ind,1,dp)-0.5*T{dp}(i,ind)*(Vth(i)-E_L(i))))/(0.5*T{dp}(i,ind)*(Vth(i)-E_L(i)));
        area_des1(ind,1) = (0.5*T{dp}(i,ind)*(Vth(i)-E_L(i)));
        area_real1(ind,1) = area_act_dep(ind,1,dp);

 		% Phase 2
        t1 = T{dp}(i,ind);
        dV1 = double(((c1)*(x1)*(l1))*exp(l1*t1)+((c2)*(x2)*(l2))*exp(l2*t1)+((c3)*(x3)*(l3))*exp(l3*t1));
        Vr_ = Vr(i)-Sp_1-norma(param(5),low(5),up(5))*csi;
        beta1 = dV1-(Vr(i)-Vth(i))/tau_m(i)-norma(param(3),low(3),up(3))/Cm(i)+norma(param(5),low(5),up(5))/Cm(i);

        clear c1 c2 c3
        
        sol2 = @(t,c1,c2,c3) Vth(i)-vpa(real(((c1*x1*exp(l1*t)+c2*x2*exp(l2*t)+c3*x3*exp(l3*t)+Sp_1))),2);

        c1 = norma(param(5),low(5),up(5));
        c2 = (l3*Vr_-beta1+norma(param(5),low(5),up(5))*l1*csi)/(l3-l2);
        c3 = (beta1-norma(param(5),low(5),up(5))*l1*csi-l2*Vr_)/(l3-l2);


        area_act_dep(ind,2,dp) = integral(@(t) 0.5*(double(sol2(t,c1,c2,c3))).*(sign(double(sol2(t,c1,c2,c3)))+1),0,T{dp}(i,ind)-tref(i));

        error_dep(ind,2,dp) = abs(area_act_dep(ind,2,dp)-0.5*(T{dp}(i,ind)-tref(i))*(Vth(i)-Vr(i)));
        %error_dep(ind,2,dp) = (abs(area_act_dep(ind,2,dp)-0.5*(T{dp}(i,ind)-tref(i))*(Vth(i)-Vr(i))))/(0.5*(T{dp}(i,ind)-tref(i))*(Vth(i)-Vr(i)));
        area_des1(ind,2) = (0.5*(T{dp}(i,ind)-tref(i))*(Vth(i)-Vr(i)));
        area_real1(ind,2) = area_act_dep(ind,2,dp);

        clear c1 c2 c3

 		% Phase 3
        tss = T{dp}(i,ind)*gainSFA(i,dp)-tref(i);
        sol3 = @(t,c1,c2,c3) Vth(i)-vpa(real(((c1*x1*exp(l1*t)+c2*x2*exp(l2*t)+c3*x3*exp(l3*t)+Sp_1))),2);

        mi = csi2*exp(l1*tss);
        eta = L2*exp(l2*tss);
        teta = L3*exp(l3*tss);
        
        c1 = norma(param(5),low(5),up(5));
        c2 = (Vr_*L3-Vr_*teta-norma(param(5),low(5),up(5))*mi-norma(param(3),low(3),up(3))+norma(param(5),low(5),up(5))*csi2)/(L3-L2+eta-teta);
        c3 = (norma(param(5),low(5),up(5))*mi+Vr_*(eta-L2)+norma(param(3),low(3),up(3))-norma(param(5),low(5),up(5))*csi2)/(L3-L2+eta-teta);
        
        area_act_dep(ind,3,dp) = integral(@(t) 0.5*(double(sol3(t,c1,c2,c3))).*(sign(double(sol3(t,c1,c2,c3)))+1),0.0,T{dp}(i,ind)*gainSFA(i,dp)-tref(i));
        
        error_dep(ind,3,dp) = abs(area_act_dep(ind,3,dp)-0.5*(T{dp}(i,ind)*gainSFA(i,dp)-tref(i))*(Vth(i)-Vr(i)));
        %error_dep(ind,3,dp) = (abs(area_act_dep(ind,3,dp)-0.5*(T{dp}(i,ind)*gainSFA(i,dp)-tref(i))*(Vth(i)-Vr(i))))/(0.5*(T{dp}(i,ind)*gainSFA(i,dp)-tref(i))*(Vth(i)-Vr(i)));
        area_des1(ind,3) = (0.5*(T{dp}(i,ind)*gainSFA(i,dp)-tref(i))*(Vth(i)-Vr(i)));
        area_real1(ind,3) = area_act_dep(ind,3,dp);
    end
end

% 2:3 per non considerare la latency del 1° spike

for p = 2:3
    error_depol1(p) = mean(error_dep(:,p,1).^2);
end

for p = 2:3
    error_depol2(p) = mean(error_dep(:,p,2).^2);
end

for p = 2:3
    error_depol3(p) = mean(error_dep(:,p,3).^2);
end

   sq_err = double(sqrt(mean([error_depol1 error_depol2 error_depol3])));

   error_all = [error_all; error_depol1 error_depol2 error_depol3];
   cf = [cf;sq_err];
end

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
    step = 0.1;     % Factor to compute the step in the t_vec of points where we search for zero crossing

    delta = (-1/tau_m(neu_ind)-norma(param(2),low(2),up(2)))^2-4*(norma(param(2),low(2),up(2))/tau_m(neu_ind)+norma(param(1),low(1),up(1))/Cm(neu_ind));     % [1/ms^2]


  %  area_Vhyp_ss = 0.5*base*(Vth(i)-Vhyp_ss);

	% Autorhythm step
	for ind = 1:length(Tton(i,:))

		% Phase 1
        sol1 = @(t,c1,c2,c3,x1,x2,x3,l1,l2,l3,Sp_1) Vth(i)-vpa(real(((c1*x1*exp(l1*t)+c2*x2*exp(l2*t)+c3*x3*exp(l3*t)+Sp_1))),2);

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
        % delta_Vm = @(t,c1,c2,c3,x1,x2,x3,l1,l2,l3,Sp_1) (c1*x1*exp(l1*t)+c2*x2*exp(l2*t)+c3*x3*exp(l3*t)+Sp_1)-Vth(i);
% 		delta_Vm = @(t) fun_subs(equs(1),[0,param(1),param(2),param(3),param(4),param(5),param(6)],c1);
%		t_act_ton(ind,1) = fzero(delta_Vm,Tton(i,ind));
        %options = optimset('PlotFcns',{@optimplotx,@optimplotfval});
        %options = optimset('PlotFcns',{@optimplotx,@optimplotfval},'Display','iter');

%        figure; syms tp; ezplot(sol1(tp,c1,c2,c3,x1,x2,x3,l1,l2,l3,Sp_1),0,1000)
% 		grid on
        area_act_ton(ind,1) = integral(@(t) 0.5*(double(sol1(t,c1,c2,c3,x1,x2,x3,l1,l2,l3,Sp_1))).*(sign(double(sol1(t,c1,c2,c3,x1,x2,x3,l1,l2,l3,Sp_1)))+1),0.0,Tton(i,ind));
%         area_act_ton(ind,1)
%          Tton(i,ind)
		error_ton(ind,1) = (abs(area_act_ton(ind,1)-0.5*Tton(i,ind)*(Vth(i)-E_L(i))))/(0.5*Tton(i,ind)*(Vth(i)-E_L(i)));
      %  0.5*Tton(i,ind)*(Vth(i)-E_L(i))

 	%	figure; syms tp; ezplot(sol1(tp,c1,c2,c3,x1,x2,x3,l1,l2,l3,Sp_1),[0:1000]); grid on; title('Vm_ton_lat');
%         vpa(sol1(tp,c1,c2,c3,x1,x2,x3,l1,l2,l3,Sp_1),2)

% 		% Phase 2
        t1 = Tton(i,ind);
%         t1 = Tton(i,ind);
        dV1 = double(((c1)*(x1)*(l1))*exp(l1*t1)+((c2)*(x2)*(l2))*exp(l2*t1)+((c3)*(x3)*(l3))*exp(l3*t1));
        Vr_ = Vr(i)-Sp_1-norma(param(5),low(5),up(5))*csi;
        beta1 = dV1-(Vr(i)-Vth(i))/tau_m(i)-norma(param(3),low(3),up(3))/Cm(i)+norma(param(5),low(5),up(5))/Cm(i);
%
         clear c1 c2 c3
%
%         syms c1 c2 c3
        sol2 = @(t,c1,c2,c3) Vth(i)-vpa(real(((c1*x1*exp(l1*t)+c2*x2*exp(l2*t)+c3*x3*exp(l3*t)+Sp_1))),2);

        c1 = norma(param(5),low(5),up(5));
        c2 = (l3*Vr_-beta1+norma(param(5),low(5),up(5))*l1*csi)/(l3-l2);
        c3 = (beta1-norma(param(5),low(5),up(5))*l1*csi-l2*Vr_)/(l3-l2);

        area_act_ton(ind,2) = integral(@(t) 0.5*(double(sol2(t,c1,c2,c3))).*(sign(double(sol2(t,c1,c2,c3)))+1),0.0,(Tton(i,ind)-tref(i)));

        error_ton(ind,2) = (abs(area_act_ton(ind,2)-0.5*(Tton(i,ind)-tref(i))*(Vth(i)-Vr(i))))/(0.5*(Tton(i,ind)-tref(i))*(Vth(i)-Vr(i)));
      %  (0.5*(Tton(i,ind)-tref(i))*(Vth(i)-Vr(i)))
   %    figure; syms tp; ezplot(sol2(tp,c1,c2,c3),[0:1000]); grid on; title('Vm_ton_1spk');
% 		clear delta_Vm


% 		% Phase 3
        clear c1 c2 c3
        tss = Tton(i,ind)-tref(i);
        sol3 = @(t,c1,c2,c3) Vth(i)-vpa(real(((c1*x1*exp(l1*t)+c2*x2*exp(l2*t)+c3*x3*exp(l3*t)+Sp_1))),2);

        mi = csi2*exp(l1*tss);
        eta = L2*exp(l2*tss);
        teta = L3*exp(l3*tss);

        c1 = norma(param(5),low(5),up(5));
        c2 = (Vr_*L3-Vr_*teta-norma(param(5),low(5),up(5))*mi-norma(param(3),low(3),up(3))+norma(param(5),low(5),up(5))*csi2)/(L3-L2+eta-teta);
        c3 = (norma(param(5),low(5),up(5))*mi+Vr_*(eta-L2)+norma(param(3),low(3),up(3))-norma(param(5),low(5),up(5))*csi2)/(L3-L2+eta-teta);

        area_act_ton(ind,3) = integral(@(t) 0.5*(double(sol3(t,c1,c2,c3))).*(sign(double(sol3(t,c1,c2,c3)))+1),0.0,Tton(i,ind)-tref(i));

        error_ton(ind,3) = (abs(area_act_ton(ind,3)-0.5*(Tton(i,ind)-tref(i))*(Vth(i)-Vr(i))))/(0.5*(Tton(i,ind)-tref(i))*(Vth(i)-Vr(i)));
    %   figure; ezplot(sol3(tp,c1,c2,c3),[0:1000]); grid on;  title('Vm_ton_ss')
      %  (0.5*(Tton(i,ind)-tref(i))*(Vth(i)-Vr(i)))
      %  vpa(sol3(tp,c1,c2,c3));
		clear c1 c2 c3
    end

    for p = 1:3
        error(p) = mean(error_ton(:,p).^2);
    end



	% Depolarization steps
	for dp = 1:3			% Number of depolarizing phases
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

			error_dep(ind,1,dp) = (abs(area_act_dep(ind,1,dp)-0.5*T{dp}(i,ind)*(Vth(i)-E_L(i))))/(0.5*T{dp}(i,ind)*(Vth(i)-E_L(i)));
			area_des1(ind,1) = (0.5*T{dp}(i,ind)*(Vth(i)-E_L(i)));
             area_real1(ind,1) = area_act_dep(ind,1,dp);
 %          figure; syms tp; ezplot((sol1(tp,c1,c2,c3,x1,x2,x3,l1,l2,l3,Sp_1)),[0 1000]); title(['Lat-dep',num2str(dp),' - Ades ',...
  %               num2str(area_des1(ind,1)),' Areal ',num2str(area_real1(ind,1))])
%
	% 		% Phase 2
			t1 = T{dp}(i,ind);
			dV1 = double(((c1)*(x1)*(l1))*exp(l1*t1)+((c2)*(x2)*(l2))*exp(l2*t1)+((c3)*(x3)*(l3))*exp(l3*t1));
            Vr_ = Vr(i)-Sp_1-norma(param(5),low(5),up(5))*csi;
            beta1 = dV1-(Vr(i)-Vth(i))/tau_m(i)-norma(param(3),low(3),up(3))/Cm(i)+norma(param(5),low(5),up(5))/Cm(i);
%
         clear c1 c2 c3

			sol2 = @(t,c1,c2,c3) Vth(i)-vpa(real(((c1*x1*exp(l1*t)+c2*x2*exp(l2*t)+c3*x3*exp(l3*t)+Sp_1))),2);

			c1 = norma(param(5),low(5),up(5));
			c2 = (l3*Vr_-beta1+norma(param(5),low(5),up(5))*l1*csi)/(l3-l2);
			c3 = (beta1-norma(param(5),low(5),up(5))*l1*csi-l2*Vr_)/(l3-l2);


            area_act_dep(ind,2,dp) = integral(@(t) 0.5*(double(sol2(t,c1,c2,c3))).*(sign(double(sol2(t,c1,c2,c3)))+1),0,T{dp}(i,ind)-tref(i));

			error_dep(ind,2,dp) = (abs(area_act_dep(ind,2,dp)-0.5*(T{dp}(i,ind)-tref(i))*(Vth(i)-Vr(i))))/(0.5*(T{dp}(i,ind)-tref(i))*(Vth(i)-Vr(i)));

             area_des1(ind,2) = (0.5*(T{dp}(i,ind)-tref(i))*(Vth(i)-Vr(i)));
             area_real1(ind,2) = area_act_dep(ind,2,dp);
        %     figure; syms tp; ezplot((sol2(tp,c1,c2,c3)),[0 1000]); title(['1spk-dep',num2str(dp),' - Ades ',...
        %         num2str(area_des1(ind,2)),' Areal ',num2str(area_real1(ind,2))])
%
			clear c1 c2 c3
	% 		clear delta_Vm
	%
	% 		% Phase 3
			tss = T{dp}(i,ind)*gainSFA(i,dp)-tref(i);
			sol3 = @(t,c1,c2,c3) Vth(i)-vpa(real(((c1*x1*exp(l1*t)+c2*x2*exp(l2*t)+c3*x3*exp(l3*t)+Sp_1))),2);

             mi = csi2*exp(l1*tss);
           eta = L2*exp(l2*tss);
           teta = L3*exp(l3*tss);

           c1 = norma(param(5),low(5),up(5));
           c2 = (Vr_*L3-Vr_*teta-norma(param(5),low(5),up(5))*mi-norma(param(3),low(3),up(3))+norma(param(5),low(5),up(5))*csi2)/(L3-L2+eta-teta);
          c3 = (norma(param(5),low(5),up(5))*mi+Vr_*(eta-L2)+norma(param(3),low(3),up(3))-norma(param(5),low(5),up(5))*csi2)/(L3-L2+eta-teta);

            area_act_dep(ind,3,dp) = integral(@(t) 0.5*(double(sol3(t,c1,c2,c3))).*(sign(double(sol3(t,c1,c2,c3)))+1),0.0,T{dp}(i,ind)*gainSFA(i,dp)-tref(i));

            error_dep(ind,3,dp) = (abs(area_act_dep(ind,3,dp)-0.5*(T{dp}(i,ind)*gainSFA(i,dp)-tref(i))*(Vth(i)-Vr(i))))/(0.5*(T{dp}(i,ind)*gainSFA(i,dp)-tref(i))*(Vth(i)-Vr(i)));
            area_des1(ind,3) = (0.5*(T{dp}(i,ind)*gainSFA(i,dp)-tref(i))*(Vth(i)-Vr(i)));
            area_real1(ind,3) = area_act_dep(ind,3,dp);
     %        figure; syms tp; ezplot((sol3(tp,c1,c2,c3)),[0 1000]); title(['Ss-dep',num2str(dp),' - Ades ',...
      %           num2str(area_des1(ind,3)),' Areal ',num2str(area_real1(ind,3))])
		end
    end


    for p = 1:3
        error_depol1(p) = mean(error_dep(:,p,1).^2);
    end
%
    for p = 1:3
        error_depol2(p) = mean(error_dep(:,p,2).^2);
    end

    for p = 1:3
        error_depol3(p) = mean(error_dep(:,p,3).^2);
    end


    % Post-hyperpolarization step
    % Phase 4A: Latency of rebound burst

    for ind = 1:length(Tlb)
        clear c1 c2 c3 Sp_1 Vr_

        Sp_1 = double(subs(Sp(1),{Ist,k_adap,k2,A2,k1,A1,Ie},{0,norma(param(1),low(1),up(1)),norma(param(2),low(2),up(2)),norma(param(3),low(3),up(3)),norma(param(4),low(4),up(4)),...
            norma(param(5),low(5),up(5)),norma(param(6),low(6),up(6))}));        % The external current is not active after hyperpol

        sol4 = @(t,c1,c2,c3) Vth(i)-vpa(real(((c1*x1*exp(l1*t)+c2*x2*exp(l2*t)+c3*x3*exp(l3*t)+Sp_1))),2);

        Vhyp_ss = double(subs(Sp(1),{Ist,k_adap,k2,A2,k1,A1,Ie},{Iinh(i),norma(param(1),low(1),up(1)),norma(param(2),low(2),up(2)),norma(param(3),low(3),up(3)),...
            norma(param(4),low(4),up(4)),norma(param(5),low(5),up(5)),norma(param(6),low(6),up(6))}));          %Iinh to check!!!
        Ihyp_ss = double(subs(Sp(2),{Ist,k_adap,k2,A2,k1,A1,Ie},{Iinh(i),norma(param(1),low(1),up(1)),norma(param(2),low(2),up(2)),norma(param(3),low(3),up(3)),...
            norma(param(4),low(4),up(4)),norma(param(5),low(5),up(5)),norma(param(6),low(6),up(6))}));

        phi = -Vhyp_ss/tau_m(i)-Ihyp_ss/Cm(i)+E_L(i)/tau_m(i)+norma(param(6),low(6),up(6))/Cm(i);
        c1 = 0;
        c2 = ((Vhyp_ss-Sp_1)*l3-phi)/(l3-l2);
        c3 = (-(Vhyp_ss-Sp_1)*l2+phi)/(l3-l2);


        area_act_lb(ind) = integral(@(t) 0.5*(double(sol4(t,c1,c2,c3))).*(sign(double(sol4(t,c1,c2,c3)))+1),0.0,Tlb(i,ind));

        error_lb(ind) = (abs(area_act_lb(ind)-0.5*Tlb(i,ind)*(Vth(i)-Vhyp_ss)))/(0.5*Tlb(i,ind)*(Vth(i)-Vhyp_ss));

       % figure; syms tp; ezplot((sol4(tp,c1,c2,c3)),[0 1000]); title(['Lat burst'])



        % Phase 4: first spike of rebound burst
        tlb = Tlb(i,ind);
        dVlb = ((c1*(x1(1))*(l1))*exp(l1*tlb)+((c2)*(x2(1))*(l2))*exp(l2*tlb)+((c3)*(x3(1))*(l3))*exp(l3*tlb));
        beta2 = dVlb-(Vr(i)-Vth(i))/tau_m(i)-(norma(param(3),low(3),up(3)))/Cm(i)+norma(param(5),low(5),up(5))/Cm(i); % Con A2 sottratto dal valore raggiunto al tempo dello spike

        clear c1 c2 c3

        sol5 = @(t,c1,c2,c3) Vth(i)-vpa(real(((c1*x1*exp(l1*t)+c2*x2*exp(l2*t)+c3*x3*exp(l3*t)+Sp_1))),2);
        Vr_ = Vr(i)-Sp_1-norma(param(5),low(5),up(5))*csi;
        c1 = norma(param(5),low(5),up(5));
        c2 = (l3*Vr_-beta2+norma(param(5),low(5),up(5))*l1*csi)/(l3-l2);
        c3 = (beta2-norma(param(5),low(5),up(5))*l1*csi-l2*Vr_)/(l3-l2);

        area_act_rb(ind)=integral(@(t) 0.5*(double(sol5(t,c1,c2,c3))).*(sign(double(sol5(t,c1,c2,c3)))+1),0.0,Tb(i,ind)-tref(i));

        error_rb(ind) = (abs(area_act_rb(ind)-0.5*(Tb(i,ind)-tref(i))*(Vth(i)-Vr(i))))/(0.5*(Tb(i,ind)-tref(i))*(Vth(i)-Vr(i)));
     %    figure; syms tp; ezplot(sol5(tp,c1,c2,c3),[0 1000]); title('Reb burst')
%         sol5(0,c1,c2,c3)
%         c1
%         c2
%         c3
    end

    error_lat_burst = mean(error_lb.^2);
    error_first_burst = mean(error_rb.^2);


    % AHP duration: att!! For oscillatory neurons, we consider the steady-state value, which
    % is the mean value of oscillations. To change for non oscillatory neurons
    % Phase 5: AHP1
%     for ind = 1:length(T_ahp)
%         clear c1 c2 c3
%         syms c1 c2 c3
%
%         t_vec = 0:(T_ahp(ind,1,i))*step:(T_ahp(ind,1,i))*fac;
%
%         sol6 = @(t,c1,c2,c3) vpa(real(((c1*x1*exp(l1*t)+c2*x2*exp(l2*t)+c3*x3*exp(l3*t)+Sp_1)-Vth(i))),2);
%
% %         Vdep1_ss = double(subs(Sp(1),{Ist,k_adap,k2,A2,k1,A1,Ie},{0,norma(param(1),low(1),up(1)),-1/tau_m(i),norma(param(2),low(2),up(2)),...
% %             norma(param(3),low(3),up(3)),norma(param(4),low(4),up(4)),norma(param(5),low(5),up(5))}));
%          Vdep1_ss = E_L(i);
%         Idep1_ss = double(subs(Sp(2),{Ist,k_adap,k2,A2,k1,A1,Ie},{Istim(1,i),norma(param(1),low(1),up(1)),-1/tau_m(i),norma(param(2),low(2),up(2)),...
%             norma(param(3),low(3),up(3)),norma(param(4),low(4),up(4)),norma(param(5),low(5),up(5))}));
%
%         phi = -Vdep1_ss/tau_m(i)-Idep1_ss/Cm(i)+E_L(i)/tau_m(i)+norma(param(5),low(5),up(5))/Cm(i);
%         c1 = 0;
%         c2 = ((Vdep1_ss-Sp_1)*l3-phi)/(l3-l2);
%         c3 = (-(Vdep1_ss-Sp_1)*l2+phi)/(l3-l2);
%
%         for ii = 1:length(t_vec)
%             yy(ii) = fzero_sym(@(t) sol6(t,c1,c2,c3),t_vec(ii),optimset('TolX',0.01));
%         end
%         zz = unique(yy);
%         zz = sort(zz);
%         z = round(1e2*zz)*1e-2; % Get rid of spurious differences
%         z = unique(z);
%
%         t_act_ahp(ind,1) = z(1);
%         clear zz z yy t_vec
%         error_ahp(ind,1) = w*abs(t_act_ahp(ind,1)-T_ahp(ind,1,i));
%
%     end
%
%
%     % Phase 6: AHP2
%     for ind = 1:length(T_ahp)
%         clear c1 c2 c3
%         syms c1 c2 c3
%
%         t_vec = 0:(T_ahp(ind,2,i))*step:(T_ahp(ind,2,i))*fac;
%
%         sol7 = @(t,c1,c2,c3) vpa(real(((c1*x1*exp(l1*t)+c2*x2*exp(l2*t)+c3*x3*exp(l3*t)+Sp_1)-Vth(i))),2);
%
% %         Vdep2_ss = double(subs(Sp(1),{Ist,k_adap,k2,A2,k1,A1,Ie},{0,norma(param(1),low(1),up(1)),-1/tau_m(i),norma(param(2),low(2),up(2)),...
% %             norma(param(3),low(3),up(3)),norma(param(4),low(4),up(4)),norma(param(5),low(5),up(5))}));
%         Vdep2_ss = E_L(i);
%         Idep2_ss = double(subs(Sp(2),{Ist,k_adap,k2,A2,k1,A1,Ie},{Istim(2,i),norma(param(1),low(1),up(1)),-1/tau_m(i),norma(param(2),low(2),up(2)),...
%             norma(param(3),low(3),up(3)),norma(param(4),low(4),up(4)),norma(param(5),low(5),up(5))}));
%
%         phi = -Vdep2_ss/tau_m(i)-Idep2_ss/Cm(i)+E_L(i)/tau_m(i)+norma(param(5),low(5),up(5))/Cm(i);
%         c1 = 0;
%         c2 = ((Vdep2_ss-Sp_1)*l3-phi)/(l3-l2);
%         c3 = (-(Vdep2_ss-Sp_1)*l2+phi)/(l3-l2);
%
%         for ii = 1:length(t_vec)
%             yy(ii) = fzero_sym(@(t) sol7(t,c1,c2,c3),t_vec(ii),optimset('TolX',0.01));
%         end
%         zz = unique(yy);
%         zz = sort(zz);
%         z = round(1e2*zz)*1e-2; % Get rid of spurious differences
%         z = unique(z);
%
%         t_act_ahp(ind,2) = z(1);
%         clear zz z yy t_vec
%         error_ahp(ind,2) = w*abs(t_act_ahp(ind,2)-T_ahp(ind,2,i));
%     end
%
%     % Phase 7: AHP3
%     for ind = 1:length(T_ahp)
%         clear c1 c2 c3
%         syms c1 c2 c3
%
%         t_vec = 0:(T_ahp(ind,3,i))*step:(T_ahp(ind,3,i))*fac;
%
%         sol8 = @(t,c1,c2,c3) vpa(real(((c1*x1*exp(l1*t)+c2*x2*exp(l2*t)+c3*x3*exp(l3*t)+Sp_1)-Vth(i))),2);
%
% %         Vdep3_ss = double(subs(Sp(1),{Ist,k_adap,k2,A2,k1,A1,Ie},{0,norma(param(1),low(1),up(1)),-1/tau_m(i),norma(param(2),low(2),up(2)),...
% %             norma(param(3),low(3),up(3)),norma(param(4),low(4),up(4)),norma(param(5),low(5),up(5))}));
%         Vdep3_ss = E_L(i);
%         Idep3_ss = double(subs(Sp(2),{Ist,k_adap,k2,A2,k1,A1,Ie},{Istim(3,i),norma(param(1),low(1),up(1)),-1/tau_m(i),norma(param(2),low(2),up(2)),...
%             norma(param(3),low(3),up(3)),norma(param(4),low(4),up(4)),norma(param(5),low(5),up(5))}));
%
%
%         phi = -Vdep3_ss/tau_m(i)-Idep3_ss/Cm(i)+E_L(i)/tau_m(i)+norma(param(5),low(5),up(5))/Cm(i);
%         c1 = 0;
%         c2 = ((Vdep3_ss-Sp_1)*l3-phi)/(l3-l2);
%         c3 = (-(Vdep3_ss-Sp_1)*l2+phi)/(l3-l2);
%
%        for ii = 1:length(t_vec)
%             yy(ii) = fzero_sym(@(t) sol8(t,c1,c2,c3),t_vec(ii),optimset('TolX',0.01));
%         end
%         zz = unique(yy);
%         zz = sort(zz);
%         z = round(1e2*zz)*1e-2; % Get rid of spurious differences
%         z = unique(z);
%
%         t_act_ahp(ind,3) = z(1);
%         error_ahp(ind,3) = w*abs(t_act_ahp(ind,3)-T_ahp(ind,3,i));
%
%     end
%
%     for p = 1:3
%
%         errorAHP(p) = mean(error_ahp(:,p).^2);
%     end
%
    sq_err = double(sqrt(mean([error error_depol1 error_depol2 error_depol3 error_lat_burst error_first_burst])));

   error_all = [error_all; error error_depol1 error_depol2 error_depol3 error_lat_burst error_first_burst];
   cf = [cf;sq_err];
end

function [ c, ceq ] = confun_eglif( param, low, up, neu_ind, Iin,  Cm, tau_m, E_L, Vth, Vinh_ss, t_ref, L2, L3, Sp, Tton, Istim, V1,V2, Vss,Tdep3,gainSFA, I_ahp1, I_ahp2, I_ahp3 )

    global con;
    
%CONFUN check the constraints (both linear and non linear)

    
    syms k_adap k2 Ist A2 k1 A1 Ie t1 tss t
 
    % Oscillation parameters
    % Frequency
    delta = (-1/tau_m(neu_ind)-norma(param(2),low(2),up(2)))^2-4*(norma(param(2),low(2),up(2))/tau_m(neu_ind)+norma(param(1),low(1),up(1))/Cm(neu_ind));     % [1/ms^2]
    Vinh_act=double(subs(Sp(1),{Ist,k_adap,k2,Ie},{Iin,norma(param(1),low(1),up(1)),norma(param(2),low(2),up(2)),norma(param(6),low(6),up(6))}));
    Sp_ton=double(subs(Sp(1),{Ist,k_adap,k2,Ie},{0,norma(param(1),low(1),up(1)),norma(param(2),low(2),up(2)),norma(param(6),low(6),up(6))}));

   % Nonlinear inequality constraints
  c = [
        

        % delta positivo!
  -delta/abs(delta)+0.000001;    
  
  
  % Steady state during hyperpolarization
   Vinh_act/(Vinh_ss-35)-1;
   Vinh_act/abs(Vinh_ss+35)+1;
   
      
 % SS Vm during tonic below Vth
  Vth(neu_ind)/Sp_ton-1;
    
];

con = [con;c'];
    % Nonlinear equality constraints
        %ceq = [param(2)+1/tau_m(neu_ind)];      % To have null damping
        ceq = [];
end

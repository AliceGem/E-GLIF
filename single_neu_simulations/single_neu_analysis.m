close all
clc
clear
n_sim = 1;

neu=1;
neu_type = 3;       % 1 = PC; 2 = MLI; 3 = DCN
neu_names = {'PC','MLI','GoC'};
PC_BPR = 0;

protocol = 1;       % 1 = protocol from literature; 2 = protocol from experiments

deltaSFA=[100 100 100 200 250];
np = 9;     % Numero fasi del protocollo di stimolazione: 9 steps di corrente

if PC_BPR
    durata = [10000 1000 1000 1000 1000 10 1000 1000 1000];
    durata2 = [10000 1000 1000 1000 1000 10 1000 1000 1000];
    ep = 2;     % Numero fasi eccitatorie
else
    durata = [10000 1000 1000 1000 1000 1000 1000 1000 1000];
    durata2 = [10000 1000 1000 1000 1000 1000 1000 1000 1000];    % Durata degli intervalli in cui calcolare la frequenza per il plot finale
    ep = 3;     % Numero fasi eccitatorie
end

rp = 9;
intSFA = 1000;      % [ms] interval where we compute SFA (f1s/fin)
delta = deltaSFA(neu); % [ms]
rebound = [8];          % intervals where we start to measure rebound properties
delta_rebound = [150 400 200 150 150];   % Window to compute reb burst for each neuron type

Istim = [{'400','800','2400'},{'12' ,'24','36'},{'142','284','426'}]

% PARAMETRI
% For each neuron we need mean frequency values in each phase
% (lasting 1 s) and in the first 2 spikes and in the ending 5 ms of each excitatory phase
% + during rebound:
% f_mean = mean(ISI)
% fin = mean(ISI(first2spikes))
% f1s = mean(ISI(last5ms))

if durata(1) == 1000
    long_phase = 6;
else
    long_phase = 1;
end
 path = '/home/alice/workspace/E-GLIF/single_neu_simulations/';


% Input current definition
ns = 5;
dur_dep = 1000;
dur_res = 30;       %[ms]
if neu_type == 2 || neu_type == 3
    if n_sim==1
        mult = load([path,'eglif_',neu_names{neu_type},'multimeter-72-0.dat']);
        spk = load([path,'eglif_',neu_names{neu_type},'spikes-71-0.gdf']);
    else
        mult = load([path,'eglif_',neu_names{neu_type},'multimeter-72-0_',neu_names{neu_type},'1.dat']);
        spk = load([path,'eglif_',neu_names{neu_type},'spikes-71-0_',neu_names{neu_type},'1.gdf']);
    end

else
    if n_sim == 1
        mult = load([path,'eglif_',neu_names{neu_type},'multimeter-76-0.dat']);
        spk = load([path,'eglif_',neu_names{neu_type},'spikes-73-0.gdf']);
    else
        mult = load([path,'eglif_',neu_names{neu_type},'multimeter-76-0_',neu_names{neu_type},'1.dat']);
        spk = load([path,'eglif_',neu_names{neu_type},'spikes-73-0_',neu_names{neu_type},'1.gdf']);
    end
end
x = [0:1:length(mult(:,1))-1];
I = mult(:,end);
step_values = [I(durata(1)+5) I(sum(durata(1:3))+5) I(sum(durata(1:5))+5)]



figure;
plot(mult(:,2),mult(:,3),'k','LineWidth',2)
xlim([0 max(mult(:,2))])
hold on
plot(mult(:,2),mult(:,4),'--','Color',[7 93 13]/255,'LineWidth',2)
hold on
plot(mult(:,2),-85+I/32,'r','LineWidth',2)
hold on
for sp = 1:length(spk)
   line([spk(sp,2) spk(sp,2)],[-20 -10],'Color','k','LineWidth',2)
   line([spk(sp,2) spk(sp,2)],[mult(1,4)-0.5 -20],'LineStyle','--','Color','k','LineWidth',1)
end
xlabel('Time [s]')
ylabel('V_m [mV]')
set(gca,'FontSize',12)


% Focus on rebound phase
REB_FOCUS = 0
if REB_FOCUS
    r = length(rebound);
    xlim([sum(durata(1:rebound(r)))-200 sum(durata)])
    ylim([-100 0])
    set(gca,'Ytick',[-180:40:0],...
         'XTick',[16800:200:18000],'XTickLabel',{'16.8','17','17.2','17.4','17.6','17.8','18'})
    line([sum(durata(1:rebound(r))) sum(durata(1:rebound(r)))],[-180 0],'Color','k','LineStyle',':')

    saveas(gcf,['Rebound',neu_names{neu_type},'.fig']);
    saveas(gcf,['Rebound',neu_names{neu_type},'.bmp']);
    saveas(gcf,['Rebound',neu_names{neu_type},'.eps']);
end

figure
plot(mult(find(mult(:,1)==neu),5),'g')
hold on
plot(mult(find(mult(:,1)==neu),6),'r')
title('Spike-triggered currents')
xlabel('t [ms]')
ylabel('I [pA]')
hold on
for k = 1:np
line([sum(durata(1:k)) sum(durata(1:k))],[-100 160],'Color',[0.8 0.8 0.8])
hold on
end

legend({'I_d_e_p','I_a_d_a_p'},'location','southeast')



for sm = 1:n_sim

if neu_type == 2 || neu_type == 3
    if n_sim==1
        mult = load([path,'eglif_',neu_names{neu_type},'multimeter-72-0.dat']);     % The number depends on the implicit index of the created node
        spk = load([path,'eglif_',neu_names{neu_type},'spikes-71-0.gdf']);
    else
        mult = load([path,'eglif_',neu_names{neu_type},'multimeter-72-0_',neu_names{neu_type},num2str(sm),'.dat']);
        spk = load([path,'eglif_',neu_names{neu_type},'spikes-71-0_',neu_names{neu_type},num2str(sm),'.gdf']);
    end
else
    if n_sim ==1
        mult = load([path,'eglif_',neu_names{neu_type},'multimeter-76-0.dat']);
        spk = load([path,'eglif_',neu_names{neu_type},'spikes-73-0.gdf']);
    else
        mult = load([path,'eglif_',neu_names{neu_type},'multimeter-72-0_',neu_names{neu_type},num2str(sm),'.dat']);
        spk = load([path,'eglif_',neu_names{neu_type},'spikes-71-0_',neu_names{neu_type},num2str(sm),'.gdf']);
    end
end
spk_neu = (spk(find(spk(:,1)==neu),2));


% Frequency
for i = 1:length(durata2)
    clear spi
    spi = spk_neu(find(spk_neu(:,1)<sum(durata(1:i))+1));
    spi = spi(find(spi(:,1)>=sum(durata(1:i))-durata(i)+1));
    if size(spi,1)>1
        ISI_neu {i}= diff(spi);
        freq_tonic_neu_ISI(i) = 1/(mean(ISI_neu{i})*0.001);    % Frequency
        % of tonic activity in each phase computed from mean ISI
        freq_tonic_neu(i) = size(spi,1)/(durata(i)/1000);
    else
        freq_tonic_neu_ISI(i) = 0;
        freq_tonic_neu(i) = 0;
    end
    if PC_BPR && i == 6
       spk_burst = spi;
    end
    if PC_BPR && i == 7
       spk_pause = spi;
    end
    if i == long_phase
        if size(spi,1)>1
             coefvar = std(ISI_neu{i})/mean(ISI_neu{i})
        else
            coefvar = NaN;
        end
    end
end


% SFA new method
SFA = [];
for i = 1:ep             % Num dep phases
    f1s(i) = 1/(mean(ISI_neu{2*i}(end-5:end))*0.001);
    fin(i) = 1/(mean(ISI_neu{2*i}(1:2))*0.001);
    SFA=[SFA f1s(i)/fin(i)];
end

if PC_BPR
    fin(i+1) = 1000/(mean(ISI_neu{2*(i+1)}));
    pause = spk_pause(1)-spk_burst(end);
end

% Rebound burst
for r = 1:length(rebound)
spk_rebound{r} = spk_neu(find(spk_neu(:,1)<sum(durata(1:rebound(r)))+delta_rebound(neu_type)));
spk_rebound{r} = spk_rebound{r}(find(spk_rebound{r}(:,1)>=sum(durata(1:rebound(r)))));
if size(spk_rebound{r},1)>1
    ISI_rebound{r} = diff(spk_rebound{r});
    freq_rebound = 1/(mean(ISI_rebound{r})*0.001);
    if length(ISI_rebound{r})<2 || isempty(ISI_rebound)
        first_reb(r) = NaN;
        lat_reb(r) = NaN;
    else
        first_reb(r) = 1/(ISI_rebound{r}(1)*0.001);
        lat_reb(r) = spk_rebound{r}(1)-sum(durata(1:rebound(r)));
    end
else
    freq_rebound = 0;
end
end


    % Save the computed paremeters for simulation sim and clear the
    % variables
    frequencies(sm,:) = freq_tonic_neu;
    frequencies_ISI(sm,:) = freq_tonic_neu_ISI;
    CVs(sm) = coefvar;
    freq_in(sm,:)=fin;
    freq_1s(sm,:)=f1s;
    SFAs(sm,:) = SFA;
%     p_fit{sm} = [pin; p1s];
    freq_rebounds(sm) = freq_rebound;
    first_rebs(sm,:) = first_reb;
    lat_rebs(sm,:) = lat_reb;


%
%     clear spi spk_neu freq_tonicGoC freq_tonicGoC_ISI ISI_GoC CV SFA f1s fin p1s pin...
%         freq_rebound spk_rebound ISI_rebound spk_res ISI_res ISI_res_step_first first_res last_res...
%         mean_res ISI_burst burst_res lat_res ISI_res_step

end
%
freq_all{protocol} = frequencies;
freq_all_ISI{protocol} = frequencies_ISI;
CV_all{protocol} = CVs;
fin_all{protocol} = freq_in;
f1s_all{protocol} = freq_1s;
SFA_all{protocol} = SFAs;
reb_all{protocol} = freq_rebounds;


%% Vm plot
plot(mult(find(mult(:,1)==neu),3),'k')
hold on
plot(mult(find(mult(:,1)==neu),4),'r--')
%xlim([0 sum(durata)])c
ylim([-150 0])
ylabel('V [mV]')
xlabel('t [ms]')
title('EGLIF neuron')
hold on
%% Plot frequencies, SFA and resonance parameters

% Frequencies

xax = [1:4]

% Plot REFERENCE values from literature, used as targets in optimization
fI_PC = 0.1;        % Masoli
mean_des_all = [60.0 ([800.0 1600.0 2400.0]*fI_PC+80).*[1 1 2.5];  ...   % PC
    8.5 30 60 90;...
    30 50 80 110];
mean_des_all = [60.0 ([400.0 800.0 2400.0]*fI_PC+80).*[1 1 2.5];  ...   % PC
    8.5 30 60 90;...
    30 50 80 110];
sd_des_all = [7 1 1 1;...
    2.7 1 5 10;...
    6 2 5 15];
mean_des_end_all = [mean_des_all(:,2:end)./[1 1.25 10;1 1 1;1.1 1.1 1.4]];              % DCN with linear adaptation
mean_des_end_all = [mean_des_all(:,2:end)./[1 1.25 10;1 1 1;1.2 1.2 1.2]];              % DCN with constant adaptation
sd_des_end_all = sd_des_all(:,2:end);
mean_des = mean_des_all(neu_type,:);
sd_des = sd_des_all(neu_type,:);
mean_des_end = mean_des_end_all(neu_type,:);
sd_des_end = sd_des_end_all(neu_type,:);
figure
plot(xax,mean_des,'go','MarkerFaceColor','g')
hold on
errorbar(xax,mean_des,sd_des,'g','LineStyle','none')
plot(xax(2:end),mean_des_end,'gsquare')
hold on
errorbar(xax(2:end),mean_des_end,sd_des_end,'g','LineStyle','none')


%Plot MODEL values
h1 = plot(xax,[mean(frequencies_ISI(:,1)) ...
    mean(freq_in(:,1)) mean(freq_in(:,2)) mean(freq_in(:,3))],'ko','MarkerFaceColor','k')
hold on
h2 = plot(xax([2:size(freq_1s,2)+1]),mean(freq_1s,1),'ksquare')
errorbar(xax([2:size(freq_1s,2)+1]),mean(freq_1s,1),std(freq_1s,0,1),'k','LineStyle','none')
errorbar(xax,[mean(frequencies_ISI(:,1)) mean(freq_in,1)],[std(frequencies_ISI(:,1)) std(freq_in,0,1)],'k','LineStyle','none')
legend([h1 h2],{'f_i_n','f_1_s'},'location','northeastoutside')

p1s_mean = polyfit(step_values(1:size(freq_1s,2)),mean(freq_1s,1),1)
pin_mean = polyfit(step_values,mean(freq_in,1),1)
plot([2:4],p1s_mean(1).*step_values+p1s_mean(2),'k:')
plot([2:4],pin_mean(1).*step_values+pin_mean(2),'k--')
set(gca,'XTick',[1:4],'XTickLabel',{'0',Istim{(neu_type-1)*3+1:(neu_type)*3}},'FontSize',14)
ylim([0 Inf])
grid on
xlim([0.5 4.5])

ylabel('f [Hz]')
xlabel('I_s_t_i_m [pA]')
title(['f-I_s_t_i_m plot ',neu_names{neu_type}])
saveas(gcf,['If_ModLit',neu_names{neu_type},'.fig']);
saveas(gcf,['If_ModLit',neu_names{neu_type},'.bmp'])
saveas(gcf,['If_ModLit',neu_names{neu_type},'.eps'])

%%
% SFA - example for GoC (not computed for MLI)
xSFA = [30 50 100 150]
SFA_des = [0.9 0.6 0.4 0.3;0.05 0.05 0.05 0.1];  %mean; SD
figure;
plot(mean(freq_in,1),mean(SFAs,1),'k^')
hold on
h=errorbar(mean(freq_in,1),mean(SFAs,1),std(SFAs,1,1),'k')
plot(xSFA,SFA_des(1,:),'g')
hold on
fill([xSFA,flip(xSFA)],[SFA_des(1,:)-1*SFA_des(2,:) flip(SFA_des(1,:)+1*SFA_des(2,:))],'r','LineStyle','none')
alpha(0.25)
%xlim([0.5 4.5])
set(gca,'XTick',[50:50:200],'FontSize',14)
ylim([0 1.2])
xlim([0 170])
xlabel('f_i_n')
ylabel('f_1_s/f_i_n')
title('SFA')
grid on
saveas(gcf,'SFA_ModLit.fig')
saveas(gcf,'SFA_ModLit.bmp')
saveas(gcf,'SFA_ModLit.eps')

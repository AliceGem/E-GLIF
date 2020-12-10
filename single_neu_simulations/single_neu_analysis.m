%close all
clc
%clear
n_sim = 10; % numero di simulazioni

neu=1; % neurone da selezionare tra quelli presenti nella prima colonna dei file
neu_type = 3;       % 1 = PC; 2 = MLI; 3 = CA1PC
neu_names = {'PC','MLI','CA1PC'};

np = 9;     % Numero fasi del protocollo di stimolazione: 9 steps di corrente
new_freq = 1;
durata = [10000 1000 1000 1000 1000 1000 1000 1000 1000];    % durata vari stage del protocollo
durata2 = [10000 1000 1000 1000 1000 1000 1000 1000 1000];    % Durata degli intervalli in cui calcolare la frequenza per il plot finale
ep = 3;     % Numero fasi eccitatorie

Istim = [{'400','800','2400'},{'12' ,'24','36'},{'100','200','300'}];

% PARAMETRI
% For each neuron we need mean frequency values in each phase
% (lasting 1 s) and in the first 2 spikes and in the ending 5 ms of each excitatory phase
% + during rebound:
% f_mean = mean(ISI)
% fin = mean(ISI(first2spikes))
% f1s = mean(ISI(last5ms))

long_phase = 1;

% Input current definition
ns = 5;
dur_dep = 1000;
dur_res = 30;       %[ms]

mult = load('eglif_CA1PC_1_multimeter-12-0.dat');
spk = load('eglif_CA1PC_1_spikes-11-0.gdf');

x = [0:1:length(mult(:,1))-1];
I = mult(:,end);
%step_values = [I(durata(1)+5) I(sum(durata(1:3))+5) I(sum(durata(1:5))+5)]

figure;
plot(mult(:,2),mult(:,3),'k','LineWidth',2)
xlim([0 max(mult(:,2))])
hold on
plot(mult(:,2),mult(:,4),'--','Color',[7 93 13]/255,'LineWidth',2)
hold on
plot(mult(:,2),-85+I/32,'r','LineWidth',2)
hold on
%Spikes
for sp = 1:length(spk)
   line([spk(sp,2) spk(sp,2)],[40 50],'Color','k','LineWidth',2)
   line([spk(sp,2) spk(sp,2)],[mult(1,4)-0.5 40],'LineStyle','--','Color','k','LineWidth',1)
end
xlabel('Time [ms]')
ylabel('V_m [mV]')
set(gca,'FontSize',12)

xlim([0 100])

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

legend({'Idep','Iadap'},'location','southeast')

xlim([0 100])

for sm = 1:n_sim
    mult = load(append('eglif_CA1PC_',num2str(sm),'_multimeter-12-0.dat'));
    spk = load(append('eglif_CA1PC_',num2str(sm),'_spikes-11-0.gdf'));
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
        if(length(ISI_neu{2*i})>2)
            if(length(ISI_neu{2*i})>5)
                f1s(i) = 1/(mean(ISI_neu{2*i}(end-5:end))*0.001);
            else
                f1s(i) = 1/(mean(ISI_neu{2*i}(end-2:end))*0.001);
            end
        else
            f1s(i) = 1/(mean(ISI_neu{2*i})*0.001);
        end
        fin(i) = 1/(mean(ISI_neu{2*i}(1))*0.001);
        SFA=[SFA f1s(i)/fin(i)];
    end

    % Save the computed paremeters for simulation sim and clear the variables
    frequencies(sm,:) = freq_tonic_neu;
    frequencies_ISI(sm,:) = freq_tonic_neu_ISI;
    frequencies_dep{sm} = ISI_neu;
    CVs(sm) = coefvar;
    freq_in(sm,:)=fin;
    freq_1s(sm,:)=f1s;
    SFAs(sm,:) = SFA;
end

frequencies_ISI
%% Plot frequencies, SFA and resonance parameters

% Frequencies
xax = [1:4];

% Plot REFERENCE values from literature, used as targets in optimization

mean_des = [0 47.7 90.8 125.2];
sd_des = [0 4.8 6.6 6.2];


mean_des_end = [23.7 30.8 41.9];
sd_des_end = [1.3 1.5 0.9];

figure
plot(xax,mean_des,'go','MarkerFaceColor','g')
hold on
errorbar(xax,mean_des,sd_des,'g','LineStyle','none')
plot(xax(2:end),mean_des_end,'gsquare')
hold on
errorbar(xax(2:end),mean_des_end,sd_des_end,'g','LineStyle','none')


%Plot MODEL values
h1 = plot(xax,[mean(frequencies_ISI(:,1)) ...
    mean(freq_in(:,1)) mean(freq_in(:,2)) mean(freq_in(:,3))],'ko','MarkerFaceColor','k');
hold on
h2 = plot(xax([2:size(freq_1s,2)+1]),mean(freq_1s,1),'ksquare');
errorbar(xax([2:size(freq_1s,2)+1]),mean(freq_1s,1),std(freq_1s,0,1),'k','LineStyle','none')
errorbar(xax,[mean(frequencies_ISI(:,1)) mean(freq_in,1)],[std(frequencies_ISI(:,1)) std(freq_in,0,1)],'k','LineStyle','none')
legend([h1 h2],{'f_i_n','f_1_s'},'location','northeastoutside')

p1s_mean = polyfit(step_values(1:size(freq_1s,2)),mean(freq_1s,1),1);
pin_mean = polyfit(step_values,mean(freq_in,1),1);
plot([2:4],p1s_mean(1).*step_values+p1s_mean(2),'k:')
plot([2:4],pin_mean(1).*step_values+pin_mean(2),'k--')
set(gca,'XTick',[1:4],'XTickLabel',{'0',Istim{(neu_type-1)*3+1:(neu_type)*3}},'FontSize',14)
ylim([0 inf])
grid on
xlim([0.5 4.5])
ylabel('f [Hz]')
xlabel('I_s_t_i_m [pA]')
title(['f-I_s_t_i_m plot ',neu_names{neu_type}])

%% TO ADD for VERIFICATION of EXP DATA
% 1) check of mean frequencies during 3 stimulation steps
% 2) check of number of spikes during 3 stimulation steps
mean_freq_des_mean = [30.4 47.5 57.3]
mean_freq_des_sd = [2.0 1.9 1.8]
spk_count_mean = [11.2 18.3 22.2]
spk_count_des = [0.8 0.9 0.7]
%%
n1 = [];
n2 = [];
n3 = [];
for sm = 1:n_sim
    freq = frequencies_dep(sm);
    n1 = [n1, length(freq{1}{2})];
    n2 = [n2, length(freq{1}{4})];
    n3 = [n3, length(freq{1}{6})];
end

N1 = 30;
N2 = 30;
N3 = 30;
N1 = min(n1);
N2 = min(n2);
N3 = min(n3);

f_matrix_1 = zeros(n_sim,N1);
f_matrix_2 = zeros(n_sim,N2);
f_matrix_3 = zeros(n_sim,N3);
for sm = 1:n_sim
    freq = frequencies_dep(sm);
    f_matrix_1(sm,:) = 1./freq{1}{2}(1:N1)*1000;
    f_matrix_2(sm,:) = 1./freq{1}{4}(1:N2)*1000;
    f_matrix_3(sm,:) = 1./freq{1}{6}(1:N3)*1000;
end
data = mean(f_matrix_3);
data = data/data(1)*100;
accuracy = round((1 - length(find(data>100))/length(data))*100,1);
data(data > 100) = [];
x = 1:length(data);
p = polyfit(x,data,4);
y = [data(1), polyval(p,x(1:end-1))];
x = 1:length(y);
p = polyfit(x,y,5);
y2 = [data(1), polyval(p,x(1:end-1))];
% figure
% plot(data,'.','MarkerSize',15)
% hold on
% plot(x,y)
% hold on
% plot(x,y2)
% legend('data','1','2')

figure
plot(data,'.','MarkerSize',15)
hold on
errorbar(data, std(f_matrix_3,0,1),'b','LineStyle','none')
hold on
plot(x,y2)
title(['SFA results with','Accuracy ',num2str(accuracy),'%'])

mean_des = [0 16 40 65];
sd_des = [0 5 6 7];
figure
plot(xax,[mean(frequencies_ISI(:,1)) ...
    mean(mean(f_matrix_1)) mean(mean(f_matrix_2)) mean(mean(f_matrix_3))],'ko','MarkerFaceColor','k','MarkerSize',8);
hold on
errorbar(xax,[mean(frequencies_ISI(:,1)) mean(mean(f_matrix_1)) mean(mean(f_matrix_2))...
    mean(mean(f_matrix_3))],[std(frequencies_ISI(:,1)) std(f_matrix_1,0,1)...
    std(mean(f_matrix_2)',0,1) std(mean(f_matrix_3)',0,1)],'k','LineStyle','none')
hold on
plot(xax,mean_des,'go','MarkerFaceColor','g','MarkerSize',6)
hold on
errorbar(xax,mean_des,sd_des,'g','LineStyle','none')
xlim([0.5 4.5])
grid on
title('Mean frequencies')

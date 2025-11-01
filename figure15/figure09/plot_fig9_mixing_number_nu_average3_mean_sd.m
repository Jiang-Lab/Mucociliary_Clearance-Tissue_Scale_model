function plot_fig9_mixing_number_nu_average3_mean_sd


clear all; clc
cd diagnosis-mixing-nu-2groups

    figure(1);clf(1)  % xy plane
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% nu=0.1
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    filename = 'mixing_nu0d1_4_uniform';
    s = load(filename);  % iter, time, mix_number
    time = s(:,2); 
    mix_number = s(:,3);  m0 = s(1,3); m1_value = mix_number/m0; 
    
    filename = 'mixing_nu0d1_5_uniform';
    s = load(filename);  % iter, time, mix_number
    mix_number = s(:,3); m0 = s(1,3);  m2_value = mix_number/m0; 

    filename = 'mixing_nu0d1_6_uniform';
    s = load(filename);  % iter, time, mix_number
    mix_number = s(:,3); m0 = s(1,3);  m3_value = mix_number/m0; 
    
    alls=[m1_value m2_value m3_value];
    avg_data = mean(alls,2); % column wise
    for k=1:length(avg_data)
        std_data(k,1) = std(alls(k,:));
    end
    % avg_data
    % std_data
    loglog(time,avg_data,'r-','linewidth',2,'markersize',8); hold on
    fill([time(2:end); flip(time(2:end))], [avg_data(2:end)-std_data(2:end); flip(avg_data(2:end)+std_data(2:end))],...
        [1 0 0],'FaceAlpha',0.1,'LineStyle','none'); hold on
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% nu=0.2
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    filename = 'mixing_nu0d2_2_uniform';
    s = load(filename);  % iter, time, mix_number
    time = s(:,2); 
    mix_number = s(:,3);  m0 = s(1,3); m1_value = mix_number/m0; 

    filename = 'mixing_nu0d2_3_uniform';
    s = load(filename);  % iter, time, mix_number
    mix_number = s(:,3); m0 = s(1,3);  m2_value =  mix_number/m0; 

    filename = 'mixing_nu0d2_4_uniform';
    s = load(filename);  % iter, time, mix_number
    mix_number = s(:,3); m0 = s(1,3);  m3_value =  mix_number/m0; 

    alls=[m1_value m2_value m3_value];
    avg_data = mean(alls,2); % column wise
    for k=1:length(avg_data)
        std_data(k,1) = std(alls(k,:));
    end

    fill([time(2:end); flip(time(2:end))], [avg_data(2:end)-std_data(2:end); flip(avg_data(2:end)+std_data(2:end))],...
        [0 0 1],'FaceAlpha',0.1,'LineStyle','none'); hold on
    loglog(time,avg_data,'b-','linewidth',2,'markersize',8); hold on
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% nu=0.4
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    filename = 'mixing_nu0d4_1_uniform';
    s = load(filename);  % iter, time, mix_number
    time = s(:,2); 
    mix_number = s(:,3);  m0 = s(1,3); m1_value = mix_number/m0; 

    filename = 'mixing_nu0d4_2_uniform';
    s = load(filename);  % iter, time, mix_number
    mix_number = s(:,3); m0 = s(1,3);  m2_value = mix_number/m0;  

    filename = 'mixing_nu0d4_3_uniform';
    s = load(filename);  % iter, time, mix_number
    mix_number = s(:,3); m0 = s(1,3);  m3_value = mix_number/m0; 

    alls=[m1_value m2_value m3_value];
    avg_data = mean(alls,2); % column wise
    for k=1:length(avg_data)
        std_data(k,1) = std(alls(k,:));
    end

    fill([time(2:end); flip(time(2:end))], [avg_data(2:end)-std_data(2:end); flip(avg_data(2:end)+std_data(2:end))],...
        [0 1 0],'FaceAlpha',0.1,'LineStyle','none'); hold on
    plot(time,avg_data,'g-','linewidth',2,'markersize',8); hold on

    % loglog(time,m_value/3,'gs-','linewidth',2); hold on


grid on
xlim([0.08,14])
ylim([0.5,1.5])
set(gca,'fontsize',30)
set(gca,'xtick',[0.1 1 10])
% set(gca,'ytick',[0.5 1])
xlabel time
ylabel {$m/m_0$}

legend('$\nu=0.1$','','','$\nu=0.2$','','$\nu=0.4$','Location','northeastoutside')

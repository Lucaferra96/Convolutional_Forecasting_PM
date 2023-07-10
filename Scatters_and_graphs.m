close all
clear
clc

%% System setup
addpath('C:\Users\ferrarilu\OneDrive - ETH Zurich\Tesi_Olgiati\Input',...
    'C:\Users\ferrarilu\OneDrive - ETH Zurich\Tesi_Olgiati\Funzioni',...
    'Output'); 

outFOLD_network1 = 'C:\Users\ferrarilu\OneDrive - ETH Zurich\Tesi_Olgiati\Output\20230613_Lombardia_CNN_IY_3In_10Out_20N_20Conv_Georef_v\';
outFOLD_network2 = 'C:\Users\ferrarilu\OneDrive - ETH Zurich\Tesi_Olgiati\Output\20230613_Australia_CNN_IY_3In_10Out_20N_20Conv_Georef_v\';
outFOLD_network3 = 'C:\Users\ferrarilu\OneDrive - ETH Zurich\Tesi_Olgiati\Output\20230613_Poland_CNN_IY_3In_10Out_20N_20Conv_Georef_v\';

net_type1 = 1;    % 0 = LSTM Network 1 = Convolutional Neural Network 2 = ConvLSTM 3 = Feed Forward
net_type2 = 1;
net_type3 = 1;

NumInputs = 3;
NumOutputs = 10; 

domain1 = 0; % 0 = Lombardia; 1 = Australia 2 = Poland
domain2 = 1;
domain3 = 2;

switch domain1
    case 0
        dm1='Lombardia_';
    case 1
        dm1='Australia_';
    case 2
        dm1='Poland_';
end

switch domain2
    case 0
        dm2='Lombardia_';
    case 1
        dm2='Australia_';
    case 2
        dm2='Poland_';
end

switch domain3
    case 0
        dm3='Lombardia_';
    case 1
        dm3='Australia_';
    case 2
        dm3='Poland_';
end

%% Load first net
path = strcat(outFOLD_network1, 'net_opt.mat');
load(path);

path = strcat(outFOLD_network1, 'Data.mat');
load(path);

%% Evaluate model adjusted
if net_type == 3
    output_test = sim(net_opt,input_test);
else
    output_test = predict(net_opt,input_test);
end

output_test = std_test.*output_test + mean_test ;

max_err_1 = max(abs(output_test - PM10_target_test), [], 2);

for i = 1:NumOutputs
    corr_test(i) = 1 - sum((output_test(i,:) - PM10_target_test(i,:)).^2)/(sum((mean(PM10_target_test(i,:)) - PM10_target_test(i,:)).^2));
    nrmse_test(i) = sqrt((sum((PM10_target_test(i,:)-output_test(i,:)).^2))/(sum((output_test(i,:)-mean(output_test(i,:))).^2))); 
end

PM10_target_test_net1 = [PM10_target_test(3,:); PM10_target_test(7,:)];
output_test_net1 = [output_test(3,:); output_test(7,:)];
corr_test_net1 = [corr_test(3); corr_test(7)];
nrmse_test_net1 = [nrmse_test(3); nrmse_test(7)];

%% Load second net
path = strcat(outFOLD_network2, 'net_opt.mat');
load(path);

path = strcat(outFOLD_network2, 'Data.mat');
load(path);

%% Evaluate model adjusted
if net_type == 3
    output_test = sim(net_opt,input_test);
else
    output_test = predict(net_opt,input_test);
end

output_test = std_test.*output_test +mean_test ;

max_err_2 = max(abs(output_test - PM10_target_test), [], 2);

for i = 1:NumOutputs
    corr_test(i) = 1 - sum((output_test(i,:) - PM10_target_test(i,:)).^2)/(sum((mean(PM10_target_test(i,:)) - PM10_target_test(i,:)).^2));
    nrmse_test(i) = sqrt((sum((PM10_target_test(i,:)-output_test(i,:)).^2))/(sum((output_test(i,:)-mean(output_test(i,:))).^2))); 
end

PM10_target_test_net2 = [PM10_target_test(3,:); PM10_target_test(7,:)];
output_test_net2 = [output_test(3,:); output_test(7,:)];
corr_test_net2 = [corr_test(3); corr_test(7)];
nrmse_test_net2 = [nrmse_test(3); nrmse_test(7)];

%% Load third net
path = strcat(outFOLD_network3, 'net_opt.mat');
load(path);

path = strcat(outFOLD_network3, 'Data.mat');
load(path);

%% Evaluate model adjusted
if net_type == 3
    output_test = sim(net_opt,input_test);
else
    output_test = predict(net_opt,input_test);
end

output_test = std_test.*output_test +mean_test ;

max_err_3 = max(abs(output_test - PM10_target_test), [], 2);

for i = 1:NumOutputs
    corr_test(i) = 1 - sum((output_test(i,:) - PM10_target_test(i,:)).^2)/(sum((mean(PM10_target_test(i,:)) - PM10_target_test(i,:)).^2));
    nrmse_test(i) = sqrt((sum((PM10_target_test(i,:)-output_test(i,:)).^2))/(sum((output_test(i,:)-mean(output_test(i,:))).^2))); 
end

PM10_target_test_net3 = [PM10_target_test(3,:); PM10_target_test(7,:)];
output_test_net3 = [output_test(3,:); output_test(7,:)];
corr_test_net3 = [corr_test(3); corr_test(7)];
nrmse_test_net3 = [nrmse_test(3); nrmse_test(7)];

%% Plots

figure;
tiledlayout(2,3);

for i = 1:2
    nexttile;

    plot(PM10_target_test_net1(i,:),output_test_net1(i,:),'r*')
    grid on
    hold on
    
    plot([0 max(max([PM10_target_test_net1(i,:) output_test_net1(i,:)]))],[0 max(max([PM10_target_test_net1(i,:) output_test_net1(i,:)]))],'b--');
    ylim([0, max(output_test_net1(i,:), [], 'all')]); 
    s1=xlabel('PM10 concentration [\mug/m^3]');
    s2=ylabel('CNN [\mug/m^3]');
    s3=text(0.02*max(max([PM10_target_test_net1(i,:) output_test_net1(i,:)])),0.95*max(output_test_net1(i,:), [], 'all'),strcat('R^2= ',num2str(corr_test_net1(i), '%10.3f')));
    s4=text(0.02*max(max([PM10_target_test_net1(i,:) output_test_net1(i,:)])),0.85*max(output_test_net1(i,:), [], 'all'),strcat('nrmse= ',num2str(nrmse_test_net1(i), '%10.3f')));
    s=[s1 s2 s3 s4];

    if i == 1
%         text(-25,max(output_test_net1(i,:), [], 'all'),'(a)');
        title('(a)');
    else
%         text(-25,max(output_test_net1(i,:), [], 'all'),'(d)');
        title('(d)');
    end

    set(s, 'fontsize', 10);
    set(gca, 'fontsize', 10);

    nexttile;
    
    plot(PM10_target_test_net2(i,:),output_test_net2(i,:),'r*')
    grid on
    hold on
    
    plot([0 max(max([PM10_target_test_net2(i,:) output_test_net2(i,:)]))],[0 max(max([PM10_target_test_net2(i,:) output_test_net2(i,:)]))],'b--');
    ylim([0, max(output_test_net2(i,:), [], 'all')]);
    s1=xlabel('PM10 concentration [\mug/m^3]');
    s2=ylabel('CNN [\mug/m^3]');
    s3=text(0.02*max(max([PM10_target_test_net2(i,:) output_test_net2(i,:)])),0.95*max(output_test_net2(i,:), [], 'all'),strcat('R^2 = ',num2str(corr_test_net2(i), '%10.3f')));
    s4=text(0.02*max(max([PM10_target_test_net2(i,:) output_test_net2(i,:)])),0.85*max(output_test_net2(i,:), [], 'all'),strcat('nrmse = ',num2str(nrmse_test_net2(i), '%10.3f')));
    s=[s1 s2 s3 s4];

    if i == 1
%         text(-25,max(output_test_net2(i,:), [], 'all'),'(b)');
        title('(b)');
    else
%         text(-25,max(output_test_net2(i,:), [], 'all'),'(e)');
        title('(e)');
    end

    set(s, 'fontsize', 10);
    set(gca, 'fontsize', 10);

    nexttile;
    
    plot(PM10_target_test_net3(i,:),output_test_net3(i,:),'r*')
    grid on
    hold on
    
    plot([0 max(max([PM10_target_test_net3(i,:) output_test_net3(i,:)]))],[0 max(max([PM10_target_test_net3(i,:) output_test_net3(i,:)]))],'b--');
    ylim([0, max(output_test_net3(i,:), [], 'all')]);
    s1=xlabel('PM10 concentration [\mug/m^3]');
    s2=ylabel('CNN [\mug/m^3]');
    s3=text(0.02*max(max([PM10_target_test_net3(i,:) output_test_net3(i,:)])),0.95*max(output_test_net3(i,:), [], 'all'),strcat('R^2 = ',num2str(corr_test_net3(i), '%10.3f')));
    s4=text(0.02*max(max([PM10_target_test_net3(i,:) output_test_net3(i,:)])),0.85*max(output_test_net3(i,:), [], 'all'),strcat('nrmse = ',num2str(nrmse_test_net3(i), '%10.3f')));
    s=[s1 s2 s3 s4];

    if i == 1
%         text(-30,max(output_test_net3(i,:), [], 'all'),'(c)');
        title('(c)');
    else
%         text(-35,max(output_test_net3(i,:), [], 'all'),'(f)');
        title('(f)');
    end

    set(s, 'fontsize', 10);
    set(gca, 'fontsize', 10);
    
    
%     path = strcat(outFOLD_network, dm, 'Day_', num2str(i), '_R2_test.png');
%     saveas(gcf, path);
%     close gcf

end

%% PM25 plot

path = strcat(outFOLD_network1, 'PM25_stats.csv');

PM25_Lombardia = readmatrix(path);

path = strcat(outFOLD_network2, 'PM25_stats.csv');

PM25_Australia = readmatrix(path);

path = strcat(outFOLD_network3, 'PM25_stats.csv');

PM25_Poland = readmatrix(path);

% path = strcat(outFOLD_network3, 'PM25_stats.csv');
% 
% PM25_Poland = readmatrix(path);

figure;
tiledlayout(2,2);
nexttile;
plot(1:NumOutputs, [PM25_Lombardia(1,:); PM25_Australia(1,:); PM25_Poland(1,:)], '.-', 'MarkerSize', 20, 'LineWidth', 2);
grid on;
ylim([-0.5, 1]);
xticks(1:NumOutputs);
xlabel('Days'); ylabel('R^2');
title('(a)');

nexttile;
plot(1:NumOutputs, [PM25_Lombardia(2,:); PM25_Australia(2,:); PM25_Poland(2,:)], '.-', 'MarkerSize', 20, 'LineWidth', 2);
grid on;
ylim([0, 1.5]);
xticks(1:NumOutputs);
xlabel('Days'); ylabel('NRMSE');
% ylim([0.2, 1.2]);
title('(b)');

nexttile;
plot(1:NumOutputs, [PM25_Lombardia(3,:); PM25_Australia(3,:); PM25_Poland(3,:)], '.-', 'MarkerSize', 20, 'LineWidth', 2);
grid on;
ylim([0, 0.6]);
xticks(1:NumOutputs);
xlabel('Days'); ylabel('NMAE');
title('(c)');

ax = nexttile;
p_ax=ax.Position;
area([NaN NaN], NaN(2, 4));
leg = legend({'Lombardy', 'GSA', 'LSOV'});
p_leg=leg.Position;
delete(ax)
ax=axes('Position',[p_ax(1:2) 0 0]);
area([NaN NaN], NaN(2, 4));
leg = legend({'Lombardy', 'GSA', 'LSOV'});
leg.Location = 'none';
leg.Interpreter = 'latex';
leg.FontSize = 16;
ax.Visible = false;
leg.Position=p_leg;

path = 'C:\Users\ferrarilu\OneDrive - ETH Zurich\Tesi_Olgiati\Output\PM25_Stats.png';
saveas(gcf, path);












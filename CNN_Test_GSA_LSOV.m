clear
close all
clc

%% System setup
addpath('C:\Users\ferrarilu\OneDrive - ETH Zurich\Tesi_Olgiati\Input',...
    'C:\Users\ferrarilu\OneDrive - ETH Zurich\Tesi_Olgiati\Funzioni',...
    'Output'); % Cartelle contenenti dati meteorologici, sugli inquinanti, stazioni, mappa e funzioni

outFOLD_network = 'C:\Users\ferrarilu\OneDrive - ETH Zurich\Tesi_Olgiati\Output\20231113_Australia_CNN\';
outFOLD_network2 = 'C:\Users\ferrarilu\OneDrive - ETH Zurich\Tesi_Olgiati\Output\20231113_Poland_CNN\';

net_type = 1;    % 0 = LSTM Network 1 = Convolutional Neural Network 2 = ConvLSTM 3 = Feed Forward
net_type2 = 1;


switch net_type
    case 0
        nt = 'LSTM';
    case 1
        nt = 'CNN';
    case 2
        nt = 'ConvLSTM';
    case 3
        nt = 'FF';
end

NumInputs = 3;
NumOutputs = 10; 

domain = 1; % 0 = Lombardia; 1 = Australia 2 = Poland
domain2 = 2;

switch domain
    case 0
        dm='Lombardia_';
    case 1
        dm='Australia_';
    case 2
        dm='Poland_';
end

% Nstaz_val = 10; % Number of random stations used in validation
% CoordFlag = true;

%% Load data
if domain == 0
    PM10=readmatrix('PM102015_2021.csv');% Concentrazione PM10 in ug/m3
    PM10 = PM10(2:end,2:end);
    coord = readmatrix('UTM_coordinates.csv');
    lat = coord(2,:);
    long = coord(3,:);
elseif domain == 1
    PM10=readmatrix('PM102015_2021_Australia.csv');% Concentrazione PM10 in ug/m3
    PM10 = PM10(:,2:end);
    coord = readmatrix('UTM_coordinates_Australia.csv');
    lat = coord(1,:);
    long = coord(2,:);
elseif domain == 2
    PM10=readmatrix('PM102015_2021_Poland.csv');% Concentrazione PM10 in ug/m3
    PM10 = PM10(:,2:end);
    coord = readmatrix('UTM_coordinates_Poland.csv');
    lat = coord(1,:);
    long = coord(2,:);
end

%% Analyze data and autocorrelation
% for i = 1:width(PM10)
%     c(i, :) = correlogram(PM10(:,i), PM10(:,i), 30);
% end
% 
% c_av = mean(c);
% 
% figure;
% lag=0:30;
% plot(lag,c_av,'b');
% grid on;
% % title(   'Correlogram'  );
% xlabel(  'lag (delay steps) [day]'  );
% ylabel(  'correlation'  );
% 
% outFOLD = 'C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Output\';
% path = strcat(outFOLD, 'Correlogram_', dm, '.png');
% saveas(gcf, path);
% close gcf

% m{1} = repmat(mean(PM10, 2, 'omitnan'), 1, width(PM10));
% s{1} = repmat(std(PM10, 0, 2, 'omitnan'), 1, width(PM10));
% X{1} = (PM10 - m{1}) ./ s{1}; 

T=365;
for i=1:width(PM10)
    % Media
    [ ~ , m{1}(:,i)] = moving_average(PM10(:,i) , T , 1 ); 
    % Deviazione standard 
    [ ~ , s2 ] = moving_average( ( PM10(:,i) - m{1}(:,i) ).^2 , T , 1 );
    s{1}(:,i) = s2 .^ 0.5;
    % Seasonal adjustment
    X{1}(:,i) = ( PM10(:,i) - m{1}(:,i)) ./ s{1}(:,i) ; 
end

% for i = 1:width(X{1})
%     c2(i, :) = correlogram(X{1}(:,i), X{1}(:,i), 30);
% end
% 
% c_av2 = mean(c2);
% 
% figure;
% lag=0:30;
% plot(lag,c_av2,'b');
% grid on;
% % title(   'Correlogram'  );
% xlabel(  'lag (delay steps) [day]'  );
% ylabel(  'correlation'  );
% 
% path = strcat(outFOLD, 'Correlogram_', dm, 'no_seasonality.png');
% saveas(gcf, path);
% close gcf

% seasonality_ratio = mean(abs(X{1}./PM10), 2);

% figure;
% plot(seasonality_ratio);


%% Preparing calibration and validation data
Data= cell(1,NumInputs);

Data{1}=PM10;% Particolato sottile (d<10um) [ug/m3]
Data{2}=lat(:,:).'; % latitudine delle stazioni [m]
Data{3}=long(:,:).'; % longitudine delle stazioni [m]


%% Persistent model
% Prepare data
[input, target, mean_norm, std_norm, PM10_target] = input_preprocess(Data, m, s, PM10,...
    NumOutputs, false, 365, 0, 0);

for i = 1:length(input{1,1})
    input_pers(i,:) = input{1,1}{1,i};
    input_pers_val(i,:) = input{1,2}{1,i};
end

for i = 1:length(target{1,1})
    target_pers(i,:) = target{1,1}{1,i};

    mean_pers(i,:) = mean_norm{1,1}{1,i};

    std_pers(i,:) = std_norm{1,1}{1,i};

    PM10_target_pers(i,:) = PM10_target{1,1}{1,i};
    PM10_target_pers_val(i,:) = PM10_target{1,2}{1,i};
end


% Results
output_pers = repmat(input_pers_val, NumOutputs, 1);

% corr = corrcoef(output_pers, PM10_target_pers);
% corr_pers_global = corr(1,2);
% nrmse_pers_global = sqrt((sum((PM10_target_pers-output_pers).^2))/(sum((output_pers-mean(output_pers)).^2))); 

for i = 1:NumOutputs
%     corr_iter = corrcoef(output_pers(i,:), PM10_target_pers(i,:));
%     corr_pers(i) = corr_iter(1,2);
    corr_pers(i) = 1 - sum((output_pers(i,:) - PM10_target_pers_val(i,:)).^2)/(sum((mean(PM10_target_pers_val(i,:)) - PM10_target_pers_val(i,:)).^2));
    nrmse_pers(i) = sqrt((sum((PM10_target_pers_val(i,:)-output_pers(i,:)).^2))/(sum((output_pers(i,:)-mean(output_pers(i,:))).^2))); 
end

err_medio_abs_norm_pers = 1/width(output_pers)*(sum(abs(PM10_target_pers_val-output_pers),2)./mean(abs(PM10_target_pers_val),2));

%% Load first net
% outFOLD = outFOLD_network;

path = strcat(outFOLD_network, 'net_opt.mat');
load(path);

path = strcat(outFOLD_network, 'Data.mat');
load(path);

%% Evaluate model adjusted
net_opt = best_net_opt;

if net_type == 3
    output_val = sim(net_opt, input_val);
else
    output_val = predict(net_opt, input_val);
end

output_val = std_val.*output_val +mean_val ;

if net_type == 3
    output_test = sim(net_opt,input_test);
else
    output_test = predict(net_opt,input_test);
end

output_test = std_test.*output_test +mean_test ;

err_medio_abs_norm = 1/NumOutputs*(sum(abs(PM10_target_test-output_test))./mean(abs(PM10_target_test)));
mean_err_medio_abs_norm = mean(err_medio_abs_norm);


for i = 1:NumOutputs
    corr_val(i) = 1 - sum((output_val(i,:) - PM10_target_val(i,:)).^2)/(sum((mean(PM10_target_val(i,:)) - PM10_target_val(i,:)).^2));
    nrmse_val(i) = sqrt((sum((PM10_target_val(i,:)-output_val(i,:)).^2))/(sum((output_val(i,:)-mean(output_val(i,:))).^2))); 
    

    corr_test(i) = 1 - sum((output_test(i,:) - PM10_target_test(i,:)).^2)/(sum((mean(PM10_target_test(i,:)) - PM10_target_test(i,:)).^2));
    nrmse_test(i) = sqrt((sum((PM10_target_test(i,:)-output_test(i,:)).^2))/(sum((output_test(i,:)-mean(output_test(i,:))).^2))); 
    
    if i == 3 || i == 7
        createScatter(PM10_target_test(i,:),output_test(i,:),'PM10 concentration [\mug/m^3]','CNN [\mug/m^3]',corr_test(i),nrmse_test(i));
        path = strcat(outFOLD_network, dm, 'Day_', num2str(i), '_R2_test.png');
        saveas(gcf, path);
        close gcf
    end
end

err_medio_abs_norm_test1 = 1/width(output_test)*(sum(abs(PM10_target_test-output_test),2)./mean(abs(PM10_target_test),2));
err_medio_abs_norm_val = 1/width(output_val)*(sum(abs(PM10_target_val-output_val),2)./mean(abs(PM10_target_val),2));


figure;
plot(1:NumOutputs, [corr_val; corr_test; corr_pers], '.-', 'MarkerSize', 20, 'LineWidth', 2);
grid on;
xticks(1:NumOutputs);
xlabel('Days'); ylabel('R^2');
legend('Validation dataset', 'Test dataset', 'Persistent model', 'Location', 'northeastoutside');

% path = strcat(outFOLD_network, dm, 'R2_val_comparison.png');
% saveas(gcf, path);
% close gcf

figure;
plot(1:NumOutputs, [nrmse_val; nrmse_test; nrmse_pers], '.-', 'MarkerSize', 20, 'LineWidth', 2);
grid on;
xticks(1:NumOutputs);
xlabel('Days'); ylabel('nrmse');
legend('Validation dataset', 'Test dataset', 'Persistent model', 'Location', 'northeastoutside');

% path = strcat(outFOLD_network, dm, 'nrmse_val_comparison.png');
% saveas(gcf, path);
% close gcf

figure;
plot(1:NumOutputs, [err_medio_abs_norm_val, err_medio_abs_norm_test1, err_medio_abs_norm_pers], '.-', 'MarkerSize', 20, 'LineWidth', 2);
grid on;
xticks(1:NumOutputs);
xlabel('Days'); ylabel('NMAE');
legend('Validation dataset', 'Test dataset', 'Persistent model', 'Location', 'northeastoutside');

% path = strcat(outFOLD_network, dm, 'nmae_val_comparison.png');
% saveas(gcf, path);
% close gcf

figure;
tiledlayout(2,2);
nexttile;
plot(1:NumOutputs, [corr_val; corr_test; corr_pers], '.-', 'MarkerSize', 20, 'LineWidth', 2);
grid on;
xticks(1:NumOutputs);
ylabel('R^2');
% text(-2.5,1,'(a)');
title('(a)');

nexttile;
plot(1:NumOutputs, [nrmse_val; nrmse_test; nrmse_pers], '.-', 'MarkerSize', 20, 'LineWidth', 2);
grid on;
xticks(1:NumOutputs);
xlabel('Days'); ylabel('NRMSE');
% text(-2.5,1.1,'(b)');
title('(b)');

nexttile;
plot(1:NumOutputs, [err_medio_abs_norm_val, err_medio_abs_norm_test1, err_medio_abs_norm_pers], '.-', 'MarkerSize', 20, 'LineWidth', 2);
grid on;
xticks(1:NumOutputs);
xlabel('Days'); ylabel('NMAE');
% text(-2.5,0.6,'(c)');
title('(c)');

ax = nexttile;
p_ax=ax.Position;
area([NaN NaN], NaN(2, 4));
leg = legend('Validation dataset', 'Test dataset', 'Persistent model');
p_leg=leg.Position;
delete(ax)
ax=axes('Position',[p_ax(1:2) 0 0]);
area([NaN NaN], NaN(2, 4));
leg = legend('Validation dataset', 'Test dataset', 'Persistent model');
leg.Location = 'none';
leg.Interpreter = 'latex';
leg.FontSize = 16;
ax.Visible = false;
leg.Position=p_leg;

% path = strcat(outFOLD_network, dm, 'metrics_comparison.png');
% saveas(gcf, path);
% close gcf

%% Load second model
path = strcat(outFOLD_network2, 'net_opt.mat');
load(path);

path = strcat(outFOLD_network2, 'Data.mat');
load(path);

%% Load data

clear input target mean_norm std_norm PM10_target input_pers input_pers_val target_pers...
    mean_pers std_pers PM10_target_pers PM10_target_pers_val


if domain2 == 0
    PM10=readmatrix('PM102015_2021.csv');% Concentrazione PM10 in ug/m3
    PM10 = PM10(2:end,2:end);
    coord = readmatrix('UTM_coordinates.csv');
    lat = coord(2,:);
    long = coord(3,:);
elseif domain2 == 1
    PM10=readmatrix('PM102015_2021_Australia.csv');% Concentrazione PM10 in ug/m3
    PM10 = PM10(:,2:end);
    coord = readmatrix('UTM_coordinates_Australia.csv');
    lat = coord(1,:);
    long = coord(2,:);
elseif domain2 == 2
    PM10=readmatrix('PM102015_2021_Poland.csv');% Concentrazione PM10 in ug/m3
    PM10 = PM10(:,2:end);
    coord = readmatrix('UTM_coordinates_Poland.csv');
    lat = coord(1,:);
    long = coord(2,:);
end

Data= cell(1,NumInputs);

Data{1}=PM10;% Particolato sottile (d<10um) [ug/m3]
Data{2}=lat(:,:).'; % latitudine delle stazioni [m]
Data{3}=long(:,:).'; % longitudine delle stazioni [m]


%% Persistent model
% Prepare data
[input, target, mean_norm, std_norm, PM10_target] = input_preprocess(Data, m, s, PM10,...
    NumOutputs, false, 365, 0, 0);

for i = 1:length(input{1,1})
    input_pers(i,:) = input{1,1}{1,i};
    input_pers_val(i,:) = input{1,2}{1,i};
end

for i = 1:length(target{1,1})
    target_pers(i,:) = target{1,1}{1,i};

    mean_pers(i,:) = mean_norm{1,1}{1,i};

    std_pers(i,:) = std_norm{1,1}{1,i};

    PM10_target_pers(i,:) = PM10_target{1,1}{1,i};
    PM10_target_pers_val(i,:) = PM10_target{1,2}{1,i};
end


% Results
output_pers = repmat(input_pers_val, NumOutputs, 1);

for i = 1:NumOutputs
    corr_pers2(i) = 1 - sum((output_pers(i,:) - PM10_target_pers_val(i,:)).^2)/(sum((mean(PM10_target_pers_val(i,:)) - PM10_target_pers_val(i,:)).^2));
    nrmse_pers2(i) = sqrt((sum((PM10_target_pers_val(i,:)-output_pers(i,:)).^2))/(sum((output_pers(i,:)-mean(output_pers(i,:))).^2))); 
end

err_medio_abs_norm_pers2 = 1/width(output_pers)*(sum(abs(PM10_target_pers_val-output_pers),2)./mean(abs(PM10_target_pers_val),2));


%% Evaluate model adjusted

clear output_val
net_opt = best_net_opt;

if net_type2 == 3
    output_val = sim(net_opt,input_val);
else
    output_val = predict(net_opt,input_val);
end

output_val = std_val.*output_val +mean_val ;

if net_type2 == 3
    output_test = sim(net_opt,input_test);
else
    output_test = predict(net_opt,input_test);
end

output_test = std_test.*output_test +mean_test ;

err_medio_abs_norm = 1/NumOutputs*(sum(abs(PM10_target_test-output_test))./mean(abs(PM10_target_test)));
mean_err_medio_abs_norm = mean(err_medio_abs_norm);

for i = 1:NumOutputs
    corr_val2(i) = 1 - sum((output_val(i,:) - PM10_target_val(i,:)).^2)/(sum((mean(PM10_target_val(i,:)) - PM10_target_val(i,:)).^2));
    nrmse_val2(i) = sqrt((sum((PM10_target_val(i,:)-output_val(i,:)).^2))/(sum((output_val(i,:)-mean(output_val(i,:))).^2))); 
   
    corr_test2(i) = 1 - sum((output_test(i,:) - PM10_target_test(i,:)).^2)/(sum((mean(PM10_target_test(i,:)) - PM10_target_test(i,:)).^2));
    nrmse_test2(i) = sqrt((sum((PM10_target_test(i,:)-output_test(i,:)).^2))/(sum((output_test(i,:)-mean(output_test(i,:))).^2))); 
end

err_medio_abs_norm_test2 = 1/width(output_test)*(sum(abs(PM10_target_test-output_test),2)./mean(abs(PM10_target_test),2));
err_medio_abs_norm_val = 1/width(output_val)*(sum(abs(PM10_target_val-output_val),2)./mean(abs(PM10_target_val),2));

% figure;
% plot(1:NumOutputs, [corr_val2; corr_test2; corr_pers], '.-', 'MarkerSize', 20, 'LineWidth', 2);
% grid on;
% xticks(1:NumOutputs);
% xlabel('Days'); ylabel('R^2');
% legend('Validation dataset', 'Test dataset', 'Persistent model', 'Location', 'northeastoutside');
% 
% path = strcat(outFOLD_network2, dm, 'R2_val_comparison.png');
% saveas(gcf, path);
% close gcf
% 
% figure;
% plot(1:NumOutputs, [nrmse_val2; nrmse_test2; nrmse_pers], '.-', 'MarkerSize', 20, 'LineWidth', 2);
% grid on;
% xticks(1:NumOutputs);
% xlabel('Days'); ylabel('nrmse');
% legend('Validation dataset', 'Test dataset', 'Persistent model', 'Location', 'northeastoutside');
% 
% path = strcat(outFOLD_network2, dm, 'nrmse_val_comparison.png');
% saveas(gcf, path);
% close gcf
% 
% figure;
% plot(1:NumOutputs, [err_medio_abs_norm_test2, err_medio_abs_norm_val, err_medio_abs_norm_pers], '.-', 'MarkerSize', 20, 'LineWidth', 2);
% grid on;
% xticks(1:NumOutputs);
% xlabel('Days'); ylabel('NMAE');
% legend('Validation dataset', 'Test dataset', 'Persistent model', 'Location', 'northeastoutside');
% 
% path = strcat(outFOLD_network2, dm, 'nmae_val_comparison.png');
% saveas(gcf, path);
% close gcf


% figure;
% tiledlayout(2,2);
% nexttile;
% plot(1:NumOutputs, [corr_val2; corr_test2; corr_pers], '.-', 'MarkerSize', 20, 'LineWidth', 2);
% grid on;
% xticks(1:NumOutputs);
% ylabel('R^2');
% title('(a)');
% 
% nexttile;
% plot(1:NumOutputs, [nrmse_val2; nrmse_test2; nrmse_pers], '.-', 'MarkerSize', 20, 'LineWidth', 2);
% grid on;
% xticks(1:NumOutputs);
% xlabel('Days'); ylabel('NRMSE');
% title('(b)');
% 
% nexttile;
% plot(1:NumOutputs, [err_medio_abs_norm_test2, err_medio_abs_norm_val, err_medio_abs_norm_pers], '.-', 'MarkerSize', 20, 'LineWidth', 2);
% grid on;
% xticks(1:NumOutputs);
% xlabel('Days'); ylabel('NMAE');
% title('(c)');
% 
% ax = nexttile;
% p_ax=ax.Position;
% area([NaN NaN], NaN(2, 4));
% leg = legend('Validation dataset', 'Test dataset', 'Persistent model');
% p_leg=leg.Position;
% delete(ax)
% ax=axes('Position',[p_ax(1:2) 0 0]);
% area([NaN NaN], NaN(2, 4));
% leg = legend('Validation dataset', 'Test dataset', 'Persistent model');
% leg.Location = 'none';
% leg.Interpreter = 'latex';
% leg.FontSize = 16;
% ax.Visible = false;
% leg.Position=p_leg;
% 
% path = strcat(outFOLD_network2, dm, 'metrics_comparison.png');
% saveas(gcf, path);
% close gcf


%% 2 models comparison
figure;
tiledlayout(2,2);
nexttile;
% plot(1:NumOutputs, [corr_test; corr_test2; corr_pers; corr_pers2], '.-', '--', '.-', '--', 'MarkerSize', 20, 'LineWidth', 2);
plot(1:NumOutputs, [corr_test], '.-b', 'MarkerSize', 20, 'LineWidth', 2);
grid on; hold on;
plot(1:NumOutputs, [corr_test2], '.-r', 'MarkerSize', 20, 'LineWidth', 2);
plot(1:NumOutputs, [corr_pers], '--b', 'MarkerSize', 20, 'LineWidth', 2);
plot(1:NumOutputs, [corr_pers2], '--r', 'MarkerSize', 20, 'LineWidth', 2);
xticks(1:NumOutputs);
xlabel('Days'); ylabel('R^2');
title('(a)');

nexttile;
% plot(1:NumOutputs, [nrmse_test; nrmse_test2; nrmse_pers; nrmse_pers2], '.-', 'MarkerSize', 20, 'LineWidth', 2);
plot(1:NumOutputs, [nrmse_test], '.-b', 'MarkerSize', 20, 'LineWidth', 2);
grid on; hold on;
plot(1:NumOutputs, [nrmse_test2], '.-r', 'MarkerSize', 20, 'LineWidth', 2);
plot(1:NumOutputs, [nrmse_pers], '--b', 'MarkerSize', 20, 'LineWidth', 2);
plot(1:NumOutputs, [nrmse_pers2], '--r', 'MarkerSize', 20, 'LineWidth', 2);
xticks(1:NumOutputs);
xlabel('Days'); ylabel('NRMSE');
% ylim([0.2, 1.2]);
title('(b)');

nexttile;
% plot(1:NumOutputs, [err_medio_abs_norm_test1, err_medio_abs_norm_test2, err_medio_abs_norm_pers, err_medio_abs_norm_pers2], '.-', 'MarkerSize', 20, 'LineWidth', 2);
plot(1:NumOutputs, [err_medio_abs_norm_test1], '.-b', 'MarkerSize', 20, 'LineWidth', 2);
grid on; hold on;
plot(1:NumOutputs, [err_medio_abs_norm_test2], '.-r', 'MarkerSize', 20, 'LineWidth', 2);
plot(1:NumOutputs, [err_medio_abs_norm_pers], '--b', 'MarkerSize', 20, 'LineWidth', 2);
plot(1:NumOutputs, [err_medio_abs_norm_pers2], '--r', 'MarkerSize', 20, 'LineWidth', 2);
xticks(1:NumOutputs);
xlabel('Days'); ylabel('NMAE');
title('(c)');

ax = nexttile;
p_ax=ax.Position;
area([NaN NaN], NaN(2, 4));
leg = legend('CNN GSA', 'CNN LSOV', 'Persistence GSA', 'Persistence LSOV');
p_leg=leg.Position;
delete(ax)
ax=axes('Position',[p_ax(1:2) 0 0]);
area([NaN NaN], NaN(2, 4));
leg = legend('CNN GSA', 'CNN LSOV', 'Persistence GSA', 'Persistence LSOV');
leg.Location = 'none';
leg.Interpreter = 'latex';
leg.FontSize = 16;
ax.Visible = false;
leg.Position=p_leg;



































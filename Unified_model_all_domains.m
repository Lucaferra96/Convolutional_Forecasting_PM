clear
close all
clc

%% System setup
addpath('C:\Users\ferrarilu\OneDrive - ETH Zurich\Tesi_Olgiati\Input',...
    'C:\Users\ferrarilu\OneDrive - ETH Zurich\Tesi_Olgiati\Funzioni',...
    'Output'); % Cartelle contenenti dati meteorologici, sugli inquinanti, stazioni, mappa e funzioni

outFOLD_network = 'C:\Users\ferrarilu\OneDrive - ETH Zurich\Tesi_Olgiati\Output\20230613_Australia_CNN_IY_3In_10Out_20N_20Conv_Georef_v\';

net_type = 1;    % 0 = LSTM Network 1 = Convolutional Neural Network 2 = ConvLSTM 3 = Feed Forward

% switch net_type
%     case 0
%         nt = 'LSTM';
%     case 1
%         nt = 'CNN';
%     case 2
%         nt = 'ConvLSTM';
%     case 3
%         nt = 'FF';
% end

NumInputs = 3;
NumOutputs = 10; 

domain = 0; % 0 = Lombardia; 1 = Australia 2 = Poland
domain2 = 2;

switch domain
    case 0
        dm1='Lombardy';
    case 1
        dm1='GSA';
    case 2
        dm1='LSOV';
end

switch domain2
    case 0
        dm2='Lombardy';
    case 1
        dm2='GSA';
    case 2
        dm2='LSOV';
end



%% Load data I case
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

%% Preparing calibration and validation data
Data= cell(1,3);

Data{1}=PM10;% Particolato sottile (d<10um) [ug/m3]
Data{2}=lat(:,:).'; % latitudine delle stazioni [m]
Data{3}=long(:,:).'; % longitudine delle stazioni [m]

%% Seasonal adjustment
%Normalizzazione di PM10, NOx, precipitazione, velocità del vento e direzione del vento
X=cell(NumInputs,1);

s=cell(NumInputs-2,1);
m=cell(NumInputs-2,1);

% Analisi della media e della varianza
T=365; % periodicità

for i=1:width(PM10)
    % Media
    [ ~ , m{1}(:,i)] = moving_average(PM10(:,i) , T , 1 ); 
    % Deviazione standard 
    [ ~ , s2 ] = moving_average( ( PM10(:,i) - m{1}(:,i) ).^2 , T , 1 );
    s{1}(:,i) = s2 .^ 0.5;
    % Seasonal adjustment
    X{1}(:,i) = ( PM10(:,i) - m{1}(:,i)) ./ s{1}(:,i) ; 
end


X{end-1} = Data{end-1};
X{end} = Data{end};

%% Prepare data CCN

[input, target, mean_norm, std_norm, PM10_target, ~] = input_preprocess(X, m, s, PM10,...
    NumOutputs, true, 365, 0, 1);

input_test = zeros(NumInputs, length(input{1,3}{1,1}));
target_test = zeros(NumOutputs, length(target{1,3}{1,1}));
mean_test = zeros(NumOutputs, length(mean_norm{1,3}{1,1}));
std_test = zeros(NumOutputs, length(std_norm{1,3}{1,1}));
PM10_target_test = zeros(NumOutputs, length(PM10_target{1,3}{1,1}));

for i = 1:length(input{1,1})
    input_test(i,:) = input{1,3}{1,i};
end

for i = 1:length(target{1,1})
    target_test(i,:) = target{1,3}{1,i};

    mean_test(i,:) = mean_norm{1,3}{1,i};

    std_test(i,:) = std_norm{1,3}{1,i};

    PM10_target_test(i,:) = PM10_target{1,3}{1,i};
end

%% Load net
path = strcat(outFOLD_network, 'net_opt.mat');
load(path);

%% Test net

if net_type == 3
    output_test = sim(net_opt,input_test);
else
    output_test = predict(net_opt,input_test);
end

output_test = std_test.*output_test +mean_test ;

% err_medio_abs_norm = 1/NumOutputs*(sum(abs(PM10_target_test-output_test))./mean(abs(PM10_target_test)));
% mean_err_medio_abs_norm = mean(err_medio_abs_norm);


for i = 1:NumOutputs
    corr_test_Australia(i) = 1 - sum((output_test(i,:) - PM10_target_test(i,:)).^2)/(sum((mean(PM10_target_test(i,:)) - PM10_target_test(i,:)).^2));
    nrmse_test_Australia(i) = sqrt((sum((PM10_target_test(i,:)-output_test(i,:)).^2))/(sum((output_test(i,:)-mean(output_test(i,:))).^2))); 
end

err_medio_abs_norm_test_Australia = 1/width(output_test)*(sum(abs(PM10_target_test-output_test),2)./mean(abs(PM10_target_test),2));

%% Load data II case
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

%% Preparing calibration and validation data
Data= cell(1,3);

Data{1}=PM10;% Particolato sottile (d<10um) [ug/m3]
Data{2}=lat(:,:).'; % latitudine delle stazioni [m]
Data{3}=long(:,:).'; % longitudine delle stazioni [m]

%% Seasonal adjustment
%Normalizzazione di PM10, NOx, precipitazione, velocità del vento e direzione del vento
X=cell(NumInputs,1);

s=cell(NumInputs-2,1);
m=cell(NumInputs-2,1);

% Analisi della media e della varianza
T=365; % periodicità

for i=1:width(PM10)
    % Media
    [ ~ , m{1}(:,i)] = moving_average(PM10(:,i) , T , 1 ); 
    % Deviazione standard 
    [ ~ , s2 ] = moving_average( ( PM10(:,i) - m{1}(:,i) ).^2 , T , 1 );
    s{1}(:,i) = s2 .^ 0.5;
    % Seasonal adjustment
    X{1}(:,i) = ( PM10(:,i) - m{1}(:,i)) ./ s{1}(:,i) ; 
end


X{end-1} = Data{end-1};
X{end} = Data{end};

%% Prepare data CCN
[input, target, mean_norm, std_norm, PM10_target, ~] = input_preprocess(X, m, s, PM10,...
    NumOutputs, true, 365, 0, 1);

input_test = zeros(NumInputs, length(input{1,3}{1,1}));
target_test = zeros(NumOutputs, length(target{1,3}{1,1}));
mean_test = zeros(NumOutputs, length(mean_norm{1,3}{1,1}));
std_test = zeros(NumOutputs, length(std_norm{1,3}{1,1}));
PM10_target_test = zeros(NumOutputs, length(PM10_target{1,3}{1,1}));

for i = 1:length(input{1,1})
    input_test(i,:) = input{1,3}{1,i};
end

for i = 1:length(target{1,1})
    target_test(i,:) = target{1,3}{1,i};

    mean_test(i,:) = mean_norm{1,3}{1,i};

    std_test(i,:) = std_norm{1,3}{1,i};

    PM10_target_test(i,:) = PM10_target{1,3}{1,i};
end

%% Test second net

if net_type == 3
    output_test = sim(net_opt,input_test);
else
    output_test = predict(net_opt,input_test);
end

output_test = std_test.*output_test +mean_test ;

% err_medio_abs_norm = 1/NumOutputs*(sum(abs(PM10_target_test-output_test))./mean(abs(PM10_target_test)));
% mean_err_medio_abs_norm = mean(err_medio_abs_norm);

for i = 1:NumOutputs
    corr_test_Poland(i) = 1 - sum((output_test(i,:) - PM10_target_test(i,:)).^2)/(sum((mean(PM10_target_test(i,:)) - PM10_target_test(i,:)).^2));
    nrmse_test_Poland(i) = sqrt((sum((PM10_target_test(i,:)-output_test(i,:)).^2))/(sum((output_test(i,:)-mean(output_test(i,:))).^2))); 
end

err_medio_abs_norm_test_Poland = 1/width(output_test)*(sum(abs(PM10_target_test-output_test),2)./mean(abs(PM10_target_test),2));

%% Domain comparison
figure;
tiledlayout(2,2);
nexttile;
plot(1:NumOutputs, [corr_test_Australia], '.-b', 'MarkerSize', 20, 'LineWidth', 2);
grid on; hold on;
plot(1:NumOutputs, [corr_test_Poland], '.-r', 'MarkerSize', 20, 'LineWidth', 2);
xticks(1:NumOutputs);
ylim([-0.5, 1]);
xlabel('Days'); ylabel('R^2');
title('(a)');

nexttile;
plot(1:NumOutputs, [nrmse_test_Australia], '.-b', 'MarkerSize', 20, 'LineWidth', 2);
grid on; hold on;
plot(1:NumOutputs, [nrmse_test_Poland], '.-r', 'MarkerSize', 20, 'LineWidth', 2);
xticks(1:NumOutputs);
ylim([0, 1.5]);
xlabel('Days'); ylabel('NRMSE');
% ylim([0.2, 1.2]);
title('(b)');

nexttile;
plot(1:NumOutputs, [err_medio_abs_norm_test_Australia], '.-b', 'MarkerSize', 20, 'LineWidth', 2);
grid on; hold on;
plot(1:NumOutputs, [err_medio_abs_norm_test_Poland], '.-r', 'MarkerSize', 20, 'LineWidth', 2);
xticks(1:NumOutputs);
ylim([0, 0.6]);
xlabel('Days'); ylabel('NMAE');
title('(c)');

ax = nexttile;
p_ax=ax.Position;
area([NaN NaN], NaN(2, 4));
leg = legend(dm1, dm2);
p_leg=leg.Position;
delete(ax)
ax=axes('Position',[p_ax(1:2) 0 0]);
area([NaN NaN], NaN(2, 4));
leg = legend(dm1, dm2);
leg.Location = 'none';
leg.Interpreter = 'latex';
leg.FontSize = 16;
ax.Visible = false;
leg.Position=p_leg;























































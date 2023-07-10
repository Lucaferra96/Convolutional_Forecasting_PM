clear
close all
clc

%% System setup
addpath('C:\Users\ferrarilu\OneDrive - ETH Zurich\Tesi_Olgiati\Input',...
    'C:\Users\ferrarilu\OneDrive - ETH Zurich\Tesi_Olgiati\Funzioni',...
    'Output'); % Cartelle contenenti dati meteorologici, sugli inquinanti, stazioni, mappa e funzioni

outFOLD = 'C:\Users\ferrarilu\OneDrive - ETH Zurich\Tesi_Olgiati\Output\20230613_Australia_CNN_IY_3In_10Out_20N_20Conv_Georef_v\';
path = strcat(outFOLD, 'net_opt.mat');
load(path);

net_type = 1;    % 0 = LSTM Network 1 = Convolutional Neural Network 2 = ConvLSTM 3 = Feed Forward

ValSelection = 1;

domain = 1; % 0 = Lombardia; 1 = Australia 2 = Poland
Nstaz_val = 0;

Tval = 365;

NumInputs = 3;
NumOutputs = 10; 

switch domain
    case 0
        dm='Lombardia_';
    case 1
        dm='Australia_';
    case 2
        dm='Poland_';
end

% Nstaz_val = 10; % Number of random stations used in validation
CoordFlag = true;

%% Load data
if domain == 0
    PM25=readmatrix('PM252015_2021.csv');% Concentrazione PM10 in ug/m3
    PM25 = PM25(2:end,2:end);
    coord = readmatrix('UTM_coordinates_PM25_Lombardia.csv');
    lat = coord(2,:);
    long = coord(3,:);
elseif domain == 1
    PM25=readmatrix('PM252015_2021_Australia.csv');% Concentrazione PM10 in ug/m3
    PM25 = PM25(:,2:end);
    coord = readmatrix('UTM_coordinates_Australia.csv');
    lat = coord(1,:);
    long = coord(2,:);
elseif domain == 2
    PM25=readmatrix('PM252015_2021_Poland.csv');% Concentrazione PM10 in ug/m3
    PM25 = PM25(:,2:end);
    coord = readmatrix('UTM_coordinates_PM25_Poland.csv');
    lat = coord(1,:);
    long = coord(2,:);
end

%% Preparing calibration and validation data

Data= cell(1,3);

Data{1}=PM25;% Particolato sottile (d<10um) [ug/m3]
Data{2}=lat(:,:).'; % latitudine delle stazioni [m]
Data{3}=long(:,:).'; % longitudine delle stazioni [m]

%% Seasonal adjustment
%Normalizzazione di PM10, NOx, precipitazione, velocità del vento e direzione del vento
X=cell(NumInputs,1);

s=cell(NumInputs-2,1);
m=cell(NumInputs-2,1);

% Analisi della media e della varianza
T=365; % periodicità

for i=1:width(PM25)
    % Media
    [ ~ , m{1}(:,i)] = moving_average(PM25(:,i) , T , 1 ); 
    % Deviazione standard 
    [ ~ , s2 ] = moving_average( ( PM25(:,i) - m{1}(:,i) ).^2 , T , 1 );
    s{1}(:,i) = s2 .^ 0.5;
    % Seasonal adjustment
    X{1}(:,i) = ( PM25(:,i) - m{1}(:,i)) ./ s{1}(:,i) ; 
end

X{end-1} = Data{end-1};
X{end} = Data{end};

%% Prepare data persistent model
[input, target, mean_norm, std_norm, PM25_target] = input_preprocess(Data, m, s, PM25,...
    NumOutputs, false, 0, 0, 0);

for i = 1:length(input{1,1})
    input_pers(i,:) = input{1,1}{1,i};
end

for i = 1:length(target{1,1})
    target_pers(i,:) = target{1,1}{1,i};

    mean_pers(i,:) = mean_norm{1,1}{1,i};

    std_pers(i,:) = std_norm{1,1}{1,i};

    PM25_target_pers(i,:) = PM25_target{1,1}{1,i};
end

clear input_id target_id mean_id std_id PM25_target_id input_val target_val mean_val std_val PM25_target_val

%% Results
output_pers = repmat(input_pers, NumOutputs, 1);

corr = corrcoef(output_pers, PM25_target_pers);
corr_pers_global = corr(1,2);
nrmse_pers_global = sqrt((sum((PM25_target_pers-output_pers).^2))/(sum((output_pers-mean(output_pers)).^2))); 

for i = 1:NumOutputs
    corr_pers(i) = 1 - sum((output_pers(i,:) - PM25_target_pers(i,:)).^2)/(sum((mean(PM25_target_pers(i,:)) - PM25_target_pers(i,:)).^2));
    nrmse_pers(i) = sqrt((sum((PM25_target_pers(i,:)-output_pers(i,:)).^2))/(sum((output_pers(i,:)-mean(output_pers(i,:))).^2))); 
end

err_medio_abs_norm_pers = 1/width(output_pers)*(sum(abs(PM25_target_pers-output_pers),2)./mean(abs(PM25_target_pers),2));

%% Prepare data CCN
[input, target, mean_norm, std_norm, PM25_target, test_stations] = input_preprocess(X, m, s, PM25,...
    NumOutputs, CoordFlag, Tval, Nstaz_val, ValSelection);


for i = 1:length(input{1,1})
    input_id(i,:) = input{1,3}{1,i};
%     input_val(i,:) = input{1,2}{1,i};
%     input_test(i,:) = input{1,3}{1,i};
end

for i = 1:length(target{1,1})
    target_id(i,:) = target{1,3}{1,i};
%     target_val(i,:) = target{1,2}{1,i};
%     target_test(i,:) = target{1,3}{1,i};

    mean_id(i,:) = mean_norm{1,3}{1,i};
%     mean_val(i,:) = mean_norm{1,2}{1,i};
%     mean_test(i,:) = mean_norm{1,3}{1,i};

    std_id(i,:) = std_norm{1,3}{1,i};
%     std_val(i,:) = std_norm{1,2}{1,i};
%     std_test(i,:) = std_norm{1,3}{1,i};

    PM25_target_id(i,:) = PM25_target{1,3}{1,i};
%     PM25_target_val(i,:) = PM25_target{1,2}{1,i};
%     PM25_target_test(i,:) = PM25_target{1,3}{1,i};
end

%% Test network
if net_type == 3
    output_id = sim(net_opt,input_id);
else
    output_id = predict(net_opt,input_id);
end

output_id = std_id.*output_id +mean_id ;

corr = 1 - sum((output_id - PM25_target_id).^2)/(sum((mean(PM25_target_id) - PM25_target_id).^2));
nrmse = sqrt((sum((PM25_target_id-output_id).^2))/(sum((output_id-mean(output_id)).^2))); 

% createScatter(PM25_target_id,output_id,'PM25 concentration [\mug/m^3]','CNN output [\mug/m^3]',corr,nrmse);

for i = 1:NumOutputs
    corr_id(i) = 1 - sum((output_id(i,:) - PM25_target_id(i,:)).^2)/(sum((mean(PM25_target_id(i,:)) - PM25_target_id(i,:)).^2));
    nrmse_id(i) = sqrt((sum((PM25_target_id(i,:)-output_id(i,:)).^2))/(sum((output_id(i,:)-mean(output_id(i,:))).^2))); 
    
%     createScatter(PM10_target_T(i,:),output_T_val(i,:),'PM25 concentration [\mug/m^3]','CNN output [\mug/m^3]',corr_T(i),nrmse_T(i));
end

err_medio_abs_norm_id = 1/width(output_id)*(sum(abs(PM25_target_id-output_id),2)./mean(abs(PM25_target_id),2));

figure;
tiledlayout(2,2);
nexttile;
plot(1:NumOutputs, [corr_id; corr_pers], '.-', 'MarkerSize', 20, 'LineWidth', 2);
grid on;
xticks(1:NumOutputs);
ylabel('R^2');
text(-2.5,1,'(a)');

nexttile;
plot(1:NumOutputs, [nrmse_id; nrmse_pers], '.-', 'MarkerSize', 20, 'LineWidth', 2);
grid on;
xticks(1:NumOutputs);
xlabel('Days'); ylabel('NRMSE');
text(-2.5,1.1,'(b)');

nexttile;
plot(1:NumOutputs, [err_medio_abs_norm_id, err_medio_abs_norm_pers], '.-', 'MarkerSize', 20, 'LineWidth', 2);
grid on;
xticks(1:NumOutputs);
xlabel('Days'); ylabel('NMAE');
text(-2.5,0.6,'(c)');

ax = nexttile;
p_ax=ax.Position;
area([NaN NaN], NaN(2, 4));
leg = legend('Test dataset', 'Persistent model');
p_leg=leg.Position;
delete(ax)
ax=axes('Position',[p_ax(1:2) 0 0]);
area([NaN NaN], NaN(2, 4));
leg = legend('Test dataset', 'Persistent model');
leg.Location = 'none';
leg.Interpreter = 'latex';
leg.FontSize = 16;
ax.Visible = false;
leg.Position=p_leg;

path = strcat(outFOLD, dm, 'PM25_metrics_comparison.png');
saveas(gcf, path);
% close gcf

path = strcat(outFOLD, 'PM25_stats.csv');
writematrix([corr_id; nrmse_id; err_medio_abs_norm_id'], path);









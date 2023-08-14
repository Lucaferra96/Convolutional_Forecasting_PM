clear
close all
clc

%% System setup
addpath('C:\Users\ferrarilu\OneDrive - ETH Zurich\Tesi_Olgiati\Input',...
    'C:\Users\ferrarilu\OneDrive - ETH Zurich\Tesi_Olgiati\Funzioni',...
    'Output'); % Cartelle contenenti dati meteorologici, sugli inquinanti, stazioni, mappa e funzioni

outFOLD_network = 'C:\Users\ferrarilu\OneDrive - ETH Zurich\Tesi_Olgiati\Output\20230613_Lombardia_CNN_IY_3In_10Out_20N_20Conv_Georef_v\';

%% System setup
%%%%%%%%%%%%%%%%%%%% TEST DATA SELECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(0=independent stations 1=independent year
ValSelection = 1;

domain = 0; % 0 = Lombardia; 1 = Australia 2 = Poland

Tval = 365;

Nstaz_val = 1;

%%%%%%%%%%%%%%%%%%%% NETWORK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define the type of neural network
NumInputs = 3; % numero di dati in input: PM10 (giorno precedente),NOx (giorno precedente), velocità del vento, precipitazione, latidudine e longitudine (UTM est)
NumOutputs = 10; % dati in uscita: PM10 giorno successivo
net_type = 1;    % 0 = LSTM Network 1 = Convolutional Neural Network 2 = ConvLSTM 3 = Feed Forward

CoordFlag = true; % parameter that determines if coordinates are included into the network inputs

%% Select individual station

ind_station = [11, 45, 54];

%% First station
% Import meteorological, emissions and monitorig stations data
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

PM10 = PM10(:,ind_station(1));
lat = lat(ind_station(1));
long = long(ind_station(1));


% Preparing calibration and validation data
Data= cell(1,3);

Data{1}=PM10;% Particolato sottile (d<10um) [ug/m3]
Data{2}=lat(:,:).'; % latitudine delle stazioni [m]
Data{3}=long(:,:).'; % longitudine delle stazioni [m]

% Seasonal adjustment
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

% figure;
% plot(m{1}(:,1));


X{end-1} = Data{end-1};
X{end} = Data{end};
% 
% X{1} = X{1}(:,1);
% X{2} = X{2}(1,1);
% X{3} = X{3}(1,1);
% 
% m{1} = m{1}(:,1);
% 
% s{1} = s{1}(:,1);
% 
% PM10 = PM10(:,1);


% Prepare data CCN
[input, target, mean_norm, std_norm, PM10_target, test_stations] = input_preprocess(X, m, s, PM10,...
    NumOutputs, CoordFlag, Tval, Nstaz_val, ValSelection);


for i = 1:length(input{1,1})
    input_id(i,:) = input{1,1}{1,i};
    input_val(i,:) = input{1,2}{1,i};
    input_test(i,:) = input{1,3}{1,i};
end

for i = 1:length(target{1,1})
    target_id(i,:) = target{1,1}{1,i};
    target_val(i,:) = target{1,2}{1,i};
    target_test(i,:) = target{1,3}{1,i};

    mean_id(i,:) = mean_norm{1,1}{1,i};
    mean_val(i,:) = mean_norm{1,2}{1,i};
    mean_test(i,:) = mean_norm{1,3}{1,i};

    std_id(i,:) = std_norm{1,1}{1,i};
    std_val(i,:) = std_norm{1,2}{1,i};
    std_test(i,:) = std_norm{1,3}{1,i};

    PM10_target_id(i,:) = PM10_target{1,1}{1,i};
    PM10_target_val(i,:) = PM10_target{1,2}{1,i};
    PM10_target_test1(i,:) = PM10_target{1,3}{1,i};
end

% Load first net
% outFOLD = outFOLD_network;

path = strcat(outFOLD_network, 'net_opt.mat');
load(path);

% Evaluate model adjusted
if net_type == 3
    output_test = sim(net_opt,input_test);
else
    output_test = predict(net_opt,input_test);
end

output_test_net1 = std_test.*output_test +mean_test ;

for i = 1:NumOutputs
    corr_test_net1(i) = 1 - sum((output_test_net1(i,:) - PM10_target_test1(i,:)).^2)/(sum((mean(PM10_target_test1(i,:)) - PM10_target_test1(i,:)).^2));
    nrmse_test_net1(i) = sqrt((sum((PM10_target_test1(i,:)-output_test_net1(i,:)).^2))/(sum((output_test_net1(i,:)-mean(output_test_net1(i,:))).^2))); 
    
end

PM10_target_test1 = [PM10_target_test1(3,:); PM10_target_test1(7,:)];
output_test_net1 = [output_test_net1(3,:); output_test_net1(7,:)];
corr_test_net1 = [corr_test_net1(3); corr_test_net1(7)];
nrmse_test_net1 = [nrmse_test_net1(3); nrmse_test_net1(7)];

[max_err_1, I] = max((output_test_net1 - PM10_target_test1), [], 2);


%% Second station
% Import meteorological, emissions and monitorig stations data
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

PM10 = PM10(:,ind_station(2));
lat = lat(ind_station(2));
long = long(ind_station(2));


% Preparing calibration and validation data
Data= cell(1,3);

Data{1}=PM10;% Particolato sottile (d<10um) [ug/m3]
Data{2}=lat(:,:).'; % latitudine delle stazioni [m]
Data{3}=long(:,:).'; % longitudine delle stazioni [m]

% Seasonal adjustment
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

% Prepare data CCN
[input, target, mean_norm, std_norm, PM10_target, test_stations] = input_preprocess(X, m, s, PM10,...
    NumOutputs, CoordFlag, Tval, Nstaz_val, ValSelection);


for i = 1:length(input{1,1})
    input_id(i,:) = input{1,1}{1,i};
    input_val(i,:) = input{1,2}{1,i};
    input_test(i,:) = input{1,3}{1,i};
end

for i = 1:length(target{1,1})
    target_id(i,:) = target{1,1}{1,i};
    target_val(i,:) = target{1,2}{1,i};
    target_test(i,:) = target{1,3}{1,i};

    mean_id(i,:) = mean_norm{1,1}{1,i};
    mean_val(i,:) = mean_norm{1,2}{1,i};
    mean_test(i,:) = mean_norm{1,3}{1,i};

    std_id(i,:) = std_norm{1,1}{1,i};
    std_val(i,:) = std_norm{1,2}{1,i};
    std_test(i,:) = std_norm{1,3}{1,i};

    PM10_target_id(i,:) = PM10_target{1,1}{1,i};
    PM10_target_val(i,:) = PM10_target{1,2}{1,i};
    PM10_target_test2(i,:) = PM10_target{1,3}{1,i};
end


% Evaluate model adjusted
if net_type == 3
    output_test = sim(net_opt,input_test);
else
    output_test = predict(net_opt,input_test);
end

output_test_net2 = std_test.*output_test +mean_test ;

max_err_2 = max(abs(output_test_net2 - PM10_target_test2), [], 2);


for i = 1:NumOutputs
    corr_test_net2(i) = 1 - sum((output_test_net2(i,:) - PM10_target_test2(i,:)).^2)/(sum((mean(PM10_target_test2(i,:)) - PM10_target_test2(i,:)).^2));
    nrmse_test_net2(i) = sqrt((sum((PM10_target_test2(i,:)-output_test_net2(i,:)).^2))/(sum((output_test_net2(i,:)-mean(output_test_net2(i,:))).^2))); 
end

PM10_target_test2 = [PM10_target_test2(3,:); PM10_target_test2(7,:)];
output_test_net2 = [output_test_net2(3,:); output_test_net2(7,:)];
corr_test_net2 = [corr_test_net2(3); corr_test_net2(7)];
nrmse_test_net2 = [nrmse_test_net2(3); nrmse_test_net2(7)];

%% Third station
% Import meteorological, emissions and monitorig stations data
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

PM10 = PM10(:,ind_station(3));
lat = lat(ind_station(3));
long = long(ind_station(3));


% Preparing calibration and validation data
Data= cell(1,3);

Data{1}=PM10;% Particolato sottile (d<10um) [ug/m3]
Data{2}=lat(:,:).'; % latitudine delle stazioni [m]
Data{3}=long(:,:).'; % longitudine delle stazioni [m]

% Seasonal adjustment
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

% Prepare data CCN
[input, target, mean_norm, std_norm, PM10_target, test_stations] = input_preprocess(X, m, s, PM10,...
    NumOutputs, CoordFlag, Tval, Nstaz_val, ValSelection);


for i = 1:length(input{1,1})
    input_id(i,:) = input{1,1}{1,i};
    input_val(i,:) = input{1,2}{1,i};
    input_test(i,:) = input{1,3}{1,i};
end

for i = 1:length(target{1,1})
    target_id(i,:) = target{1,1}{1,i};
    target_val(i,:) = target{1,2}{1,i};
    target_test(i,:) = target{1,3}{1,i};

    mean_id(i,:) = mean_norm{1,1}{1,i};
    mean_val(i,:) = mean_norm{1,2}{1,i};
    mean_test(i,:) = mean_norm{1,3}{1,i};

    std_id(i,:) = std_norm{1,1}{1,i};
    std_val(i,:) = std_norm{1,2}{1,i};
    std_test(i,:) = std_norm{1,3}{1,i};

    PM10_target_id(i,:) = PM10_target{1,1}{1,i};
    PM10_target_val(i,:) = PM10_target{1,2}{1,i};
    PM10_target_test3(i,:) = PM10_target{1,3}{1,i};
end


% Evaluate model adjusted
if net_type == 3
    output_test = sim(net_opt,input_test);
else
    output_test = predict(net_opt,input_test);
end

output_test_net3 = std_test.*output_test +mean_test ;

max_err_3 = max(abs(output_test_net3 - PM10_target_test3), [], 2);

for i = 1:NumOutputs
    corr_test_net3(i) = 1 - sum((output_test_net3(i,:) - PM10_target_test3(i,:)).^2)/(sum((mean(PM10_target_test3(i,:)) - PM10_target_test3(i,:)).^2));
    nrmse_test_net3(i) = sqrt((sum((PM10_target_test3(i,:)-output_test_net3(i,:)).^2))/(sum((output_test_net3(i,:)-mean(output_test_net3(i,:))).^2))); 
end

PM10_target_test3 = [PM10_target_test3(3,:); PM10_target_test3(7,:)];
output_test_net3 = [output_test_net3(3,:); output_test_net3(7,:)];
corr_test_net3 = [corr_test_net3(3); corr_test_net3(7)];
nrmse_test_net3 = [nrmse_test_net3(3); nrmse_test_net3(7)];

%% Plots

figure;
tiledlayout(2,3);

for i = 1:2
    nexttile;

    plot(PM10_target_test1(i,:),output_test_net1(i,:),'r*')
    grid on
    hold on
    
    plot([0 max(max([PM10_target_test1(i,:) output_test_net1(i,:)]))],[0 max(max([PM10_target_test1(i,:) output_test_net1(i,:)]))],'b--');
    ylim([0, max(output_test_net1(i,:), [], 'all')]); 
    s1=xlabel('PM10 concentration [\mug/m^3]');
    s2=ylabel('CNN [\mug/m^3]');
    s3=text(0.02*max(max([PM10_target_test1(i,:) output_test_net1(i,:)])),0.95*max(output_test_net1(i,:), [], 'all'),strcat('R^2= ',num2str(corr_test_net1(i), '%10.3f')));
    s4=text(0.02*max(max([PM10_target_test1(i,:) output_test_net1(i,:)])),0.85*max(output_test_net1(i,:), [], 'all'),strcat('nrmse= ',num2str(nrmse_test_net1(i), '%10.3f')));
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
    
    plot(PM10_target_test2(i,:),output_test_net2(i,:),'r*')
    grid on
    hold on
    
    plot([0 max(max([PM10_target_test2(i,:) output_test_net2(i,:)]))],[0 max(max([PM10_target_test2(i,:) output_test_net2(i,:)]))],'b--');
    ylim([0, max(output_test_net2(i,:), [], 'all')]);
    s1=xlabel('PM10 concentration [\mug/m^3]');
    s2=ylabel('CNN [\mug/m^3]');
    s3=text(0.02*max(max([PM10_target_test2(i,:) output_test_net2(i,:)])),0.95*max(output_test_net2(i,:), [], 'all'),strcat('R^2 = ',num2str(corr_test_net2(i), '%10.3f')));
    s4=text(0.02*max(max([PM10_target_test2(i,:) output_test_net2(i,:)])),0.85*max(output_test_net2(i,:), [], 'all'),strcat('nrmse = ',num2str(nrmse_test_net2(i), '%10.3f')));
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
    
    plot(PM10_target_test3(i,:),output_test_net3(i,:),'r*')
    grid on
    hold on
    
    plot([0 max(max([PM10_target_test3(i,:) output_test_net3(i,:)]))],[0 max(max([PM10_target_test3(i,:) output_test_net3(i,:)]))],'b--');
    ylim([0, max(output_test_net3(i,:), [], 'all')]);
    s1=xlabel('PM10 concentration [\mug/m^3]');
    s2=ylabel('CNN [\mug/m^3]');
    s3=text(0.02*max(max([PM10_target_test3(i,:) output_test_net3(i,:)])),0.95*max(output_test_net3(i,:), [], 'all'),strcat('R^2 = ',num2str(corr_test_net3(i), '%10.3f')));
    s4=text(0.02*max(max([PM10_target_test3(i,:) output_test_net3(i,:)])),0.85*max(output_test_net3(i,:), [], 'all'),strcat('nrmse = ',num2str(nrmse_test_net3(i), '%10.3f')));
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

%% Errors
% Import meteorological, emissions and monitorig stations data
if domain == 0
    PM10_original=readmatrix('PM102015_2021.csv');% Concentrazione PM10 in ug/m3
    PM10_original = PM10_original(2:end,2:end);
    coord_original = readmatrix('UTM_coordinates.csv');
    lat_original = coord_original(2,:);
    long_original = coord_original(3,:);
elseif domain == 1
    PM10_original=readmatrix('PM102015_2021_Australia.csv');% Concentrazione PM10 in ug/m3
    PM10_original = PM10_original(:,2:end);
    coord_original = readmatrix('UTM_coordinates_Australia.csv');
    lat_original = coord_original(1,:);
    long_original = coord_original(2,:);
elseif domain == 2
    PM10_original=readmatrix('PM102015_2021_Poland.csv');% Concentrazione PM10 in ug/m3
    PM10_original = PM10_original(:,2:end);
    coord_original = readmatrix('UTM_coordinates_Poland.csv');
    lat_original = coord_original(1,:);
    long_original = coord_original(2,:);
end


for NS = 1:width(PM10_original)

    PM10 = PM10_original(:,NS);
    lat = lat_original(NS);
    long = long_original(NS);
    
    
    
    % Preparing calibration and validation data
    Data= cell(1,3);
    
    Data{1}=PM10;% Particolato sottile (d<10um) [ug/m3]
    Data{2}=lat(:,:).'; % latitudine delle stazioni [m]
    Data{3}=long(:,:).'; % longitudine delle stazioni [m]
    
    % Seasonal adjustment
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
    
    % figure;
    % plot(m{1}(:,1));
    
    
    X{end-1} = Data{end-1};
    X{end} = Data{end};
    % 
    % X{1} = X{1}(:,1);
    % X{2} = X{2}(1,1);
    % X{3} = X{3}(1,1);
    % 
    % m{1} = m{1}(:,1);
    % 
    % s{1} = s{1}(:,1);
    % 
    % PM10 = PM10(:,1);
    
    
    % Prepare data CCN
    [input, target, mean_norm, std_norm, PM10_target, test_stations] = input_preprocess(X, m, s, PM10,...
        NumOutputs, CoordFlag, Tval, Nstaz_val, ValSelection);
    
    
    for i = 1:length(input{1,1})
        input_id(i,:) = input{1,1}{1,i};
        input_val(i,:) = input{1,2}{1,i};
        input_test(i,:) = input{1,3}{1,i};
    end
    
    for i = 1:length(target{1,1})
        target_id(i,:) = target{1,1}{1,i};
        target_val(i,:) = target{1,2}{1,i};
        target_test(i,:) = target{1,3}{1,i};
    
        mean_id(i,:) = mean_norm{1,1}{1,i};
        mean_val(i,:) = mean_norm{1,2}{1,i};
        mean_test(i,:) = mean_norm{1,3}{1,i};
    
        std_id(i,:) = std_norm{1,1}{1,i};
        std_val(i,:) = std_norm{1,2}{1,i};
        std_test(i,:) = std_norm{1,3}{1,i};
    
        PM10_target_id(i,:) = PM10_target{1,1}{1,i};
        PM10_target_val(i,:) = PM10_target{1,2}{1,i};
        PM10_target_test(i,:) = PM10_target{1,3}{1,i};
    end
    
    path = strcat(outFOLD_network, 'net_opt.mat');
    load(path);
    
    % Evaluate model adjusted
    if net_type == 3
        output_test = sim(net_opt,input_test);
    else
        output_test = predict(net_opt,input_test);
    end
    
    output_test_net = std_test.*output_test +mean_test ;
    
    for i = 1:NumOutputs
        corr_test_net(i) = 1 - sum((output_test_net(i,:) - PM10_target_test(i,:)).^2)/(sum((mean(PM10_target_test(i,:)) - PM10_target_test(i,:)).^2));
        nrmse_test_net(i,NS) = sqrt((sum((PM10_target_test(i,:)-output_test_net(i,:)).^2))/(sum((output_test_net(i,:)-mean(output_test_net(i,:))).^2))); 
        
    end

    err_medio_abs_norm_pers(:,NS) = 1/width(output_test_net)*(sum(abs(PM10_target_test-output_test_net),2)./mean(abs(PM10_target_test),2));

    
    max_err(:,NS) = max(abs(output_test_net - PM10_target_test), [], 2);
    
    clear input target mean_norm std_norm PM10_target test_stations
end


NMAE_m1 = mean(err_medio_abs_norm_pers);
NMAE_m2 = mean(err_medio_abs_norm_pers(1:3,:));

h = 2 * iqr(NMAE_m1) * width(NMAE_m1)^(-1/3);
bins = round((max(NMAE_m1) - min(NMAE_m1))/h);

figure;
h1 = histogram(NMAE_m1, bins);
xlabel('NMAE');
ylabel('Instances');
hold on;
h2 = histogram(NMAE_m2, bins);
legend('10 days horizon', '3 days horizon');
h1.FaceAlpha = 0.4;
h2.FaceAlpha = 0.4;


% figure;
% boxchart(err_medio_abs_norm_pers);
% xlabel('Station number');
% ylabel('NRMSE');
% grid on;
















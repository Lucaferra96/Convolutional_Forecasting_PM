close all
clear
clc

global_start = tic;

%% System setup
addpath('C:\Users\ferrarilu\OneDrive - ETH Zurich\Tesi_Olgiati\Input',...
    'C:\Users\ferrarilu\OneDrive - ETH Zurich\Tesi_Olgiati\Funzioni',...
    'Output'); % Cartelle contenenti dati meteorologici, sugli inquinanti, stazioni, mappa e funzioni

%%%%%%%%%%%%%%%%%%%% TEST DATA SELECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(0=independent stations 1=independent year
ValSelection = 1;

domain = 0; % 0 = Lombardia; 1 = Australia 2 = Poland

Tval = 365;

switch domain
    case 0
        Nstaz_val = 10; % Number of random stations used in validation
    case 1
        Nstaz_val = 6; % Number of random stations used in validation
    case 2
        Nstaz_val = 4; % Number of random stations used in validation
end

%%%%%%%%%%%%%%%%%%%% NETWORK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define the type of neural network
NumInputs = 3; % numero di dati in input: PM10 (giorno precedente),NOx (giorno precedente), velocità del vento, precipitazione, latidudine e longitudine (UTM est)
NumOutputs = 10; % dati in uscita: PM10 giorno successivo
net_type = 3;    % 0 = LSTM Network 1 = Convolutional Neural Network 2 = ConvLSTM 3 = Feed Forward

if net_type == 3
    neurons = [64, 32]; % Number of neurons for the FFNN, 
else
    neurons = 20;
end

filterSize = 20;

numHiddenUnits=15; % parametro che rappresenta il numero di collegamenti tra le variabili (all'aumentare del parametro aumenta la rappresentazione dei picchi di concentrazione)

CoordFlag = true; % parameter that determines if coordinates are included into the network inputs

%%%%%%%%%%%%%%%%%%%% EPOCHS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define the pre-ides epochs
PreEpochs=30; %default 30
%Define training epochs
Epochs=270;%270 %default 300Als 150EMR

%% Output fil ename
switch ValSelection
    case 0
        ss='';
    case 1
        ss='IY_';
end
switch net_type
    case 0
        nt = 'LSTM_';
    case 1
        nt = 'CNN_';
    case 2
        nt = 'ConvLSTM_';
    case 3
        nt = 'FF_';
end
switch domain
    case 0
        dm='Lombardia_';
    case 1
        dm='Australia_';
    case 2
        dm='Poland_';
end

in = num2str(NumInputs);
out = num2str(NumOutputs);
neurons_str = num2str(neurons(1));
conv_str = num2str(filterSize);

switch CoordFlag
    case true
        coord_str = 'Georef_';
    case false
        coord_str = 'No_coord_';
end

outFOLD=strcat('C:\Users\ferrarilu\OneDrive - ETH Zurich\Tesi_Olgiati\Output\',...
    datestr(now,'yyyymmdd'),'_',dm,nt,ss,in,'In_',out,'Out_',neurons_str,'N_',...
    conv_str,'Conv_',coord_str, 'v','\');

mkdir(outFOLD);

%% Import meteorological, emissions and monitorig stations data
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

% figure;
% plot(m{1}(:,1));


X{end-1} = Data{end-1};
X{end} = Data{end};

%% Prepare data persistent model
[input, target, mean_norm, std_norm, PM10_target] = input_preprocess(Data, m, s, PM10,...
    NumOutputs, false, 0, 0, 0);

for i = 1:length(input{1,1})
    input_pers(i,:) = input{1,1}{1,i};
end

for i = 1:length(target{1,1})
    target_pers(i,:) = target{1,1}{1,i};

    mean_pers(i,:) = mean_norm{1,1}{1,i};

    std_pers(i,:) = std_norm{1,1}{1,i};

    PM10_target_pers(i,:) = PM10_target{1,1}{1,i};
end

clear input_id target_id mean_id std_id PM10_target_id input_val target_val mean_val std_val PM10_target_val

%% Results
output_pers = repmat(input_pers, NumOutputs, 1);

corr = corrcoef(output_pers, PM10_target_pers);
corr_pers_global = corr(1,2);
nrmse_pers_global = sqrt((sum((PM10_target_pers-output_pers).^2))/(sum((output_pers-mean(output_pers)).^2))); 

for i = 1:NumOutputs
    corr_pers(i) = 1 - sum((output_pers(i,:) - PM10_target_pers(i,:)).^2)/(sum((mean(PM10_target_pers(i,:)) - PM10_target_pers(i,:)).^2));
    nrmse_pers(i) = sqrt((sum((PM10_target_pers(i,:)-output_pers(i,:)).^2))/(sum((output_pers(i,:)-mean(output_pers(i,:))).^2))); 
end

% path = strcat(outFOLD, 'Data_persistent.mat');
% save(path, 'input_pers','output_pers','target_pers','mean_pers','std_pers','PM10_target_pers',...
%     'corr_pers_global','nrmse_pers_global','corr_pers','nrmse_pers');

figure;
plot(1:NumOutputs, [corr_pers], '.-', 'MarkerSize', 20, 'LineWidth', 2);
grid on;
xticks(1:NumOutputs);
xlabel('Days'); ylabel('R^2');
% path = strcat(outFOLD, 'R2_persistent.png');
% saveas(gcf, path);
close gcf;

figure;
plot(1:NumOutputs, [nrmse_pers], '.-', 'MarkerSize', 20, 'LineWidth', 2);
grid on;
xticks(1:NumOutputs);
xlabel('Days'); ylabel('nrmse');
% path = strcat(outFOLD, 'nrmse_persistent.png');
% saveas(gcf, path);
close gcf;

%% Prepare data CCN
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

%%
% Nstaz = width(X{1,1});
% Data_id{1} = X{1}(1:end-NumOutputs,:);
% Tval = Tval*2;
% N = (length(Data_id{1})-Tval)*Nstaz;
% 
% idx_row_id = false(length(Data_id{1}), 1);
% idx_row_id(1:end-Tval) = 1;
% test_stations = [];
% 
% input_id{1} = reshape(Data_id{1}(idx_row_id,:), N, 1);
% %%
% input_id{2} = repmat(X{end-1}(:), (length(Data_id{1})-Tval), 1);
% input_id{3} = repmat(X{end}(:), (length(Data_id{1})-Tval), 1);
% input_id{2} = reshape(input_id{2}, 63, 1815);
% input_id{2} = input_id{2}';
% input_id{2} = reshape(input_id{2}, N, 1);
% input_id{3} = reshape(input_id{3}, 63, 1815);
% input_id{3} = input_id{3}';
% input_id{3} = reshape(input_id{3}, N, 1);

%%

% if CoordFlag == true
%     input_id{2} = repmat(Data{end-1}(:), (length(Data_id{1})-Tval), 1);
%     input_id{3} = repmat(Data{end}(:), (length(Data_id{1})-Tval), 1);
%     input_id{2} = reshape(input_id{2}, N, 1);
%     input_id{3} = reshape(input_id{3}, N, 1);
% end




%% Struttura della rete a grafo (GNN)
if CoordFlag == false
    NumInputs = NumInputs-2;
end

if net_type == 3
    net = feedforwardnet(neurons);
    net.trainParam.epochs=Epochs;
    net.divideFcn = 'divideind';
    net.divideParam.trainInd = 1:length(input_id);
    net.divideParam.valInd = length(input_id)+1:length(input_id)+Tval;
else
    [CNN] = build_network(NumInputs, NumOutputs, net_type, neurons, filterSize, numHiddenUnits); % creo la rete neurale
    analyzeNetwork(CNN);
    options = trainingOptions('adam', ...
    'MaxEpochs', Epochs, ...
    'Verbose',1, ...
    'ValidationData',{input_val, target_val}, ...
    'ValidationFrequency', 30, ...
    'OutputNetwork', 'best-validation-loss' ...
    );
end


%% Calibration and validation CNN 

for i = 1:PreEpochs
    if net_type == 3
        training_start = tic;
        net=train(net,[input_id, input_val], [target_id, target_val], 'useParallel','yes','useGPU','yes');
        training_time(i) = toc(training_start);
        running_start = tic;
        output_id = sim(net,input_id);
        running_time(i) = toc(running_start);
    else
        training_start = tic;
        net = trainNetwork(input_id,target_id,CNN,options); 
        training_time(i) = toc(training_start);
        running_start = tic;
        output_id = predict(net,input_id);
        running_time(i) = toc(running_start);
    end
    
%     output_id = ((max(target_id, [], 2) - min(target_id, [], 2)) .* output_id) + min(target_id, [], 2);

    output_id = std_id.*output_id +mean_id ;

    % Validation set
%     if net_type == 3
%         output_val = sim(net,input_val);
%     else
%         output_val = predict(net,input_val);
%     end
%     
%     output_val = std_val.*output_val +mean_val ;
% 
%     R2_val(i) = 1 - sum((output_val(i,:) - PM10_target_val(i,:)).^2)/(sum((mean(PM10_target_val(i,:)) - PM10_target_val(i,:)).^2));
% 
%     nrmse_val(i) = sqrt((sum((PM10_target_val-output_val).^2))/(sum((output_val-mean(output_val)).^2))); 

    % Test set
    if net_type == 3
        output_test = sim(net,input_test);
    else
        output_test = predict(net,input_test);
    end

%     output_test = ((max(target_test, [], 2) - min(target_test, [], 2)) .* output_test) + min(target_test, [], 2);

    output_test = std_test.*output_test +mean_test ;

%     R2_test(i) = 1 - sum((output_test(i,:) - PM10_target_test(i,:)).^2)/(sum((mean(PM10_target_test(i,:)) - PM10_target_test(i,:)).^2));

    nrmse_test(i) = sqrt((sum((PM10_target_test-output_test).^2))/(sum((output_test-mean(output_test)).^2))); 
    
    if nrmse_test(i) <= min(nrmse_test)
        net_opt = net ; % trovo la GNN che fornisce l'R2 migliore
        opt_time = [training_time(i); running_time(i)];
    end

end 

fprintf('\nThe optimal network was trained in: %f seconds and run in %f seconds\n', opt_time(1), opt_time(2))

path = strcat(outFOLD, 'net_opt.mat');
save(path, 'net_opt');
% load(path);

path = strcat(outFOLD, 'Data.mat');
save(path, 'input_id','input_val','input_test','target_id','target_val','target_test',...
    'mean_id','mean_val','mean_test','std_id','std_val','std_test','PM10_target_id',...
    'PM10_target_val','PM10_target_test', 'test_stations', 'net_type', 'NumOutputs');


%% Create figures

if net_type == 3
    output_val = sim(net_opt,input_val);
else
    output_val = predict(net_opt,input_val);
end

% output_val = ((max(target_val, [], 2) - min(target_val, [], 2)) .* output_val) + min(target_val, [], 2);
output_val = std_val.*output_val +mean_val ;

corr = 1 - sum((output_val - PM10_target_val).^2)/(sum((mean(PM10_target_val) - PM10_target_val).^2));
nrmse = sqrt((sum((PM10_target_val-output_val).^2))/(sum((output_val-mean(output_val)).^2))); 

createScatter(PM10_target_val,output_val,'PM10 concentration [\mug/m^3]','CNN validation [\mug/m^3]',corr,nrmse);
% saveas(gcf, 'C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Output\correlogram_temporal_val.png');

if net_type == 3
    output_test = sim(net_opt,input_test);
else
    output_test = predict(net_opt,input_test);
end

% output_test = ((max(target_test, [], 2) - min(target_test, [], 2)) .* output_test) + min(target_test, [], 2);
output_test = std_test.*output_test +mean_test ;


corr = 1 - sum((output_test - PM10_target_test).^2)/(sum((mean(PM10_target_test) - PM10_target_test).^2));
nrmse = sqrt((sum((PM10_target_test-output_test).^2))/(sum((output_test-mean(output_test)).^2))); 

createScatter(PM10_target_test,output_test,'PM10 concentration [\mug/m^3]','CNN test [\mug/m^3]',corr,nrmse);
% saveas(gcf, 'C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Output\correlogram_spatial_val.png');


for i = 1:NumOutputs
    corr_val_step(i) = 1 - sum((output_val(i,:) - PM10_target_val(i,:)).^2)/(sum((mean(PM10_target_val(i,:)) - PM10_target_val(i,:)).^2));
    nrmse_val_step(i) = sqrt((sum((PM10_target_val(i,:)-output_val(i,:)).^2))/(sum((output_val(i,:)-mean(output_val(i,:))).^2))); 
    
%     createScatter(PM10_target_val(i,:),output_val(i,:),'PM10 concentration [\mug/m^3]','CNN validation [\mug/m^3]',corr_val_step(i),nrmse_val_step(i));

    corr_test_step(i) = 1 - sum((output_test(i,:) - PM10_target_test(i,:)).^2)/(sum((mean(PM10_target_test(i,:)) - PM10_target_test(i,:)).^2));
    nrmse_test_step(i) = sqrt((sum((PM10_target_test(i,:)-output_test(i,:)).^2))/(sum((output_test(i,:)-mean(output_test(i,:))).^2))); 
    
%     createScatter(PM10_target_test(i,:),output_test(i,:),'PM10 concentration [\mug/m^3]','CNN test [\mug/m^3]',corr_test_step(i),nrmse_test_step(i));
end

figure;
plot(1:NumOutputs, [corr_val_step; corr_test_step; corr_pers], '.-', 'MarkerSize', 20, 'LineWidth', 2);
grid on;
xticks(1:NumOutputs);
xlabel('Days'); ylabel('R^2');
legend('Validation dataset', 'Test dataset', 'Persistent model', 'Location', 'northeastoutside');

path = strcat(outFOLD, dm, 'R2_val_comparison.png');
saveas(gcf, path);

figure;
plot(1:NumOutputs, [nrmse_val_step; nrmse_test_step; nrmse_pers], '.-', 'MarkerSize', 20, 'LineWidth', 2);
grid on;
xticks(1:NumOutputs);
xlabel('Days'); ylabel('nrmse');
legend('Validation dataset', 'Test dataset', 'Persistent model', 'Location', 'northeastoutside');

path = strcat(outFOLD, dm, 'nrmse_val_comparison.png');
saveas(gcf, path);


%% End and save times
global_end = toc(global_start);

path = strcat(outFOLD, 'times.txt');

fid =fopen(path, 'w' );
fprintf(fid, 'The total time was: %f seconds\n', global_end);
fprintf(fid, '\nThe optimal network was trained in: %f seconds and run in %f seconds\n', opt_time(1), opt_time(2));
fclose(fid);

fprintf('\nThe total time was: %f seconds\n', global_end)
% load handel
% sound(y,Fs);


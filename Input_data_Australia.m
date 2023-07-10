clear
close all
clc

addpath('C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Input',...
    'C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Dati',...
    'C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Dati\Data_Australia',...
'C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Funzioni');

%% Load stations
Stazioni_aria = readtable("air-quality-monitoring-sites-summary.xlsx");

idx_staz = contains(Stazioni_aria{:,7}, 'Decommissioned');
Stazioni_aria(idx_staz,:) = [];

%% Load data
Air_data = readtable("NSW daylyPM10 2015-21.xls");
Air_data = Air_data(1:end-31,:);
year_2021 = readtable("Sydney Daily PM10 2021.xls");
Air_data = vertcat(Air_data, year_2021);

Air_data(:,[false;idx_staz]) = [];

%% Adjust leap years
idx = strcmp(Air_data{:,1}, '29/02/2016');
Air_data(idx,:) = [];
idx = strcmp(Air_data{:,1}, '29/02/2020');
Air_data(idx,:) = [];

%% Remove Alexandra station
Stazioni_aria(2,:) = [];
Air_data(:,3) = [];

%% Coordinate of selected stations
[Coord(2,:),Coord(1,:)] = deg2utm(Stazioni_aria{:,3},Stazioni_aria{:,4});

writematrix(Coord, 'C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Input\UTM_coordinates_Australia.csv');

%% Coord for maps
coord_maps = table(Stazioni_aria{:,4}, Stazioni_aria{:,3}, 'VariableNames', {'X', 'Y'});
writetable(coord_maps, 'C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Mappe\Stazioni_PM10_Australia.csv');

%% Substitute unavailable data
% Gaps in the data of 5 days or less are filled through linear
% interpolation of the available data from the same station. Longer gaps
% are filled by applying annverse distance weighting interpolation,
% considering the KNN stations.

PM10_fill = Air_data{:,2:end};

PM10_fill(PM10_fill < 0)=NaN; 

% T=365;
% for i = 1:width(PM10_fill)
%     [~, ym(:,i)] = moving_average(PM10_fill(:,i), T, 5);
%     [~, s2] = moving_average((PM10_fill(:,i) - ym(:,i)).^2, T, 5);
%     s(:,i) = s2.^0.5;
%     PM10_fill(:,i) = (PM10_fill(:,i) - ym(:,i))./s(:,i); 
% end

PM10_fill = fillmissing(PM10_fill, 'linear', 'MaxGap', 1);

% PM10_fill = PM102015_2019(2:end,:);

Coord = Coord';

K = round(sqrt(length(Coord)));

for i = 1:length(Coord)

    idx = isnan(PM10_fill(:,i));

    if sum(idx) > 0
        % Find KNN stations and distance
        X = Coord(:,1:2);
        X(i,:) = [];
        Y = Coord(i,1:2);
        [mIdx,mD] = knnsearch(X,Y,'K',K,'Distance','seuclidean');
        dati_temp = [PM10_fill(:,mIdx)];

        % Weights
        W = 1./(mD.^2);   

        % Interpolation
        PM10_fill(idx,i) =  sum(W.*dati_temp(idx,:), 2, 'omitnan')/sum(W);
    end
end

% PM10_fill = s.*PM10_fill+ym;

Air_data{:,2:end} = PM10_fill;

%% Create table with data
writetable(Air_data, 'C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Input\PM102015_2021_Australia.csv');



%% PM2.5
clear

%% Load stations
Stazioni_aria = readtable("air-quality-monitoring-sites-summary.xlsx");

idx_staz = contains(Stazioni_aria{:,7}, 'Decommissioned');
Stazioni_aria(idx_staz,:) = [];

%% Load data
Air_data = readtable("Sydney dayly PM25 2015-21.xls");
Air_data(:,[false;idx_staz]) = [];

%% Adjust leap years
idx = strcmp(Air_data{:,1}, '29/02/2016');
Air_data(idx,:) = [];
idx = strcmp(Air_data{:,1}, '29/02/2020');
Air_data(idx,:) = [];

%% Remove Alexandra and empty station
Air_data(:,3) = [];
Stazioni_aria(2,:) = [];

%% Coordinate of selected stations
[Coord(2,:),Coord(1,:)] = deg2utm(Stazioni_aria{:,3},Stazioni_aria{:,4});

writematrix(Coord, 'C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Input\UTM_coordinates_PM25_Australia.csv');

%% Coord for maps
coord_maps = table(Stazioni_aria{:,4}, Stazioni_aria{:,3}, 'VariableNames', {'X', 'Y'});
writetable(coord_maps, 'C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Mappe\Stazioni_PM25_Australia.csv');

%% Substitute unavailable data
% Gaps in the data of 5 days or less are filled through linear
% interpolation of the available data from the same station. Longer gaps
% are filled by applying annverse distance weighting interpolation,
% considering the KNN stations.

PM25_fill = Air_data{:,2:end};

PM25_fill(PM25_fill < 0)=NaN; 

PM25_fill = fillmissing(PM25_fill, 'linear', 'MaxGap', 1);

% PM10_fill = PM102015_2019(2:end,:);

Coord = Coord';

K = round(sqrt(length(Coord)));

for i = 1:length(Coord)

    idx = isnan(PM25_fill(:,i));

    if sum(idx) > 0
        % Find KNN stations and distance
        X = Coord(:,1:2);
        X(i,:) = [];
        Y = Coord(i,1:2);
        [mIdx,mD] = knnsearch(X,Y,'K',K,'Distance','seuclidean');
        dati_temp = [PM25_fill(:,mIdx)];

        % Weights
        W = 1./(mD.^2);   

        % Interpolation
        PM25_fill(idx,i) =  sum(W.*dati_temp(idx,:), 2, 'omitnan')/sum(W);
    end
end

Air_data{:,2:end} = PM25_fill;

%% Create table with data
writetable(Air_data, 'C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Input\PM252015_2021_Australia.csv');















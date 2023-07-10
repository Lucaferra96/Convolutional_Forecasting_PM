clear
close all
clc

addpath('C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Input',...
    'C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Dati',...
    'C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Dati\Data_Poland',...
'C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Funzioni');

%% Load stations
Stazioni_aria = readtable("./Dati/Data_Poland/Stations/Stations.csv");
Stazioni_aria = sortrows(Stazioni_aria, 4);
ID_PM10 = string(Stazioni_aria{:,4});

%% Load data
opts = detectImportOptions('2015_PM10_24g.xlsx');
stations = intersect(opts.SelectedVariableNames, ID_PM10, 'stable');
stations = ['KodStacji'; stations];
opts = detectImportOptions('2016_PM10_24g.xlsx');
stations = intersect(opts.SelectedVariableNames, stations, 'stable');

opts = detectImportOptions('2015_PM10_24g.xlsx');
opts.VariableNamesRange = 'A1';
opts.DataRange = 'A4';
opts.VariableNames = stations;


Dati_2015aria = readtable("2015_PM10_24g.xlsx", opts);

opts = detectImportOptions('2016_PM10_24g.xlsx');
opts.VariableNamesRange = 'A1';
opts.DataRange = 'A6';
opts.VariableNames = stations;


Dati_2016aria = readtable("2016_PM10_24g.xlsx", opts);
Dati_2017aria = readtable("2017_PM10_24g.xlsx", opts);
Dati_2018aria = readtable("2018_PM10_24g.xlsx", opts);
Dati_2019aria = readtable("2019_PM10_24g.xlsx", opts);
Dati_2020aria = readtable("2020_PM10_24g.xlsx", opts);
Dati_2021aria = readtable("2021_PM10_24g.xlsx", opts);

%% Adjust leap years
idx = Dati_2016aria{:,1} == '29/02/2016';
Dati_2016aria(idx,:) = [];
idx = Dati_2020aria{:,1} == '29/02/2020';
Dati_2020aria(idx,:) = [];

%% Rearrange data 
PM102015_2019 = vertcat(Dati_2015aria, Dati_2016aria, Dati_2017aria, Dati_2018aria, Dati_2019aria, Dati_2020aria, Dati_2021aria);

PM10 = PM102015_2019{:,2:end};

%% Coordinate of selected stations
idx = ismember(Stazioni_aria{:,4}, stations);
Stazioni_aria = Stazioni_aria(idx,:);

[Coord(2,:),Coord(1,:)] = deg2utm(Stazioni_aria{:,2},Stazioni_aria{:,1});

writematrix(Coord, 'C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Input\UTM_coordinates_Poland.csv');

%% Coord for maps
coord_maps = table(Stazioni_aria{:,1}, Stazioni_aria{:,2}, 'VariableNames', {'X', 'Y'});
writetable(coord_maps, 'C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Mappe\Stazioni_PM10_Polonia.csv');

%% Substitute unavailable data
% Gaps in the data of 1 days  are filled through linear
% interpolation of the available data from the same station. Longer gaps
% are filled by applying annverse distance weighting interpolation,
% considering the KNN stations.

PM10_fill = PM102015_2019{:,2:end};

PM10_fill(PM10_fill < 0)=NaN; 

% T=365;
% for i = 1:width(PM10_fill)
%     [~, ym(:,i)] = moving_average(PM10_fill(:,i), T, 5);
%     [~, s2] = moving_average((PM10_fill(:,i) - ym(:,i)).^2, T, 5);
%     s(:,i) = s2.^0.5;
%     PM10_fill(:,i) = (PM10_fill(:,i) - ym(:,i))./s(:,i); 
% end

PM10_fill = fillmissing(PM10_fill, 'linear', 'MaxGap', 1);

Coord = Coord';

K = round(sqrt(length(Coord)));

for i = 1:length(Coord)

    idx = isnan(PM10_fill(:,i));

    if sum(idx) > 0
        % Find KNN stations and distance
        X = Coord;
        X(i,:) = [];
        Y = Coord(i,:);
        [mIdx,mD] = knnsearch(X,Y,'K',K,'Distance','seuclidean');
        dati_temp = [PM10_fill(:,mIdx)];

        % Weights
        W = 1./(mD.^2);   

        % Interpolation
        PM10_fill(idx,i) =  sum(W.*dati_temp(idx,:), 2, 'omitnan')/sum(W);
    end
end

% PM10_fill = s.*PM10_fill+ym;

PM102015_2019{:,2:end} = PM10_fill;

%% Create table with data
writetable(PM102015_2019, 'C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Input\PM102015_2021_Poland.csv');

%% PM2.5
% clear
clear Coord
%% Load stations
Stazioni_aria = readtable("./Dati/Data_Poland/Stations/Stations.csv");
Stazioni_aria = sortrows(Stazioni_aria, 4);
ID_PM10 = string(Stazioni_aria{:,4});

%% Load data
opts = detectImportOptions('2015_PM25_24g.xlsx');
stations = intersect(opts.SelectedVariableNames, ID_PM10, 'stable');
stations = ['KodStacji'; stations];
opts = detectImportOptions('2016_PM25_24g.xlsx');
stations = intersect(opts.SelectedVariableNames, stations, 'stable');

opts = detectImportOptions('2015_PM25_24g.xlsx');
opts.VariableNamesRange = 'A1';
opts.DataRange = 'A4';
opts.VariableNames = stations;


Dati_2015aria = readtable("2015_PM25_24g.xlsx", opts);

opts = detectImportOptions('2016_PM25_24g.xlsx');
opts.VariableNamesRange = 'A1';
opts.DataRange = 'A6';
opts.VariableNames = stations;


Dati_2016aria = readtable("2016_PM25_24g.xlsx", opts);
Dati_2017aria = readtable("2017_PM25_24g.xlsx", opts);
Dati_2018aria = readtable("2018_PM25_24g.xlsx", opts);
Dati_2019aria = readtable("2019_PM25_24g.xlsx", opts);
Dati_2020aria = readtable("2020_PM25_24g.xlsx", opts);
Dati_2021aria = readtable("2021_PM25_24g.xlsx", opts);

%% Adjust leap years
idx = Dati_2016aria{:,1} == '29/02/2016';
Dati_2016aria(idx,:) = [];
idx = Dati_2020aria{:,1} == '29/02/2020';
Dati_2020aria(idx,:) = [];

%% Rearrange data 
PM252015_2019 = vertcat(Dati_2015aria, Dati_2016aria, Dati_2017aria, Dati_2018aria, Dati_2019aria, Dati_2020aria, Dati_2021aria);

PM25 = PM252015_2019{:,2:end};

%% Coordinate of selected stations
idx = ismember(Stazioni_aria{:,4}, stations);
Stazioni_aria = Stazioni_aria(idx,:);

[Coord(2,:),Coord(1,:)] = deg2utm(Stazioni_aria{:,2},Stazioni_aria{:,1});

writematrix(Coord, 'C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Input\UTM_coordinates_PM25_Poland.csv');

%% Coord for maps
coord_maps = table(Stazioni_aria{:,1}, Stazioni_aria{:,2}, 'VariableNames', {'X', 'Y'});
writetable(coord_maps, 'C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Mappe\Stazioni_PM25_Polonia.csv');

%% Substitute unavailable data
% Gaps in the data of 1 days  are filled through linear
% interpolation of the available data from the same station. Longer gaps
% are filled by applying annverse distance weighting interpolation,
% considering the KNN stations.

PM25_fill = PM252015_2019{:,2:end};

PM25_fill(PM25_fill < 0)=NaN; 

PM25_fill = fillmissing(PM25_fill, 'linear', 'MaxGap', 1);

Coord = Coord';

K = round(sqrt(length(Coord)));

for i = 1:length(Coord)

    idx = isnan(PM25_fill(:,i));

    if sum(idx) > 0
        % Find KNN stations and distance
        X = Coord;
        X(i,:) = [];
        Y = Coord(i,:);
        [mIdx,mD] = knnsearch(X,Y,'K',K,'Distance','seuclidean');
        dati_temp = [PM25_fill(:,mIdx)];

        % Weights
        W = 1./(mD.^2);   

        % Interpolation
        PM25_fill(idx,i) =  sum(W.*dati_temp(idx,:), 2, 'omitnan')/sum(W);
    end
end

PM252015_2019{:,2:end} = PM25_fill;

%% Create table with data
writetable(PM252015_2019, 'C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Input\PM252015_2021_Poland.csv');






















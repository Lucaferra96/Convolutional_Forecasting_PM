clear
close all
clc

addpath('C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Input',...
    'C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Dati',...
'C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Funzioni');

%% Load stations
Stazioni_aria = readtable("Stazioni_qualit__dell_aria.csv");

idx = isnat(Stazioni_aria{:,11});
Stazioni_aria = Stazioni_aria(idx,:);

idx_staz = strcmp('PM10 (SM2005)', Stazioni_aria{:,2});
Stazioni_PM10 = Stazioni_aria(idx_staz,:);

ID_PM10 = Stazioni_PM10(:,1);

% % NO2
% 
% idx_staz = strcmp('Ossidi di Azoto', Stazioni_aria{:,2});
% Stazioni_NO2 = Stazioni_aria(idx_staz,:);
% 
% ID_NO2 = Stazioni_NO2(:,4);

%% Load data
opts = delimitedTextImportOptions("NumVariables", 5);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["IdSensore", "Data", "Valore", "Stato", "idOperatore"];
opts.VariableTypes = ["double", "datetime", "double", "categorical", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "Stato", "EmptyFieldRule", "auto");
opts = setvaropts(opts, "Data", "InputFormat", "dd/MM/yyyy HH:mm:ss");

% Import the data
% Dati_2015aria = readtable("Air_2015.csv", opts);
% Dati_2016aria = readtable("Air_2016.csv", opts);
% Dati_2017aria = readtable("Air_2017.csv", opts);
% Dati_2018aria = readtable("Air_2018.csv", opts);
% Dati_2019aria = readtable("Air_2019.csv", opts);
% Dati_2020aria = readtable("Air_2020.csv", opts);
% Dati_2021aria = readtable("Air_2021.csv", opts);

load('Dati_Aria.mat')

clear opts

%% Adjust leap years

idx = Dati_2016aria{:,2} == '29/02/2016 00:00:00';
Dati_2016aria(idx,:) = [];
idx = Dati_2020aria{:,2} == '29/02/2020 00:00:00';
Dati_2020aria(idx,:) = [];

%% Select common stations
ID_PM10 = intersect(Dati_2015aria(:,1), ID_PM10);
ID_PM10 = intersect(Dati_2016aria(:,1), ID_PM10);
ID_PM10 = intersect(Dati_2017aria(:,1), ID_PM10);
ID_PM10 = intersect(Dati_2018aria(:,1), ID_PM10);
ID_PM10 = intersect(Dati_2019aria(:,1), ID_PM10);
ID_PM10 = intersect(Dati_2020aria(:,1), ID_PM10);
ID_PM10 = intersect(Dati_2021aria(:,1), ID_PM10);

idx = ismember(Stazioni_PM10{:,1}, ID_PM10{:,:});
Stazioni_PM10 = Stazioni_PM10(idx,:);

% Select stations common to NO2 and PM10
% ID_PM10 = Stazioni_PM10(:,4);
% ID_NO2 = intersect(ID_PM10, ID_NO2);
% ID_PM10 = ID_NO2;
% 
% idx = ismember(Stazioni_NO2{:,4}, ID_NO2{:,:});
% Stazioni_NO2 = Stazioni_NO2(idx,:);
% 
% idx = ismember(Stazioni_PM10{:,4}, ID_PM10{:,:});
% Stazioni_PM10 = Stazioni_PM10(idx,:);
% 
% ID_PM10 = Stazioni_PM10(:,1);
% ID_NO2 = Stazioni_NO2(:,1);


%% Coordinate of selected stations
Coord = Stazioni_PM10{:, [1,12,13]};

Coord = sortrows(Coord,1);
Coord = Coord';

writematrix(Coord, 'C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Input\UTM_coordinates.csv');

%% Coord for maps
coord_maps = table(Stazioni_PM10{:,15}, Stazioni_PM10{:,14}, 'VariableNames', {'X', 'Y'});
writetable(coord_maps, 'C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Mappe\Stazioni_PM10_Lombardia.csv');

%% Extract PM10 of selected stations
% PM10 2015
idx = ismember(Dati_2015aria(:,1), ID_PM10); %mette 1 alle righe che contengono gli ID_PM10
PM10_r=Dati_2015aria(idx,1:3); %seleziona tutti i dati degli PM10
PM10_s=sortrows(PM10_r,1); %ordina i dati in base all'IDstazione
data= convertTo(PM10_s.Data,'excel'); %converte le date come valore excel
PM10_t=[data,PM10_s.IdSensore,PM10_s.Valore]; %converte la tabella in una matrice
PM10_p=pivottable(PM10_t,1,2,3,@sum); %tabella pivot creata date sulle righe e IdSensore sulle colonne
Vuoti=cellfun('isempty',PM10_p); %seleziona le celle vuote
PM10_p(Vuoti)={NaN};%sostituisce le celle vpute con NaN
PM10_2015=cell2mat(PM10_p(:,2:end));

% PM10 2016
idx=ismember(Dati_2016aria(:,1),ID_PM10); %mette 1 alle righe che contengono gli ID_PM10
PM10_r=Dati_2016aria(idx,1:3); %seleziona tutti i dati degli PM10
PM10_s=sortrows(PM10_r,1); %ordina i dati in base all'IDstazione
data= convertTo(PM10_s.Data,'excel'); %converte le date come valore excel
PM10_t=[data,PM10_s.IdSensore,PM10_s.Valore]; %converte la tabella in una matrice
PM10_p=pivottable(PM10_t,1,2,3,@sum); %tabella pivot creata date sulle righe e IdSensore sulle colonne
Vuoti=cellfun('isempty',PM10_p); %seleziona le celle vuote
PM10_p(Vuoti)={NaN};%sostituisce le celle vpute con NaN
PM10_2016=cell2mat(PM10_p(:,2:end));

% PM10 2017
idx=ismember(Dati_2017aria(:,1),ID_PM10); %mette 1 alle righe che contengono gli ID_PM10
PM10_r=Dati_2017aria(idx,1:3); %seleziona tutti i dati degli PM10
PM10_s=sortrows(PM10_r,1); %ordina i dati in base all'IDstazione
data= convertTo(PM10_s.Data,'excel'); %converte le date come valore excel
PM10_t=[data,PM10_s.IdSensore,PM10_s.Valore]; %converte la tabella in una matrice
PM10_p=pivottable(PM10_t,1,2,3,@sum); %tabella pivot creata date sulle righe e IdSensore sulle colonne
Vuoti=cellfun('isempty',PM10_p); %seleziona le celle vuote
PM10_p(Vuoti)={NaN};%sostituisce le celle vpute con NaN
PM10_2017=cell2mat(PM10_p(:,2:end));

% PM10 2018
idx=ismember(Dati_2018aria(:,1),ID_PM10); %mette 1 alle righe che contengono gli ID_PM10
PM10_r=Dati_2018aria(idx,1:3); %seleziona tutti i dati degli PM10
PM10_s=sortrows(PM10_r,1); %ordina i dati in base all'IDstazione
data= convertTo(PM10_s.Data,'excel'); %converte le date come valore excel
PM10_t=[data,PM10_s.IdSensore,PM10_s.Valore]; %converte la tabella in una matrice
PM10_p=pivottable(PM10_t,1,2,3,@sum); %tabella pivot creata date sulle righe e IdSensore sulle colonne
Vuoti=cellfun('isempty',PM10_p); %seleziona le celle vuote
PM10_p(Vuoti)={NaN};%sostituisce le celle vpute con NaN
PM10_2018=cell2mat(PM10_p(:,2:end));

% PM10 2019
idx=ismember(Dati_2019aria(:,1),ID_PM10); %mette 1 alle righe che contengono gli ID_PM10
PM10_r=Dati_2019aria(idx,1:3); %seleziona tutti i dati degli PM10
PM10_s=sortrows(PM10_r,1); %ordina i dati in base all'IDstazione
data= convertTo(PM10_s.Data,'excel'); %converte le date come valore excel
PM10_t=[data,PM10_s.IdSensore,PM10_s.Valore]; %converte la tabella in una matrice
PM10_p=pivottable(PM10_t,1,2,3,@sum); %tabella pivot creata date sulle righe e IdSensore sulle colonne
Vuoti=cellfun('isempty',PM10_p); %seleziona le celle vuote
PM10_p(Vuoti)={NaN};%sostituisce le celle vpute con NaN
PM10_2019=cell2mat(PM10_p(:,2:end));

% PM10 2020
idx=ismember(Dati_2020aria(:,1),ID_PM10); %mette 1 alle righe che contengono gli ID_PM10
PM10_r=Dati_2020aria(idx,1:3); %seleziona tutti i dati degli PM10
PM10_s=sortrows(PM10_r,1); %ordina i dati in base all'IDstazione
data= convertTo(PM10_s.Data,'excel'); %converte le date come valore excel
PM10_t=[data,PM10_s.IdSensore,PM10_s.Valore]; %converte la tabella in una matrice
PM10_p=pivottable(PM10_t,1,2,3,@sum); %tabella pivot creata date sulle righe e IdSensore sulle colonne
Vuoti=cellfun('isempty',PM10_p); %seleziona le celle vuote
PM10_p(Vuoti)={NaN};%sostituisce le celle vpute con NaN
PM10_2020=cell2mat(PM10_p(:,2:end));

% PM10 2021
idx=ismember(Dati_2021aria(:,1),ID_PM10); %mette 1 alle righe che contengono gli ID_PM10
PM10_r=Dati_2021aria(idx,1:3); %seleziona tutti i dati degli PM10
PM10_s=sortrows(PM10_r,1); %ordina i dati in base all'IDstazione
data= convertTo(PM10_s.Data,'excel'); %converte le date come valore excel
PM10_t=[data,PM10_s.IdSensore,PM10_s.Valore]; %converte la tabella in una matrice
PM10_p=pivottable(PM10_t,1,2,3,@sum); %tabella pivot creata date sulle righe e IdSensore sulle colonne
Vuoti=cellfun('isempty',PM10_p); %seleziona le celle vuote
PM10_p(Vuoti)={NaN};%sostituisce le celle vpute con NaN
PM10_2021=cell2mat(PM10_p(:,2:end));

clear PM10_r PM10_s data PM10_t PM10_p idx Vuoti

%% Add 01/01/2015

PM10_2015 = [PM10_2015(1,:); NaN(1,width(PM10_2015)); PM10_2015(2:end,:)];

%% Join the data together

PM102015_2019=[PM10_2015; PM10_2016(2:end,:); PM10_2017(2:end,:); PM10_2018(2:end,:); PM10_2019(2:end,:); PM10_2020(2:end,:); PM10_2021(2:end,:)]; 

%% Substitute unavailable data
% Gaps in the data of 5 days or less are filled through linear
% interpolation of the available data from the same station. Longer gaps
% are filled by applying annverse distance weighting interpolation,
% considering the KNN stations.

PM102015_2019(PM102015_2019 < 0)=NaN; 

PM10_fill = PM102015_2019(2:end,:);

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
        X = Coord(:,2:3);
        X(i,:) = [];
        Y = Coord(i,2:3);
        [mIdx,mD] = knnsearch(X,Y,'K',K,'Distance','seuclidean');
        dati_temp = [PM10_fill(:,mIdx)];

        % Weights
        W = 1./(mD.^2);   

        % Interpolation
        PM10_fill(idx,i) =  sum(W.*dati_temp(idx,:), 2, 'omitnan')/sum(W);
    end
end

% PM10_fill = s.*PM10_fill+ym;

PM102015_2019(2:end,:) = PM10_fill;

%% Create table with data

t1= datetime(2015,1,1,0,0,0);
t2=datetime(2021,12,31,23,0,0);
dn=(t1:days(1):t2).';
dn(dn=='29-Feb-2016') = [];
dn(dn=='29-Feb-2020') = [];
dn = [NaT; dn];

PM102015_2019_T = table(dn, PM102015_2019);
PM102015_2019_T = splitvars(PM102015_2019_T);

writetable(PM102015_2019_T, 'C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Input\PM102015_2021.csv');

%% PM2.5

idx_staz = strcmp('Particelle sospese PM2.5', Stazioni_aria{:,2});
Stazioni_PM25 = Stazioni_aria(idx_staz,:);
ID_PM25 = Stazioni_PM25(:,1);

%% Select common stations
ID_PM25 = intersect(Dati_2015aria(:,1), ID_PM25);
ID_PM25 = intersect(Dati_2016aria(:,1), ID_PM25);
ID_PM25 = intersect(Dati_2017aria(:,1), ID_PM25);
ID_PM25 = intersect(Dati_2018aria(:,1), ID_PM25);
ID_PM25 = intersect(Dati_2019aria(:,1), ID_PM25);
ID_PM25 = intersect(Dati_2020aria(:,1), ID_PM25);
ID_PM25 = intersect(Dati_2021aria(:,1), ID_PM25);

idx = ismember(Stazioni_PM25{:,1}, ID_PM25{:,:});
Stazioni_PM25 = Stazioni_PM25(idx,:);

%% Coordinate of selected stations
Coord = Stazioni_PM25{:, [1,12:13]};

Coord = sortrows(Coord,1);
Coord = Coord';

writematrix(Coord, 'C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Input\UTM_coordinates_PM25_Lombardia.csv');

%% Coord for maps
coord_maps = table(Stazioni_PM25{:,15}, Stazioni_PM25{:,14}, 'VariableNames', {'X', 'Y'});
writetable(coord_maps, 'C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Mappe\Stazioni_PM25_Lombardia.csv');

%% Extract PM25 of selected stations
% PM25 2015
idx = ismember(Dati_2015aria(:,1), ID_PM25); %mette 1 alle righe che contengono gli ID_PM25
PM25_r=Dati_2015aria(idx,1:3); %seleziona tutti i dati degli PM25
PM25_s=sortrows(PM25_r,1); %ordina i dati in base all'IDstazione
data= convertTo(PM25_s.Data,'excel'); %converte le date come valore excel
PM25_t=[data,PM25_s.IdSensore,PM25_s.Valore]; %converte la tabella in una matrice
PM25_p=pivottable(PM25_t,1,2,3,@sum); %tabella pivot creata date sulle righe e IdSensore sulle colonne
Vuoti=cellfun('isempty',PM25_p); %seleziona le celle vuote
PM25_p(Vuoti)={NaN};%sostituisce le celle vpute con NaN
PM25_2015=cell2mat(PM25_p(:,2:end));

% PM25 2016
idx=ismember(Dati_2016aria(:,1),ID_PM25); %mette 1 alle righe che contengono gli ID_PM25
PM25_r=Dati_2016aria(idx,1:3); %seleziona tutti i dati degli PM25
PM25_s=sortrows(PM25_r,1); %ordina i dati in base all'IDstazione
data= convertTo(PM25_s.Data,'excel'); %converte le date come valore excel
PM25_t=[data,PM25_s.IdSensore,PM25_s.Valore]; %converte la tabella in una matrice
PM25_p=pivottable(PM25_t,1,2,3,@sum); %tabella pivot creata date sulle righe e IdSensore sulle colonne
Vuoti=cellfun('isempty',PM25_p); %seleziona le celle vuote
PM25_p(Vuoti)={NaN};%sostituisce le celle vpute con NaN
PM25_2016=cell2mat(PM25_p(:,2:end));

% PM25 2017
idx=ismember(Dati_2017aria(:,1),ID_PM25); %mette 1 alle righe che contengono gli ID_PM25
PM25_r=Dati_2017aria(idx,1:3); %seleziona tutti i dati degli PM25
PM25_s=sortrows(PM25_r,1); %ordina i dati in base all'IDstazione
data= convertTo(PM25_s.Data,'excel'); %converte le date come valore excel
PM25_t=[data,PM25_s.IdSensore,PM25_s.Valore]; %converte la tabella in una matrice
PM25_p=pivottable(PM25_t,1,2,3,@sum); %tabella pivot creata date sulle righe e IdSensore sulle colonne
Vuoti=cellfun('isempty',PM25_p); %seleziona le celle vuote
PM25_p(Vuoti)={NaN};%sostituisce le celle vpute con NaN
PM25_2017=cell2mat(PM25_p(:,2:end));

% PM25 2018
idx=ismember(Dati_2018aria(:,1),ID_PM25); %mette 1 alle righe che contengono gli ID_PM25
PM25_r=Dati_2018aria(idx,1:3); %seleziona tutti i dati degli PM25
PM25_s=sortrows(PM25_r,1); %ordina i dati in base all'IDstazione
data= convertTo(PM25_s.Data,'excel'); %converte le date come valore excel
PM25_t=[data,PM25_s.IdSensore,PM25_s.Valore]; %converte la tabella in una matrice
PM25_p=pivottable(PM25_t,1,2,3,@sum); %tabella pivot creata date sulle righe e IdSensore sulle colonne
Vuoti=cellfun('isempty',PM25_p); %seleziona le celle vuote
PM25_p(Vuoti)={NaN};%sostituisce le celle vpute con NaN
PM25_2018=cell2mat(PM25_p(:,2:end));

% PM25 2019
idx=ismember(Dati_2019aria(:,1),ID_PM25); %mette 1 alle righe che contengono gli ID_PM25
PM25_r=Dati_2019aria(idx,1:3); %seleziona tutti i dati degli PM25
PM25_s=sortrows(PM25_r,1); %ordina i dati in base all'IDstazione
data= convertTo(PM25_s.Data,'excel'); %converte le date come valore excel
PM25_t=[data,PM25_s.IdSensore,PM25_s.Valore]; %converte la tabella in una matrice
PM25_p=pivottable(PM25_t,1,2,3,@sum); %tabella pivot creata date sulle righe e IdSensore sulle colonne
Vuoti=cellfun('isempty',PM25_p); %seleziona le celle vuote
PM25_p(Vuoti)={NaN};%sostituisce le celle vpute con NaN
PM25_2019=cell2mat(PM25_p(:,2:end));

% PM25 2020
idx=ismember(Dati_2020aria(:,1),ID_PM25); %mette 1 alle righe che contengono gli ID_PM25
PM25_r=Dati_2020aria(idx,1:3); %seleziona tutti i dati degli PM25
PM25_s=sortrows(PM25_r,1); %ordina i dati in base all'IDstazione
data= convertTo(PM25_s.Data,'excel'); %converte le date come valore excel
PM25_t=[data,PM25_s.IdSensore,PM25_s.Valore]; %converte la tabella in una matrice
PM25_p=pivottable(PM25_t,1,2,3,@sum); %tabella pivot creata date sulle righe e IdSensore sulle colonne
Vuoti=cellfun('isempty',PM25_p); %seleziona le celle vuote
PM25_p(Vuoti)={NaN};%sostituisce le celle vpute con NaN
PM25_2020=cell2mat(PM25_p(:,2:end));

% PM25 2021
idx=ismember(Dati_2021aria(:,1),ID_PM25); %mette 1 alle righe che contengono gli ID_PM25
PM25_r=Dati_2021aria(idx,1:3); %seleziona tutti i dati degli PM25
PM25_s=sortrows(PM25_r,1); %ordina i dati in base all'IDstazione
data= convertTo(PM25_s.Data,'excel'); %converte le date come valore excel
PM25_t=[data,PM25_s.IdSensore,PM25_s.Valore]; %converte la tabella in una matrice
PM25_p=pivottable(PM25_t,1,2,3,@sum); %tabella pivot creata date sulle righe e IdSensore sulle colonne
Vuoti=cellfun('isempty',PM25_p); %seleziona le celle vuote
PM25_p(Vuoti)={NaN};%sostituisce le celle vpute con NaN
PM25_2021=cell2mat(PM25_p(:,2:end));

clear PM25_r PM25_s data PM25_t PM25_p idx Vuoti

%% Add 01/01/2015

PM25_2015 = [PM25_2015(1,:); NaN(1,width(PM25_2015)); PM25_2015(2:end,:)];

%% Join the data together

PM252015_2019=[PM25_2015; PM25_2016(2:end,:); PM25_2017(2:end,:); PM25_2018(2:end,:); PM25_2019(2:end,:); PM25_2020(2:end,:); PM25_2021(2:end,:)]; 

%% Substitute unavailable data
% Gaps in the data of 5 days or less are filled through linear
% interpolation of the available data from the same station. Longer gaps
% are filled by applying annverse distance weighting interpolation,
% considering the KNN stations.

PM252015_2019(PM252015_2019 < 0)=NaN; 

PM25_fill = PM252015_2019(2:end,:);

PM25_fill = fillmissing(PM25_fill, 'linear', 'MaxGap', 1);

PM25_fill(PM25_fill < 0)=NaN; 

Coord = Coord';

K = round(sqrt(length(Coord)));

for i = 1:length(Coord)

    idx = isnan(PM25_fill(:,i));

    if sum(idx) > 0
        % Find KNN stations and distance
        X = Coord(:,2:3);
        X(i,:) = [];
        Y = Coord(i,2:3);
        [mIdx,mD] = knnsearch(X,Y,'K',K,'Distance','seuclidean');
        dati_temp = [PM25_fill(:,mIdx)];

        % Weights
        W = 1./(mD.^2);   

        % Interpolation
        PM25_fill(idx,i) =  sum(W.*dati_temp(idx,:), 2, 'omitnan')/sum(W);
    end
end

PM252015_2019(2:end,:) = PM25_fill;

%% Create table with data

t1= datetime(2015,1,1,0,0,0);
t2=datetime(2021,12,31,23,0,0);
dn=(t1:days(1):t2).';
dn(dn=='29-Feb-2016') = [];
dn(dn=='29-Feb-2020') = [];
dn = [NaT; dn];

PM252015_2019_T = table(dn, PM252015_2019);
PM252015_2019_T = splitvars(PM252015_2019_T);

writetable(PM252015_2019_T, 'C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Input\PM252015_2021.csv');














% %% Extract NO2 of selected stations
% % NO2 2015
% idx = ismember(Dati_2015aria(:,1), ID_NO2); %mette 1 alle righe che contengono gli ID_NO2
% NO2_r=Dati_2015aria(idx,1:3); %seleziona tutti i dati degli NO2
% NO2_s=sortrows(NO2_r,1); %ordina i dati in base all'IDstazione
% data= convertTo(NO2_s.Data,'excel'); %converte le date come valore excel
% NO2_t=[data,NO2_s.IdSensore,NO2_s.Valore]; %converte la tabella in una matrice
% NO2_p=pivottable(NO2_t,1,2,3,@sum); %tabella pivot creata date sulle righe e IdSensore sulle colonne
% Vuoti=cellfun('isempty',NO2_p); %seleziona le celle vuote
% NO2_p(Vuoti)={NaN};%sostituisce le celle vpute con NaN
% NO2_2015=cell2mat(NO2_p(:,2:end));
% 
% % NO2 2016
% idx=ismember(Dati_2016aria(:,1),ID_NO2); %mette 1 alle righe che contengono gli ID_NO2
% NO2_r=Dati_2016aria(idx,1:3); %seleziona tutti i dati degli NO2
% NO2_s=sortrows(NO2_r,1); %ordina i dati in base all'IDstazione
% data= convertTo(NO2_s.Data,'excel'); %converte le date come valore excel
% NO2_t=[data,NO2_s.IdSensore,NO2_s.Valore]; %converte la tabella in una matrice
% NO2_p=pivottable(NO2_t,1,2,3,@sum); %tabella pivot creata date sulle righe e IdSensore sulle colonne
% Vuoti=cellfun('isempty',NO2_p); %seleziona le celle vuote
% NO2_p(Vuoti)={NaN};%sostituisce le celle vpute con NaN
% NO2_2016=cell2mat(NO2_p(:,2:end));
% 
% % NO2 2017
% idx=ismember(Dati_2017aria(:,1),ID_NO2); %mette 1 alle righe che contengono gli ID_NO2
% NO2_r=Dati_2017aria(idx,1:3); %seleziona tutti i dati degli NO2
% NO2_s=sortrows(NO2_r,1); %ordina i dati in base all'IDstazione
% data= convertTo(NO2_s.Data,'excel'); %converte le date come valore excel
% NO2_t=[data,NO2_s.IdSensore,NO2_s.Valore]; %converte la tabella in una matrice
% NO2_p=pivottable(NO2_t,1,2,3,@sum); %tabella pivot creata date sulle righe e IdSensore sulle colonne
% Vuoti=cellfun('isempty',NO2_p); %seleziona le celle vuote
% NO2_p(Vuoti)={NaN};%sostituisce le celle vpute con NaN
% NO2_2017=cell2mat(NO2_p(:,2:end));
% 
% % NO2 2018
% idx=ismember(Dati_2018aria(:,1),ID_NO2); %mette 1 alle righe che contengono gli ID_NO2
% NO2_r=Dati_2018aria(idx,1:3); %seleziona tutti i dati degli NO2
% NO2_s=sortrows(NO2_r,1); %ordina i dati in base all'IDstazione
% data= convertTo(NO2_s.Data,'excel'); %converte le date come valore excel
% NO2_t=[data,NO2_s.IdSensore,NO2_s.Valore]; %converte la tabella in una matrice
% NO2_p=pivottable(NO2_t,1,2,3,@sum); %tabella pivot creata date sulle righe e IdSensore sulle colonne
% Vuoti=cellfun('isempty',NO2_p); %seleziona le celle vuote
% NO2_p(Vuoti)={NaN};%sostituisce le celle vpute con NaN
% NO2_2018=cell2mat(NO2_p(:,2:end));
% 
% % NO2 2019
% idx=ismember(Dati_2019aria(:,1),ID_NO2); %mette 1 alle righe che contengono gli ID_NO2
% NO2_r=Dati_2019aria(idx,1:3); %seleziona tutti i dati degli NO2
% NO2_s=sortrows(NO2_r,1); %ordina i dati in base all'IDstazione
% data= convertTo(NO2_s.Data,'excel'); %converte le date come valore excel
% NO2_t=[data,NO2_s.IdSensore,NO2_s.Valore]; %converte la tabella in una matrice
% NO2_p=pivottable(NO2_t,1,2,3,@sum); %tabella pivot creata date sulle righe e IdSensore sulle colonne
% Vuoti=cellfun('isempty',NO2_p); %seleziona le celle vuote
% NO2_p(Vuoti)={NaN};%sostituisce le celle vpute con NaN
% NO2_2019=cell2mat(NO2_p(:,2:end));
% 
% % NO2 2020
% idx=ismember(Dati_2020aria(:,1),ID_NO2); %mette 1 alle righe che contengono gli ID_NO2
% NO2_r=Dati_2020aria(idx,1:3); %seleziona tutti i dati degli NO2
% NO2_s=sortrows(NO2_r,1); %ordina i dati in base all'IDstazione
% data= convertTo(NO2_s.Data,'excel'); %converte le date come valore excel
% NO2_t=[data,NO2_s.IdSensore,NO2_s.Valore]; %converte la tabella in una matrice
% NO2_p=pivottable(NO2_t,1,2,3,@sum); %tabella pivot creata date sulle righe e IdSensore sulle colonne
% Vuoti=cellfun('isempty',NO2_p); %seleziona le celle vuote
% NO2_p(Vuoti)={NaN};%sostituisce le celle vpute con NaN
% NO2_2020=cell2mat(NO2_p(:,2:end));
% 
% % NO2 2021
% idx=ismember(Dati_2021aria(:,1),ID_NO2); %mette 1 alle righe che contengono gli ID_NO2
% NO2_r=Dati_2021aria(idx,1:3); %seleziona tutti i dati degli NO2
% NO2_s=sortrows(NO2_r,1); %ordina i dati in base all'IDstazione
% data= convertTo(NO2_s.Data,'excel'); %converte le date come valore excel
% NO2_t=[data,NO2_s.IdSensore,NO2_s.Valore]; %converte la tabella in una matrice
% NO2_p=pivottable(NO2_t,1,2,3,@sum); %tabella pivot creata date sulle righe e IdSensore sulle colonne
% Vuoti=cellfun('isempty',NO2_p); %seleziona le celle vuote
% NO2_p(Vuoti)={NaN};%sostituisce le celle vpute con NaN
% NO2_2021=cell2mat(NO2_p(:,2:end));
% 
% clear NO2_r NO2_s data NO2_t NO2_p idx Vuoti
% 
% %% Add 01/01/2015
% 
% NO2_2015 = [NO2_2015(1,:); NaN(1,width(NO2_2015)); NO2_2015(2:end,:)];
% 
% %% Join the data together
% 
% NO22015_2019=[NO2_2015; NO2_2016(2:end,:); NO2_2017(2:end,:); NO2_2018(2:end,:); NO2_2019(2:end,:); NO2_2020(2:end,:); NO2_2021(2:end,:)]; 
% 
% %% Substitute unavailable data
% % Gaps in the data of 5 days or less are filled through linear
% % interpolation of the available data from the same station. Longer gaps
% % are filled by applying annverse distance weighting interpolation,
% % considering the KNN stations.
% 
% NO22015_2019(NO22015_2019 < 0)=NaN; 
% 
% NO2_fill = fillmissing(NO22015_2019(2:end,:), 'linear', 'MaxGap', 1);
% 
% % NO2_fill = NO22015_2019(2:end,:);
% 
% Coord = Coord';
% 
% K = round(sqrt(length(Coord)));
% 
% for i = 1:length(Coord)
% 
%     idx = isnan(NO2_fill(:,i));
% 
%     if sum(idx) > 0
%         % Find KNN stations and distance
%         X = Coord(:,2:3);
%         X(i,:) = [];
%         Y = Coord(i,2:3);
%         [mIdx,mD] = knnsearch(X,Y,'K',K,'Distance','seuclidean');
%         dati_temp = [NO2_fill(:,mIdx)];
% 
%         % Weights
%         W = 1./(mD.^2);   
% 
%         % Interpolation
%         NO2_fill(idx,i) =  sum(W.*dati_temp(idx,:), 2, 'omitnan')/sum(W);
%     end
% end
% 
% NO22015_2019(2:end,:) = NO2_fill;
% 
% %% Create table with data
% 
% t1= datetime(2015,1,1,0,0,0);
% t2=datetime(2021,12,31,23,0,0);
% dn=(t1:days(1):t2).';
% dn(dn=='29-Feb-2016') = [];
% dn(dn=='29-Feb-2020') = [];
% dn = [NaT; dn];
% 
% NO22015_2019_T = table(dn, NO22015_2019);
% NO22015_2019_T = splitvars(NO22015_2019_T);
% 
% writetable(NO22015_2019_T, 'C:\Users\Luca\OneDrive - Politecnico di Milano\AgriAir\Tesi_Olgiati\Input\NO22015_2019_New.csv');
































%% Meteorology

% 
% Stazioni_meteo = readtable("Stazioni_Meteorologiche.csv");
% 
% idx = isnat(Stazioni_meteo{:,9});
% Stazioni_meteo = Stazioni_meteo(idx,:);
% 
% 
% idx_staz = strcmp('Precipitazione', Stazioni_meteo{:,2});
% Stazioni_Prec = Stazioni_meteo(idx_staz,:);
% idx_staz = strcmp('Velocità Vento', Stazioni_meteo{:,2});
% Stazioni_Wind = Stazioni_meteo(idx_staz,:);

%% Find common stations that measure all the data
% Stazioni = intersect(Stazioni_PM10{:,4}, Stazioni_NO2{:,4});
% Stazioni = intersect(Stazioni_PM10{:,4}, Stazioni_Prec{:,4});
% Stazioni = intersect(Stazioni_PM10{:,4}, Stazioni_Wind{:,4});

% Stazioni_inq = intersect(Stazioni_PM10{:,4}, Stazioni_NO2{:,4});
% Stazioni_m = intersect(Stazioni_Wind{:,4}, Stazioni_Prec{:,4});

% idx_staz = ismember(Stazioni_PM10{:,4}, Stazioni);
% Stazioni = Stazioni_PM10(idx_staz,:);




% dati_15 = readtable("Air_2015.csv");
% 
% idx_staz = strcmp('PM10 (SM2005)', Stazioni_aria{:,2});
% Stazioni_PM10 = Stazioni(idx_staz,[1,12,13]);
% 
% pm1015 = ismember(dati_15{:,1}, Stazioni_PM10{:,1});
% 
% pm1015 = dati_15(pm1015,:);








































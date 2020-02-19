%Ethan Green
%February 14th, 2020
%Pull charge data from BiGGmetData
%using parsed uncompartmentalized mets
%from LrGG_Model

%% Loading in relevant data
load MetParsed.mat 
load BiGGmetData_01_04_2017.mat
load LrGG_Model.mat

%% Pullling data and exporting to LrGG Model
[yn, id] = ismember(mets, BiGGmetData(:, 2));
model.metCharges(yn, 1) = str2double(BiGGmetData(id(yn), 7));
[yn, id] = ismember(mets, BiGGmetData(:, 2));
model.metFormulas(yn, 1) = BiGGmetData(id(yn), 6);
[yn, id] = ismember(mets, BiGGmetData(:, 2));
model.mets(yn, 1) = BiGGmetData(id(yn), 1);

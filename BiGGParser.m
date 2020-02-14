%Ethan Green
%February 14th, 2020
%Pull charge data from BiGGmetData
%using parsed uncompartmentalized mets
%from LrGG_Model

%% Loading in relevant data
load MetParsed.mat 
load BiGGmetData_01_04_2017.mat

%% Pullling charge data
metCharges = zeros(length(mets),1);
for i = 1:length(mets)
    ind = find(ismember(BiGGmetData(:,2),mets(i))==1);
    chargestr = (BiGGmetData(ind(1),7));
    metCharges(i) = str2double(chargestr{1});
end

%% Exporting data to LrGG Model
load LrGG_Model.mat
model.metCharges = metCharges;
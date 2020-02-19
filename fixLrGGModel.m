%Ethan Green
%February 14th, 2020
%Attempting to fix the LrGG Model

%% Loading in relevant data
clear
load MetParsed.mat 
load BiGGmetData_01_04_2017.mat
load LrGG_Model.mat

%% General model statistics
%Number of reactions
disp("The LrGG Model has "+length(model.rxns)+" reactions,")
%Number of metabolites 
disp(length(model.mets)+" metabolites,")
%Number of genes 
disp(length(model.genes)+" genes,")
%Number of transport reactions 
%disp(length(model.
%Number of extracellular metabolites

%% Fixing metCharges, metFormulas, and mets
[yn, id] = ismember(mets, BiGGmetData(:, 2));
model.metCharges(yn, 1) = str2double(BiGGmetData(id(yn), 7));
[yn, id] = ismember(mets, BiGGmetData(:, 2));
model.metFormulas(yn, 1) = BiGGmetData(id(yn), 6);
[yn, id] = ismember(mets, BiGGmetData(:, 2));
model.mets(yn, 1) = BiGGmetData(id(yn), 1);
disp("The model's Charges, metabolite names and formulas were standardized using the BiGG Database.")

%% Classifying R and X group containing formulae as unknowns
%getElementalComposition creates to matrices, elements lists all of the
%elements in the metabolite(s) examined and metEle lists the number of
%elements seen in the metabaolite(s) referencing the elements vector
[metEle, elements] = getElementalComposition(model.metFormulas);
%Using pulled data to identify any problematic 'R' or 'X'groups in a
%compound
genEle = ismember(elements, {'R', 'X'});
%Removing all metFormulas identified to have R
l = length(model.metFormulas(any(metEle(:, genEle), 2)));
disp(l+" metabolites with R or X groups had their formulae removed from the model.")
model.metFormulas(any(metEle(:, genEle), 2)) = {''};